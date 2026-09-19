"""
One rate-limited, retrying HTTP transport shared by every PROVESID web client.

Before this module each web-API client carried its own copy of "pause, request,
decide what the status code meant, maybe give up" --- four copies that had
drifted apart, none of which honoured ``Retry-After`` and only one of which
retried at all. :class:`HTTPClient` is the single place that decides *when to
ask again*; the clients keep everything that encodes their upstream's contract:
the URLs they build, the bodies they parse and the exceptions they raise.

The one thing services genuinely disagree about is what a response *means*.
PubChem answers a momentary overload with ``PUGVIEW.ServerBusy`` behind a 404,
so its status code alone is a lie; the NCI resolver returns plain text and a
bare 404 for absence. That disagreement is the ``classify`` callback --- a
function from a response to an :class:`Outcome` --- and it is the only part of
the policy a caller is expected to supply.

Example:
    >>> client = HTTPClient(min_interval=0.2, timeout=30)
    >>> client.get_text("https://example.org/thing")   # doctest: +SKIP
    'an answer'
"""

import email.utils
import logging
import threading
import time
from enum import Enum
from typing import Any, Callable, Dict, Optional, Type
from urllib.parse import urlsplit

import requests

logger = logging.getLogger(__name__)


class ServiceError(Exception):
    """
    Base for every failure this package reports from a web service.

    Each client raises its own subclass --- :class:`~provesid.pubchemview.PubChemViewError`,
    :class:`~provesid.resolver.NCIResolverError` and so on --- so a caller can
    catch one service or, through this base, all of them. A raw ``requests``
    exception never escapes :class:`HTTPClient`.

    The message is the whole of what most callers want, so it stays the single
    positional argument and ``str(exc)`` is unchanged. The response detail is
    keyword-only and defaults to None, which is what lets a client raise one of
    these by hand --- ``raise PubChemError("bad request")`` --- exactly as
    before. It exists because two services distinguish their failures by
    status: CAS Common Chemistry reports a rejected key as 401 and an unknown
    CAS number as 404, and both have to become different entries in the dict it
    returns.

    Args:
        message: The human-readable description.
        status_code: The HTTP status that caused the failure, when there was a
            response at all. None for a timeout or a connection error.
        url: The URL that was requested.
        response: The raw response, for a caller that needs to read the body.
            None when no response arrived.

    Attributes:
        status_code: As above.
        url: As above.
        response: As above.

    Example:
        >>> issubclass(NotFoundError, ServiceError)
        True
        >>> exc = ServiceError("nope", status_code=404)
        >>> str(exc), exc.status_code
        ('nope', 404)
        >>> ServiceError("by hand").status_code is None
        True
    """

    def __init__(self, message: str = "", *, status_code: Optional[int] = None,
                 url: Optional[str] = None,
                 response: Optional[requests.Response] = None):
        super().__init__(message)
        self.status_code = status_code
        self.url = url
        self.response = response


class NotFoundError(ServiceError):
    """
    The service answered, and its answer was that the record does not exist.

    This is a statement about the data, not about the request. It is never
    retried, because asking again cannot change it.

    Example:
        >>> issubclass(NotFoundError, ServiceError)
        True
    """
    pass


class RateLimitError(ServiceError):
    """
    The service kept throttling the request until the retry budget ran out.

    Raised only when every attempt was refused with a rate-limit response; a
    429 that clears on a later attempt is invisible to the caller.

    Example:
        >>> issubclass(RateLimitError, ServiceError)
        True
    """
    pass


class ServiceTimeoutError(ServiceError):
    """
    Every attempt timed out or the connection could not be made.

    Example:
        >>> issubclass(ServiceTimeoutError, ServiceError)
        True
    """
    pass


class Outcome(Enum):
    """
    What a response means, once a service-specific classifier has read it.

    Attributes:
        OK: Use the response.
        ABSENT: The record does not exist --- raise the client's not-found
            exception without retrying.
        RETRY: A transient condition; ask again after backing off.
        FATAL: A permanent error that is not absence --- a malformed request,
            a rejected key. Retrying cannot help.

    Example:
        >>> Outcome.RETRY.name
        'RETRY'
    """
    OK = "ok"
    ABSENT = "absent"
    RETRY = "retry"
    FATAL = "fatal"


#: Status codes that mean "the service is momentarily unwilling", whatever
#: else the body says. 429 is explicit throttling; 5xx is the service failing
#: on its own side.
RETRYABLE_STATUS = frozenset({429, 500, 502, 503, 504, 507, 509})


def default_classify(response: requests.Response) -> Outcome:
    """
    Classify a response by its HTTP status alone.

    The right reading for any service that uses status codes honestly: 2xx is
    an answer, 404 is absence, 429 and 5xx are worth retrying, and any other
    4xx is a permanent error in the request itself. Services that describe
    errors in the body --- PubChem does --- supply their own classifier
    instead.

    Args:
        response: The response to classify.

    Returns:
        The :class:`Outcome` for this response.

    Example:
        >>> class R: status_code = 404
        >>> default_classify(R()).name
        'ABSENT'
        >>> class R: status_code = 503
        >>> default_classify(R()).name
        'RETRY'
    """
    status = response.status_code
    if 200 <= status < 300:
        return Outcome.OK
    if status in RETRYABLE_STATUS:
        return Outcome.RETRY
    if status == 404:
        return Outcome.ABSENT
    if 400 <= status < 500:
        return Outcome.FATAL
    return Outcome.RETRY


def retry_after_seconds(response: requests.Response) -> Optional[float]:
    """
    Read a ``Retry-After`` header, in either form the standard allows.

    RFC 9110 permits a number of seconds or an HTTP date. A service that
    troubles itself to say when to come back knows better than any backoff
    curve, so :class:`HTTPClient` prefers this over its own schedule.

    Args:
        response: The response whose headers to read.

    Returns:
        The wait in seconds, or None when the header is absent, unparseable
        or in the past.

    Example:
        >>> class R: headers = {"Retry-After": "12"}
        >>> retry_after_seconds(R())
        12.0
        >>> class R: headers = {}
        >>> retry_after_seconds(R()) is None
        True
    """
    try:
        raw = response.headers.get("Retry-After")
    except AttributeError:
        return None
    if not raw:
        return None

    raw = str(raw).strip()
    try:
        return max(0.0, float(raw))
    except ValueError:
        pass

    try:
        when = email.utils.parsedate_to_datetime(raw)
    except (TypeError, ValueError):
        return None
    if when is None:
        return None
    delay = when.timestamp() - time.time()
    return delay if delay > 0 else None


class RateLimiter:
    """
    The clock a host's requests are paced against.

    A client's ``min_interval`` is its own promise about how fast it will ask.
    The clock it measures that promise against belongs to the *host*, because
    the limit being respected does too: PubChem publishes five requests per
    second **per IP**, not per Python object. Two clients aimed at PubChem in
    one process --- a :class:`~provesid.pubchem.PubChemAPI` and a
    :class:`~provesid.pubchemview.PubChemView`, which is the ordinary way to
    use this package --- each kept their own clock before this class existed,
    so each could believe it was pacing correctly while together they asked
    twice as fast as PubChem allows.

    Sharing the clock means a request waits ``min_interval`` after the last
    request *anyone* made to that host. The lock is held across the sleep, so
    threads queue rather than all waking at once.

    Attributes:
        last_request_time: When any client last asked this host, as a Unix
            timestamp; 0.0 before the first request.

    Example:
        >>> limiter = RateLimiter()
        >>> limiter.last_request_time
        0.0
        >>> _ = limiter.wait(0.0)       # pacing off: returns at once
        >>> limiter.last_request_time > 0
        True
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self.last_request_time = 0.0

    def wait(self, min_interval: float) -> float:
        """
        Sleep until ``min_interval`` has passed since this host was last asked.

        Args:
            min_interval: Seconds the caller promises to leave between
                requests. 0 or less disables the wait but still records the
                request, because the request is happening either way and the
                next caller needs to know when.

        Returns:
            The time at which the caller may proceed, as a Unix timestamp.

        Example:
            >>> limiter = RateLimiter()
            >>> first = limiter.wait(0.05)
            >>> limiter.wait(0.05) - first >= 0.045
            True
        """
        with self._lock:
            if min_interval > 0:
                elapsed = time.time() - self.last_request_time
                if elapsed < min_interval:
                    time.sleep(min_interval - elapsed)
            self.last_request_time = time.time()
            return self.last_request_time


#: Every host a limiter has been asked for, so that clients aimed at the same
#: service find the same clock. Keyed by lower-cased ``host:port``. Entries are
#: never removed --- there are a handful of them and each is two floats.
_host_limiters: Dict[str, RateLimiter] = {}
_host_limiters_lock = threading.Lock()


def host_limiter(url_or_host: str) -> RateLimiter:
    """
    Return the process-wide :class:`RateLimiter` for one host, creating it once.

    Args:
        url_or_host: A full URL, whose host is used, or a bare host. Accepting
            either is deliberate: a client already holds its service's base
            URL, so it can pass that and needs no second piece of
            configuration to get its pacing shared correctly.

    Returns:
        The limiter for that host. Two calls naming the same host return the
        same object.

    Example:
        >>> a = host_limiter("https://pubchem.ncbi.nlm.nih.gov/rest/pug")
        >>> b = host_limiter("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view")
        >>> a is b
        True
        >>> a is host_limiter("https://www.ebi.ac.uk/chebi")
        False
    """
    parsed = urlsplit(url_or_host)
    host = (parsed.netloc or parsed.path or url_or_host).strip().lower()
    with _host_limiters_lock:
        limiter = _host_limiters.get(host)
        if limiter is None:
            limiter = _host_limiters[host] = RateLimiter()
        return limiter


class HTTPClient:
    """
    Rate-limited HTTP client with retry and back-off, shared by every
    PROVESID web-API module.

    One instance belongs to one client object and paces that client's own
    requests; it holds no state beyond the time of its last request and the
    policy it was given.

    Args:
        min_interval: Minimum seconds between two requests from this client.
            0 disables pacing.
        timeout: Default per-request timeout in seconds.
        max_retries: Retries *after* the first attempt, so a request is made
            at most ``max_retries + 1`` times.
        backoff: Base for the exponential wait, ``backoff * 2 ** attempt``
            seconds. 0 retries immediately, which is what tests want.
        max_backoff: Ceiling on any single wait, including one a
            ``Retry-After`` header asks for.
        max_elapsed: Ceiling on the *total* time spent waiting between
            attempts. Retrying stops once the next wait would take the sum past
            it, whatever ``max_retries`` allows. None means ``max_retries`` is
            the only bound.

            This exists because a service that says ``Retry-After: 30`` --- as
            PubChem does when it throttles an IP --- turns three retries into a
            ninety-second call. Dropping the retry would be worse: a caller
            resolving ten thousand compounds loses one to every transient 503.
            One wait the service itself asked for recovers most of them; the
            budget is what stops the rest of the curve from being charged to a
            caller who is waiting at a prompt.
        headers: Headers sent with every request. Omitted entirely when None,
            so a stub that accepts only ``(url, timeout=...)`` still works.
        classify: Maps a response to an :class:`Outcome`. Defaults to
            :func:`default_classify`.
        error_cls: Raised for a fatal response --- a malformed request, a
            rejected key --- and, unless ``retry_exhausted_cls`` says
            otherwise, for an exhausted retry budget.
        not_found_cls: Raised for :attr:`Outcome.ABSENT`.
        timeout_cls: Raised when every attempt timed out or could not connect.
            Defaults to ``error_cls``.
        rate_limit_cls: Raised when every attempt was throttled. Defaults to
            ``error_cls``.
        retry_exhausted_cls: Raised when a transient condition that was neither
            a throttle nor a timeout outlived the retry budget --- a service
            that stayed busy. Defaults to ``error_cls``. PubChem passes its
            ``PubChemServerError`` here, because "the service kept failing" and
            "the request was wrong" are different things to its callers.
        session: A ``requests.Session`` to make the calls through, for
            connection pooling and persistent headers. When None the module
            functions ``requests.get`` / ``requests.post`` are called, which is
            what lets a test stub them.
        pace_host: The service this client shares its pacing clock with, as a
            host or as the base URL it already holds. Given, the client waits
            ``min_interval`` after the last request *any* client made to that
            host --- the only way to honour a limit expressed per IP, such as
            PubChem's five per second. Omitted, the client paces alone, which
            is right for a stub and for a service with no shared budget.
        logger: Logger for the DEBUG line per request and the WARNING per
            retry. Defaults to this module's logger.

    Example:
        >>> client = HTTPClient(min_interval=0.2, timeout=10, max_retries=2)
        >>> client.get_json("https://example.org/data.json")   # doctest: +SKIP
        {'ok': True}
    """

    def __init__(self, *, min_interval: float = 0.0, timeout: float = 30,
                 max_retries: int = 3, backoff: float = 1.0,
                 max_backoff: float = 60.0,
                 max_elapsed: Optional[float] = None,
                 headers: Optional[Dict[str, str]] = None,
                 classify: Callable[[requests.Response], Outcome] = default_classify,
                 error_cls: Type[Exception] = ServiceError,
                 not_found_cls: Type[Exception] = NotFoundError,
                 timeout_cls: Optional[Type[Exception]] = None,
                 rate_limit_cls: Optional[Type[Exception]] = None,
                 retry_exhausted_cls: Optional[Type[Exception]] = None,
                 session: Optional[requests.Session] = None,
                 pace_host: Optional[str] = None,
                 logger: Optional[logging.Logger] = None):
        self.min_interval = min_interval
        self.timeout = timeout
        self.max_retries = max_retries
        self.backoff = backoff
        self.max_backoff = max_backoff
        self.max_elapsed = max_elapsed
        self.headers = headers
        self.classify = classify
        self.error_cls = error_cls
        self.not_found_cls = not_found_cls
        self.timeout_cls = timeout_cls or error_cls
        self.rate_limit_cls = rate_limit_cls or error_cls
        self.retry_exhausted_cls = retry_exhausted_cls or error_cls
        self.session = session
        self.logger = logger or logging.getLogger(__name__)
        self.last_request_time = 0.0
        self.limiter = host_limiter(pace_host) if pace_host else RateLimiter()

    def rate_limit(self) -> None:
        """
        Sleep, if needed, so requests to this service stay ``min_interval``
        apart.

        Called by :meth:`request` before every attempt --- including retries,
        which is the point: a service that is shedding load should not be
        asked again faster than a service that is not.

        The interval is this client's own; the clock it is measured against
        belongs to :attr:`limiter`, which is shared with every other client
        aimed at the same host when ``pace_host`` was given. So the wait is
        ``min_interval`` since *anybody* last asked that service, not since
        this object did.

        Example:
            >>> client = HTTPClient(min_interval=0.0)
            >>> client.rate_limit()     # returns at once when pacing is off
        """
        self.last_request_time = self.limiter.wait(self.min_interval)

    def _send(self, method: str, url: str, *, params=None, data=None,
              json=None, headers=None, timeout=None,
              stream: bool = False) -> requests.Response:
        """
        Make one HTTP call, passing only the arguments that were actually
        given.

        Every optional argument is omitted from the call when it is None, so
        the request reads ``requests.get(url, timeout=30)`` in the common
        case. Several test suites stub ``requests.get`` with exactly that
        signature; a client that always passed ``params=None, headers=None``
        would break them for no gain.

        The call goes through :attr:`session` when the client was given one,
        and otherwise through the ``requests`` module functions. Which of the
        two is visible from outside: a session's ``get`` carries the session's
        persistent headers and pooled connection, and is patched as
        ``requests.Session.get``.

        Args:
            method: ``"GET"`` or ``"POST"``.
            url: The full URL.
            params: Query-string parameters.
            data: Form or raw body for a POST.
            json: JSON body for a POST.
            headers: Per-request headers, merged over the client's own.
            timeout: Overrides the client's timeout.
            stream: Passed through to ``requests`` for large downloads.

        Returns:
            The raw response.

        Raises:
            ValueError: If ``method`` is neither GET nor POST.
        """
        kwargs: Dict[str, Any] = {"timeout": self.timeout if timeout is None else timeout}
        merged = {**self.headers, **headers} if self.headers and headers else (headers or self.headers)
        if merged:
            kwargs["headers"] = merged
        if params is not None:
            kwargs["params"] = params
        if data is not None:
            kwargs["data"] = data
        if json is not None:
            kwargs["json"] = json
        if stream:
            kwargs["stream"] = True

        caller = self.session if self.session is not None else requests
        verb = method.upper()
        if verb == "GET":
            return caller.get(url, **kwargs)
        if verb == "POST":
            return caller.post(url, **kwargs)
        raise ValueError(f"Unsupported HTTP method: {method}")

    def _wait(self, attempt: int, response: Optional[requests.Response]) -> float:
        """
        Decide how long to wait before the next attempt.

        A ``Retry-After`` the service sent wins over the exponential curve,
        because it is the service's own estimate. Either way the wait is
        capped at ``max_backoff``.

        Args:
            attempt: Zero-based index of the attempt that just failed.
            response: The response that failed, or None for a timeout or
                connection error.

        Returns:
            Seconds to sleep.
        """
        if response is not None:
            asked = retry_after_seconds(response)
            if asked is not None:
                return min(asked, self.max_backoff)
        return min(self.backoff * (2 ** attempt), self.max_backoff)

    def request(self, method: str, url: str, **kwargs: Any) -> requests.Response:
        """
        Make a request, retrying transient failures, and return the response.

        Args:
            method: ``"GET"`` or ``"POST"``.
            url: The full URL.
            **kwargs: ``params``, ``data``, ``json``, ``headers``, ``timeout``
                and ``stream``, all optional --- see :meth:`_send`.

        Returns:
            The response, already classified as :attr:`Outcome.OK`.

        Raises:
            not_found_cls: The service reported the record as absent.
            error_cls: A permanent error.
            timeout_cls: Every attempt timed out or failed to connect.
            rate_limit_cls: Every attempt was throttled.
            retry_exhausted_cls: A transient condition outlived the retry
                budget, either in attempts or in ``max_elapsed`` seconds.
                Defaults to ``error_cls``.

        Every one of those carries the status, the URL and the response on the
        exception when there was a response --- see :class:`ServiceError` --- so
        a client can tell a rejected key from an unknown record without
        re-reading the wire.

        Example:
            >>> client = HTTPClient()
            >>> client.request("GET", "https://example.org/x").text   # doctest: +SKIP
            'ok'
        """
        last_error: Optional[str] = None
        last_status: Optional[int] = None
        timed_out = False
        throttled = False
        waited = 0.0
        out_of_time = False

        for attempt in range(self.max_retries + 1):
            self.rate_limit()
            response: Optional[requests.Response] = None

            try:
                self.logger.debug(f"{method.upper()} {url}")
                response = self._send(method, url, **kwargs)
            except requests.Timeout as exc:
                timed_out = True
                last_error = f"request timed out: {exc}"
            except requests.ConnectionError as exc:
                timed_out = True
                last_error = f"connection failed: {exc}"
            except requests.RequestException as exc:
                # Anything else requests can raise --- a malformed URL, a
                # broken redirect chain. Not worth a retry, and it must not
                # reach the caller as a requests exception.
                raise self._fail(self.error_cls, f"Request to {url} failed: {exc}",
                                 url=url) from exc

            if response is not None:
                verdict = self.classify(response)
                last_status = response.status_code

                if verdict is Outcome.OK:
                    return response

                if verdict is Outcome.ABSENT:
                    self.logger.debug(f"absent: {last_status} for {url}")
                    raise self._fail(
                        self.not_found_cls,
                        f"No data for {url} (HTTP {last_status})",
                        status_code=last_status, url=url, response=response,
                    )

                if verdict is Outcome.FATAL:
                    raise self._fail(
                        self.error_cls,
                        f"HTTP {last_status} for {url}: {self._body_excerpt(response)}",
                        status_code=last_status, url=url, response=response,
                    )

                throttled = last_status == 429
                last_error = f"HTTP {last_status}"

            if attempt == self.max_retries:
                break

            wait = self._wait(attempt, response)
            if self.max_elapsed is not None and waited + wait > self.max_elapsed:
                # The service is willing to be asked again, just not soon
                # enough to be worth the caller's time.
                out_of_time = True
                self.logger.warning(
                    f"{last_error} for {url}; giving up rather than waiting "
                    f"another {wait:.1f}s on top of {waited:.1f}s "
                    f"(max_elapsed={self.max_elapsed:g}s)"
                )
                break

            self.logger.warning(
                f"{last_error} for {url}; retrying in {wait:.1f}s "
                f"(attempt {attempt + 2} of {self.max_retries + 1})"
            )
            if wait > 0:
                time.sleep(wait)
                waited += wait

        if out_of_time:
            # "Spent" would be a lie when the budget stopped the very first
            # wait, which is the usual case against a service asking for more
            # than the whole budget --- so say how much was actually used.
            message = (f"Request to {url} failed; stopped retrying after "
                       f"{waited:.0f}s of its {self.max_elapsed:g}s retry "
                       f"budget: {last_error}")
        else:
            attempts = self.max_retries + 1
            message = (f"Request to {url} failed after {attempts} attempt(s): "
                       f"{last_error}")
        if throttled:
            failure = self.rate_limit_cls
        elif timed_out:
            failure = self.timeout_cls
        else:
            failure = self.retry_exhausted_cls
        raise self._fail(failure, message, status_code=last_status, url=url,
                         response=response)

    @staticmethod
    def _fail(cls: Type[Exception], message: str, *,
              status_code: Optional[int] = None, url: Optional[str] = None,
              response: Optional[requests.Response] = None) -> Exception:
        """
        Build the exception to raise, with the response detail when it fits.

        Every exception class in this package descends from
        :class:`ServiceError` and so accepts the keyword detail. A caller is
        free to pass a plain ``Exception`` subclass as ``error_cls``, though,
        and handing that one keywords it never declared would turn a service
        failure into a ``TypeError``. So the detail is attached only when the
        class is known to take it.

        Args:
            cls: The exception class to instantiate.
            message: The message.
            status_code: The HTTP status, when there was a response.
            url: The URL requested.
            response: The raw response, when one arrived.

        Returns:
            The exception, not yet raised.

        Example:
            >>> HTTPClient._fail(ServiceError, "busy", status_code=503).status_code
            503
            >>> isinstance(HTTPClient._fail(ValueError, "busy", status_code=503), ValueError)
            True
        """
        if issubclass(cls, ServiceError):
            return cls(message, status_code=status_code, url=url,
                       response=response)
        return cls(message)

    @staticmethod
    def _body_excerpt(response: requests.Response, limit: int = 200) -> str:
        """
        Return the beginning of a response body for an error message.

        Args:
            response: The response to read.
            limit: Maximum characters to return.

        Returns:
            The first ``limit`` characters of the body, or ``''`` when it
            cannot be read.
        """
        try:
            return response.text[:limit]
        except Exception:
            return ""

    def get(self, url: str, **kwargs: Any) -> requests.Response:
        """
        GET a URL and return the raw response.

        Args:
            url: The full URL.
            **kwargs: See :meth:`request`.

        Returns:
            The response.

        Example:
            >>> HTTPClient().get("https://example.org/img.png").content   # doctest: +SKIP
            b'\\x89PNG...'
        """
        return self.request("GET", url, **kwargs)

    def post(self, url: str, **kwargs: Any) -> requests.Response:
        """
        POST to a URL and return the raw response.

        Args:
            url: The full URL.
            **kwargs: See :meth:`request`.

        Returns:
            The response.

        Example:
            >>> HTTPClient().post("https://example.org/q", data={"cid": 2244})   # doctest: +SKIP
            <Response [200]>
        """
        return self.request("POST", url, **kwargs)

    def get_text(self, url: str, **kwargs: Any) -> str:
        """
        GET a URL and return its body as stripped text.

        Args:
            url: The full URL.
            **kwargs: See :meth:`request`.

        Returns:
            The response body, with surrounding whitespace removed.

        Example:
            >>> HTTPClient().get_text("https://example.org/smiles")   # doctest: +SKIP
            'CCO'
        """
        return self.request("GET", url, **kwargs).text.strip()

    def get_json(self, url: str, **kwargs: Any) -> Any:
        """
        GET a URL and return its parsed JSON body.

        Args:
            url: The full URL.
            **kwargs: See :meth:`request`.

        Returns:
            The decoded JSON.

        Raises:
            error_cls: The body is not valid JSON. A service that answers 200
                with something unparseable has failed as surely as one that
                answers 500, and the caller should hear about it in the same
                way.

        Example:
            >>> HTTPClient().get_json("https://example.org/data.json")   # doctest: +SKIP
            {'ok': True}
        """
        return self.decode_json(self.request("GET", url, **kwargs), url)

    def post_json(self, url: str, **kwargs: Any) -> Any:
        """
        POST to a URL and return its parsed JSON body.

        Args:
            url: The full URL.
            **kwargs: See :meth:`request`.

        Returns:
            The decoded JSON.

        Raises:
            error_cls: The body is not valid JSON.

        Example:
            >>> HTTPClient().post_json("https://example.org/q", json={"n": 1})   # doctest: +SKIP
            {'ok': True}
        """
        return self.decode_json(self.request("POST", url, **kwargs), url)

    def decode_json(self, response: requests.Response,
                    url: Optional[str] = None) -> Any:
        """
        Parse a response body as JSON, reporting failure as a service error.

        Public because a client that decides for itself whether a body is JSON
        --- ChEBI reads the ``Content-Type``, because it serves molfiles and SVG
        from the same API as its records --- needs the same failure reported the
        same way.

        Args:
            response: The response to decode.
            url: The URL, for the error message. Falls back to the response's
                own when it has one.

        Returns:
            The decoded JSON.

        Raises:
            error_cls: The body is not valid JSON.

        Example:
            >>> class R:
            ...     text = 'not json'
            ...     def json(self): raise ValueError("nope")
            >>> HTTPClient().decode_json(R(), "https://example.org/x")
            Traceback (most recent call last):
                ...
            provesid.http.ServiceError: Response from https://example.org/x is not JSON: not json
        """
        if url is None:
            url = getattr(response, "url", "") or "the service"
        try:
            return response.json()
        except ValueError as exc:
            raise self._fail(
                self.error_cls,
                f"Response from {url} is not JSON: {self._body_excerpt(response)}",
                status_code=getattr(response, "status_code", None), url=url,
                response=response,
            ) from exc
