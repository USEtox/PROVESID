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
import time
from enum import Enum
from typing import Any, Callable, Dict, Optional, Type

import requests

logger = logging.getLogger(__name__)


class ServiceError(Exception):
    """
    Base for every failure this package reports from a web service.

    Each client raises its own subclass --- :class:`~provesid.pubchemview.PubChemViewError`,
    :class:`~provesid.resolver.NCIResolverError` and so on --- so a caller can
    catch one service or, through this base, all of them. A raw ``requests``
    exception never escapes :class:`HTTPClient`.

    Example:
        >>> issubclass(NotFoundError, ServiceError)
        True
    """
    pass


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
        headers: Headers sent with every request. Omitted entirely when None,
            so a stub that accepts only ``(url, timeout=...)`` still works.
        classify: Maps a response to an :class:`Outcome`. Defaults to
            :func:`default_classify`.
        error_cls: Raised for a fatal response and for an exhausted retry
            budget.
        not_found_cls: Raised for :attr:`Outcome.ABSENT`.
        timeout_cls: Raised when every attempt timed out or could not connect.
            Defaults to ``error_cls``.
        rate_limit_cls: Raised when every attempt was throttled. Defaults to
            ``error_cls``.
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
                 headers: Optional[Dict[str, str]] = None,
                 classify: Callable[[requests.Response], Outcome] = default_classify,
                 error_cls: Type[Exception] = ServiceError,
                 not_found_cls: Type[Exception] = NotFoundError,
                 timeout_cls: Optional[Type[Exception]] = None,
                 rate_limit_cls: Optional[Type[Exception]] = None,
                 logger: Optional[logging.Logger] = None):
        self.min_interval = min_interval
        self.timeout = timeout
        self.max_retries = max_retries
        self.backoff = backoff
        self.max_backoff = max_backoff
        self.headers = headers
        self.classify = classify
        self.error_cls = error_cls
        self.not_found_cls = not_found_cls
        self.timeout_cls = timeout_cls or error_cls
        self.rate_limit_cls = rate_limit_cls or error_cls
        self.logger = logger or logging.getLogger(__name__)
        self.last_request_time = 0.0

    def rate_limit(self) -> None:
        """
        Sleep, if needed, so this client's requests stay ``min_interval``
        apart.

        Called by :meth:`request` before every attempt --- including retries,
        which is the point: a service that is shedding load should not be
        asked again faster than a service that is not.

        Example:
            >>> client = HTTPClient(min_interval=0.0)
            >>> client.rate_limit()     # returns at once when pacing is off
        """
        if self.min_interval <= 0:
            self.last_request_time = time.time()
            return
        elapsed = time.time() - self.last_request_time
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self.last_request_time = time.time()

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

        verb = method.upper()
        if verb == "GET":
            return requests.get(url, **kwargs)
        if verb == "POST":
            return requests.post(url, **kwargs)
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
            error_cls: A permanent error, or a transient one that outlived the
                retry budget.
            timeout_cls: Every attempt timed out or failed to connect.
            rate_limit_cls: Every attempt was throttled.

        Example:
            >>> client = HTTPClient()
            >>> client.request("GET", "https://example.org/x").text   # doctest: +SKIP
            'ok'
        """
        last_error: Optional[str] = None
        last_status: Optional[int] = None
        timed_out = False
        throttled = False

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
                raise self.error_cls(f"Request to {url} failed: {exc}") from exc

            if response is not None:
                verdict = self.classify(response)
                last_status = response.status_code

                if verdict is Outcome.OK:
                    return response

                if verdict is Outcome.ABSENT:
                    self.logger.debug(f"absent: {last_status} for {url}")
                    raise self.not_found_cls(
                        f"No data for {url} (HTTP {last_status})"
                    )

                if verdict is Outcome.FATAL:
                    raise self.error_cls(
                        f"HTTP {last_status} for {url}: {self._body_excerpt(response)}"
                    )

                throttled = last_status == 429
                last_error = f"HTTP {last_status}"

            if attempt == self.max_retries:
                break

            wait = self._wait(attempt, response)
            self.logger.warning(
                f"{last_error} for {url}; retrying in {wait:.1f}s "
                f"(attempt {attempt + 2} of {self.max_retries + 1})"
            )
            if wait > 0:
                time.sleep(wait)

        attempts = self.max_retries + 1
        message = f"Request to {url} failed after {attempts} attempt(s): {last_error}"
        if throttled:
            raise self.rate_limit_cls(message)
        if timed_out:
            raise self.timeout_cls(message)
        raise self.error_cls(message)

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
        return self._decode(self.request("GET", url, **kwargs), url)

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
        return self._decode(self.request("POST", url, **kwargs), url)

    def _decode(self, response: requests.Response, url: str) -> Any:
        """
        Parse a response body as JSON, reporting failure as a service error.

        Args:
            response: The response to decode.
            url: The URL, for the error message.

        Returns:
            The decoded JSON.

        Raises:
            error_cls: The body is not valid JSON.
        """
        try:
            return response.json()
        except ValueError as exc:
            raise self.error_cls(
                f"Response from {url} is not JSON: {self._body_excerpt(response)}"
            ) from exc
