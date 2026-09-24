"""
The ClassyFire web service (Wishart lab):
[`ClassyFireAPI`][provesid.classyfire.ClassyFireAPI].

.. warning::
   The service has classified nothing new since February 2023. It still
   serves queries submitted before then, but a new submission is never
   processed. For ChEBI classes computed offline use
   [`provesid.taxonomy.ChebifierClassifier`][provesid.taxonomy.ChebifierClassifier].

The three raw calls still use ``requests`` directly, with no pacing and no
retry, and return `requests.Response` objects rather than parsed data. Plan
step 7 replaces them; until then they are documented as they are.
[`get_classification`][provesid.classyfire.ClassyFireAPI.get_classification]
is the parsed wrapper, on the shared transport.

The server cuts every response off after about 108 KB, while still
answering HTTP 200: a whole query's JSON or CSV never arrives complete, and
``curl`` reports "transfer closed with outstanding read data remaining".
Ask for it a page at a time (``get_query(..., page=, per_page=)``, or
``get_classification``).

Examples:
    >>> from provesid.classyfire import ClassyFireAPI
    >>> result = ClassyFireAPI.get_classification(1)           # doctest: +SKIP
    >>> result["classification_status"], len(result["entities"])  # doctest: +SKIP
    ('Done', 655)
"""

import logging

import requests
from .cache import cached, clear_cache, get_cache_info
from .http import HTTPClient, NotFoundError, ServiceError

logger = logging.getLogger(__name__)


class ClassyFireError(ServiceError):
    """
    A ClassyFire request failed, including a response the server cut off.

    Raised by
    [`get_classification`][provesid.classyfire.ClassyFireAPI.get_classification];
    the three raw calls return None instead.

    Examples:
        >>> issubclass(ClassyFireError, ServiceError)
        True
    """


class ClassyFireNotFoundError(ClassyFireError, NotFoundError):
    """
    ClassyFire has no query with this ID.

    Examples:
        >>> issubclass(ClassyFireNotFoundError, NotFoundError)
        True
    """


class ClassyFireAPI:
    """
    Class to interact with the ClassyFire API. The class is converted from the original Ruby code
    provided at https://bitbucket.org/wishartlab/classyfire_api/src/master/lib/classyfire_api.rb

    Every method is static. Submit a structure with
    [`submit_query`][provesid.classyfire.ClassyFireAPI.submit_query], poll
    [`query_status`][provesid.classyfire.ClassyFireAPI.query_status], then
    fetch the result with
    [`get_query`][provesid.classyfire.ClassyFireAPI.get_query]. See the module
    warning: new submissions are no longer processed.

    Examples:
        >>> response = ClassyFireAPI.submit_query(
        ...     "Example Query", "C1=CC(=CC=C1[N+](=O)[O-])Cl")    # doctest: +SKIP
        >>> query_id = response.json()["id"]                       # doctest: +SKIP
        >>> ClassyFireAPI.query_status(query_id).text             # doctest: +SKIP
        'Done'
    """
    URL = 'http://classyfire.wishartlab.com'

    TIMEOUT = 60
    """Seconds the raw calls wait for an answer. They had no timeout, so a
    stalled connection waited forever."""

    PAGE_SIZE = 10
    """Entities per page for
    [`get_classification`][provesid.classyfire.ClassyFireAPI.get_classification].
    The server cuts responses off after about 108 KB. On query 1 (2026-09-24)
    10 entities made 48 KB and 25 made 95 KB, while 50 were cut off."""

    MIN_INTERVAL = 2.5
    """Seconds between two requests from
    [`get_classification`][provesid.classyfire.ClassyFireAPI.get_classification].
    Measured on 2026-09-24: three requests in quick succession pass and the
    fourth gets HTTP 429, with no ``Retry-After``; ten requests 2 s apart all
    passed. At 0.5 s, fetching query 1 drew 52 429s."""

    _http = HTTPClient(
        min_interval=MIN_INTERVAL,
        timeout=TIMEOUT,
        pace_host=URL,
        error_cls=ClassyFireError,
        not_found_cls=ClassyFireNotFoundError,
        headers={"Accept": "application/json"},
        logger=logger,
    )

    @staticmethod
    def clear_cache():
        """
        Delete every cached ClassyFire response, in memory and on disk.

        Examples:
            >>> ClassyFireAPI.clear_cache()
            >>> ClassyFireAPI.get_cache_info()['file_count']
            0
        """
        clear_cache(service='classyfire')

    @staticmethod
    def get_cache_info():
        """
        Size and location of the ClassyFire cache.

        Returns:
            The statistics
            [`provesid.cache.get_cache_info`][provesid.cache.get_cache_info]
            reports for
                ``service='classyfire'``.

        Examples:
            >>> ClassyFireAPI.get_cache_info()['cache_directory'].endswith('classyfire')
            True
        """
        return get_cache_info(service='classyfire')

    @staticmethod
    @cached(service='classyfire')
    def submit_query(label, input, type='STRUCTURE', use_cache=True):
        """
        Submit a structure for classification.

        Args:
            label: A name for the query, kept by ClassyFire.
            input: The structure, as SMILES (or several, one per line).
            type: ClassyFire's ``query_type``; ``"STRUCTURE"`` by default.
            use_cache: When False, do not read the cache.

        Returns:
            (requests.Response): The service's answer, whose JSON carries the
            query ``id`` on success. An HTTP error comes back as its response,
            not raised; a connection failure comes back as None.

        Note:
            New submissions have not been processed since February 2023.

        Examples:
            >>> response = ClassyFireAPI.submit_query("test", "CCO")  # doctest: +SKIP
            >>> response.json()["id"]                                  # doctest: +SKIP
        """
        try:
            response = requests.post(
                f"{ClassyFireAPI.URL}/queries",
                json={"label": label, "query_input": input, "query_type": type},
                headers={"Accept": "application/json", "Content-Type": "application/json"},
                timeout=ClassyFireAPI.TIMEOUT,
            )
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            return e.response
        except requests.exceptions.RequestException as e:
            logger.warning("ClassyFire submit_query failed: %s", e)
            return e.response
        return response

    @staticmethod
    @cached(service='classyfire')
    def query_status(query_id, use_cache=True):
        """
        Retrieves the status of a query.

        Args:
            query_id: The ID
                [`submit_query`][provesid.classyfire.ClassyFireAPI.submit_query]
                returned.
            use_cache: When False, do not read the cache.

        Returns:
            (requests.Response): The status response, or None when the request
            failed for any reason (logged at WARNING). Its body is plain
            text such as ``Done``, not JSON: read ``.text``.

        Note:
            The status is cached like everything else here, so polling with
            the cache on sees the first answer forever; poll with
            ``use_cache=False``.

        Examples:
            >>> ClassyFireAPI.query_status(1, use_cache=False).status_code  # doctest: +SKIP
            200
        """
        try:
            response = requests.get(
                f"{ClassyFireAPI.URL}/queries/{query_id}/status.json",
                headers={"Accept": "application/json", "Content-Type": "application/json"},
                timeout=ClassyFireAPI.TIMEOUT,
            )
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as e:
            logger.warning("ClassyFire query_status(%s) failed: %s", query_id, e)
            return None

    @staticmethod
    @cached(service='classyfire')
    def get_query(query_id, format="json", page=None, per_page=None, use_cache=True):
        """
        Retrieve a query's classification.

        Args:
            query_id: The query ID.
            format: ``"json"``, ``"sdf"`` or ``"csv"``.
            page: Which page of entities, from 1. None asks for the whole
                query, which the server cuts off beyond about 108 KB.
            per_page: Entities per page. The server's default is 100, which
                is too many: see
                [`PAGE_SIZE`][provesid.classyfire.ClassyFireAPI.PAGE_SIZE].
            use_cache: When False, do not read the cache.

        Returns:
            (requests.Response): The result in the requested format. An HTTP
            error comes back as its response, not raised. A connection
            failure, a timeout or a response cut off mid-body comes back as
            None, logged at WARNING.

        Raises:
            ValueError: If ``format`` is not one of the three.

        Examples:
            >>> page = ClassyFireAPI.get_query(1, page=1, per_page=10).json()  # doctest: +SKIP
            >>> page["classification_status"], page["number_of_pages"]         # doctest: +SKIP
            ('Done', 66)
        """
        accept = {
            "json": "application/json",
            "sdf": "chemical/x-mdl-sdfile",
            "csv": "text/csv",
        }
        if format not in accept:
            raise ValueError("Invalid format. Must be one of json, sdf, or csv.")
        params = {"page": page, "per_page": per_page}
        try:
            response = requests.get(
                f"{ClassyFireAPI.URL}/queries/{query_id}.{format}",
                headers={"Accept": accept[format]},
                params={key: value for key, value in params.items() if value is not None},
                timeout=ClassyFireAPI.TIMEOUT,
            )
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            return e.response
        except requests.exceptions.RequestException as e:
            # Query 1 in full: ChunkedEncodingError, "Connection broken:
            # IncompleteRead", because the server stops sending after ~108 KB.
            logger.warning("ClassyFire get_query(%s) failed: %s", query_id, e)
            return e.response
        return response

    @staticmethod
    @cached(service='classyfire')
    def get_classification(query_id, per_page=None, use_cache=True):
        """
        Fetch a finished query's classification, a page at a time, as a dict.

        A convenience layer over ClassyFire's ``/queries/{id}.json``. The
        server cuts off any response beyond about 108 KB, so asking for a
        query whole fails; this asks for pages small enough to arrive
        complete and joins them. Requests go through the shared transport,
        with pacing, retries and this module's exceptions.

        Args:
            query_id: The query ID.
            per_page: Entities per page;
                [`PAGE_SIZE`][provesid.classyfire.ClassyFireAPI.PAGE_SIZE]
                when None.
            use_cache: When False, do not read the cache.

        Returns:
            (dict): The first page's fields (``id``, ``label``,
            ``classification_status``, ``number_of_elements``, ...) with
            ``entities`` and ``invalid_entities`` joined across every page.
            ``number_of_pages`` is the count for ``per_page``.

        Raises:
            ClassyFireNotFoundError: There is no query with this ID.
            ClassyFireError: A page failed, or came back cut off.

        Examples:
            >>> result = ClassyFireAPI.get_classification(1)         # doctest: +SKIP
            >>> result["entities"][0]["kingdom"]["name"]             # doctest: +SKIP
            'Organic compounds'
        """
        per_page = per_page or ClassyFireAPI.PAGE_SIZE
        url = f"{ClassyFireAPI.URL}/queries/{query_id}.json"

        def fetch(page):
            return ClassyFireAPI._http.get_json(url, params={"page": page, "per_page": per_page})

        result = fetch(1)
        result["entities"] = list(result.get("entities") or [])
        result["invalid_entities"] = _as_list(result.get("invalid_entities"))
        for page in range(2, (result.get("number_of_pages") or 1) + 1):
            data = fetch(page)
            result["entities"].extend(data.get("entities") or [])
            result["invalid_entities"].extend(_as_list(data.get("invalid_entities")))
        return result


def _as_list(value):
    """
    ``invalid_entities`` as a list: a whole query gives a list, a page gives 0.

    Args:
        value: The field as ClassyFire returned it.

    Returns:
        The list, or ``[]`` for anything that is not one.

    Examples:
        >>> _as_list(0), _as_list([{"id": 1}])
        ([], [{'id': 1}])
    """
    return list(value) if isinstance(value, list) else []
