"""
The ClassyFire web service (Wishart lab): :class:`ClassyFireAPI`.

.. warning::
   The service has classified nothing new since February 2023. It still
   serves queries submitted before then, but a new submission is never
   processed. For ChEBI classes computed offline use
   :class:`provesid.taxonomy.ChebifierClassifier`.

This is the last module still calling ``requests`` directly, with no pacing
and no retry, and its methods return :class:`requests.Response` objects rather
than parsed data. Plan step 7 replaces it; until then it is documented as it
is.

Examples:
    >>> from provesid.classyfire import ClassyFireAPI
    >>> response = ClassyFireAPI.get_query(1)                  # doctest: +SKIP
    >>> response.json()["classification_status"]               # doctest: +SKIP
    'Done'
"""

import requests
from .cache import cached, clear_cache, get_cache_info

class ClassyFireAPI:
    """
    Class to interact with the ClassyFire API. The class is converted from the original Ruby code
    provided at https://bitbucket.org/wishartlab/classyfire_api/src/master/lib/classyfire_api.rb

    Every method is static. Submit a structure with :meth:`submit_query`, poll
    :meth:`query_status`, then fetch the result with :meth:`get_query`. See
    the module warning: new submissions are no longer processed.

    Examples:
        >>> response = ClassyFireAPI.submit_query(
        ...     "Example Query", "C1=CC(=CC=C1[N+](=O)[O-])Cl")    # doctest: +SKIP
        >>> query_id = response.json()["id"]                       # doctest: +SKIP
        >>> ClassyFireAPI.query_status(query_id).json()            # doctest: +SKIP
    """
    URL = 'http://classyfire.wishartlab.com'

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
            The statistics :func:`provesid.cache.get_cache_info` reports for
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
            requests.Response: The service's answer, whose JSON carries the
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
                headers={"Accept": "application/json", "Content-Type": "application/json"}
            )
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            return e.response
        except requests.exceptions.RequestException as e:
            return e.response
        return response

    @staticmethod
    @cached(service='classyfire')
    def query_status(query_id, use_cache=True):
        """
        Retrieves the status of a query.

        Args:
            query_id: The ID :meth:`submit_query` returned.
            use_cache: When False, do not read the cache.

        Returns:
            requests.Response: The status response, or None when the request
            failed for any reason.

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
                headers={"Accept": "application/json", "Content-Type": "application/json"}
            )
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as e:
            return None

    @staticmethod
    @cached(service='classyfire')
    def get_query(query_id, format="json", use_cache=True):
        """
        Retrieve a query's classification.

        Args:
            query_id: The query ID.
            format: ``"json"``, ``"sdf"`` or ``"csv"``.
            use_cache: When False, do not read the cache.

        Returns:
            requests.Response: The result in the requested format. An HTTP
            error comes back as its response, not raised; a connection failure
            comes back as None.

        Raises:
            ValueError: If ``format`` is not one of the three.

        Examples:
            >>> result = ClassyFireAPI.get_query(1).json()          # doctest: +SKIP
            >>> result["classification_status"], result["number_of_elements"]  # doctest: +SKIP
            ('Done', 655)
        """
        try:
            if format == "json":
                response = requests.get(
                    f"{ClassyFireAPI.URL}/queries/{query_id}.json",
                    headers={"Accept": "application/json"}
                )
            elif format == "sdf":
                response = requests.get(
                    f"{ClassyFireAPI.URL}/queries/{query_id}.sdf",
                    headers={"Accept": "chemical/x-mdl-sdfile"}
                )
            elif format == "csv":
                response = requests.get(
                    f"{ClassyFireAPI.URL}/queries/{query_id}.csv",
                    headers={"Accept": "text/csv"}
                )
            else:
                raise ValueError("Invalid format. Must be one of json, sdf, or csv.")
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            return e.response
        except requests.exceptions.RequestException as e:
            return e.response
        return response
