import logging
import time
from typing import Any
import requests
from .cache import cached, clear_cache, get_cache_info
from .http import (
    HTTPClient,
    NotFoundError,
    Outcome,
    ServiceError,
    ServiceTimeoutError,
    default_classify,
)
from py2opsin import py2opsin

#: Seconds between two requests to OPSIN. The service is a name parser rather
#: than a database lookup, and EBI publishes no per-IP figure for it, so this is
#: politeness: ten requests a second. Note that OPSIN and ChEBI are both served
#: from ``www.ebi.ac.uk``, so they share a pacing clock and ten per second is
#: the budget for the two of them together --- which is the point of keying the
#: clock by host.
OPSIN_MIN_INTERVAL = 0.1


class OPSINError(ServiceError):
    """
    Base exception for OPSIN web-service failures.

    The methods below do not let it reach the caller: they report a failure in
    the ``status`` key of the dict they return, which is the contract they have
    always had.

    Example:
        >>> issubclass(OPSINTimeoutError, OPSINError)
        True
    """
    pass


class OPSINNotFoundError(OPSINError, NotFoundError):
    """
    Raised for an absence OPSIN reports with no usable body.

    Rare in practice: OPSIN explains an unparseable name in the body of its
    404, and :func:`opsin_classify` treats that body as the answer rather than
    as absence.

    Example:
        >>> issubclass(OPSINNotFoundError, NotFoundError)
        True
    """
    pass


class OPSINTimeoutError(OPSINError, ServiceTimeoutError):
    """
    Every attempt timed out or the connection could not be made.

    Example:
        >>> issubclass(OPSINTimeoutError, ServiceTimeoutError)
        True
    """
    pass


def opsin_classify(response: requests.Response) -> Outcome:
    """
    Decide what an OPSIN response means, treating its 404 as an answer.

    OPSIN answers a name it cannot parse with HTTP 404 and a complete JSON
    body: ``{"status": "FAILURE", "message": "<name> was uninterpretable due
    to the following section of the name: ..."}``. The 404 is an envelope, not
    the answer --- the answer is in the body, and it is the only place the
    *reason* exists. Verified live on 2026-09-19 against
    ``notachemical12345``.

    Treating that as :attr:`~provesid.http.Outcome.ABSENT` would raise before
    the body was read, which is what the module did by hand before: it mapped
    the status to the string ``"FAILURE"`` and returned, so ``message`` was
    empty for every failure the module ever reported --- including in the
    WARNING that ``get_id_from_list`` logs, which has therefore always said
    "Failed to get ID for x: ".

    Args:
        response: The response to classify.

    Returns:
        :attr:`~provesid.http.Outcome.OK` for a 404, whose body is the answer;
        otherwise whatever :func:`~provesid.http.default_classify` says.

    Example:
        >>> class R: status_code = 404
        >>> opsin_classify(R()).name
        'OK'
        >>> class R: status_code = 503
        >>> opsin_classify(R()).name
        'RETRY'
    """
    if response.status_code == 404:
        return Outcome.OK
    return default_classify(response)


def _not_resolved(result: Any) -> bool:
    """
    Report whether an OPSIN lookup failed to produce a structure.

    A name OPSIN could not parse today it cannot parse tomorrow either, so the
    absence itself would be safe to cache --- but a timeout is reported the
    same way, and that one must not become permanent. The two cannot be told
    apart from the returned dict, so neither is stored; re-asking a name parser
    is cheap.

    Args:
        result: The value returned by ``get_id``.

    Returns:
        True unless the result reports ``status`` as ``"SUCCESS"``.

    Example:
        >>> _not_resolved({"status": "FAILURE", "message": "uninterpretable"})
        True
        >>> _not_resolved({"status": "SUCCESS", "smiles": "CCO"})
        False
    """
    return not (isinstance(result, dict) and result.get("status") == "SUCCESS")


def _any_not_resolved(results: Any) -> bool:
    """
    Report whether a batch of OPSIN lookups contains anything unresolved.

    ``get_id_from_list`` is a loop over ``get_id``, which caches each name on
    its own, so the aggregate adds nothing worth cementing a failure for: a
    batch of a thousand names where one request timed out would otherwise
    remember that timeout every time the same batch was asked for.

    Args:
        results: The list ``get_id_from_list`` returned.

    Returns:
        True when the list is empty or any entry is not a success.

    Example:
        >>> _any_not_resolved([{"status": "SUCCESS"}, {"status": "FAILURE"}])
        True
        >>> _any_not_resolved([{"status": "SUCCESS"}])
        False
    """
    if not results:
        return True
    return any(_not_resolved(result) for result in results)


class OPSIN:
    """
    A Python interface to the OPSIN name-to-structure web service.

    OPSIN turns a systematic IUPAC name into a structure. This class calls the
    hosted service; :class:`PYOPSIN` runs the same parser locally through
    ``py2opsin``, which is faster and works offline, and is the better choice
    for anything but a handful of names.

    Attributes:
        base_url: The OPSIN web service endpoint.
        use_cache: Whether lookups are served from the cache.

    Example:
        >>> OPSIN().base_url
        'https://www.ebi.ac.uk/opsin/ws/'
    """

    def __init__(self, use_cache: bool = True):
        """
        Initialize the OPSIN web-service client.

        Args:
            use_cache: Whether to use the cache for lookups (default True).
                When False, skips the cache lookup but still stores results.
        """
        # This module used to name the Cambridge host, opsin.ch.cam.ac.uk,
        # which answers every request with a 301 to here --- so asking here
        # directly saves a redirect per name. Verified live on 2026-09-19.
        self.base_url = "https://www.ebi.ac.uk/opsin/ws/"
        self.use_cache = use_cache
        self.logger = logging.getLogger(__name__)

        # One shared transport. Before this the module had no pacing, no retry
        # and no timeout handling: a momentary 503 reached the caller as a raw
        # KeyError out of the status-code table.
        self._http = HTTPClient(
            min_interval=OPSIN_MIN_INTERVAL,
            classify=opsin_classify,
            error_cls=OPSINError,
            not_found_cls=OPSINNotFoundError,
            timeout_cls=OPSINTimeoutError,
            pace_host=self.base_url,
            logger=self.logger,
        )

    #: Bumped whenever an entry written by an earlier version must not be
    #: served. Version 2 is this change: the old code cached *failures*, so a
    #: name that hit one momentary 503 is on disk as a permanent ``"FAILURE"``
    #: for a name OPSIN parses perfectly well --- and every cached failure
    #: carries an empty ``message``, because the old code never read the body
    #: that explains it. Neither is fixable in place, so the old entries are
    #: made unreachable instead.
    CACHE_SCHEMA_VERSION = 2

    def __cache_key__(self) -> tuple:
        """
        Identify this client for cache-key purposes.

        Without this the key is the class name alone, which is stable across
        processes --- what :func:`~provesid.cache.stable_key_part` is for --- but
        carries no way to retire an entry whose content is now known to be
        wrong. The schema version is that way.

        Returns:
            Tuple of the class path, the configured base URL and
            :attr:`CACHE_SCHEMA_VERSION`.

        Example:
            >>> OPSIN().__cache_key__()
            ('provesid.opsin.OPSIN', 'https://www.ebi.ac.uk/opsin/ws/', 2)
        """
        return ("provesid.opsin.OPSIN", self.base_url,
                self.CACHE_SCHEMA_VERSION)

    def clear_cache(self):
        """Clear the cache for all OPSIN methods"""
        clear_cache(service='opsin')
    
    def get_cache_info(self):
        """Get information about the current cache state"""
        return get_cache_info(service='opsin')

    @cached(service='opsin', skip_if=_not_resolved)
    def get_id(self, iupac_name: str, timeout=30):
        """
        Return the structure OPSIN parses out of an IUPAC name.

        Adapted from the IUPAC WorldFAIR book:
        https://iupac.github.io/WFChemCookbook/tools/opsin_api_jn.html

        Args:
            iupac_name: The systematic name to parse.
            timeout: Request timeout in seconds.

        Returns:
            A dict in the shape of :meth:`_empty_res`. ``status`` is OPSIN's
            own, ``"SUCCESS"`` or ``"FAILURE"``, and on a failure ``message``
            carries OPSIN's explanation of which part of the name it could not
            read. On a transport failure ``status`` is ``"FAILURE"`` and
            ``message`` describes the failure; nothing is raised, and nothing
            is cached.

        Example:
            >>> OPSIN().get_id("ethanol")["smiles"]        # doctest: +SKIP
            'CCO'
            >>> OPSIN().get_id("notachemical")["status"]   # doctest: +SKIP
            'FAILURE'
        """
        apiurl = self.base_url + iupac_name + '.json'
        res = self._empty_res()
        res["iupac_name"] = iupac_name

        try:
            jsondata = self._http.get_json(apiurl, timeout=timeout)
        except OPSINError as e:
            # A parse failure is not an error --- opsin_classify reads that out
            # of the 404 body. Reaching here means the service itself failed.
            res["status"] = "FAILURE"
            res["message"] = str(e)
            self.logger.warning(f"OPSIN request failed for {iupac_name}: {e}")
            return res

        res["status"] = jsondata.get("status", "FAILURE")
        res["message"] = jsondata.get("message", "")
        if res["status"] != "SUCCESS":
            return res

        res["smiles"] = jsondata["smiles"]
        res["stdinchi"] = jsondata["stdinchi"]
        res["stdinchikey"] = jsondata["stdinchikey"]
        res["inchi"] = jsondata["inchi"]
        return res

    @staticmethod
    def _empty_res():
        """
        create an empty response dictionary of the following format:
        {
            "status": "SUCCESS",
            "message": "",
            "inchi": "InChI=1/C2H2Cl4/c3-1(4)2(5)6/h1-2H",
            "stdinchi": "InChI=1S/C2H2Cl4/c3-1(4)2(5)6/h1-2H",
            "stdinchikey": "QPFMBZIOSGYJDE-UHFFFAOYSA-N",
            "smiles": "ClC(C(Cl)Cl)Cl"
        }
        """
        return {
            "iupac_name": "",
            "status": "",
            "message": "",
            "inchi": "",
            "stdinchi": "",
            "stdinchikey": "",
            "smiles": ""
        }
    
    @cached(service='opsin', skip_if=_any_not_resolved)
    def get_id_from_list(self, iupac_names: list, timeout=30, pause_time=0.5):
        """
        Parse a list of IUPAC names, one request each.

        Args:
            iupac_names: The systematic names to parse.
            timeout: Per-request timeout in seconds.
            pause_time: Seconds to sleep after each name, on top of the pacing
                the transport already applies.

        Returns:
            One dict per input name, in order, each in the shape
            :meth:`get_id` returns. A name OPSIN could not parse is a
            ``"FAILURE"`` entry carrying OPSIN's explanation in ``message``,
            not a missing row.

        Example:
            >>> OPSIN().get_id_from_list(["ethanol", "benzene"])   # doctest: +SKIP
            [{'iupac_name': 'ethanol', ...}, {'iupac_name': 'benzene', ...}]
        """
        results = []
        for iupac_name in iupac_names:
            res = self.get_id(iupac_name, timeout)
            if res["status"] != "SUCCESS":
                self.logger.warning(
                    f"Failed to get ID for {iupac_name}: {res['message']}"
                )
            results.append(res)
            time.sleep(pause_time)
        return results
    
class PYOPSIN:
    """
    This class uses py2opsin (https://github.com/JacksonBurns/py2opsin) package 
    that can be installed via pip:
    pip install py2opsin
    It provides all functionalities as the OPSIN class above -and some more- but works offline 
    and faster. You need to have Java installed on your system.
    Note:
    output_format (str, optional) – One of “SMILES”, “ExtendedSMILES”, 
    “CML”, “InChI”, “StdInChI”, or “StdInChIKey”. Defaults to “SMILES”.
    """
    def __init__(self, jar_fpath = "default"):
        self.jar_fpath = jar_fpath

    def get_smiles(self, iupac_name: str):
        smiles = py2opsin(iupac_name, output_format="SMILES", jar_fpath=self.jar_fpath)
        return smiles
    
    def get_extended_smiles(self, iupac_name: str):
        ext_smiles = py2opsin(iupac_name, output_format="ExtendedSMILES", jar_fpath=self.jar_fpath)
        return ext_smiles
    
    def get_inchi(self, iupac_name: str):
        inchi = py2opsin(iupac_name, output_format="InChI", jar_fpath=self.jar_fpath)
        return inchi
    
    def get_std_inchi(self, iupac_name: str):
        std_inchi = py2opsin(iupac_name, output_format="StdInChI", jar_fpath=self.jar_fpath)
        return std_inchi
    
    def get_std_inchikey(self, iupac_name: str):
        std_inchikey = py2opsin(iupac_name, output_format="StdInChIKey", jar_fpath=self.jar_fpath)
        return std_inchikey
    
    def get_CML(self, iupac_name: str):
        cml = py2opsin(iupac_name, output_format="CML", jar_fpath=self.jar_fpath)
        return cml

    def get_id(self, iupac_name: str):
        res = self._empty_res()
        res["iupac_name"] = iupac_name
        res["smiles"] = self.get_smiles(iupac_name)
        res["extended_smiles"] = self.get_extended_smiles(iupac_name)
        res["inchi"] = self.get_inchi(iupac_name)
        res["stdinchi"] = self.get_std_inchi(iupac_name)
        res["stdinchikey"] = self.get_std_inchikey(iupac_name)
        res["cml"] = self.get_CML(iupac_name)
        if res["smiles"] == "":
            res["status"] = "FAILURE"
        else:
            res["status"] = "SUCCESS"                
        return res

    def get_id_from_list(self, iupac_names: list):
        """
        This function does not need a loop since the py2opsin package handles lists internally.
        """
        res = self.get_id(iupac_names)
        # convert to list of dicts
        results = []
        for i, name in enumerate(iupac_names):
            single_res = self._empty_res()
            single_res["iupac_name"] = name
            single_res["smiles"] = res["smiles"][i]
            single_res["extended_smiles"] = res["extended_smiles"][i]
            single_res["inchi"] = res["inchi"][i]
            single_res["stdinchi"] = res["stdinchi"][i]
            single_res["stdinchikey"] = res["stdinchikey"][i]
            single_res["cml"] = res["cml"][i]
            single_res["status"] = res["status"][i]
            results.append(single_res)
        return results

    @staticmethod
    def _empty_res():
        """
        “SMILES”, “ExtendedSMILES”, 
    “CML”, “InChI”, “StdInChI”, or “StdInChIKey”. Defaults to “SMILES”.
        create an empty response dictionary of the following format:
        {
            "status": "SUCCESS",
            "message": "",
            "inchi": "InChI=1/C2H2Cl4/c3-1(4)2(5)6/h1-2H",
            "stdinchi": "InChI=1S/C2H2Cl4/c3-1(4)2(5)6/h1-2H",
            "stdinchikey": "QPFMBZIOSGYJDE-UHFFFAOYSA-N",
            "smiles": "ClC(C(Cl)Cl)Cl",
            "extended_smiles": "",
            "cml": ""
        }
        """
        return {
            "status": "",
            "message": "",
            "inchi": "",
            "stdinchi": "",
            "stdinchikey": "",
            "smiles": "",
            "extended_smiles": "",
            "cml": ""
        }

