"""
PROVESID: chemical identifiers and properties, offline first.

:class:`Search` resolves CAS numbers, names, SMILES, InChIs, InChIKeys,
DTXSIDs and formulas against the offline databases that are installed ---
PubChem's CAS-bearing compounds (:class:`PubChemID`), EPA CompTox
(:class:`CompToxID`), ChEBI (:class:`ChebiSDF`) and ChEMBL (:class:`CheMBL`),
with ZeroPM (:class:`ZeroPM`) on request --- and scores each answer by how
many of them agree. :mod:`provesid.datasets` installs the databases by name;
nothing is downloaded on your behalf.

The online clients --- :class:`PubChemAPI`, :class:`PubChemView`,
:class:`NCIChemicalIdentifierResolver`, :class:`ChEBI`, :class:`OPSIN` and
:class:`CASCommonChem` --- share one transport (:mod:`provesid.http`) that
paces, retries and caches. :class:`Search` asks two of them only when
``online_fallback=True`` and no offline source answered.

Examples:
    >>> import provesid
    >>> provesid.datasets.missing()                 # what Search would lack  # doctest: +SKIP
    []
    >>> df = provesid.Search("cas", show_progress=False).search("50-78-2")
    >>> df.loc[0, "InChIKey"]
    'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
"""

__version__ = "0.7.0"

from .http import (
    HTTPClient,
    ServiceError,
    NotFoundError,
    RateLimitError,
    ServiceTimeoutError,
    release_holds,
)
from . import datasets
from .datasets import (
    DATASETS,
    Dataset,
    DownloadError,
    MissingDatasetError,
    download_file,
    human_bytes,
    md5_of_file,
    read_checksum,
)
from .cascommonchem import (
    CASCommonChem,
    CASCommonChemError,
    CASCommonChemNotFoundError,
    CASCommonChemTimeoutError,
)
from .chebi import (
    ChEBI,
    ChEBIError,
    ChEBINotFoundError,
    ChEBITimeoutError,
    get_chebi_entity,
    search_chebi,
)
from .chebi_sdf import ChebiSDF
from .chembl import CheMBL, ChEMBLError
from .classyfire import ClassyFireAPI
from .opsin import OPSIN, OPSINError, OPSINNotFoundError, OPSINTimeoutError, PYOPSIN
from .pubchem import (
    PubChemAPI,
    CompoundProperties,
    PubChemNotFoundError,
    PubChemError,
    PubChemServerError,
    PubChemTimeoutError,
    Domain,
)
from .pubchem_id import PubChemID
from .comptox import CompToxID
from .sqlite_client import SQLiteClient, DatabaseClosedError
from .pubchemview import (
    PubChemView,
    PropertyData,
    PubChemViewError,
    PubChemViewNotFoundError,
    get_experimental_property,
    get_all_experimental_properties,
    get_property_values_only,
    get_property_table,
)
from .pubchemview_parse import ParsedValue, parse_value
from .config import set_cas_api_key, get_cas_api_key, remove_cas_api_key, show_config
from .cache import (
    CACHE_KEY_VERSION,
    CACHE_SERVICES,
    clear_cache,
    get_cache_info,
    get_all_cache_info,
    export_cache,
    import_cache,
    get_cache_size,
    get_service_cache,
    set_cache_warning_threshold,
    enable_cache_warnings,
)
from .resolver import (
    NCIChemicalIdentifierResolver,
    NCIResolverError,
    NCIResolverNotFoundError,
    nci_cas_to_mol,
    nci_id_to_mol,
    nci_resolver,
    nci_smiles_to_names,
    nci_name_to_smiles,
    nci_inchi_to_smiles,
    nci_cas_to_inchi,
    nci_get_molecular_weight,
    nci_get_formula,
)
from .utils import check_CASRN
from .zeropm import ZeroPM
from .reach import REACHDossierID
from .search import (
    Search,
    normalize_structure,
    strip_salts,
    resolve_cascade,
    mw_within,
    OUTPUT_COLUMNS,
)
from .taxonomy import (
    ChebifierClassifier,
    ChebifierError,
    ChebifierMissingError,
    classify_chebifier,
    chebifier_available,
    default_ensemble_available,
    missing_ensemble_modules,
    ensure_v244_indices,
    ensure_element_class_mappings,
    TAXONOMY_COLUMNS,
)
