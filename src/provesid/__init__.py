# chemid __init__.py
__version__ = "0.2.0"

from .cascommonchem import CASCommonChem
from .chebi import ChEBI, ChEBIError, get_chebi_entity, search_chebi
from .classyfire import ClassyFireAPI
from .opsin import OPSIN
from .pubchem import PubChemAPI, CompoundProperties, PubChemNotFoundError, PubChemError, Domain
from .pubchemview import (
    PubChemView, 
    PropertyData, 
    PubChemViewError, 
    PubChemViewNotFoundError,
    get_experimental_property,
    get_all_experimental_properties,
    get_property_values_only,
    get_property_table
)
from .cache import (
    clear_cache,
    get_cache_info,
    export_cache,
    import_cache,
    get_cache_size,
    set_cache_warning_threshold,
    enable_cache_warnings
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
    nci_get_formula
)
from .utils import check_CASRN