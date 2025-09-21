from .cascommonchem import CASCommonChem
from .pubchem import PubChemAPI
from typing import Union, List, Optional
import pandas as pd
from tqdm import tqdm
import logging

# Optional RDKit import
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    Chem = None
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available. Install with: pip install rdkit-pypi")

def smiles_to_canonical(smiles: str) -> Optional[str]:
    """
    Convert SMILES to canonical SMILES using RDKit.
    
    Args:
        smiles: Input SMILES string
        
    Returns:
        Canonical SMILES string or None if conversion fails
    """
    if not RDKIT_AVAILABLE or Chem is None:
        logging.warning("RDKit not available. Returning original SMILES.")
        return smiles
    
    if pd.isna(smiles) or smiles == "" or smiles == "nan":
        return None
        
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception as e:
        logging.warning(f"Failed to canonicalize SMILES '{smiles}': {e}")
        return None

def casrn_to_compounds(cas_rn: Union[str, List[str]], 
                      show_progress: bool = True) -> pd.DataFrame:
    """
    Process a list of CAS RN and return a merged DataFrame with compound IDs and properties.
    
    First looks for compounds in CAS Common Chemistry, then searches PubChem for missing compounds.
    Finally cleans up and merges the DataFrames.
    
    Args:
        cas_rn: Single CAS RN string or list of CAS RN strings
        show_progress: Whether to show progress bars
        
    Returns:
        DataFrame with columns: CASRN, name, SMILES, canonical_smiles, InChI, 
        InChIKey, molecular_formula, molecular_mass, foundby, source
    """
    # Initialize APIs
    ccc_api = CASCommonChem()
    pubchem_api = PubChemAPI()
    
    # Convert input to list
    if isinstance(cas_rn, str):
        casrn_list = [cas_rn]
    else:
        casrn_list = cas_rn
    
    # Step 1: Query CAS Common Chemistry
    ccc_results = []
    iterator = tqdm(casrn_list, desc="Querying CAS Common Chemistry") if show_progress else casrn_list
    
    for casrn in iterator:
        try:
            res = ccc_api.cas_to_detail(casrn)
            res["cas_rn"] = casrn
            ccc_results.append(res)
        except Exception as e:
            logging.warning(f"Failed to query CAS {casrn}: {e}")
            # Add empty result to maintain structure
            ccc_results.append({
                "cas_rn": casrn,
                "canonicalSmile": "",
                "smile": "",
                "status": "Error"
            })
    
    df_ccc = pd.DataFrame(ccc_results)
    
    # Find missing compounds (those without canonical SMILES)
    missing_casrns = df_ccc[
        (df_ccc["canonicalSmile"] == "") | 
        (df_ccc["canonicalSmile"].isna())
    ]["cas_rn"].tolist()
    
    # Step 2: Query PubChem for missing compounds
    pubchem_results = []
    if missing_casrns:
        iterator = tqdm(missing_casrns, desc="Querying PubChem for missing compounds") if show_progress else missing_casrns
    
    cid_list = []
    for casrn in iterator:
        res = pubchem_api.find_cids_comprehensive(casrn)
        cid_list.append(res)
    df_cid_missing = pd.DataFrame(cid_list)

    pubchem_results = []
    cid_iterator = tqdm(df_cid_missing.iterrows(), total=len(df_cid_missing), desc="Fetching compound details from PubChem") if show_progress else df_cid_missing.iterrows()
    for i, row in cid_iterator:
        for cid in row["total_unique_cids"]:
            res = pubchem_api.get_all_compound_info(cid)
            res["casrn"] = row["query"]
            res["synonyms"] = pubchem_api.get_compound_synonyms(cid)
            pubchem_results.append(res)
    
    # Step 3: Process and standardize CAS Common Chemistry results
    df_ccc_clean = df_ccc[df_ccc["smile"].notna() & (df_ccc["smile"] != "")].copy()
    
    if not df_ccc_clean.empty:
        df_ccc_clean.rename(columns={
            "cas_rn": "CASRN", 
            "canonicalSmile": "canonical_smiles", 
            "inchi": "InChI", 
            "inchiKey": "InChIKey",
            "molecularFormula": "molecular_formula",
            "molecularMass": "molecular_mass", 
            "name": "name", 
            "smile": "SMILES"
        }, inplace=True)
        
        df_ccc_light = df_ccc_clean[[
            "CASRN", "name", "SMILES", "canonical_smiles", "InChI", 
            "InChIKey", "molecular_formula", "molecular_mass"
        ]].copy()
        df_ccc_light["foundby"] = "CASRN"
        df_ccc_light["source"] = "CAS Common Chemistry"
    else:
        df_ccc_light = pd.DataFrame()
    
    # Step 4: Process and standardize PubChem results
    if pubchem_results:
        df_pubchem = pd.DataFrame(pubchem_results)
        df_pubchem = df_pubchem[df_pubchem["CID"].notna()].copy()
        
        if not df_pubchem.empty:
            df_cid_missing.rename(columns={"query": "casrn", "recommended_domain": "source"}, inplace=True)
            df_cid_missing["source"] = df_cid_missing["source"]+" PubChem"
            df_pubchem = df_pubchem.merge(df_cid_missing[["casrn", "source"]], on="casrn", how="left")
            df_pubchem["foundby"] = "CASRN"
            df_pubchem.rename(columns={
                "MolecularFormula": "molecular_formula", 
                "MolecularWeight": "molecular_mass",
                "Title": "name",
                "casrn": "CASRN"
            }, inplace=True)
            
            df_pubchem["foundby"] = "CASRN"
            df_pubchem["SMILES"] = df_pubchem["SMILES"].astype(str)
            
            # Generate canonical SMILES
            if show_progress:
                tqdm.pandas(desc="Generating canonical SMILES")
                df_pubchem["canonical_smiles"] = df_pubchem["SMILES"].progress_apply(smiles_to_canonical)
            else:
                df_pubchem["canonical_smiles"] = df_pubchem["SMILES"].apply(smiles_to_canonical)
            
            df_pubchem_light = df_pubchem[[
                "CASRN", "name", "SMILES", "canonical_smiles", "InChI", 
                "InChIKey", "molecular_formula", "molecular_mass", "foundby", "source"
            ]].copy()
        else:
            df_pubchem_light = pd.DataFrame()
    else:
        df_pubchem_light = pd.DataFrame()
    
    # Step 5: Combine results
    dataframes_to_concat = []
    if not df_ccc_light.empty:
        dataframes_to_concat.append(df_ccc_light)
    if not df_pubchem_light.empty:
        dataframes_to_concat.append(df_pubchem_light)
    
    if dataframes_to_concat:
        df_all_compounds = pd.concat(dataframes_to_concat, ignore_index=True)
    else:
        # Return empty DataFrame with expected columns
        df_all_compounds = pd.DataFrame(columns=[
            "CASRN", "name", "SMILES", "canonical_smiles", "InChI", 
            "InChIKey", "molecular_formula", "molecular_mass", "foundby", "source"
        ])
    
    return df_all_compounds