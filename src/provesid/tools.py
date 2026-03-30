from .pubchem import PubChemID
from .opsin import PYOPSIN
from .comptox import CompToxID
from .chebi import ChebiSDF
from .zeropm import ZeroPM
from .chembl import CheMBL
from typing import Union, List, Optional, Dict, Any, Tuple
import pandas as pd
from tqdm import tqdm
import logging
import re
from difflib import SequenceMatcher

# Optional RDKit import
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    Chem = None
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available. Install with: pip install rdkit-pypi")

def iupac_name_to_id(iupac_name_list: list):
    """
    converts a list of iupac_names to a dataframe that contains the iupac_name, inchi, inchikey, and smiles for each name. It uses the 
    PYOPSIN class to call the py2opsin package. Make sure that java is installed on your system.
    Note that at the moment, the column extended_smiles seems to have an encoding issue and the column cml is not strictly necessary for most users.
    """
    opsin = PYOPSIN()
    results = opsin.get_id_from_list(iupac_name_list)
    return pd.DataFrame(results)
    
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


def _is_missing(value: Any) -> bool:
    """Return True for None, NaN, or empty-string-like values."""
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip() == "" or value.strip().lower() == "nan"
    try:
        return bool(pd.isna(value))
    except Exception:
        return False


def _pick_first(*values: Any) -> Any:
    """Return the first non-missing value from values."""
    for value in values:
        if not _is_missing(value):
            return value
    return None


def _normalize_synonyms(value: Any) -> Optional[str]:
    """Convert synonym values to a normalized string representation."""
    if _is_missing(value):
        return None

    if isinstance(value, (list, tuple, set)):
        cleaned = [str(v).strip() for v in value if not _is_missing(v)]
        return "; ".join(cleaned) if cleaned else None

    text = str(value).strip()
    return text if text else None


_CAS_PATTERN = re.compile(r"\b\d{2,7}-\d{2}-\d\b")


def _to_float(value: Any) -> Optional[float]:
    """Safely convert value to float."""
    if _is_missing(value):
        return None
    try:
        return float(value)
    except Exception:
        return None


def _text_similarity(a: Optional[str], b: Optional[str]) -> float:
    """Return normalized text similarity in [0, 1]."""
    if _is_missing(a) or _is_missing(b):
        return 0.0
    a_text = str(a).strip().lower()
    b_text = str(b).strip().lower()
    if a_text == b_text:
        return 1.0
    return SequenceMatcher(None, a_text, b_text).ratio()


def _extract_cas_values(value: Any) -> List[str]:
    """Extract CAS numbers from strings, collections, and mapping values."""
    found: List[str] = []

    if value is None:
        return found

    if isinstance(value, dict):
        for dict_value in value.values():
            found.extend(_extract_cas_values(dict_value))
    elif isinstance(value, (list, tuple, set)):
        for item in value:
            found.extend(_extract_cas_values(item))
    else:
        text = str(value)
        found.extend(_CAS_PATTERN.findall(text))

    deduped = sorted(set(found))
    return deduped


def _inchi_to_smiles(inchi: Optional[str]) -> Optional[str]:
    """Convert InChI to SMILES when RDKit is available."""
    if _is_missing(inchi) or not RDKIT_AVAILABLE or Chem is None:
        return None
    try:
        mol = Chem.MolFromInchi(str(inchi))
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def _inchikey_from_smiles(smiles: Optional[str]) -> Optional[str]:
    """Convert SMILES to InChIKey when RDKit is available."""
    if _is_missing(smiles) or not RDKIT_AVAILABLE or Chem is None:
        return None
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return None
        inchi = Chem.MolToInchi(mol)
        if _is_missing(inchi):
            return None
        return Chem.InchiToInchiKey(inchi)
    except Exception:
        return None


def _first_cas(cas_values: List[str]) -> Optional[str]:
    """Return a deterministic first CAS value."""
    return cas_values[0] if cas_values else None


def _make_candidate(
    source: str,
    *,
    name: Optional[str] = None,
    iupac_name: Optional[str] = None,
    molecular_formula: Optional[str] = None,
    smiles: Optional[str] = None,
    inchi: Optional[str] = None,
    inchikey: Optional[str] = None,
    dtxsid: Optional[str] = None,
    molecular_mass: Optional[float] = None,
    synonyms: Optional[str] = None,
    cas_candidates: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Create normalized candidate record for cross-dataset consensus."""
    canonical_smiles, rdkit_mass = _smiles_to_canonical_and_mass(smiles)
    return {
        "source": source,
        "name": name,
        "IUPAC_name": iupac_name,
        "molecular_formula": molecular_formula,
        "SMILES": smiles,
        "canonical_smiles": canonical_smiles,
        "InChI": inchi,
        "InChIKey": inchikey,
        "DTXSID": dtxsid,
        "molecular_mass": _pick_first(_to_float(molecular_mass), rdkit_mass),
        "Synonyms": synonyms,
        "CAS_candidates": sorted(set(cas_candidates or [])),
    }


def _candidate_similarity(left: Dict[str, Any], right: Dict[str, Any]) -> float:
    """Score agreement between two normalized candidate records in [0, 1]."""
    if left is None or right is None:
        return 0.0

    score = 0.0
    weight = 0.0

    left_cas = set(left.get("CAS_candidates") or [])
    right_cas = set(right.get("CAS_candidates") or [])
    if left_cas and right_cas:
        weight += 3.0
        if left_cas.intersection(right_cas):
            score += 3.0

    left_smiles = left.get("canonical_smiles")
    right_smiles = right.get("canonical_smiles")
    if not _is_missing(left_smiles) and not _is_missing(right_smiles):
        weight += 4.0
        if left_smiles == right_smiles:
            score += 4.0

    left_ik = left.get("InChIKey")
    right_ik = right.get("InChIKey")
    if not _is_missing(left_ik) and not _is_missing(right_ik):
        weight += 3.0
        if str(left_ik) == str(right_ik):
            score += 3.0

    left_inchi = left.get("InChI")
    right_inchi = right.get("InChI")
    if not _is_missing(left_inchi) and not _is_missing(right_inchi):
        weight += 2.0
        if str(left_inchi) == str(right_inchi):
            score += 2.0

    left_formula = left.get("molecular_formula")
    right_formula = right.get("molecular_formula")
    if not _is_missing(left_formula) and not _is_missing(right_formula):
        weight += 1.0
        if str(left_formula) == str(right_formula):
            score += 1.0

    left_mass = _to_float(left.get("molecular_mass"))
    right_mass = _to_float(right.get("molecular_mass"))
    if left_mass is not None and right_mass is not None:
        weight += 2.0
        diff = abs(left_mass - right_mass)
        if diff <= 0.2:
            score += 2.0
        elif diff <= 1.0:
            score += 1.0

    name_sim = _text_similarity(left.get("name"), right.get("name"))
    if name_sim > 0.0:
        weight += 1.0
        if name_sim >= 0.9:
            score += 1.0
        elif name_sim >= 0.7:
            score += 0.5

    if weight == 0.0:
        return 0.0
    return score / weight


def _best_candidate_by_name(rows: List[Dict[str, Any]], query_name: str, name_key: str) -> Optional[Dict[str, Any]]:
    """Choose the most likely row for a name query."""
    if not rows:
        return None

    ranked = []
    for row in rows:
        sim = _text_similarity(query_name, row.get(name_key))
        completeness = sum(
            [
                1 for k in ("SMILES", "INCHI", "INCHIKEY", "MOLECULAR_FORMULA", "mf", "CASRN")
                if not _is_missing(row.get(k))
            ]
        )
        ranked.append((sim, completeness, row))

    ranked.sort(key=lambda item: (item[0], item[1]), reverse=True)
    return ranked[0][2]


def _candidate_compatible_with_consensus(
    candidate: Optional[Dict[str, Any]],
    consensus: Optional[Dict[str, Any]],
    threshold: float = 0.35,
) -> bool:
    """Return True when a candidate is compatible with the consensus anchor."""
    if candidate is None:
        return False
    if consensus is None:
        return True
    if candidate.get("source") == consensus.get("source"):
        return True
    return _candidate_similarity(candidate, consensus) >= threshold


def _apply_candidate_to_result(result: Dict[str, Any], candidate: Optional[Dict[str, Any]]) -> None:
    """Populate result fields from one candidate while preserving existing values."""
    if candidate is None:
        return

    result["CASRN"] = _pick_first(result.get("CASRN"), _first_cas(candidate.get("CAS_candidates") or []))
    result["name"] = _pick_first(result.get("name"), candidate.get("name"))
    result["IUPAC_name"] = _pick_first(result.get("IUPAC_name"), candidate.get("IUPAC_name"))
    result["molecular_formula"] = _pick_first(result.get("molecular_formula"), candidate.get("molecular_formula"))
    result["SMILES"] = _pick_first(result.get("SMILES"), candidate.get("SMILES"))
    result["InChI"] = _pick_first(result.get("InChI"), candidate.get("InChI"))
    result["InChIKey"] = _pick_first(result.get("InChIKey"), candidate.get("InChIKey"))
    result["DTXSID"] = _pick_first(result.get("DTXSID"), candidate.get("DTXSID"))
    result["molecular_mass"] = _pick_first(result.get("molecular_mass"), candidate.get("molecular_mass"))
    result["Synonyms"] = _pick_first(result.get("Synonyms"), _normalize_synonyms(candidate.get("Synonyms")))

    if _is_missing(result.get("source")) and not _is_missing(candidate.get("SMILES")):
        result["source"] = candidate.get("source")


def _compute_consensus(candidates: Dict[str, Optional[Dict[str, Any]]]) -> Tuple[Optional[str], Dict[str, float], float]:
    """Compute consensus source and per-source agreement scores."""
    valid = {k: v for k, v in candidates.items() if v is not None}
    if not valid:
        return None, {}, 0.0

    support: Dict[str, float] = {}
    sources = list(valid.keys())
    for source in sources:
        others = [other for other in sources if other != source]
        if not others:
            support[source] = 1.0
            continue
        sims = [_candidate_similarity(valid[source], valid[other]) for other in others]
        support[source] = sum(sims) / len(sims)

    priority = ["chebi", "comptox", "pubchem", "zeropm", "chembl"]
    consensus_source = sorted(
        support.keys(),
        key=lambda src: (support[src], -priority.index(src) if src in priority else -99),
        reverse=True,
    )[0]

    consensus_candidate = valid[consensus_source]
    source_match_scores = {
        src: (1.0 if src == consensus_source else _candidate_similarity(consensus_candidate, valid[src]))
        for src in sources
    }
    overall = sum(source_match_scores.values()) / len(source_match_scores)
    return consensus_source, source_match_scores, overall


def _candidate_from_chebi_row(row: Dict[str, Any]) -> Dict[str, Any]:
    """Normalize one ChEBI row into a candidate record."""
    return _make_candidate(
        "ChEBI",
        name=row.get("ChEBI NAME"),
        iupac_name=row.get("ChEBI NAME"),
        molecular_formula=row.get("FORMULA"),
        smiles=row.get("SMILES"),
        inchi=row.get("INCHI"),
        inchikey=row.get("INCHIKEY"),
        synonyms=_normalize_synonyms(row.get("SYNONYM")),
        cas_candidates=_extract_cas_values(row),
    )


def _candidate_from_comptox_row(row: Dict[str, Any]) -> Dict[str, Any]:
    """Normalize one CompTox row into a candidate record."""
    return _make_candidate(
        "CompTox",
        name=row.get("PREFERRED_NAME"),
        iupac_name=row.get("IUPAC_NAME"),
        molecular_formula=row.get("MOLECULAR_FORMULA"),
        smiles=row.get("SMILES"),
        inchi=row.get("INCHI"),
        inchikey=row.get("INCHIKEY"),
        dtxsid=row.get("DTXSID"),
        molecular_mass=_pick_first(row.get("AVERAGE_MASS"), row.get("MONOISOTOPIC_MASS")),
        synonyms=_normalize_synonyms(row.get("identifiers")),
        cas_candidates=_extract_cas_values([row.get("CASRN"), row.get("identifiers")]),
    )


def _candidate_from_pubchem_row(row: Dict[str, Any]) -> Dict[str, Any]:
    """Normalize one PubChemID row into a candidate record."""
    return _make_candidate(
        "PubChemID",
        name=row.get("cmpdname"),
        iupac_name=row.get("iupacname"),
        molecular_formula=row.get("mf"),
        smiles=row.get("smiles"),
        inchi=row.get("inchi"),
        inchikey=row.get("inchikey"),
        molecular_mass=row.get("mw"),
        synonyms=_normalize_synonyms(row.get("synonyms")),
        cas_candidates=_extract_cas_values(row.get("cas_numbers")),
    )


def _candidate_from_zeropm_name_table(name: str, table: pd.DataFrame) -> Optional[Dict[str, Any]]:
    """Normalize ZeroPM name lookup table into one candidate record."""
    if table is None or table.empty:
        return None

    working = table.copy()
    if "rank" in working.columns:
        working = working.sort_values(by="rank", ascending=True)

    first = working.iloc[0]
    inchi = first.get("inchi")
    smiles = _inchi_to_smiles(inchi)

    cas_values = []
    if "cas" in working.columns:
        cas_values = [str(v) for v in working["cas"].dropna().astype(str).tolist()]

    synonyms = None
    if "name" in working.columns and not working["name"].dropna().empty:
        synonyms = _normalize_synonyms(working["name"].dropna().astype(str).unique().tolist())

    return _make_candidate(
        "ZeroPM",
        name=name,
        iupac_name=name,
        smiles=smiles,
        inchi=inchi,
        inchikey=first.get("inchikey"),
        synonyms=synonyms,
        cas_candidates=_extract_cas_values(cas_values),
    )


def _candidate_from_zeropm_smiles(smiles_query: str, zeropm: ZeroPM) -> Optional[Dict[str, Any]]:
    """Normalize ZeroPM SMILES lookup path into one candidate record."""
    cas_result = zeropm.get_cas_from_smiles(smiles_query)
    cas_values = _extract_cas_values(cas_result)
    if not cas_values:
        return None

    tables = []
    for cas in cas_values[:5]:
        table = zeropm.get_id_table_from_cas(cas)
        if table is not None and not table.empty:
            tables.append(table)

    if not tables:
        return _make_candidate("ZeroPM", smiles=smiles_query, cas_candidates=cas_values)

    combined = pd.concat(tables, ignore_index=True)
    first = combined.iloc[0]
    inchi = first.get("inchi")
    smiles = _pick_first(smiles_query, _inchi_to_smiles(inchi))

    synonyms = None
    if "synonyms" in combined.columns and not combined["synonyms"].dropna().empty:
        synonyms = _normalize_synonyms(combined["synonyms"].dropna().astype(str).unique().tolist())

    return _make_candidate(
        "ZeroPM",
        smiles=smiles,
        inchi=inchi,
        inchikey=first.get("inchikey"),
        synonyms=synonyms,
        cas_candidates=_extract_cas_values(cas_values),
    )


def _candidate_from_chembl_row(row: Dict[str, Any], chembl: Optional[CheMBL] = None) -> Dict[str, Any]:
    """Normalize one ChEMBL row into a candidate record."""
    props = None
    molregno = row.get("molregno")
    if chembl is not None and not _is_missing(molregno):
        try:
            props = chembl.get_properties(int(molregno))
        except Exception:
            props = None

    return _make_candidate(
        "ChEMBL",
        name=row.get("pref_name"),
        molecular_formula=None,
        smiles=row.get("canonical_smiles"),
        inchi=row.get("standard_inchi"),
        inchikey=row.get("standard_inchi_key"),
        molecular_mass=(props or {}).get("mw_freebase"),
        synonyms=_normalize_synonyms(row.get("synonyms")),
        cas_candidates=_extract_cas_values(row.get("synonyms")),
    )


def ids_from_name(
    name: str,
    chebi: Optional[ChebiSDF] = None,
    comptox: Optional[CompToxID] = None,
    pubchem: Optional[PubChemID] = None,
    zeropm: Optional[ZeroPM] = None,
    chembl: Optional[CheMBL] = None,
) -> Dict[str, Any]:
    """
    Resolve one chemical name into a unified identifier record from offline sources.

    The function queries the offline datasets (ChEBI, CompTox, PubChemID, ZeroPM),
    computes a cross-source consensus score based on structure and identifier
    agreement, and then populates output fields using the same source priority as
    ids_from_CAS.

    Args:
        name: Chemical name input.
        chebi: Optional initialized ChebiSDF client.
        comptox: Optional initialized CompToxID client.
        pubchem: Optional initialized PubChemID client.
        zeropm: Optional initialized ZeroPM client.
        chembl: Optional initialized CheMBL client.

    Returns:
        Dictionary similar to ids_from_CAS with additional consensus fields:
        match_score, consensus_source, and source_match_scores.
    """
    result: Dict[str, Any] = {
        "query": name,
        "CASRN": None,
        "name": None,
        "IUPAC_name": None,
        "molecular_formula": None,
        "SMILES": None,
        "canonical_smiles": None,
        "InChI": None,
        "InChIKey": None,
        "DTXSID": None,
        "molecular_mass": None,
        "Synonyms": None,
        "foundby": "name",
        "source": None,
        "match_score": 0.0,
        "consensus_source": None,
        "source_match_scores": {},
    }

    candidates: Dict[str, Optional[Dict[str, Any]]] = {
        "chebi": None,
        "comptox": None,
        "pubchem": None,
        "zeropm": None,
        "chembl": None,
    }

    if chebi is not None:
        try:
            chebi_rows = chebi.search_by_name(name, exact=True)
            if not chebi_rows:
                chebi_rows = chebi.search_by_name(name, exact=False)
            if not chebi_rows:
                chebi_rows = chebi.search_by_synonym(name, exact=False)
            if chebi_rows:
                best_chebi = _best_candidate_by_name(chebi_rows, name, "ChEBI NAME")
                if best_chebi is not None:
                    candidates["chebi"] = _candidate_from_chebi_row(best_chebi)
        except Exception as e:
            logging.warning(f"ChEBI name lookup failed for '{name}': {e}")

    if comptox is not None:
        try:
            comptox_rows: List[Dict[str, Any]] = []
            exact_row = comptox.get_by_name(name)
            if exact_row:
                comptox_rows.append(exact_row)
            comptox_rows.extend(comptox.search_by_name(name, exact=False, limit=5))
            if comptox_rows:
                best = _best_candidate_by_name(comptox_rows, name, "PREFERRED_NAME")
                if best is not None:
                    candidates["comptox"] = _candidate_from_comptox_row(best)
        except Exception as e:
            logging.warning(f"CompTox name lookup failed for '{name}': {e}")

    if pubchem is not None:
        try:
            pubchem_rows = pubchem.search_by_name(name, exact=True, limit=5)
            if not pubchem_rows:
                pubchem_rows = pubchem.search_by_name(name, exact=False, limit=5)
            if pubchem_rows:
                best = _best_candidate_by_name(pubchem_rows, name, "cmpdname")
                if best is not None:
                    candidates["pubchem"] = _candidate_from_pubchem_row(best)
        except Exception as e:
            logging.warning(f"PubChemID name lookup failed for '{name}': {e}")

    if zeropm is not None:
        try:
            zeropm_table = zeropm.get_id_table_from_name(name)
            candidates["zeropm"] = _candidate_from_zeropm_name_table(name, zeropm_table)
        except Exception as e:
            logging.warning(f"ZeroPM name lookup failed for '{name}': {e}")

    if chembl is not None:
        try:
            chembl_rows = chembl.search_by_name(name, limit=5)
            if chembl_rows:
                best = _best_candidate_by_name(chembl_rows, name, "pref_name")
                if best is not None:
                    candidates["chembl"] = _candidate_from_chembl_row(best, chembl=chembl)
        except Exception as e:
            logging.warning(f"ChEMBL name lookup failed for '{name}': {e}")

    consensus_source, source_match_scores, match_score = _compute_consensus(candidates)
    consensus_candidate = candidates.get(consensus_source) if consensus_source else None

    for source_key in ["chebi", "comptox", "pubchem", "zeropm"]:
        candidate = candidates.get(source_key)
        if _candidate_compatible_with_consensus(candidate, consensus_candidate):
            _apply_candidate_to_result(result, candidate)

    # Keep ChEMBL as a local enrichment layer after primary source population.
    chembl_candidate = candidates.get("chembl")
    if _candidate_compatible_with_consensus(chembl_candidate, consensus_candidate):
        _apply_candidate_to_result(result, chembl_candidate)

    canonical_smiles, rdkit_mass = _smiles_to_canonical_and_mass(result["SMILES"])
    result["canonical_smiles"] = canonical_smiles
    result["molecular_mass"] = _pick_first(result["molecular_mass"], rdkit_mass)
    result["name"] = _pick_first(result["name"], result["IUPAC_name"])
    result["IUPAC_name"] = _pick_first(result["IUPAC_name"], result["name"])
    result["source"] = _pick_first(result["source"], consensus_candidate.get("source") if consensus_candidate else None)

    result["consensus_source"] = _pick_first(
        consensus_candidate.get("source") if consensus_candidate else None,
        None,
    )
    result["source_match_scores"] = {
        (
            candidates[src].get("source") if candidates.get(src) else src
        ): round(score, 4)
        for src, score in source_match_scores.items()
    }
    result["match_score"] = round(match_score, 4)

    return result


def ids_from_SMILES(
    smiles: str,
    chebi: Optional[ChebiSDF] = None,
    comptox: Optional[CompToxID] = None,
    pubchem: Optional[PubChemID] = None,
    zeropm: Optional[ZeroPM] = None,
    chembl: Optional[CheMBL] = None,
) -> Dict[str, Any]:
    """
    Resolve one SMILES input into a unified identifier record from offline sources.

    The function queries the offline datasets (ChEBI, CompTox, PubChemID, ZeroPM),
    computes a cross-source consensus score based on structure and identifier
    agreement, and then populates output fields using the same source priority as
    ids_from_CAS.

    Args:
        smiles: SMILES input string.
        chebi: Optional initialized ChebiSDF client.
        comptox: Optional initialized CompToxID client.
        pubchem: Optional initialized PubChemID client.
        zeropm: Optional initialized ZeroPM client.
        chembl: Optional initialized CheMBL client.

    Returns:
        Dictionary similar to ids_from_CAS with additional consensus fields:
        match_score, consensus_source, and source_match_scores.
    """
    result: Dict[str, Any] = {
        "query": smiles,
        "CASRN": None,
        "name": None,
        "IUPAC_name": None,
        "molecular_formula": None,
        "SMILES": smiles,
        "canonical_smiles": None,
        "InChI": None,
        "InChIKey": None,
        "DTXSID": None,
        "molecular_mass": None,
        "Synonyms": None,
        "foundby": "SMILES",
        "source": None,
        "match_score": 0.0,
        "consensus_source": None,
        "source_match_scores": {},
    }

    candidates: Dict[str, Optional[Dict[str, Any]]] = {
        "chebi": None,
        "comptox": None,
        "pubchem": None,
        "zeropm": None,
        "chembl": None,
    }

    canonical_input, _ = _smiles_to_canonical_and_mass(smiles)
    inchikey_input = _inchikey_from_smiles(smiles)

    if chebi is not None and not _is_missing(inchikey_input):
        try:
            chebi_row = chebi.search_by_inchikey(str(inchikey_input))
            if chebi_row:
                candidates["chebi"] = _candidate_from_chebi_row(chebi_row)
        except Exception as e:
            logging.warning(f"ChEBI SMILES lookup failed for '{smiles}': {e}")

    if comptox is not None:
        try:
            comptox_row = comptox.get_by_smiles(smiles)
            if comptox_row is None and not _is_missing(canonical_input):
                comptox_row = comptox.get_by_smiles(str(canonical_input))
            if comptox_row:
                candidates["comptox"] = _candidate_from_comptox_row(comptox_row)
        except Exception as e:
            logging.warning(f"CompTox SMILES lookup failed for '{smiles}': {e}")

    if pubchem is not None:
        try:
            pubchem_row = pubchem.get_by_smiles(smiles)
            if pubchem_row is None and not _is_missing(canonical_input):
                pubchem_row = pubchem.get_by_smiles(str(canonical_input))
            if pubchem_row:
                candidates["pubchem"] = _candidate_from_pubchem_row(pubchem_row)
        except Exception as e:
            logging.warning(f"PubChemID SMILES lookup failed for '{smiles}': {e}")

    if zeropm is not None:
        try:
            candidates["zeropm"] = _candidate_from_zeropm_smiles(smiles, zeropm)
        except Exception as e:
            logging.warning(f"ZeroPM SMILES lookup failed for '{smiles}': {e}")

    if chembl is not None:
        try:
            chembl_row = chembl.search_by_smiles(smiles)
            if chembl_row is None and not _is_missing(canonical_input):
                chembl_row = chembl.search_by_smiles(str(canonical_input))
            if chembl_row:
                candidates["chembl"] = _candidate_from_chembl_row(chembl_row, chembl=chembl)
        except Exception as e:
            logging.warning(f"ChEMBL SMILES lookup failed for '{smiles}': {e}")

    consensus_source, source_match_scores, match_score = _compute_consensus(candidates)
    consensus_candidate = candidates.get(consensus_source) if consensus_source else None

    for source_key in ["chebi", "comptox", "pubchem", "zeropm"]:
        candidate = candidates.get(source_key)
        if _candidate_compatible_with_consensus(candidate, consensus_candidate):
            _apply_candidate_to_result(result, candidate)

    # Keep ChEMBL as a local enrichment layer after primary source population.
    chembl_candidate = candidates.get("chembl")
    if _candidate_compatible_with_consensus(chembl_candidate, consensus_candidate):
        _apply_candidate_to_result(result, chembl_candidate)

    canonical_smiles, rdkit_mass = _smiles_to_canonical_and_mass(result["SMILES"])
    result["canonical_smiles"] = canonical_smiles
    result["molecular_mass"] = _pick_first(result["molecular_mass"], rdkit_mass)
    result["name"] = _pick_first(result["name"], result["IUPAC_name"])
    result["IUPAC_name"] = _pick_first(result["IUPAC_name"], result["name"])
    result["source"] = _pick_first(result["source"], consensus_candidate.get("source") if consensus_candidate else None)

    result["consensus_source"] = _pick_first(
        consensus_candidate.get("source") if consensus_candidate else None,
        None,
    )
    result["source_match_scores"] = {
        (
            candidates[src].get("source") if candidates.get(src) else src
        ): round(score, 4)
        for src, score in source_match_scores.items()
    }
    result["match_score"] = round(match_score, 4)

    return result


def _smiles_to_canonical_and_mass(smiles: Optional[str]) -> Tuple[Optional[str], Optional[float]]:
    """Generate canonical SMILES and molecular mass from one RDKit parse."""
    if _is_missing(smiles):
        return None, None

    if not RDKIT_AVAILABLE or Chem is None:
        return str(smiles), None

    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return None, None
        canonical = Chem.MolToSmiles(mol, canonical=True)
        mass = float(Descriptors.MolWt(mol)) if Descriptors is not None else None
        return canonical, mass
    except Exception as e:
        logging.warning(f"Failed to parse SMILES '{smiles}': {e}")
        return None, None


def ids_from_CAS(
    cas: str,
    chebi: Optional[ChebiSDF] = None,
    comptox: Optional[CompToxID] = None,
    pubchem: Optional[PubChemID] = None,
    zeropm: Optional[ZeroPM] = None,
    chembl: Optional[CheMBL] = None,
) -> Dict[str, Any]:
    """
    Build a unified identifier record for one CAS number using offline databases.

    The lookup strategy is offline-first and follows the source priority:
    ChEBI -> CompTox -> PubChemID -> ZeroPM. ChEMBL is used as a local
    enrichment step when a SMILES string is available.

    Args:
        cas: CAS Registry Number.
        chebi: Optional initialized ChebiSDF client.
        comptox: Optional initialized CompToxID client.
        pubchem: Optional initialized PubChemID client.
        zeropm: Optional initialized ZeroPM client.
        chembl: Optional initialized CheMBL client.

    Returns:
        Dictionary with standardized fields:
        CASRN, name, IUPAC_name, molecular_formula, SMILES, canonical_smiles,
        InChI, InChIKey, DTXSID, molecular_mass, Synonyms, foundby, source.
    """
    result: Dict[str, Any] = {
        "CASRN": cas,
        "name": None,
        "IUPAC_name": None,
        "molecular_formula": None,
        "SMILES": None,
        "canonical_smiles": None,
        "InChI": None,
        "InChIKey": None,
        "DTXSID": None,
        "molecular_mass": None,
        "Synonyms": None,
        "foundby": "CASRN",
        "source": None,
    }

    # 1) ChEBI
    if chebi is not None:
        try:
            chebi_matches = chebi.search_by_cas(cas)
            if chebi_matches:
                row = chebi_matches[0]
                result["name"] = _pick_first(result["name"], row.get("ChEBI NAME"))
                result["IUPAC_name"] = _pick_first(result["IUPAC_name"], row.get("ChEBI NAME"))
                result["molecular_formula"] = _pick_first(result["molecular_formula"], row.get("FORMULA"))
                result["SMILES"] = _pick_first(result["SMILES"], row.get("SMILES"))
                result["InChI"] = _pick_first(result["InChI"], row.get("INCHI"))
                result["InChIKey"] = _pick_first(result["InChIKey"], row.get("INCHIKEY"))
                result["Synonyms"] = _pick_first(result["Synonyms"], _normalize_synonyms(row.get("SYNONYM")))
                if not _is_missing(result["SMILES"]):
                    result["source"] = "ChEBI"
        except Exception as e:
            logging.warning(f"ChEBI lookup failed for CAS {cas}: {e}")

    # 2) CompTox
    if _is_missing(result["SMILES"]) and comptox is not None:
        try:
            comptox_row = comptox.get_by_casrn(cas)
            if comptox_row:
                result["name"] = _pick_first(result["name"], comptox_row.get("PREFERRED_NAME"))
                result["IUPAC_name"] = _pick_first(result["IUPAC_name"], comptox_row.get("IUPAC_NAME"))
                result["molecular_formula"] = _pick_first(result["molecular_formula"], comptox_row.get("MOLECULAR_FORMULA"))
                result["SMILES"] = _pick_first(result["SMILES"], comptox_row.get("SMILES"))
                result["InChI"] = _pick_first(result["InChI"], comptox_row.get("INCHI"))
                result["InChIKey"] = _pick_first(result["InChIKey"], comptox_row.get("INCHIKEY"))
                result["DTXSID"] = _pick_first(result["DTXSID"], comptox_row.get("DTXSID"))
                result["Synonyms"] = _pick_first(result["Synonyms"], _normalize_synonyms(comptox_row.get("identifiers")))
                result["molecular_mass"] = _pick_first(
                    result["molecular_mass"],
                    comptox_row.get("AVERAGE_MASS"),
                    comptox_row.get("MONOISOTOPIC_MASS"),
                )
                if not _is_missing(result["SMILES"]):
                    result["source"] = "CompTox"
        except Exception as e:
            logging.warning(f"CompTox lookup failed for CAS {cas}: {e}")

    # 3) PubChemID (offline SQLite)
    if _is_missing(result["SMILES"]) and pubchem is not None:
        try:
            pubchem_row = pubchem.get_by_cas(cas)
            if pubchem_row:
                result["name"] = _pick_first(result["name"], pubchem_row.get("cmpdname"))
                result["IUPAC_name"] = _pick_first(result["IUPAC_name"], pubchem_row.get("iupacname"))
                result["molecular_formula"] = _pick_first(result["molecular_formula"], pubchem_row.get("mf"))
                result["SMILES"] = _pick_first(result["SMILES"], pubchem_row.get("smiles"))
                result["InChI"] = _pick_first(result["InChI"], pubchem_row.get("inchi"))
                result["InChIKey"] = _pick_first(result["InChIKey"], pubchem_row.get("inchikey"))
                result["Synonyms"] = _pick_first(result["Synonyms"], _normalize_synonyms(pubchem_row.get("synonyms")))
                result["molecular_mass"] = _pick_first(result["molecular_mass"], pubchem_row.get("mw"))
                if not _is_missing(result["SMILES"]):
                    result["source"] = "PubChemID"
        except Exception as e:
            logging.warning(f"PubChemID lookup failed for CAS {cas}: {e}")

    # 4) ZeroPM (convert first InChI to SMILES only if needed)
    if _is_missing(result["SMILES"]) and zeropm is not None:
        try:
            zeropm_table = zeropm.get_id_table_from_cas(cas)
            if zeropm_table is not None and not zeropm_table.empty:
                row = zeropm_table.iloc[0]
                inchi = _pick_first(result["InChI"], row.get("inchi"))
                inchikey = _pick_first(result["InChIKey"], row.get("inchikey"))
                synonyms = _pick_first(result["Synonyms"], _normalize_synonyms(row.get("synonyms")))

                smiles_from_inchi = None
                if not _is_missing(inchi) and RDKIT_AVAILABLE and Chem is not None:
                    mol = Chem.MolFromInchi(str(inchi))
                    if mol is not None:
                        smiles_from_inchi = Chem.MolToSmiles(mol)

                result["InChI"] = inchi
                result["InChIKey"] = inchikey
                result["Synonyms"] = synonyms
                result["SMILES"] = _pick_first(result["SMILES"], smiles_from_inchi)
                if not _is_missing(result["SMILES"]):
                    result["source"] = "ZeroPM"
        except Exception as e:
            logging.warning(f"ZeroPM lookup failed for CAS {cas}: {e}")

    # 5) ChEMBL enrichment from SMILES (still offline)
    if not _is_missing(result["SMILES"]) and chembl is not None:
        try:
            chembl_row = chembl.search_by_smiles(str(result["SMILES"]))
            if chembl_row:
                result["name"] = _pick_first(result["name"], chembl_row.get("pref_name"))
                result["InChI"] = _pick_first(result["InChI"], chembl_row.get("standard_inchi"))
                result["InChIKey"] = _pick_first(result["InChIKey"], chembl_row.get("standard_inchi_key"))
                result["Synonyms"] = _pick_first(result["Synonyms"], _normalize_synonyms(chembl_row.get("synonyms")))

                molregno = chembl_row.get("molregno")
                if not _is_missing(molregno):
                    props = chembl.get_properties(int(molregno))
                    if props:
                        result["molecular_mass"] = _pick_first(result["molecular_mass"], props.get("mw_freebase"))
        except Exception as e:
            logging.warning(f"ChEMBL enrichment failed for CAS {cas}: {e}")

    canonical_smiles, rdkit_mass = _smiles_to_canonical_and_mass(result["SMILES"])
    result["canonical_smiles"] = canonical_smiles
    result["molecular_mass"] = _pick_first(result["molecular_mass"], rdkit_mass)

    # Keep IUPAC_name and name aligned when only one is available.
    result["name"] = _pick_first(result["name"], result["IUPAC_name"])
    result["IUPAC_name"] = _pick_first(result["IUPAC_name"], result["name"])

    return result

def casrn_to_compounds(cas_rn: Union[str, List[str]], ccc_api_key: Optional[str] = None,
                      show_progress: bool = True) -> pd.DataFrame:
    """
    Resolve CAS numbers into compound identifiers using offline databases only.

    This function uses the shared offline lookup pipeline implemented in
    ids_from_CAS with source priority ChEBI -> CompTox -> PubChemID -> ZeroPM,
    and then optionally enriches with ChEMBL using SMILES.
    
    Args:
        cas_rn: Single CAS RN string or list of CAS RN strings.
        ccc_api_key: Unused legacy parameter retained in the signature.
        show_progress: Whether to show progress bars.
        
    Returns:
        DataFrame with columns: CASRN, name, IUPAC_name, molecular_formula,
        SMILES, canonical_smiles, InChI, InChIKey, DTXSID, molecular_mass,
        Synonyms, foundby, source
    """
    if ccc_api_key is not None:
        logging.warning("casrn_to_compounds ignores ccc_api_key and uses offline databases only.")

    # Convert input to list
    if isinstance(cas_rn, str):
        casrn_list = [cas_rn]
    else:
        casrn_list = cas_rn

    # Initialize offline databases once to avoid repeated connection overhead.
    clients: Dict[str, Any] = {
        "chebi": None,
        "comptox": None,
        "pubchem": None,
        "zeropm": None,
        "chembl": None,
    }

    for name, factory in [
        ("chebi", ChebiSDF),
        ("comptox", CompToxID),
        ("pubchem", PubChemID),
        ("zeropm", ZeroPM),
        ("chembl", CheMBL),
    ]:
        try:
            clients[name] = factory()
        except Exception as e:
            logging.warning(f"Could not initialize offline source {name}: {e}")

    iterator = tqdm(casrn_list, desc="Resolving CASRN from offline databases") if show_progress else casrn_list
    rows = [
        ids_from_CAS(
            casrn,
            chebi=clients["chebi"],
            comptox=clients["comptox"],
            pubchem=clients["pubchem"],
            zeropm=clients["zeropm"],
            chembl=clients["chembl"],
        )
        for casrn in iterator
    ]

    if not rows:
        return pd.DataFrame(
            columns=[
                "CASRN",
                "name",
                "IUPAC_name",
                "molecular_formula",
                "SMILES",
                "canonical_smiles",
                "InChI",
                "InChIKey",
                "DTXSID",
                "molecular_mass",
                "Synonyms",
                "foundby",
                "source",
            ]
        )

    return pd.DataFrame(rows)