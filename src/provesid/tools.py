"""Candidate records and the consensus vote behind :class:`~provesid.Search`.

A *candidate* is one source's answer about one compound, normalised into a plain
dict so that a ChEBI row, a CompTox row and a ZeroPM row can be compared without
caring where each came from. :func:`make_candidate` builds one; the
``candidate_from_*`` adapters build one from a particular source's row shape.

:func:`compute_consensus` is the vote: it scores every candidate against every
other and returns the source whose answer the others corroborate best, together
with per-source agreement scores. That is what :class:`~provesid.Search` turns
into the ``confidence`` column and what ``min_source_support`` filters on.

The rest are the small predicates and converters those two need — missing-value
handling, CAS extraction, RDKit round-trips. They are public because
:mod:`provesid.search` imports them across the module boundary, not because
callers are expected to reach for them directly.
"""

from typing import List, Optional, Dict, Any, Tuple
import pandas as pd
import logging
import re
from difflib import SequenceMatcher

from .zeropm import ZeroPM
from .chembl import CheMBL

# Optional RDKit import
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    Chem = None
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available. Install with: pip install rdkit-pypi")



def is_missing(value: Any) -> bool:
    """Report whether a value carries no information.

    Sources disagree about how to say "nothing": ``None``, ``float('nan')``,
    an empty string, and the literal string ``"nan"`` all turn up in rows read
    from SQLite and from pandas. This treats all of them the same.

    Args:
        value: Any value read from a source row.

    Returns:
        True when the value is None, NaN, blank, or the string ``"nan"``.

    Examples:
        >>> is_missing(None), is_missing("nan"), is_missing("  ")
        (True, True, True)
        >>> is_missing(0)
        False
    """
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip() == "" or value.strip().lower() == "nan"
    try:
        return bool(pd.isna(value))
    except Exception:
        return False


def pick_first(*values: Any) -> Any:
    """Return the first argument that carries information.

    Used to fill a field from a preferred source, falling back through less
    preferred ones, without a chain of conditionals.

    Args:
        *values: Candidate values, most preferred first.

    Returns:
        The first value for which :func:`is_missing` is False, or None when
        every argument is missing.

    Examples:
        >>> pick_first(None, float("nan"), "aspirin", "ASA")
        'aspirin'
    """
    for value in values:
        if not is_missing(value):
            return value
    return None


def normalize_synonyms(value: Any) -> Optional[str]:
    """Render synonyms as one semicolon-separated string.

    Sources hand back synonyms as a list, a set, or an already-joined string.
    Candidate records store one string, so every shape collapses to the same
    representation before comparison.

    Args:
        value: A synonym collection or a string of synonyms.

    Returns:
        The synonyms joined by ``"; "``, or None when there are none.

    Examples:
        >>> normalize_synonyms(["aspirin", "ASA", None])
        'aspirin; ASA'
    """
    if is_missing(value):
        return None

    if isinstance(value, (list, tuple, set)):
        cleaned = [str(v).strip() for v in value if not is_missing(v)]
        return "; ".join(cleaned) if cleaned else None

    text = str(value).strip()
    return text if text else None


_CAS_PATTERN = re.compile(r"\b\d{2,7}-\d{2}-\d\b")


def to_float(value: Any) -> Optional[float]:
    """Convert a value to float, or to None when it will not convert.

    Molecular masses arrive as floats, as strings, and as NaN, sometimes in the
    same column. Comparisons need a float or nothing, never an exception.

    Args:
        value: The value to convert.

    Returns:
        The value as a float, or None when it is missing or unparseable.

    Examples:
        >>> to_float("180.16"), to_float("n/a")
        (180.16, None)
    """
    if is_missing(value):
        return None
    try:
        return float(value)
    except Exception:
        return None


def text_similarity(a: Optional[str], b: Optional[str]) -> float:
    """Score how alike two names are, ignoring case and surrounding space.

    A cheap ``difflib`` ratio, used only as a weak signal in
    :func:`candidate_similarity`: names corroborate a match but never decide
    one, because two sources routinely use different names for the same
    structure.

    Args:
        a: One name, or None.
        b: The other name, or None.

    Returns:
        1.0 for an exact match after normalisation, 0.0 when either side is
        missing, otherwise the ``SequenceMatcher`` ratio in [0, 1].

    Examples:
        >>> text_similarity("Aspirin", "aspirin ")
        1.0
        >>> round(text_similarity("aspirin", "asprin"), 2)
        0.92
    """
    if is_missing(a) or is_missing(b):
        return 0.0
    a_text = str(a).strip().lower()
    b_text = str(b).strip().lower()
    if a_text == b_text:
        return 1.0
    return SequenceMatcher(None, a_text, b_text).ratio()


def extract_cas_values(value: Any) -> List[str]:
    r"""Find every CAS Registry Number anywhere inside a value.

    Walks dicts, lists, tuples and sets recursively and pattern-matches the
    text of everything else, so an entire source row can be handed over
    without knowing which of its columns holds a CAS.

    The pattern is structural (``\d{2,7}-\d{2}-\d``) and does **not** verify
    the check digit, so it can pick up a number-shaped string that is not a
    registered CAS.

    Args:
        value: A row, a collection, or a single value of any type.

    Returns:
        The distinct CAS-shaped strings found, sorted, or an empty list.

    Examples:
        >>> extract_cas_values({"CASRN": "50-78-2", "syn": ["ASA", "50-78-2"]})
        ['50-78-2']
    """
    found: List[str] = []

    if value is None:
        return found

    if isinstance(value, dict):
        for dict_value in value.values():
            found.extend(extract_cas_values(dict_value))
    elif isinstance(value, (list, tuple, set)):
        for item in value:
            found.extend(extract_cas_values(item))
    else:
        text = str(value)
        found.extend(_CAS_PATTERN.findall(text))

    deduped = sorted(set(found))
    return deduped


def inchi_to_smiles(inchi: Optional[str]) -> Optional[str]:
    """Convert an InChI string to SMILES.

    Args:
        inchi: The InChI string, or None.

    Returns:
        The SMILES string, or None when the input is missing, RDKit is not
        installed, or RDKit cannot parse the InChI.

    Examples:
        >>> inchi_to_smiles("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
        'CCO'
        >>> inchi_to_smiles(None) is None
        True
    """
    if is_missing(inchi) or not RDKIT_AVAILABLE or Chem is None:
        return None
    try:
        mol = Chem.MolFromInchi(str(inchi))
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def inchikey_from_smiles(smiles: Optional[str]) -> Optional[str]:
    """Derive an InChIKey from a SMILES string.

    Lets a source that publishes a structure but no InChIKey still be matched
    against one that publishes the key, which is how most cross-source
    agreement is actually established.

    Args:
        smiles: The SMILES string, or None.

    Returns:
        The InChIKey, or None when the input is missing, RDKit is not
        installed, or RDKit cannot parse the SMILES.

    Examples:
        >>> inchikey_from_smiles("OCC")
        'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
    """
    if is_missing(smiles) or not RDKIT_AVAILABLE or Chem is None:
        return None
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return None
        inchi = Chem.MolToInchi(mol)
        if is_missing(inchi):
            return None
        return Chem.InchiToInchiKey(inchi)
    except Exception:
        return None


def first_cas(cas_values: List[str]) -> Optional[str]:
    """Pick one CAS number out of a candidate's list.

    The list from :func:`extract_cas_values` is sorted, so this is stable
    across runs rather than dependent on row order.

    Args:
        cas_values: CAS numbers, as returned by :func:`extract_cas_values`.

    Returns:
        The first CAS number, or None when the list is empty.

    Note:
        "First" is first as a string, which is not the best number: aspirin's
        CompTox row sorts the retired ``11126-35-5`` ahead of ``50-78-2``.

    Examples:
        >>> first_cas(extract_cas_values("50-78-2 | 11126-35-5"))
        '11126-35-5'
        >>> first_cas([]) is None
        True
    """
    return cas_values[0] if cas_values else None


def make_candidate(
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
    """Build one source's answer in the shape every comparison expects.

    A candidate is a plain dict with a fixed set of keys, so a ChEBI row and a
    ZeroPM row can be scored against each other without either side knowing
    where the other came from. The SMILES is canonicalised on the way in, and
    the molecular mass is taken from the source when it gives one and computed
    from the structure when it does not — both so that two sources stating the
    same compound differently still compare equal.

    Args:
        source: Display name of the source, e.g. ``"ChEBI"``. This is what
            appears in ``source_details`` and in the consensus report.
        name: The source's preferred name for the compound.
        iupac_name: The IUPAC name, where the source distinguishes it.
        molecular_formula: The molecular formula as the source states it.
        smiles: The structure as SMILES.
        inchi: The structure as InChI.
        inchikey: The InChIKey.
        dtxsid: The DSSTox identifier, for sources that carry one.
        molecular_mass: The mass the source states; falls back to the mass
            RDKit computes from ``smiles``.
        synonyms: Synonyms, already flattened by :func:`normalize_synonyms`.
        cas_candidates: Every CAS the row mentions, deduplicated and sorted.

    Returns:
        The candidate record: a dict with the keys ``source``, ``name``,
        ``IUPAC_name``, ``molecular_formula``, ``SMILES``,
        ``canonical_smiles``, ``InChI``, ``InChIKey``, ``DTXSID``,
        ``molecular_mass``, ``Synonyms`` and ``CAS_candidates``.

    Examples:
        >>> cand = make_candidate("ChEBI", name="aspirin", smiles="CC(=O)Oc1ccccc1C(=O)O")
        >>> cand["canonical_smiles"]
        'CC(=O)Oc1ccccc1C(=O)O'
        >>> round(cand["molecular_mass"], 2)
        180.16
    """
    canonical_smiles, rdkit_mass = smiles_to_canonical_and_mass(smiles)
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
        "molecular_mass": pick_first(to_float(molecular_mass), rdkit_mass),
        "Synonyms": synonyms,
        "CAS_candidates": sorted(set(cas_candidates or [])),
    }


def candidate_similarity(left: Dict[str, Any], right: Dict[str, Any]) -> float:
    """Score how strongly two candidates agree that they describe one compound.

    Each field the two candidates *both* carry contributes its weight to the
    denominator and, when the values match, to the numerator. Fields only one
    side has are ignored entirely, so a sparse source is neither rewarded nor
    punished for its silence — it simply has less to say.

    The weights rank the evidence: canonical SMILES (4) above CAS overlap and
    InChIKey (3 each), above InChI (2) and mass agreement (2), above formula
    (1) and name similarity (1). Mass and name score partially — a mass within
    0.2 scores full, within 1.0 scores half; a name similarity of 0.9 scores
    full, 0.7 scores half.

    Args:
        left: One candidate record.
        right: The other candidate record.

    Returns:
        Weighted agreement in [0, 1]. Returns 0.0 when either side is None or
        when the two share no comparable field at all — note that "no shared
        evidence" and "shared evidence that disagrees" both come back as 0.0.

    Examples:
        >>> a = make_candidate("ChEBI", smiles="CC(=O)Oc1ccccc1C(=O)O")
        >>> b = make_candidate("CompTox", smiles="CC(=O)Oc1ccccc1C(=O)O")
        >>> candidate_similarity(a, b)
        1.0
    """
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
    if not is_missing(left_smiles) and not is_missing(right_smiles):
        weight += 4.0
        if left_smiles == right_smiles:
            score += 4.0

    left_ik = left.get("InChIKey")
    right_ik = right.get("InChIKey")
    if not is_missing(left_ik) and not is_missing(right_ik):
        weight += 3.0
        if str(left_ik) == str(right_ik):
            score += 3.0

    left_inchi = left.get("InChI")
    right_inchi = right.get("InChI")
    if not is_missing(left_inchi) and not is_missing(right_inchi):
        weight += 2.0
        if str(left_inchi) == str(right_inchi):
            score += 2.0

    left_formula = left.get("molecular_formula")
    right_formula = right.get("molecular_formula")
    if not is_missing(left_formula) and not is_missing(right_formula):
        weight += 1.0
        if str(left_formula) == str(right_formula):
            score += 1.0

    left_mass = to_float(left.get("molecular_mass"))
    right_mass = to_float(right.get("molecular_mass"))
    if left_mass is not None and right_mass is not None:
        weight += 2.0
        diff = abs(left_mass - right_mass)
        if diff <= 0.2:
            score += 2.0
        elif diff <= 1.0:
            score += 1.0

    name_sim = text_similarity(left.get("name"), right.get("name"))
    if name_sim > 0.0:
        weight += 1.0
        if name_sim >= 0.9:
            score += 1.0
        elif name_sim >= 0.7:
            score += 0.5

    if weight == 0.0:
        return 0.0
    return score / weight




def candidate_compatible_with_consensus(
    candidate: Optional[Dict[str, Any]],
    consensus: Optional[Dict[str, Any]],
    threshold: float = 0.35,
) -> bool:
    """Decide whether a candidate may contribute to a result the consensus anchors.

    A source that disagrees with the consensus is describing a different
    compound, and letting it fill empty fields would assemble one record out
    of two substances. This is the gate that keeps that from happening.

    Args:
        candidate: The candidate under consideration, or None.
        consensus: The consensus candidate to measure against, or None when
            no consensus was reached.
        threshold: Minimum :func:`candidate_similarity` required. The default
            of 0.35 is permissive by design: it rejects a different compound
            without rejecting a sparse source that agrees on what little it
            states.

    Returns:
        True when the candidate agrees with the consensus closely enough, when
        it *is* the consensus source, or when there is no consensus to
        contradict. False when the candidate is None.

    Examples:
        >>> aspirin = make_candidate("ChEBI", smiles="CC(=O)Oc1ccccc1C(=O)O")
        >>> also_aspirin = make_candidate("CompTox", smiles="CC(=O)OC1=C(C=CC=C1)C(O)=O")
        >>> ethanol = make_candidate("ZeroPM", smiles="CCO")
        >>> candidate_compatible_with_consensus(also_aspirin, aspirin)
        True
        >>> candidate_compatible_with_consensus(ethanol, aspirin)
        False
    """
    if candidate is None:
        return False
    if consensus is None:
        return True
    if candidate.get("source") == consensus.get("source"):
        return True
    return candidate_similarity(candidate, consensus) >= threshold


def apply_candidate_to_result(result: Dict[str, Any], candidate: Optional[Dict[str, Any]]) -> None:
    """Fill a result's empty fields from a candidate, in place.

    Never overwrites: a field already carrying a value is left alone, so
    applying candidates in priority order means the most trusted source that
    had something to say wins each field independently. A result can therefore
    take its structure from one source and its CAS from another.

    Args:
        result: The result dict to fill, modified in place.
        candidate: The candidate to read from. None is a no-op.

    Returns:
        None. The mutation is the point.

    Examples:
        >>> result = {"CASRN": "50-78-2", "name": None}
        >>> apply_candidate_to_result(result, make_candidate(
        ...     "CompTox", name="Aspirin", smiles="CC(=O)OC1=C(C=CC=C1)C(O)=O",
        ...     cas_candidates=["11126-35-5"]))
        >>> result["CASRN"], result["name"], result["source"]
        ('50-78-2', 'Aspirin', 'CompTox')
    """
    if candidate is None:
        return

    result["CASRN"] = pick_first(result.get("CASRN"), first_cas(candidate.get("CAS_candidates") or []))
    result["name"] = pick_first(result.get("name"), candidate.get("name"))
    result["IUPAC_name"] = pick_first(result.get("IUPAC_name"), candidate.get("IUPAC_name"))
    result["molecular_formula"] = pick_first(result.get("molecular_formula"), candidate.get("molecular_formula"))
    result["SMILES"] = pick_first(result.get("SMILES"), candidate.get("SMILES"))
    result["InChI"] = pick_first(result.get("InChI"), candidate.get("InChI"))
    result["InChIKey"] = pick_first(result.get("InChIKey"), candidate.get("InChIKey"))
    result["DTXSID"] = pick_first(result.get("DTXSID"), candidate.get("DTXSID"))
    result["molecular_mass"] = pick_first(result.get("molecular_mass"), candidate.get("molecular_mass"))
    result["Synonyms"] = pick_first(result.get("Synonyms"), normalize_synonyms(candidate.get("Synonyms")))

    if is_missing(result.get("source")) and not is_missing(candidate.get("SMILES")):
        result["source"] = candidate.get("source")


def compute_consensus(candidates: Dict[str, Optional[Dict[str, Any]]]) -> Tuple[Optional[str], Dict[str, float], float]:
    """Hold the vote: which source's answer do the others corroborate?

    Every candidate is scored against every other with
    :func:`candidate_similarity` and given the mean of those scores as its
    support. The winner is the best-supported source — but among sources
    within 0.05 of the top score, the more reputable one wins instead. That
    tie-break matters because support is an average over *comparable* fields:
    a source stating almost nothing can agree perfectly on that little and
    score higher than a richer source that agrees about far more.

    Reputation order is ChEBI, CompTox, PubChemID, ZeroPM, ChEMBL; a source
    not on that list sorts last.

    Args:
        candidates: Source key to candidate record. Entries whose value is
            None are ignored, so a source that found nothing does not vote.

    Returns:
        A tuple of:

        - the winning source key, or None when no source found anything;
        - per-source agreement with the winner, the winner itself scoring 1.0;
        - the mean of those scores, which is what becomes ``confidence``.

    Examples:
        >>> a = make_candidate("ChEBI", smiles="CC(=O)Oc1ccccc1C(=O)O")
        >>> b = make_candidate("CompTox", smiles="CC(=O)Oc1ccccc1C(=O)O")
        >>> source, scores, overall = compute_consensus({"chebi": a, "comptox": b})
        >>> source, overall
        ('chebi', 1.0)
    """
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
        sims = [candidate_similarity(valid[source], valid[other]) for other in others]
        support[source] = sum(sims) / len(sims)

    priority = ["chebi", "comptox", "pubchem", "zeropm", "chembl"]
    # Among sources within _PRIORITY_TOLERANCE of the top support score, prefer
    # by reputation (priority list).  This prevents a source missing structural
    # fields (e.g. no SMILES) from artificially inflating its support score
    # and winning over a more reputable source with essentially the same agreement.
    _PRIORITY_TOLERANCE = 0.05
    max_support = max(support.values())
    top_sources = [
        src for src in support if max_support - support[src] <= _PRIORITY_TOLERANCE
    ]
    top_sources.sort(
        key=lambda src: priority.index(src) if src in priority else len(priority)
    )
    consensus_source = top_sources[0]

    consensus_candidate = valid[consensus_source]
    source_match_scores = {
        src: (1.0 if src == consensus_source else candidate_similarity(consensus_candidate, valid[src]))
        for src in sources
    }
    overall = sum(source_match_scores.values()) / len(source_match_scores)
    return consensus_source, source_match_scores, overall


def candidate_from_chebi_row(row: Dict[str, Any]) -> Dict[str, Any]:
    """Adapt one ChEBI SDF row into a candidate record.

    Args:
        row: A row as :class:`~provesid.ChebiSDF` returns it.

    Returns:
        The candidate record. ChEBI states no mass, so the mass comes from
        RDKit via :func:`make_candidate`.

    Examples:
        >>> from provesid import ChebiSDF
        >>> row = ChebiSDF().get_compound_by_id("CHEBI:15365")   # doctest: +SKIP
        >>> cand = candidate_from_chebi_row(row)                 # doctest: +SKIP
        >>> cand["name"], cand["CAS_candidates"], round(cand["molecular_mass"], 2)  # doctest: +SKIP
        ('acetylsalicylic acid', ['50-78-2'], 180.16)
    """
    return make_candidate(
        "ChEBI",
        name=row.get("ChEBI NAME"),
        iupac_name=row.get("ChEBI NAME"),
        molecular_formula=row.get("FORMULA"),
        smiles=row.get("SMILES"),
        inchi=row.get("INCHI"),
        inchikey=row.get("INCHIKEY"),
        synonyms=normalize_synonyms(row.get("SYNONYM")),
        cas_candidates=extract_cas_values(row),
    )


def candidate_from_comptox_row(row: Dict[str, Any]) -> Dict[str, Any]:
    """Adapt one CompTox row into a candidate record.

    Args:
        row: A row as :class:`~provesid.CompToxID` returns it.

    Returns:
        The candidate record, carrying the DTXSID and preferring the average
        mass over the monoisotopic one.

    Examples:
        >>> from provesid import CompToxID
        >>> cand = candidate_from_comptox_row(CompToxID().get_by_casrn("50-78-2"))  # doctest: +SKIP
        >>> cand["DTXSID"], cand["molecular_mass"]               # doctest: +SKIP
        ('DTXSID5020108', 180.159)
    """
    return make_candidate(
        "CompTox",
        name=row.get("PREFERRED_NAME"),
        iupac_name=row.get("IUPAC_NAME"),
        molecular_formula=row.get("MOLECULAR_FORMULA"),
        smiles=row.get("SMILES"),
        inchi=row.get("INCHI"),
        inchikey=row.get("INCHIKEY"),
        dtxsid=row.get("DTXSID"),
        molecular_mass=pick_first(row.get("AVERAGE_MASS"), row.get("MONOISOTOPIC_MASS")),
        synonyms=normalize_synonyms(row.get("identifiers")),
        cas_candidates=extract_cas_values([row.get("CASRN"), row.get("identifiers")]),
    )


def candidate_from_pubchem_row(row: Dict[str, Any]) -> Dict[str, Any]:
    """Adapt one PubChem row into a candidate record.

    Args:
        row: A row as :class:`~provesid.PubChemID` returns it.

    Returns:
        The candidate record.

    Examples:
        >>> from provesid import PubChemID
        >>> cand = candidate_from_pubchem_row(PubChemID().get_by_cid(2244))  # doctest: +SKIP
        >>> cand["name"], cand["molecular_mass"], cand["CAS_candidates"]     # doctest: +SKIP
        ('Aspirin', 180.16, ['50-78-2'])
    """
    return make_candidate(
        "PubChemID",
        name=row.get("cmpdname"),
        iupac_name=row.get("iupacname"),
        molecular_formula=row.get("mf"),
        smiles=row.get("smiles"),
        inchi=row.get("inchi"),
        inchikey=row.get("inchikey"),
        molecular_mass=row.get("mw"),
        synonyms=normalize_synonyms(row.get("synonyms")),
        cas_candidates=extract_cas_values(row.get("cas_numbers")),
    )


def candidate_from_zeropm_name_table(name: str, table: pd.DataFrame) -> Optional[Dict[str, Any]]:
    """Adapt a ZeroPM name-lookup table into a single candidate record.

    ZeroPM answers a name with a ranked table rather than a row. The
    best-ranked entry supplies the structure; every name in the table becomes
    a synonym and every CAS a candidate CAS, which is what makes ZeroPM a
    useful corroborator of *identifiers* even where its structures are thin.
    ZeroPM publishes InChI but not SMILES, so the SMILES is derived.

    Args:
        name: The queried name, kept as the candidate's name.
        table: The lookup table, sorted by ``rank`` when that column exists.

    Returns:
        The candidate record, or None when the table is empty or None.

    Examples:
        >>> table = pd.DataFrame({"rank": [2, 1],
        ...                       "inchi": ["InChI=1S/CH4/h1H4", "InChI=1S/CH2O/c1-2/h1H2"],
        ...                       "inchikey": ["VNWKTOKETHGBQD-UHFFFAOYSA-N", "WSFSSNUMVMOOMR-UHFFFAOYSA-N"],
        ...                       "cas": ["74-82-8", "50-00-0"]})
        >>> cand = candidate_from_zeropm_name_table("Formaldehyde", table)
        >>> cand["SMILES"], cand["InChIKey"], cand["CAS_candidates"]
        ('C=O', 'WSFSSNUMVMOOMR-UHFFFAOYSA-N', ['50-00-0', '74-82-8'])
    """
    if table is None or table.empty:
        return None

    working = table.copy()
    if "rank" in working.columns:
        working = working.sort_values(by="rank", ascending=True)

    first = working.iloc[0]
    inchi = first.get("inchi")
    smiles = inchi_to_smiles(inchi)

    cas_values = []
    if "cas" in working.columns:
        cas_values = [str(v) for v in working["cas"].dropna().astype(str).tolist()]

    synonyms = None
    if "name" in working.columns and not working["name"].dropna().empty:
        synonyms = normalize_synonyms(working["name"].dropna().astype(str).unique().tolist())

    return make_candidate(
        "ZeroPM",
        name=name,
        iupac_name=name,
        smiles=smiles,
        inchi=inchi,
        inchikey=first.get("inchikey"),
        synonyms=synonyms,
        cas_candidates=extract_cas_values(cas_values),
    )


def candidate_from_zeropm_smiles(smiles_query: str, zeropm: ZeroPM) -> Optional[Dict[str, Any]]:
    """Adapt a ZeroPM structure lookup into a single candidate record.

    ZeroPM cannot be queried by structure directly. The SMILES is resolved to
    CAS numbers first, and up to five of those are looked up and pooled — a
    cap, because a structure that matches many registry entries would
    otherwise cost one query each for no added agreement.

    Args:
        smiles_query: The structure to look up, as SMILES.
        zeropm: An initialised :class:`~provesid.ZeroPM` client.

    Returns:
        The candidate record. When the CAS numbers resolve to no rows, a
        minimal candidate carrying just the query structure and those CAS
        numbers is returned instead — they are still evidence. None when the
        structure resolves to no CAS at all.

        The structure is taken from a row whose InChIKey is the query's, when
        there is one: the pooled CAS numbers include relatives, and for
        ``"CCO"`` the first is 13C-labelled ethanol.

    Examples:
        >>> from provesid import ZeroPM
        >>> cand = candidate_from_zeropm_smiles("CCO", ZeroPM())  # doctest: +SKIP
        >>> cand["InChIKey"], "64-17-5" in cand["CAS_candidates"]  # doctest: +SKIP
        ('LFQSCWFLJHTTHZ-UHFFFAOYSA-N', True)
    """
    cas_result = zeropm.get_cas_from_smiles(smiles_query)
    cas_values = extract_cas_values(cas_result)
    if not cas_values:
        return None

    tables = []
    for cas in cas_values[:5]:
        table = zeropm.get_id_table_from_cas(cas)
        if table is not None and not table.empty:
            tables.append(table)

    if not tables:
        return make_candidate("ZeroPM", smiles=smiles_query, cas_candidates=cas_values)

    combined = pd.concat(tables, ignore_index=True)
    # The CAS numbers come sorted, not ranked, and the first can belong to a
    # relative of the query: for "CCO" it is 14742-23-5, 13C-labelled ethanol.
    # Take the structure from a row that is the query, when one is.
    query_key = inchikey_from_smiles(smiles_query)
    same_structure = combined[combined["inchikey"] == query_key] if query_key else combined.iloc[0:0]
    first = same_structure.iloc[0] if not same_structure.empty else combined.iloc[0]
    inchi = first.get("inchi")
    smiles = pick_first(smiles_query, inchi_to_smiles(inchi))

    synonyms = None
    if "synonyms" in combined.columns and not combined["synonyms"].dropna().empty:
        synonyms = normalize_synonyms(combined["synonyms"].dropna().astype(str).unique().tolist())

    return make_candidate(
        "ZeroPM",
        smiles=smiles,
        inchi=inchi,
        inchikey=first.get("inchikey"),
        synonyms=synonyms,
        cas_candidates=extract_cas_values(cas_values),
    )


def candidate_from_chembl_row(row: Dict[str, Any], chembl: Optional[CheMBL] = None) -> Dict[str, Any]:
    """Adapt one ChEMBL row into a candidate record.

    Args:
        row: A row as :class:`~provesid.CheMBL` returns it.
        chembl: An optional client, used to fetch the molecular mass, which
            lives in a properties table rather than in the row. Without it the
            mass falls back to RDKit. A failed fetch is swallowed: the
            candidate is worth having without its mass.

    Returns:
        The candidate record. ChEMBL states no formula.

    Examples:
        >>> row = {"pref_name": "ASPIRIN", "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
        ...        "standard_inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        ...        "synonyms": ["Aspirin", "50-78-2"]}
        >>> cand = candidate_from_chembl_row(row)
        >>> cand["name"], cand["CAS_candidates"], round(cand["molecular_mass"], 2)
        ('ASPIRIN', ['50-78-2'], 180.16)
    """
    props = None
    molregno = row.get("molregno")
    if chembl is not None and not is_missing(molregno):
        try:
            props = chembl.get_properties(int(molregno))
        except Exception:
            props = None

    return make_candidate(
        "ChEMBL",
        name=row.get("pref_name"),
        molecular_formula=None,
        smiles=row.get("canonical_smiles"),
        inchi=row.get("standard_inchi"),
        inchikey=row.get("standard_inchi_key"),
        molecular_mass=(props or {}).get("mw_freebase"),
        synonyms=normalize_synonyms(row.get("synonyms")),
        cas_candidates=extract_cas_values(row.get("synonyms")),
    )


def candidate_from_pubchem_online(
    row: Dict[str, Any], synonyms: Optional[List[str]] = None
) -> Dict[str, Any]:
    """Adapt one PUG-REST property row into a candidate record.

    The online counterpart of :func:`candidate_from_pubchem_row`. It is kept
    apart from it, under its own source name, so that a result the network
    supplied can never be mistaken for one the local database did.

    Args:
        row: One row of :meth:`~provesid.PubChemAPI.get_properties_for_cids`,
            asked for ``Title``, ``IUPACName``, ``MolecularFormula``,
            ``SMILES``, ``InChI``, ``InChIKey`` and ``MolecularWeight``.
        synonyms: The compound's synonyms from
            :meth:`~provesid.PubChemAPI.get_compound_synonyms`, which is where
            PubChem keeps its CAS numbers.

    Returns:
        The candidate record, with source ``"PubChem (online)"``.

    Examples:
        >>> cand = candidate_from_pubchem_online(
        ...     {"CID": 2244, "Title": "Aspirin", "SMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
        ...      "MolecularWeight": "180.16"},
        ...     ["aspirin", "50-78-2"])
        >>> cand["source"], cand["CAS_candidates"], cand["molecular_mass"]
        ('PubChem (online)', ['50-78-2'], 180.16)
    """
    return make_candidate(
        "PubChem (online)",
        name=row.get("Title"),
        iupac_name=row.get("IUPACName"),
        molecular_formula=row.get("MolecularFormula"),
        smiles=row.get("SMILES"),
        inchi=row.get("InChI"),
        inchikey=row.get("InChIKey"),
        molecular_mass=row.get("MolecularWeight"),
        synonyms=normalize_synonyms(synonyms),
        cas_candidates=extract_cas_values(synonyms),
    )


def candidate_from_cactus(smiles: str, names: Optional[List[str]] = None) -> Dict[str, Any]:
    """Adapt an NCI/CADD Chemical Identifier Resolver answer into a candidate.

    CACTUS answers one representation per request, so the caller asks for the
    two that matter --- the structure and the name list --- and everything
    else is derived here: the InChIKey by RDKit, the CAS numbers from the
    names, among which CACTUS lists them.

    Args:
        smiles: The SMILES CACTUS resolved the identifier to.
        names: The ``names`` representation, one name per entry, most
            preferred first.

    Returns:
        The candidate record, with source ``"CACTUS"``.

    Examples:
        >>> cand = candidate_from_cactus("CC(=O)Oc1ccccc1C(O)=O", ["Aspirin", "50-78-2"])
        >>> cand["source"], cand["name"], cand["CAS_candidates"]
        ('CACTUS', 'Aspirin', ['50-78-2'])
    """
    names = [name for name in (names or []) if not is_missing(name)]
    return make_candidate(
        "CACTUS",
        name=names[0] if names else None,
        smiles=smiles,
        inchikey=inchikey_from_smiles(smiles),
        synonyms=normalize_synonyms(names),
        cas_candidates=extract_cas_values(names),
    )





def smiles_to_canonical_and_mass(smiles: Optional[str]) -> Tuple[Optional[str], Optional[float]]:
    """Canonicalise a SMILES string and weigh it in a single RDKit parse.

    Both are needed for every candidate, and parsing is the expensive part, so
    they are produced together.

    Args:
        smiles: The SMILES string, or None.

    Returns:
        A ``(canonical_smiles, molecular_mass)`` tuple. Both are None when the
        input is missing or RDKit cannot parse it. Without RDKit installed,
        the SMILES is passed through unchanged and the mass is None — an
        uncanonicalised structure still matches an identical string from
        another source.

    Examples:
        >>> smiles, mass = smiles_to_canonical_and_mass("OCC")
        >>> smiles, round(mass, 3)
        ('CCO', 46.069)
        >>> smiles_to_canonical_and_mass("not a smiles")
        (None, None)
    """
    if is_missing(smiles):
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


