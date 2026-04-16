"""PROVESID Search module — unified chemical identifier resolver.

Provides the :class:`Search` class for resolving chemical identifiers across multiple
offline databases (ChEBI, CompTox, PubChemID, ZeroPM, ChEMBL) with structure-aware
matching, confidence scoring, fuzzy name search, Tanimoto similarity search,
InChIKey-skeleton matching, and salt/solvent stripping.

Supported identifier types:

- ``"cas"``     — CAS Registry Number
- ``"name"``    — Chemical name (common or IUPAC)
- ``"smiles"``  — SMILES string
- ``"inchi"``   — InChI string
- ``"inchikey"``— InChIKey
- ``"dtxsid"``  — CompTox DTXSID
- ``"formula"`` — Molecular formula

Example usage::

    from provesid import Search

    # Resolve a list of CAS numbers
    s = Search("cas")
    df = s.search(["50-00-0", "64-17-5"])

    # Fuzzy name search (handles typos)
    s_name = Search("name", fuzzy=True)
    df = s_name.search(["asprin", "caffiene"])

    # SMILES with salt stripping and structure similarity
    s_smiles = Search("smiles", strip_salts=True, similarity_threshold=0.8)
    df = s_smiles.search("CC(=O)Oc1ccccc1C(=O)O")

    # InChIKey skeleton matching (same connectivity, any stereochemistry)
    s_ik = Search("inchikey", inchikey_skeleton=True)
    df = s_ik.search("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import pandas as pd
from tqdm import tqdm

from .chebi import ChebiSDF
from .chembl import CheMBL
from .comptox import CompToxID
from .pubchem import PubChemID
from .zeropm import ZeroPM
from .tools import (
    _apply_candidate_to_result,
    _candidate_compatible_with_consensus,
    _candidate_from_chebi_row,
    _candidate_from_chembl_row,
    _candidate_from_comptox_row,
    _candidate_from_pubchem_row,
    _candidate_from_zeropm_name_table,
    _candidate_from_zeropm_smiles,
    _compute_consensus,
    _extract_cas_values,
    _first_cas,
    _inchi_to_smiles,
    _inchikey_from_smiles,
    _is_missing,
    _make_candidate,
    _normalize_synonyms,
    _pick_first,
    _smiles_to_canonical_and_mass,
    _text_similarity,
    _to_float,
)

# ── Optional RDKit ─────────────────────────────────────────────────────────────
try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem, Descriptors
    from rdkit.Chem.SaltRemover import SaltRemover as _SaltRemover

    RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover
    Chem = None  # type: ignore[assignment]
    DataStructs = None  # type: ignore[assignment]
    AllChem = None  # type: ignore[assignment]
    Descriptors = None  # type: ignore[assignment]
    _SaltRemover = None  # type: ignore[assignment]
    RDKIT_AVAILABLE = False

# ── Optional rapidfuzz ─────────────────────────────────────────────────────────
try:
    from rapidfuzz import fuzz as _fuzz
    from rapidfuzz import process as _rfprocess

    RAPIDFUZZ_AVAILABLE = True
except ImportError:  # pragma: no cover
    _fuzz = None  # type: ignore[assignment]
    _rfprocess = None  # type: ignore[assignment]
    RAPIDFUZZ_AVAILABLE = False

# ── Patterns & constants ───────────────────────────────────────────────────────
_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_INCHI_PREFIX = "InChI="

# Base confidence scores per match method.
_BASE_CONFIDENCE: Dict[str, float] = {
    "exact_inchikey": 1.0,
    "exact_smiles": 0.95,
    "inchi": 0.95,
    "exact_cas": 0.90,
    "dtxsid": 0.90,
    "exact_name": 0.80,
    "inchikey_skeleton": 0.75,
    "tanimoto": 0.0,   # filled dynamically from Tanimoto score
    "fuzzy_name": 0.0,  # filled dynamically from rapidfuzz score
    "formula": 0.30,
    "unknown": 0.50,
}

# Canonical column order for the output DataFrame.
OUTPUT_COLUMNS: List[str] = [
    "query",
    "CASRN",
    "name",
    "IUPAC_name",
    "molecular_formula",
    "SMILES",
    "canonical_smiles",
    "kekulized_smiles",
    "InChI",
    "InChIKey",
    "DTXSID",
    "molecular_mass",
    "Synonyms",
    "parent_smiles",
    "parent_inchikey",
    "foundby",
    "source",
    "source_details",
    "confidence",
    "match_method",
    "match_score",
    "consensus_source",
    "source_match_scores",
]

# Name-normalization: prefixes to strip before fuzzy matching.
_NAME_PREFIXES = re.compile(
    r"^(?:"
    r"\(\u00b1\)-|"   # (±)-
    r"\(\+\)-|\(-\)-|"
    r"rac-|dl-|d-|l-|"
    r"\(r\)-|\(s\)-|\(rs\)-|"
    r"\(e\)-|\(z\)-"
    r")",
    re.IGNORECASE,
)

# Small built-in abbreviation map used by _normalize_name.
_ABBREVIATIONS: Dict[str, str] = {
    "mek": "methyl ethyl ketone",
    "mibk": "methyl isobutyl ketone",
    "dmf": "dimethylformamide",
    "dmso": "dimethyl sulfoxide",
    "thf": "tetrahydrofuran",
    "dcm": "dichloromethane",
    "etoh": "ethanol",
    "meoh": "methanol",
    "acn": "acetonitrile",
    "egme": "ethylene glycol monomethyl ether",
}

log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Module-level structure utility (used inside & outside the class)
# ─────────────────────────────────────────────────────────────────────────────

def normalize_structure(smiles: Optional[str]) -> Dict[str, Any]:
    """Convert a SMILES string into a normalized structure record.

    Runs a single RDKit parse and derives canonical SMILES, Kekulized SMILES,
    InChI, InChIKey, and molecular weight from it.  All fields are ``None``
    when RDKit is unavailable or the SMILES is invalid.

    Args:
        smiles: Input SMILES string.

    Returns:
        Dictionary with keys:
        ``canonical_smiles``, ``kekulized_smiles``, ``inchi``, ``inchikey``,
        ``mol_weight``, and ``mol`` (the RDKit Mol object; not serialized).

    Example::

        rec = normalize_structure("c1ccccc1")
        rec["canonical_smiles"]  # "c1ccccc1"
        rec["kekulized_smiles"]  # "C1=CC=CC=C1"
    """
    empty: Dict[str, Any] = {
        "canonical_smiles": None,
        "kekulized_smiles": None,
        "inchi": None,
        "inchikey": None,
        "mol_weight": None,
        "mol": None,
    }
    if _is_missing(smiles) or not RDKIT_AVAILABLE or Chem is None:
        return empty

    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return empty

        canonical = Chem.MolToSmiles(mol, canonical=True)

        # Kekulize on a copy so the original mol is unmodified
        try:
            mol_kek = Chem.RWMol(mol)
            Chem.Kekulize(mol_kek, clearAromaticFlags=False)
            kekulized = Chem.MolToSmiles(mol_kek, kekuleSmiles=True)
        except Exception:
            kekulized = None

        try:
            inchi = Chem.MolToInchi(mol)
            inchikey = Chem.InchiToInchiKey(inchi) if inchi else None
        except Exception:
            inchi = None
            inchikey = None

        mol_weight = float(Descriptors.MolWt(mol)) if Descriptors is not None else None

        return {
            "canonical_smiles": canonical,
            "kekulized_smiles": kekulized,
            "inchi": inchi,
            "inchikey": inchikey,
            "mol_weight": mol_weight,
            "mol": mol,
        }
    except Exception as exc:
        log.warning("normalize_structure failed for SMILES %r: %s", smiles, exc)
        return empty


def strip_salts(
    smiles: Optional[str],
    extra_smarts: Optional[List[str]] = None,
) -> Optional[str]:
    """Remove salt/solvent fragments from a SMILES and return the parent SMILES.

    Uses RDKit's ``SaltRemover`` with its default salt definitions, then picks
    the largest fragment by heavy-atom count when multiple fragments remain.

    Args:
        smiles: Input SMILES (may contain ``.``-separated fragments).
        extra_smarts: Optional list of additional SMARTS patterns to strip.

    Returns:
        SMILES of the parent (desalted) molecule, or ``None`` when RDKit is
        unavailable or the input is invalid.  Returns the original SMILES
        unchanged when no fragments are removed.

    Example::

        strip_salts("[Na+].[Cl-].CC(=O)O")  # "CC(=O)O"
    """
    if _is_missing(smiles) or not RDKIT_AVAILABLE or Chem is None or _SaltRemover is None:
        return smiles  # type: ignore[return-value]

    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return None

        # Build remover with optional extra patterns
        if extra_smarts:
            smarts_block = "\n".join(f"[{s}]" if not s.startswith("[") else s for s in extra_smarts)
            remover = _SaltRemover(defnData=smarts_block)
        else:
            remover = _SaltRemover()

        stripped = remover.StripMol(mol)
        if stripped is None:
            stripped = mol

        # Pick the largest fragment if still multi-component
        frags = Chem.rdmolops.GetMolFrags(stripped, asMols=True)
        if not frags:
            # SaltRemover stripped everything (all fragments are known salts).
            # Fall back to the largest fragment of the original molecule.
            frags = Chem.rdmolops.GetMolFrags(mol, asMols=True)
        if len(frags) > 1:
            stripped = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        elif len(frags) == 1:
            stripped = frags[0]

        result = Chem.MolToSmiles(stripped, canonical=True)
        return result if result else None
    except Exception as exc:
        log.warning("strip_salts failed for SMILES %r: %s", smiles, exc)
        return smiles  # type: ignore[return-value]


# ─────────────────────────────────────────────────────────────────────────────
# Search class
# ─────────────────────────────────────────────────────────────────────────────

class Search:
    """Unified chemical identifier resolver using offline databases.

    Accepts any single identifier type — CAS, name, SMILES, InChI, InChIKey,
    DTXSID, or molecular formula — and queries ChEBI, CompTox, PubChemID,
    ZeroPM, and ChEMBL to build a harmonised result.

    Features:

    - **Structure-aware matching**: canonicalisation and kekulisation via RDKit;
      InChIKey always derived and reported.
    - **Confidence scoring**: each result carries a ``confidence`` score in
      [0, 1] based on the match method and cross-source consensus.
    - **Fuzzy name matching**: rapidfuzz WRatio scorer with configurable
      cut-off (enabled with ``fuzzy=True``).
    - **Tanimoto similarity search**: Morgan-fingerprint-based fallback when
      ``similarity_threshold > 0``.
    - **InChIKey skeleton matching**: 14-character connectivity-layer prefix
      search (enabled with ``inchikey_skeleton=True``).
    - **Salt/solvent stripping**: RDKit SaltRemover + largest-fragment picker
      (enabled with ``strip_salts=True``); ``parent_smiles`` and
      ``parent_inchikey`` populated in results.
    - **Traceability**: ``source_details`` field records which sources were
      queried, whether they matched, and which output fields they contributed.

    Attributes:
        identifier_type (str): Input identifier type used for all queries.
        strip_salts (bool): Strip salts/solvents and report parent molecule.
        fuzzy (bool): Enable fuzzy name matching via rapidfuzz.
        similarity_threshold (float): Minimum Tanimoto similarity for
            structure-based fallback search (0.0 disables it).
        inchikey_skeleton (bool): Enable InChIKey 14-char skeleton matching.
        show_progress (bool): Display tqdm progress bar during batch queries.
        salt_smarts (list[str]): Additional SMARTS patterns to remove during
            salt stripping.

    Example::

        from provesid import Search

        s = Search("cas")
        df = s.search(["50-00-0", "64-17-5"])
        print(df[["CASRN", "name", "canonical_smiles", "confidence"]])

        s_fuzzy = Search("name", fuzzy=True)
        df = s_fuzzy.search(["asprin", "paracetamol"])
    """

    SUPPORTED_TYPES: frozenset = frozenset(
        ["cas", "name", "smiles", "inchi", "inchikey", "dtxsid", "formula"]
    )

    _SOURCE_KEYS: List[str] = ["chebi", "comptox", "pubchem", "zeropm", "chembl"]
    _SOURCE_DISPLAY: Dict[str, str] = {
        "chebi": "ChEBI",
        "comptox": "CompTox",
        "pubchem": "PubChemID",
        "zeropm": "ZeroPM",
        "chembl": "ChEMBL",
    }

    def __init__(
        self,
        identifier_type: str = "cas",
        *,
        strip_salts: bool = False,
        fuzzy: bool = False,
        similarity_threshold: float = 0.0,
        inchikey_skeleton: bool = False,
        show_progress: bool = True,
        salt_smarts: Optional[List[str]] = None,
        chebi: Optional[ChebiSDF] = None,
        comptox: Optional[CompToxID] = None,
        pubchem: Optional[PubChemID] = None,
        zeropm: Optional[ZeroPM] = None,
        chembl: Optional[CheMBL] = None,
    ) -> None:
        """Initialise a Search resolver.

        Args:
            identifier_type: Type of identifier to resolve.  One of ``"cas"``,
                ``"name"``, ``"smiles"``, ``"inchi"``, ``"inchikey"``,
                ``"dtxsid"``, ``"formula"``.  Defaults to ``"cas"``.
            strip_salts: Strip salt/solvent fragments and populate
                ``parent_smiles`` / ``parent_inchikey`` columns.
            fuzzy: Enable fuzzy name matching when an exact name match fails.
                Requires rapidfuzz.
            similarity_threshold: Tanimoto similarity threshold in [0, 1].
                When > 0 a Morgan-fingerprint similarity search is run as a
                fallback for SMILES queries with no exact match.  0.0 disables
                the search entirely.
            inchikey_skeleton: When True, fall back to 14-character InChIKey
                prefix matching when an exact InChIKey match fails.
            show_progress: Display a tqdm progress bar during batch queries.
            salt_smarts: Additional SMARTS patterns passed to
                :func:`strip_salts` when ``strip_salts=True``.
            chebi: Pre-initialised :class:`~provesid.ChebiSDF` client.  When
                ``None`` the client is created lazily on first use.
            comptox: Pre-initialised :class:`~provesid.CompToxID` client.
            pubchem: Pre-initialised :class:`~provesid.PubChemID` client.
            zeropm: Pre-initialised :class:`~provesid.ZeroPM` client.
            chembl: Pre-initialised :class:`~provesid.CheMBL` client.

        Raises:
            ValueError: If ``identifier_type`` is not one of the supported
                values.
        """
        if identifier_type not in self.SUPPORTED_TYPES:
            raise ValueError(
                f"identifier_type must be one of {sorted(self.SUPPORTED_TYPES)}, "
                f"got {identifier_type!r}"
            )

        self.identifier_type = identifier_type
        self.strip_salts = strip_salts
        self.fuzzy = fuzzy
        self.similarity_threshold = float(similarity_threshold)
        self.inchikey_skeleton = inchikey_skeleton
        self.show_progress = show_progress
        self.salt_smarts: List[str] = list(salt_smarts or [])

        # Client references — may be None until _ensure_clients() is called.
        self._chebi = chebi
        self._comptox = comptox
        self._pubchem = pubchem
        self._zeropm = zeropm
        self._chembl = chembl

        # Track whether automatic client init has been attempted.
        self._clients_initialized: bool = any(
            c is not None for c in [chebi, comptox, pubchem, zeropm, chembl]
        )

    # ── Client lifecycle ──────────────────────────────────────────────────────

    def _ensure_clients(self) -> None:
        """Lazily initialise all offline source clients.

        This method is idempotent — it only runs once per Search instance.
        Individual clients that fail to initialise are set to ``None`` and a
        warning is logged; the search continues with the remaining sources.
        """
        if self._clients_initialized:
            return

        for attr, factory in [
            ("_chebi", ChebiSDF),
            ("_comptox", CompToxID),
            ("_pubchem", PubChemID),
            ("_zeropm", ZeroPM),
            ("_chembl", CheMBL),
        ]:
            if getattr(self, attr) is None:
                try:
                    setattr(self, attr, factory())
                except Exception as exc:
                    log.warning("Could not initialise offline source %s: %s", attr[1:], exc)

        self._clients_initialized = True

    # ── Public entry point ────────────────────────────────────────────────────

    def search(
        self,
        queries: Union[str, List[str], pd.DataFrame, Path],
        *,
        column: Optional[str] = None,
    ) -> pd.DataFrame:
        """Resolve one or more chemical identifiers and return a DataFrame.

        Args:
            queries: Input identifiers in any of the following forms:

                - A single string — returns a one-row DataFrame.
                - A list of strings — one row per query.
                - A :class:`pandas.DataFrame` — the column given by ``column``
                  is used as the query list.  All other columns are preserved
                  in the output (left-joined on the original index).
                - A file path (:class:`pathlib.Path` or string ending in
                  ``.csv`` / ``.parquet``) — read into a DataFrame first;
                  ``column`` must be provided.

            column: Column name to read from a DataFrame or file input.
                Required when ``queries`` is a DataFrame or file path.

        Returns:
            DataFrame with columns defined in :data:`OUTPUT_COLUMNS`.  Rows
            correspond to input queries in the same order.

        Raises:
            ValueError: If a DataFrame/file input is given but ``column`` is
                not specified.
            FileNotFoundError: If the given file path does not exist.

        Example::

            s = Search("cas")
            df = s.search(["50-00-0", "64-17-5"])
            df = s.search(Path("compounds.csv"), column="CAS")
        """
        self._ensure_clients()

        query_list, extra_df = self._coerce_queries(queries, column)

        iterator = (
            tqdm(query_list, desc=f"Resolving {self.identifier_type.upper()}")
            if self.show_progress
            else query_list
        )

        rows = [self._resolve_single(q) for q in iterator]

        result_df = pd.DataFrame(rows)
        # Ensure all output columns are present (fill missing with None)
        for col in OUTPUT_COLUMNS:
            if col not in result_df.columns:
                result_df[col] = None
        result_df = result_df[OUTPUT_COLUMNS]

        # Merge extra columns from the original DataFrame if provided
        if extra_df is not None:
            extra_cols = [c for c in extra_df.columns if c not in result_df.columns]
            if extra_cols:
                result_df = pd.concat(
                    [result_df, extra_df[extra_cols].reset_index(drop=True)],
                    axis=1,
                )

        return result_df

    # ── Input normalisation ───────────────────────────────────────────────────

    def _coerce_queries(
        self,
        queries: Union[str, List[str], pd.DataFrame, Path],
        column: Optional[str],
    ) -> Tuple[List[str], Optional[pd.DataFrame]]:
        """Convert the ``queries`` argument to a plain list of strings.

        Args:
            queries: Raw input from :meth:`search`.
            column: Column name for DataFrame/file inputs.

        Returns:
            Tuple of (query_list, optional extra DataFrame for merge).

        Raises:
            ValueError: If a DataFrame/file is given without a column name.
        """
        # File path
        if isinstance(queries, (str, Path)):
            p = Path(queries)
            if p.exists() and p.suffix in {".csv", ".parquet"}:
                if column is None:
                    raise ValueError(
                        "Provide column= when passing a file path as queries."
                    )
                if p.suffix == ".parquet":
                    df = pd.read_parquet(p)
                else:
                    df = pd.read_csv(p)
                return df[column].astype(str).tolist(), df

            # Treat as a bare string query
            return [str(queries)], None

        # DataFrame
        if isinstance(queries, pd.DataFrame):
            if column is None:
                raise ValueError(
                    "Provide column= when passing a DataFrame as queries."
                )
            return queries[column].astype(str).tolist(), queries

        # List of strings
        if isinstance(queries, list):
            return [str(q) for q in queries], None

        return [str(queries)], None

    # ── Single-query dispatcher ───────────────────────────────────────────────

    def _resolve_single(self, query: str) -> Dict[str, Any]:
        """Dispatch one query to the appropriate resolver and return a result dict.

        Args:
            query: A single identifier string.

        Returns:
            Fully populated result dict matching :data:`OUTPUT_COLUMNS`.
        """
        dispatch = {
            "cas": self._resolve_cas,
            "name": self._resolve_name,
            "smiles": self._resolve_smiles,
            "inchi": self._resolve_inchi,
            "inchikey": self._resolve_inchikey,
            "dtxsid": self._resolve_dtxsid,
            "formula": self._resolve_formula,
        }
        return dispatch[self.identifier_type](query)

    # ── Empty result template ─────────────────────────────────────────────────

    def _empty_result(self, query: str, foundby: str) -> Dict[str, Any]:
        """Return a result dict with all fields initialised to None/defaults.

        Args:
            query: The original query string.
            foundby: The identifier type used for the search.

        Returns:
            Dict with all :data:`OUTPUT_COLUMNS` keys present.
        """
        return {
            "query": query,
            "CASRN": None,
            "name": None,
            "IUPAC_name": None,
            "molecular_formula": None,
            "SMILES": None,
            "canonical_smiles": None,
            "kekulized_smiles": None,
            "InChI": None,
            "InChIKey": None,
            "DTXSID": None,
            "molecular_mass": None,
            "Synonyms": None,
            "parent_smiles": None,
            "parent_inchikey": None,
            "foundby": foundby,
            "source": None,
            "source_details": {},
            "confidence": 0.0,
            "match_method": "unknown",
            "match_score": 0.0,
            "consensus_source": None,
            "source_match_scores": {},
        }

    # ── CAS resolver ─────────────────────────────────────────────────────────

    def _resolve_cas(self, cas: str) -> Dict[str, Any]:
        """Resolve a CAS Registry Number into a unified identifier record.

        Queries ChEBI → CompTox → PubChemID → ZeroPM with waterfall priority,
        then enriches via ChEMBL.

        Args:
            cas: CAS Registry Number string.

        Returns:
            Fully populated result dict.
        """
        result = self._empty_result(cas, "CASRN")
        result["CASRN"] = cas

        candidates: Dict[str, Optional[Dict[str, Any]]] = {k: None for k in self._SOURCE_KEYS}

        # ChEBI
        if self._chebi is not None:
            try:
                rows = self._chebi.search_by_cas(cas)
                if rows:
                    candidates["chebi"] = _candidate_from_chebi_row(rows[0])
            except Exception as exc:
                log.warning("ChEBI CAS lookup failed for %r: %s", cas, exc)

        # CompTox
        if self._comptox is not None:
            try:
                row = self._comptox.get_by_casrn(cas)
                if row:
                    candidates["comptox"] = _candidate_from_comptox_row(row)
            except Exception as exc:
                log.warning("CompTox CAS lookup failed for %r: %s", cas, exc)

        # PubChemID
        if self._pubchem is not None:
            try:
                row = self._pubchem.get_by_cas(cas)
                if row:
                    candidates["pubchem"] = _candidate_from_pubchem_row(row)
            except Exception as exc:
                log.warning("PubChemID CAS lookup failed for %r: %s", cas, exc)

        # ZeroPM
        if self._zeropm is not None:
            try:
                table = self._zeropm.get_id_table_from_cas(cas)
                candidates["zeropm"] = _candidate_from_zeropm_name_table(cas, table)
            except Exception as exc:
                log.warning("ZeroPM CAS lookup failed for %r: %s", cas, exc)

        # ChEMBL — enrichment via SMILES after primary sources
        smiles_so_far = _first_smiles_from_candidates(candidates)
        if self._chembl is not None and not _is_missing(smiles_so_far):
            try:
                row = self._chembl.search_by_smiles(str(smiles_so_far))
                if row:
                    candidates["chembl"] = _candidate_from_chembl_row(row, self._chembl)
            except Exception as exc:
                log.warning("ChEMBL CAS enrichment failed for %r: %s", cas, exc)

        source_details = self._build_source_details(candidates)
        return self._finalise_result(result, candidates, "exact_cas", source_details)

    # ── Name resolver ─────────────────────────────────────────────────────────

    def _resolve_name(self, name: str) -> Dict[str, Any]:
        """Resolve a chemical name into a unified identifier record.

        First attempts exact-name matching across all sources, then falls back
        to fuzzy matching when ``self.fuzzy`` is True and no exact match is
        found.

        Args:
            name: Chemical name string (common or IUPAC).

        Returns:
            Fully populated result dict.
        """
        result = self._empty_result(name, "name")
        candidates = self._candidates_from_name(name, exact=True)

        fuzzy_score: Optional[float] = None
        match_method = "exact_name"

        # Fuzzy fallback
        if not _any_candidate(candidates) and self.fuzzy:
            candidates, fuzzy_score = self._fuzzy_name_candidates(name)
            match_method = "fuzzy_name"

        source_details = self._build_source_details(candidates)
        result = self._finalise_result(result, candidates, match_method, source_details)
        if fuzzy_score is not None:
            result["confidence"] = self._compute_confidence(
                "fuzzy_name", result["match_score"], fuzzy_score=fuzzy_score
            )
        return result

    def _candidates_from_name(
        self, name: str, exact: bool = True
    ) -> Dict[str, Optional[Dict[str, Any]]]:
        """Build candidates dict from a name query across all sources.

        Args:
            name: Chemical name to search.
            exact: Whether to use exact matching.

        Returns:
            Dict mapping source keys to candidate records.
        """
        candidates: Dict[str, Optional[Dict[str, Any]]] = {k: None for k in self._SOURCE_KEYS}

        if self._chebi is not None:
            try:
                rows = self._chebi.search_by_name(name, exact=exact)
                if not rows:
                    rows = self._chebi.search_by_synonym(name, exact=exact)
                if rows:
                    candidates["chebi"] = _candidate_from_chebi_row(rows[0])
            except Exception as exc:
                log.warning("ChEBI name lookup failed for %r: %s", name, exc)

        if self._comptox is not None:
            try:
                row = self._comptox.get_by_name(name)
                if row is None:
                    matches = self._comptox.search_by_name(name, exact=False, limit=5)
                    row = matches[0] if matches else None
                if row:
                    candidates["comptox"] = _candidate_from_comptox_row(row)
            except Exception as exc:
                log.warning("CompTox name lookup failed for %r: %s", name, exc)

        if self._pubchem is not None:
            try:
                rows = self._pubchem.search_by_name(name, exact=exact, limit=5)
                if rows:
                    candidates["pubchem"] = _candidate_from_pubchem_row(rows[0])
            except Exception as exc:
                log.warning("PubChemID name lookup failed for %r: %s", name, exc)

        if self._zeropm is not None:
            try:
                table = self._zeropm.get_id_table_from_name(name)
                candidates["zeropm"] = _candidate_from_zeropm_name_table(name, table)
            except Exception as exc:
                log.warning("ZeroPM name lookup failed for %r: %s", name, exc)

        if self._chembl is not None:
            try:
                rows = self._chembl.search_by_name(name, limit=5)
                if rows:
                    candidates["chembl"] = _candidate_from_chembl_row(rows[0], self._chembl)
            except Exception as exc:
                log.warning("ChEMBL name lookup failed for %r: %s", name, exc)

        return candidates

    # ── SMILES resolver ───────────────────────────────────────────────────────

    def _resolve_smiles(self, smiles: str) -> Dict[str, Any]:
        """Resolve a SMILES string into a unified identifier record.

        Canonicalises the input, derives an InChIKey, and queries sources by
        InChIKey (ChEBI) or canonical SMILES.  Falls back to Tanimoto
        similarity search when ``self.similarity_threshold > 0`` and no exact
        match is found.

        Args:
            smiles: SMILES string.

        Returns:
            Fully populated result dict.
        """
        result = self._empty_result(smiles, "SMILES")
        result["SMILES"] = smiles

        norm = normalize_structure(smiles)
        canonical = norm["canonical_smiles"] or smiles
        inchikey = norm["inchikey"] or _inchikey_from_smiles(smiles)

        candidates: Dict[str, Optional[Dict[str, Any]]] = {k: None for k in self._SOURCE_KEYS}
        match_method = "exact_smiles"

        # ChEBI — lookup by InChIKey
        if self._chebi is not None and not _is_missing(inchikey):
            try:
                row = self._chebi.search_by_inchikey(str(inchikey))
                if row:
                    candidates["chebi"] = _candidate_from_chebi_row(row)
            except Exception as exc:
                log.warning("ChEBI SMILES lookup failed for %r: %s", smiles, exc)

        # CompTox
        if self._comptox is not None:
            try:
                row = self._comptox.get_by_smiles(smiles)
                if row is None and not _is_missing(canonical):
                    row = self._comptox.get_by_smiles(canonical)
                if row:
                    candidates["comptox"] = _candidate_from_comptox_row(row)
            except Exception as exc:
                log.warning("CompTox SMILES lookup failed for %r: %s", smiles, exc)

        # PubChemID
        if self._pubchem is not None:
            try:
                row = self._pubchem.get_by_smiles(smiles)
                if row is None and not _is_missing(canonical):
                    row = self._pubchem.get_by_smiles(canonical)
                if row:
                    candidates["pubchem"] = _candidate_from_pubchem_row(row)
            except Exception as exc:
                log.warning("PubChemID SMILES lookup failed for %r: %s", smiles, exc)

        # ZeroPM
        if self._zeropm is not None:
            try:
                candidates["zeropm"] = _candidate_from_zeropm_smiles(smiles, self._zeropm)
            except Exception as exc:
                log.warning("ZeroPM SMILES lookup failed for %r: %s", smiles, exc)

        # ChEMBL
        if self._chembl is not None:
            try:
                row = self._chembl.search_by_smiles(smiles)
                if row:
                    candidates["chembl"] = _candidate_from_chembl_row(row, self._chembl)
            except Exception as exc:
                log.warning("ChEMBL SMILES lookup failed for %r: %s", smiles, exc)

        # Tanimoto similarity fallback
        if not _any_candidate(candidates) and self.similarity_threshold > 0:
            sim_candidates, tanimoto_score = self._tanimoto_candidates(smiles)
            if _any_candidate(sim_candidates):
                candidates = sim_candidates
                match_method = "tanimoto"
                source_details = self._build_source_details(candidates)
                result = self._finalise_result(result, candidates, match_method, source_details)
                result["confidence"] = self._compute_confidence(
                    "tanimoto", result["match_score"], tanimoto=tanimoto_score
                )
                return result

        source_details = self._build_source_details(candidates)
        return self._finalise_result(result, candidates, match_method, source_details)

    # ── InChI resolver ────────────────────────────────────────────────────────

    def _resolve_inchi(self, inchi: str) -> Dict[str, Any]:
        """Resolve an InChI string into a unified identifier record.

        Converts the InChI to InChIKey via RDKit and delegates to
        :meth:`_resolve_inchikey`.  Also queries sources that store InChI
        directly (ChEBI, CompTox, PubChemID).

        Args:
            inchi: InChI string (must start with ``"InChI="``).

        Returns:
            Fully populated result dict.
        """
        result = self._empty_result(inchi, "InChI")
        result["InChI"] = inchi

        # Derive InChIKey and SMILES via RDKit
        inchikey: Optional[str] = None
        smiles: Optional[str] = None
        if RDKIT_AVAILABLE and Chem is not None and inchi.startswith(_INCHI_PREFIX):
            try:
                mol = Chem.MolFromInchi(str(inchi))
                if mol is not None:
                    inchikey = Chem.InchiToInchiKey(inchi)
                    smiles = Chem.MolToSmiles(mol)
            except Exception:
                pass

        # Pre-populate InChIKey so _finalise_result can use it even without a source match
        if not _is_missing(inchikey):
            result["InChIKey"] = inchikey
        if not _is_missing(smiles):
            result["SMILES"] = smiles

        candidates: Dict[str, Optional[Dict[str, Any]]] = {k: None for k in self._SOURCE_KEYS}

        # ChEBI
        if self._chebi is not None:
            try:
                row = self._chebi.search_by_inchi(inchi)
                if row:
                    candidates["chebi"] = _candidate_from_chebi_row(row)
            except Exception as exc:
                log.warning("ChEBI InChI lookup failed for %r: %s", inchi[:40], exc)

        # CompTox — lookup by InChIKey if derived
        if self._comptox is not None and not _is_missing(inchikey):
            try:
                row = self._comptox.get_by_inchikey(str(inchikey))
                if row:
                    candidates["comptox"] = _candidate_from_comptox_row(row)
            except Exception as exc:
                log.warning("CompTox InChI lookup failed: %s", exc)

        # PubChemID
        if self._pubchem is not None:
            try:
                row = self._pubchem.get_by_inchi(inchi)
                if row:
                    candidates["pubchem"] = _candidate_from_pubchem_row(row)
            except Exception as exc:
                log.warning("PubChemID InChI lookup failed: %s", exc)

        # ZeroPM
        if self._zeropm is not None:
            try:
                table = self._zeropm.get_id_table_from_inchi(inchi)
                candidates["zeropm"] = _candidate_from_zeropm_name_table(inchi, table)
            except Exception as exc:
                log.warning("ZeroPM InChI lookup failed: %s", exc)

        # ChEMBL — via SMILES
        if self._chembl is not None and not _is_missing(smiles):
            try:
                row = self._chembl.search_by_smiles(str(smiles))
                if row:
                    candidates["chembl"] = _candidate_from_chembl_row(row, self._chembl)
            except Exception as exc:
                log.warning("ChEMBL InChI lookup failed: %s", exc)

        source_details = self._build_source_details(candidates)
        return self._finalise_result(result, candidates, "inchi", source_details)

    # ── InChIKey resolver ─────────────────────────────────────────────────────

    def _resolve_inchikey(self, inchikey: str) -> Dict[str, Any]:
        """Resolve an InChIKey into a unified identifier record.

        Queries all offline sources by InChIKey.  Falls back to 14-character
        skeleton matching when ``self.inchikey_skeleton`` is True and no exact
        match is found.

        Args:
            inchikey: Full 27-character InChIKey
                (``XXXXXXXXXXXXXX-XXXXXXXXXX-X``).

        Returns:
            Fully populated result dict.
        """
        result = self._empty_result(inchikey, "InChIKey")
        result["InChIKey"] = inchikey

        candidates: Dict[str, Optional[Dict[str, Any]]] = {k: None for k in self._SOURCE_KEYS}
        match_method = "exact_inchikey"

        # ChEBI
        if self._chebi is not None:
            try:
                row = self._chebi.search_by_inchikey(inchikey)
                if row:
                    candidates["chebi"] = _candidate_from_chebi_row(row)
            except Exception as exc:
                log.warning("ChEBI InChIKey lookup failed for %r: %s", inchikey, exc)

        # CompTox
        if self._comptox is not None:
            try:
                row = self._comptox.get_by_inchikey(inchikey)
                if row:
                    candidates["comptox"] = _candidate_from_comptox_row(row)
            except Exception as exc:
                log.warning("CompTox InChIKey lookup failed for %r: %s", inchikey, exc)

        # PubChemID
        if self._pubchem is not None:
            try:
                row = self._pubchem.get_by_inchikey(inchikey)
                if row:
                    candidates["pubchem"] = _candidate_from_pubchem_row(row)
            except Exception as exc:
                log.warning("PubChemID InChIKey lookup failed for %r: %s", inchikey, exc)

        # ZeroPM
        if self._zeropm is not None:
            try:
                table = self._zeropm.get_id_table_from_inchikey(inchikey)
                candidates["zeropm"] = _candidate_from_zeropm_name_table(inchikey, table)
            except Exception as exc:
                log.warning("ZeroPM InChIKey lookup failed for %r: %s", inchikey, exc)

        # ChEMBL
        if self._chembl is not None:
            try:
                row = self._chembl.search_by_inchikey(inchikey)
                if row:
                    candidates["chembl"] = _candidate_from_chembl_row(row, self._chembl)
            except Exception as exc:
                log.warning("ChEMBL InChIKey lookup failed for %r: %s", inchikey, exc)

        # InChIKey skeleton fallback
        if not _any_candidate(candidates) and self.inchikey_skeleton:
            skel_candidates, skeleton = self._skeleton_candidates(inchikey)
            if _any_candidate(skel_candidates):
                candidates = skel_candidates
                match_method = "inchikey_skeleton"

        source_details = self._build_source_details(candidates)
        return self._finalise_result(result, candidates, match_method, source_details)

    # ── DTXSID resolver ───────────────────────────────────────────────────────

    def _resolve_dtxsid(self, dtxsid: str) -> Dict[str, Any]:
        """Resolve a CompTox DTXSID into a unified identifier record.

        Queries CompTox as the primary source, then cross-references other
        sources using the InChIKey derived from the CompTox result.

        Args:
            dtxsid: CompTox DTXSID string (e.g., ``"DTXSID7020182"``).

        Returns:
            Fully populated result dict.
        """
        result = self._empty_result(dtxsid, "DTXSID")
        result["DTXSID"] = dtxsid

        candidates: Dict[str, Optional[Dict[str, Any]]] = {k: None for k in self._SOURCE_KEYS}

        # CompTox primary
        if self._comptox is not None:
            try:
                row = self._comptox.get_by_dtxsid(dtxsid)
                if row:
                    candidates["comptox"] = _candidate_from_comptox_row(row)
            except Exception as exc:
                log.warning("CompTox DTXSID lookup failed for %r: %s", dtxsid, exc)

        # Cross-reference other sources by InChIKey
        comptox_cand = candidates.get("comptox")
        inchikey = comptox_cand.get("InChIKey") if comptox_cand else None

        if not _is_missing(inchikey):
            # ChEBI
            if self._chebi is not None:
                try:
                    row = self._chebi.search_by_inchikey(str(inchikey))
                    if row:
                        candidates["chebi"] = _candidate_from_chebi_row(row)
                except Exception as exc:
                    log.warning("ChEBI DTXSID cross-ref failed: %s", exc)

            # PubChemID
            if self._pubchem is not None:
                try:
                    row = self._pubchem.get_by_inchikey(str(inchikey))
                    if row:
                        candidates["pubchem"] = _candidate_from_pubchem_row(row)
                except Exception as exc:
                    log.warning("PubChemID DTXSID cross-ref failed: %s", exc)

            # ZeroPM
            if self._zeropm is not None:
                try:
                    table = self._zeropm.get_id_table_from_inchikey(str(inchikey))
                    candidates["zeropm"] = _candidate_from_zeropm_name_table(dtxsid, table)
                except Exception as exc:
                    log.warning("ZeroPM DTXSID cross-ref failed: %s", exc)

            # ChEMBL
            if self._chembl is not None:
                try:
                    row = self._chembl.search_by_inchikey(str(inchikey))
                    if row:
                        candidates["chembl"] = _candidate_from_chembl_row(row, self._chembl)
                except Exception as exc:
                    log.warning("ChEMBL DTXSID cross-ref failed: %s", exc)

        source_details = self._build_source_details(candidates)
        return self._finalise_result(result, candidates, "dtxsid", source_details)

    # ── Formula resolver ──────────────────────────────────────────────────────

    def _resolve_formula(self, formula: str) -> Dict[str, Any]:
        """Resolve a molecular formula into a unified identifier record.

        Returns the most complete match across CompTox, PubChemID, and ChEBI.
        Formulas are not unique identifiers, so confidence is capped at 0.30.

        Args:
            formula: Molecular formula string (e.g., ``"C9H8O4"``).

        Returns:
            Fully populated result dict for the best (most complete) match.
        """
        result = self._empty_result(formula, "formula")
        result["molecular_formula"] = formula

        candidates: Dict[str, Optional[Dict[str, Any]]] = {k: None for k in self._SOURCE_KEYS}

        if self._chebi is not None:
            try:
                rows = self._chebi.search_by_formula(formula)
                if rows:
                    best = _most_complete_row(rows)
                    candidates["chebi"] = _candidate_from_chebi_row(best)
            except Exception as exc:
                log.warning("ChEBI formula lookup failed for %r: %s", formula, exc)

        if self._comptox is not None:
            try:
                rows = self._comptox.search_by_formula(formula)
                if rows:
                    best = _most_complete_row(rows)
                    candidates["comptox"] = _candidate_from_comptox_row(best)
            except Exception as exc:
                log.warning("CompTox formula lookup failed for %r: %s", formula, exc)

        if self._pubchem is not None:
            try:
                rows = self._pubchem.search_by_formula(formula)
                if rows:
                    best = _most_complete_row(rows)
                    candidates["pubchem"] = _candidate_from_pubchem_row(best)
            except Exception as exc:
                log.warning("PubChemID formula lookup failed for %r: %s", formula, exc)

        source_details = self._build_source_details(candidates)
        return self._finalise_result(result, candidates, "formula", source_details)

    # ── Fuzzy name search ─────────────────────────────────────────────────────

    def _fuzzy_name_candidates(
        self, name: str
    ) -> Tuple[Dict[str, Optional[Dict[str, Any]]], Optional[float]]:
        """Search for a chemical name using fuzzy matching via rapidfuzz.

        Normalises the query name, queries each source for fuzzy name matches,
        and returns the best candidate per source plus the overall fuzzy score.

        Args:
            name: Chemical name to search (may contain typos or variations).

        Returns:
            Tuple of:
            - Dict mapping source keys to the best fuzzy-matched candidate.
            - Best fuzzy score in [0, 1], or ``None`` if rapidfuzz is
              unavailable.
        """
        candidates: Dict[str, Optional[Dict[str, Any]]] = {k: None for k in self._SOURCE_KEYS}

        if not RAPIDFUZZ_AVAILABLE or _rfprocess is None:
            log.warning("rapidfuzz not available; fuzzy name matching skipped.")
            return candidates, None

        norm_name = self._normalize_name(name)
        best_score: float = 0.0

        # ZeroPM has a built-in similar-name method
        if self._zeropm is not None:
            try:
                results = self._zeropm.query_similar_name(norm_name)
                if results is not None and not (
                    isinstance(results, pd.DataFrame) and results.empty
                ):
                    # query_similar_name may return a list or DataFrame
                    if isinstance(results, pd.DataFrame) and not results.empty:
                        table = results
                    else:
                        table = None
                    if table is not None:
                        cand = _candidate_from_zeropm_name_table(name, table)
                        if cand:
                            candidates["zeropm"] = cand
                            best_score = max(best_score, 0.7)
            except Exception as exc:
                log.warning("ZeroPM fuzzy name search failed for %r: %s", name, exc)

        # For other sources we use rapidfuzz directly against their search methods
        # (they accept fuzzy/partial inputs via exact=False)
        fuzzy_candidates = self._candidates_from_name(norm_name, exact=False)
        for key, cand in fuzzy_candidates.items():
            if cand is not None and candidates.get(key) is None:
                candidates[key] = cand

        # Compute best name similarity score across all found candidates
        for cand in candidates.values():
            if cand is None:
                continue
            sim = _text_similarity(name, cand.get("name"))
            best_score = max(best_score, sim)

        return candidates, best_score if best_score > 0 else None

    def _normalize_name(self, name: str) -> str:
        """Normalise a chemical name for fuzzy matching.

        Lowercases, strips whitespace, removes common stereochemistry prefixes,
        collapses multiple spaces, and expands known abbreviations.

        Args:
            name: Raw chemical name.

        Returns:
            Normalised name suitable for fuzzy comparison.

        Example::

            Search._normalize_name("D-Aspirin")  # "aspirin"
            Search._normalize_name("MEK")        # "methyl ethyl ketone"
        """
        n = name.strip().lower()
        n = _NAME_PREFIXES.sub("", n)
        n = re.sub(r"\s+", " ", n).strip()
        return _ABBREVIATIONS.get(n, n)

    # ── Tanimoto similarity search ────────────────────────────────────────────

    def _tanimoto_candidates(
        self, query_smiles: str
    ) -> Tuple[Dict[str, Optional[Dict[str, Any]]], Optional[float]]:
        """Find structurally similar compounds using Tanimoto similarity.

        Computes a Morgan fingerprint for ``query_smiles`` and queries each
        source with its similarity search capabilities.  Returns candidates
        that meet ``self.similarity_threshold``.

        Args:
            query_smiles: Query SMILES string.

        Returns:
            Tuple of:
            - Candidates dict (best match per source at or above threshold).
            - Best Tanimoto score observed, or ``None`` if RDKit is unavailable.

        Note:
            This is an initial implementation that uses per-source lookup; a
            future Parquet + vectorised fingerprint approach will be faster for
            large datasets.
        """
        candidates: Dict[str, Optional[Dict[str, Any]]] = {k: None for k in self._SOURCE_KEYS}

        if not RDKIT_AVAILABLE or Chem is None or DataStructs is None or AllChem is None:
            log.warning("RDKit not available; Tanimoto search skipped.")
            return candidates, None

        try:
            mol = Chem.MolFromSmiles(query_smiles)
            if mol is None:
                return candidates, None
            query_fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        except Exception as exc:
            log.warning("Could not compute fingerprint for %r: %s", query_smiles, exc)
            return candidates, None

        best_tanimoto: float = 0.0

        def _tanimoto_from_smiles(smiles: Optional[str]) -> float:
            if _is_missing(smiles) or Chem is None:
                return 0.0
            try:
                m = Chem.MolFromSmiles(str(smiles))
                if m is None:
                    return 0.0
                fp = AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=2048)
                return DataStructs.TanimotoSimilarity(query_fp, fp)
            except Exception:
                return 0.0

        # ChEMBL provides a native similarity search
        if self._chembl is not None:
            try:
                row = self._chembl.search_by_smiles(query_smiles)
                if row:
                    t = _tanimoto_from_smiles(row.get("canonical_smiles"))
                    if t >= self.similarity_threshold:
                        candidates["chembl"] = _candidate_from_chembl_row(row, self._chembl)
                        best_tanimoto = max(best_tanimoto, t)
            except Exception as exc:
                log.warning("ChEMBL Tanimoto search failed: %s", exc)

        # PubChemID — try canonical SMILES lookup as a proxy
        if self._pubchem is not None:
            try:
                norm = normalize_structure(query_smiles)
                if not _is_missing(norm["canonical_smiles"]):
                    row = self._pubchem.get_by_smiles(norm["canonical_smiles"])
                    if row:
                        t = _tanimoto_from_smiles(row.get("smiles") or row.get("canonical_smiles"))
                        if t >= self.similarity_threshold:
                            candidates["pubchem"] = _candidate_from_pubchem_row(row)
                            best_tanimoto = max(best_tanimoto, t)
            except Exception as exc:
                log.warning("PubChemID Tanimoto search failed: %s", exc)

        return candidates, best_tanimoto if best_tanimoto > 0 else None

    # ── InChIKey skeleton search ──────────────────────────────────────────────

    def _skeleton_candidates(
        self, inchikey: str
    ) -> Tuple[Dict[str, Optional[Dict[str, Any]]], str]:
        """Search by the 14-character InChIKey skeleton (connectivity layer).

        The first block of an InChIKey encodes the molecular skeleton.
        Matching on this prefix finds compounds with the same connectivity
        regardless of stereochemistry, isotopes, or charge.

        Args:
            inchikey: Full 27-character InChIKey.

        Returns:
            Tuple of:
            - Candidates dict populated from skeleton matches.
            - The 14-character skeleton prefix used.

        Example::

            candidates, skeleton = s._skeleton_candidates(
                "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
            )
            # skeleton == "BSYNRYMUTXBXSQ"
        """
        skeleton = inchikey[:14]
        candidates: Dict[str, Optional[Dict[str, Any]]] = {k: None for k in self._SOURCE_KEYS}

        # CompTox — SQL LIKE query on inchikey column
        if self._comptox is not None:
            try:
                row = self._comptox.get_by_inchikey(inchikey)
                if row is None:
                    # Partial-match: try all keys that start with skeleton
                    rows = _comptox_skeleton_search(self._comptox, skeleton)
                    row = rows[0] if rows else None
                if row:
                    candidates["comptox"] = _candidate_from_comptox_row(row)
            except Exception as exc:
                log.warning("CompTox skeleton search failed for %r: %s", skeleton, exc)

        # PubChemID
        if self._pubchem is not None:
            try:
                row = self._pubchem.get_by_inchikey(inchikey)
                if row is None:
                    rows = _pubchem_skeleton_search(self._pubchem, skeleton)
                    row = rows[0] if rows else None
                if row:
                    candidates["pubchem"] = _candidate_from_pubchem_row(row)
            except Exception as exc:
                log.warning("PubChemID skeleton search failed for %r: %s", skeleton, exc)

        # ChEBI — index-based prefix scan
        if self._chebi is not None:
            try:
                rows = _chebi_skeleton_search(self._chebi, skeleton)
                if rows:
                    candidates["chebi"] = _candidate_from_chebi_row(rows[0])
            except Exception as exc:
                log.warning("ChEBI skeleton search failed for %r: %s", skeleton, exc)

        return candidates, skeleton

    # ── Source details ────────────────────────────────────────────────────────

    def _build_source_details(
        self, candidates: Dict[str, Optional[Dict[str, Any]]]
    ) -> Dict[str, Dict[str, Any]]:
        """Build a per-source traceability record from the candidates dict.

        For each source, records whether it was found and which output fields
        it has non-null values for.

        Args:
            candidates: Mapping of source key → candidate record.

        Returns:
            Dict mapping display source name to
            ``{"found": bool, "fields": [str, ...]}``.

        Example::

            {
                "ChEBI": {"found": True, "fields": ["name", "SMILES", "InChIKey"]},
                "CompTox": {"found": False, "fields": []},
                ...
            }
        """
        _FIELD_MAP = {
            "name": "name",
            "IUPAC_name": "IUPAC_name",
            "molecular_formula": "molecular_formula",
            "SMILES": "SMILES",
            "InChI": "InChI",
            "InChIKey": "InChIKey",
            "DTXSID": "DTXSID",
            "molecular_mass": "molecular_mass",
            "Synonyms": "Synonyms",
        }

        details: Dict[str, Dict[str, Any]] = {}
        for key in self._SOURCE_KEYS:
            display = self._SOURCE_DISPLAY[key]
            cand = candidates.get(key)
            if cand is None:
                details[display] = {"found": False, "fields": []}
            else:
                fields: List[str] = []
                for cand_field, out_field in _FIELD_MAP.items():
                    val = cand.get(cand_field)
                    if not _is_missing(val):
                        fields.append(out_field)
                # CAS
                cas_vals = cand.get("CAS_candidates") or []
                if cas_vals:
                    fields.append("CASRN")
                details[display] = {"found": True, "fields": sorted(set(fields))}

        return details

    # ── Confidence scoring ────────────────────────────────────────────────────

    def _compute_confidence(
        self,
        match_method: str,
        consensus_score: float,
        *,
        fuzzy_score: Optional[float] = None,
        tanimoto: Optional[float] = None,
    ) -> float:
        """Compute the final confidence score for a result.

        The base confidence depends on the match method.  For fuzzy and
        Tanimoto methods, the raw score is used as the base.  The base is then
        modulated by the cross-source consensus score so that multi-source
        agreement boosts confidence.

        Formula::

            final = base × (0.5 + 0.5 × consensus_score)

        Args:
            match_method: One of the keys in :data:`_BASE_CONFIDENCE`.
            consensus_score: Cross-source consensus agreement in [0, 1].
            fuzzy_score: rapidfuzz similarity in [0, 1]; used when
                ``match_method == "fuzzy_name"``.
            tanimoto: Tanimoto similarity in [0, 1]; used when
                ``match_method == "tanimoto"``.

        Returns:
            Confidence value in [0, 1].
        """
        base = _BASE_CONFIDENCE.get(match_method, 0.5)

        if match_method == "fuzzy_name":
            base = fuzzy_score if fuzzy_score is not None else 0.5
        elif match_method == "tanimoto":
            base = (tanimoto * 0.85) if tanimoto is not None else 0.5

        modulated = base * (0.5 + 0.5 * max(0.0, min(1.0, consensus_score)))
        return round(min(1.0, max(0.0, modulated)), 4)

    # ── Result finalisation ───────────────────────────────────────────────────

    def _finalise_result(
        self,
        result: Dict[str, Any],
        candidates: Dict[str, Optional[Dict[str, Any]]],
        match_method: str,
        source_details: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Apply candidates to result, run consensus, and fill derived fields.

        This method:

        1. Runs :func:`~provesid.tools._compute_consensus` to pick the anchor.
        2. Applies each compatible candidate to ``result`` via
           :func:`~provesid.tools._apply_candidate_to_result`.
        3. Runs :func:`normalize_structure` on the final SMILES to populate
           ``canonical_smiles``, ``kekulized_smiles``, ``InChIKey``, etc.
        4. Validates the source-provided InChIKey against the RDKit-derived one.
        5. Runs salt stripping if ``self.strip_salts`` is True.
        6. Computes and stores the ``confidence`` score.

        Args:
            result: Partially populated result dict.
            candidates: Per-source candidate records.
            match_method: How the match was obtained.
            source_details: Pre-built source traceability dict.

        Returns:
            Fully populated result dict.
        """
        consensus_source, source_match_scores, match_score = _compute_consensus(candidates)
        consensus_candidate = candidates.get(consensus_source) if consensus_source else None

        for source_key in ["chebi", "comptox", "pubchem", "zeropm"]:
            candidate = candidates.get(source_key)
            if _candidate_compatible_with_consensus(candidate, consensus_candidate):
                _apply_candidate_to_result(result, candidate)

        # ChEMBL enrichment
        chembl_cand = candidates.get("chembl")
        if _candidate_compatible_with_consensus(chembl_cand, consensus_candidate):
            _apply_candidate_to_result(result, chembl_cand)

        # Structure normalisation
        norm = normalize_structure(result.get("SMILES"))
        result["canonical_smiles"] = norm["canonical_smiles"]
        result["kekulized_smiles"] = norm["kekulized_smiles"]
        result["molecular_mass"] = _pick_first(result.get("molecular_mass"), norm["mol_weight"])

        # InChI / InChIKey — prefer RDKit-derived values; warn on mismatch
        rdkit_inchi = norm["inchi"]
        rdkit_ik = norm["inchikey"]

        if not _is_missing(result.get("InChIKey")) and not _is_missing(rdkit_ik):
            if result["InChIKey"] != rdkit_ik:
                log.debug(
                    "InChIKey mismatch for query %r: source=%r rdkit=%r",
                    result["query"],
                    result["InChIKey"],
                    rdkit_ik,
                )
        result["InChI"] = _pick_first(result.get("InChI"), rdkit_inchi)
        result["InChIKey"] = _pick_first(result.get("InChIKey"), rdkit_ik)

        # Align name / IUPAC_name
        result["name"] = _pick_first(result.get("name"), result.get("IUPAC_name"))
        result["IUPAC_name"] = _pick_first(result.get("IUPAC_name"), result.get("name"))

        # Source
        result["source"] = _pick_first(
            result.get("source"),
            consensus_candidate.get("source") if consensus_candidate else None,
        )

        # Consensus & score metadata
        result["consensus_source"] = (
            consensus_candidate.get("source") if consensus_candidate else None
        )
        result["source_match_scores"] = {
            (candidates[src].get("source") if candidates.get(src) else src): round(score, 4)
            for src, score in source_match_scores.items()
        }
        result["match_score"] = round(match_score, 4)
        result["source_details"] = source_details
        result["match_method"] = match_method
        result["confidence"] = self._compute_confidence(match_method, match_score)

        # Salt stripping
        if self.strip_salts and not _is_missing(result.get("SMILES")):
            parent = strip_salts(result["SMILES"], self.salt_smarts or None)
            canonical = result.get("canonical_smiles")
            if not _is_missing(parent) and parent != canonical:
                result["parent_smiles"] = parent
                parent_norm = normalize_structure(parent)
                result["parent_inchikey"] = parent_norm["inchikey"]

        return result


# ─────────────────────────────────────────────────────────────────────────────
# Module-level private helpers
# ─────────────────────────────────────────────────────────────────────────────

def _any_candidate(candidates: Dict[str, Optional[Dict[str, Any]]]) -> bool:
    """Return True if at least one candidate is non-None.

    Args:
        candidates: Dict mapping source keys to candidate records.

    Returns:
        True when at least one value is not None.
    """
    return any(v is not None for v in candidates.values())


def _first_smiles_from_candidates(
    candidates: Dict[str, Optional[Dict[str, Any]]]
) -> Optional[str]:
    """Return the first non-missing SMILES found among the candidates.

    Priority order: chebi, comptox, pubchem, zeropm, chembl.

    Args:
        candidates: Dict mapping source keys to candidate records.

    Returns:
        SMILES string or None.
    """
    for key in ["chebi", "comptox", "pubchem", "zeropm", "chembl"]:
        cand = candidates.get(key)
        if cand is None:
            continue
        smiles = cand.get("SMILES") or cand.get("canonical_smiles")
        if not _is_missing(smiles):
            return smiles
    return None


def _most_complete_row(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Select the most data-complete row from a list of source records.

    Completeness is measured as the number of non-null values in the row.

    Args:
        rows: List of source record dicts.

    Returns:
        The row with the most non-null fields, or the first row if the list
        has only one element.
    """
    if not rows:
        return {}
    if len(rows) == 1:
        return rows[0]
    return max(rows, key=lambda r: sum(1 for v in r.values() if not _is_missing(v)))


def _comptox_skeleton_search(
    comptox: CompToxID, skeleton: str
) -> List[Dict[str, Any]]:
    """Search CompTox for InChIKeys sharing the same 14-character skeleton.

    This function queries the CompTox SQLite database with a LIKE predicate on
    the inchikey column.

    Args:
        comptox: Initialised :class:`~provesid.CompToxID` client.
        skeleton: 14-character InChIKey connectivity prefix.

    Returns:
        List of matching rows (may be empty).
    """
    try:
        import sqlite3

        conn = comptox._conn  # type: ignore[attr-defined]
        cur = conn.execute(
            "SELECT * FROM chemicals WHERE INCHIKEY LIKE ? LIMIT 20",
            (f"{skeleton}%",),
        )
        cols = [d[0] for d in cur.description]
        return [dict(zip(cols, row)) for row in cur.fetchall()]
    except Exception as exc:
        log.warning("CompTox skeleton search (SQL) failed: %s", exc)
        return []


def _pubchem_skeleton_search(
    pubchem: PubChemID, skeleton: str
) -> List[Dict[str, Any]]:
    """Search PubChemID SQLite for InChIKeys sharing the same skeleton.

    Args:
        pubchem: Initialised :class:`~provesid.PubChemID` client.
        skeleton: 14-character InChIKey connectivity prefix.

    Returns:
        List of matching rows (may be empty).
    """
    try:
        conn = pubchem._conn  # type: ignore[attr-defined]
        cur = conn.execute(
            "SELECT * FROM compounds WHERE inchikey LIKE ? LIMIT 20",
            (f"{skeleton}%",),
        )
        cols = [d[0] for d in cur.description]
        return [dict(zip(cols, row)) for row in cur.fetchall()]
    except Exception as exc:
        log.warning("PubChemID skeleton search (SQL) failed: %s", exc)
        return []


def _chebi_skeleton_search(
    chebi: ChebiSDF, skeleton: str
) -> List[Dict[str, Any]]:
    """Search the ChebiSDF in-memory index for skeleton-matching InChIKeys.

    Args:
        chebi: Initialised :class:`~provesid.ChebiSDF` client.
        skeleton: 14-character InChIKey connectivity prefix.

    Returns:
        List of matching compound dicts (may be empty).
    """
    try:
        results = []
        ik_index: Dict[str, Any] = chebi.index.get("inchikey_to_id", {})  # type: ignore[attr-defined]
        for ik, chebi_id in ik_index.items():
            if ik.startswith(skeleton):
                compound = chebi.get_compound_by_id(chebi_id)
                if compound:
                    results.append(compound)
                if len(results) >= 20:
                    break
        return results
    except Exception as exc:
        log.warning("ChEBI skeleton search failed: %s", exc)
        return []
