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
from .opsin import PYOPSIN
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
    "opsin": 0.97,
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
    "hit_rank",
    "n_source_support",
    "opsin_smiles",
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
    - **Fuzzy name matching**: rapidfuzz ``ratio`` scorer with configurable
      cut-off (enabled with ``fuzzy=True``).
    - **Tanimoto similarity search**: Morgan-fingerprint-based fallback when
      ``similarity_threshold > 0``.
    - **InChIKey skeleton matching**: 14-character connectivity-layer prefix
      search (enabled with ``inchikey_skeleton=True``).
    - **Salt/solvent stripping**: RDKit SaltRemover + largest-fragment picker
      (enabled with ``strip_salts=True``); ``parent_smiles`` and
      ``parent_inchikey`` populated in results.
    - **Candidate pooling + structure clustering**: the top-``k`` candidates from
      every source are pooled and clustered into distinct compounds by InChIKey
      (skeleton-aware), so a wrong top hit from one source no longer dominates.
    - **Query-aware ranking**: clusters are ranked by a confidence that combines
      the method base, how well each candidate matches the *query itself*
      (name similarity / Tanimoto), and cross-source support.
    - **Multi-hit output**: ``n_hits`` returns the best ``N`` (or ``"all"``)
      distinct compounds per query, ranked with a ``hit_rank`` column.  Default
      is one row per query.
    - **PYOPSIN structure anchoring** (opt-in, ``use_opsin=True``): IUPAC names
      are converted to SMILES offline and used as a high-confidence anchor;
      requires a Java runtime.
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
        n_hits (int | str): Default hits to return per query (int or ``"all"``).
        min_confidence (float): Confidence floor applied before truncation.
        use_opsin (bool): Enable PYOPSIN IUPAC→structure anchoring (needs Java).
        top_k_per_source (int): Candidates pulled per source before pooling.
        cluster_by_skeleton (bool): Merge stereo/charge variants when clustering.
        fuzzy_score_cutoff (float): Fuzzy score cut-off in [0, 100].
        fuzzy_scorer (str): rapidfuzz scorer name.
        consensus_compat_threshold (float): Min similarity to merge with anchor.
        query_weight (float): Weight of query agreement in the confidence score.
        return_alternatives (bool): Attach runner-up summaries when ``n_hits=1``.

    Example::

        from provesid import Search

        s = Search("cas")
        df = s.search(["50-00-0", "64-17-5"])
        print(df[["CASRN", "name", "canonical_smiles", "confidence"]])

        s_fuzzy = Search("name", fuzzy=True)
        df = s_fuzzy.search(["asprin", "paracetamol"])

        # Inspect every plausible interpretation of an ambiguous name
        df = Search("name").search("xylene", n_hits="all")
        print(df[["hit_rank", "name", "InChIKey", "confidence"]])

        # Anchor IUPAC names to a real structure via OPSIN (needs Java)
        df = Search("name", use_opsin=True).search("2-(acetyloxy)benzoic acid")
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

    # rapidfuzz scorer whitelist (name -> scorer callable resolved lazily).
    _FUZZY_SCORERS: frozenset = frozenset(
        ["WRatio", "ratio", "partial_ratio", "token_sort_ratio",
         "token_set_ratio", "QRatio"]
    )

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
        n_hits: Union[int, str] = 1,
        min_confidence: float = 0.0,
        use_opsin: bool = False,
        opsin_jar_fpath: str = "default",
        top_k_per_source: int = 5,
        cluster_by_skeleton: bool = True,
        fuzzy_score_cutoff: float = 80.0,
        fuzzy_scorer: str = "ratio",
        consensus_compat_threshold: float = 0.35,
        query_weight: float = 0.5,
        return_alternatives: bool = False,
        data_dir: Optional[Union[str, Path]] = None,
        redownload: bool = False,
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
            n_hits: Default number of ranked hits to return per query.  Either a
                positive integer or the literal ``"all"``.  Defaults to ``1``
                (one row per query).  Can be overridden per-call in
                :meth:`search`.
            min_confidence: Drop hits whose confidence is below this value
                before truncating to ``n_hits``.  Defaults to ``0.0``.
            use_opsin: Enable PYOPSIN IUPAC-name → structure anchoring for name
                queries.  Requires a Java runtime; falls back to plain name
                matching (with a one-time warning) when unavailable.  Defaults
                to ``False``.
            opsin_jar_fpath: ``jar_fpath`` passed to :class:`~provesid.PYOPSIN`.
            top_k_per_source: Number of candidate rows pulled from each source
                before pooling / clustering.  Defaults to ``5``.
            cluster_by_skeleton: Merge stereo/charge/isotope variants when
                clustering candidates by structure (14-char InChIKey skeleton).
                Defaults to ``True``.
            fuzzy_score_cutoff: rapidfuzz / ZeroPM fuzzy score cut-off in
                [0, 100].  Defaults to ``80.0``.
            fuzzy_scorer: rapidfuzz scorer name; one of ``WRatio``, ``ratio``,
                ``partial_ratio``, ``token_sort_ratio``, ``token_set_ratio``,
                ``QRatio``.  Defaults to ``"ratio"``.  Avoid ``WRatio`` and
                ``partial_ratio``: their partial-ratio term scores a short
                name highly whenever it appears anywhere inside the query, so
                ``fuzzy_score_cutoff`` stops discriminating (see
                :meth:`_name_score`).
            consensus_compat_threshold: Minimum candidate similarity for a
                candidate to be merged with the consensus anchor.  Defaults to
                ``0.35``.
            query_weight: Weight (in [0, 1]) of the query-agreement term versus
                the method base in the confidence formula.  Defaults to ``0.5``.
            return_alternatives: When ``n_hits == 1``, attach compact runner-up
                summaries in an ``alternatives`` column.  Defaults to ``False``.
            data_dir: Optional shared data root used when lazily initialising
                source clients.
            redownload: If True, lazily initialised source clients force a
                fresh dataset download.
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

        # Multi-hit / tuning attributes
        self.n_hits = self._validate_n_hits(n_hits)
        self.min_confidence = float(min_confidence)
        self.use_opsin = bool(use_opsin)
        self.opsin_jar_fpath = opsin_jar_fpath
        self.top_k_per_source = max(1, int(top_k_per_source))
        self.cluster_by_skeleton = bool(cluster_by_skeleton)
        self.fuzzy_score_cutoff = float(fuzzy_score_cutoff)
        if fuzzy_scorer not in self._FUZZY_SCORERS:
            raise ValueError(
                f"fuzzy_scorer must be one of {sorted(self._FUZZY_SCORERS)}, "
                f"got {fuzzy_scorer!r}"
            )
        self.fuzzy_scorer = fuzzy_scorer
        self.consensus_compat_threshold = float(consensus_compat_threshold)
        self.query_weight = float(query_weight)
        self.return_alternatives = bool(return_alternatives)

        self.data_dir = str(data_dir) if data_dir is not None else None
        self.redownload = redownload

        # OPSIN client — created lazily; disabled for the session on failure.
        self._opsin: Optional[PYOPSIN] = None
        self._opsin_available: bool = use_opsin

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
                    setattr(
                        self,
                        attr,
                        factory(data_dir=self.data_dir, redownload=self.redownload),
                    )
                except Exception as exc:
                    log.warning("Could not initialise offline source %s: %s", attr[1:], exc)

        self._clients_initialized = True

    @staticmethod
    def _validate_n_hits(n_hits: Union[int, str]) -> Union[int, str]:
        """Validate and normalise the ``n_hits`` argument.

        Args:
            n_hits: Either a positive integer or the literal ``"all"``.

        Returns:
            ``"all"`` or a positive ``int``.

        Raises:
            ValueError: If ``n_hits`` is neither ``"all"`` nor a positive int.
        """
        if isinstance(n_hits, str):
            if n_hits.lower() == "all":
                return "all"
            raise ValueError(f"n_hits string must be 'all', got {n_hits!r}")
        if isinstance(n_hits, bool) or not isinstance(n_hits, int) or n_hits < 1:
            raise ValueError(f"n_hits must be a positive int or 'all', got {n_hits!r}")
        return n_hits

    def _get_opsin(self) -> Optional[PYOPSIN]:
        """Lazily create the PYOPSIN client; disable for the session on failure.

        Returns:
            A :class:`~provesid.PYOPSIN` instance, or ``None`` when OPSIN is
            disabled or unavailable (e.g. no Java runtime).
        """
        if not self._opsin_available:
            return None
        if self._opsin is None:
            try:
                self._opsin = PYOPSIN(jar_fpath=self.opsin_jar_fpath)
            except Exception as exc:  # pragma: no cover - environment dependent
                log.warning(
                    "PYOPSIN unavailable (%s); disabling OPSIN anchoring for this "
                    "session.", exc,
                )
                self._opsin_available = False
                return None
        return self._opsin

    def _opsin_anchor(self, name: str) -> Optional[Dict[str, Any]]:
        """Convert an IUPAC name to a normalised structure anchor via PYOPSIN.

        Args:
            name: Chemical (IUPAC) name.

        Returns:
            Dict with keys ``smiles``, ``canonical_smiles``, ``inchikey`` when
            OPSIN parsed the name, else ``None``.
        """
        opsin = self._get_opsin()
        if opsin is None:
            return None
        try:
            smiles = opsin.get_smiles(name)
        except Exception as exc:  # pragma: no cover - environment dependent
            log.warning(
                "PYOPSIN parse failed for %r (%s); disabling OPSIN for session.",
                name, exc,
            )
            self._opsin_available = False
            return None
        if _is_missing(smiles) or not str(smiles).strip():
            return None
        norm = normalize_structure(str(smiles))
        return {
            "smiles": str(smiles),
            "canonical_smiles": norm["canonical_smiles"] or str(smiles),
            "inchikey": norm["inchikey"],
        }

    # ── Public entry point ────────────────────────────────────────────────────

    def search(
        self,
        queries: Union[str, List[str], pd.DataFrame, Path],
        *,
        column: Optional[str] = None,
        n_hits: Optional[Union[int, str]] = None,
        min_confidence: Optional[float] = None,
    ) -> pd.DataFrame:
        """Resolve one or more chemical identifiers and return a DataFrame.

        Args:
            queries: Input identifiers in any of the following forms:

                - A single string — returns a one-row DataFrame.
                - A list of strings — one row per query.
                - A :class:`pandas.DataFrame` — the column given by ``column``
                  is used as the query list.  All other columns are preserved
                  in the output (broadcast across the hit rows of each query).
                - A file path (:class:`pathlib.Path` or string ending in
                  ``.csv`` / ``.parquet``) — read into a DataFrame first;
                  ``column`` must be provided.

            column: Column name to read from a DataFrame or file input.
                Required when ``queries`` is a DataFrame or file path.
            n_hits: Per-call override of the instance ``n_hits`` (positive int
                or ``"all"``).  When ``None`` the instance default is used.
            min_confidence: Per-call override of the instance
                ``min_confidence``.  When ``None`` the instance default is used.

        Returns:
            DataFrame with columns defined in :data:`OUTPUT_COLUMNS`.  When
            ``n_hits == 1`` (the default) there is one row per query; otherwise
            up to ``n_hits`` ranked rows per query, ordered by descending
            confidence with a ``hit_rank`` column (0 = best).

        Raises:
            ValueError: If a DataFrame/file input is given but ``column`` is
                not specified, or if ``n_hits`` is invalid.
            FileNotFoundError: If the given file path does not exist.

        Example::

            s = Search("cas")
            df = s.search(["50-00-0", "64-17-5"])
            df = s.search(Path("compounds.csv"), column="CAS")

            # Return every plausible interpretation of an ambiguous name
            s_name = Search("name")
            df = s_name.search("xylene", n_hits="all")
        """
        self._ensure_clients()

        effective_n_hits = (
            self.n_hits if n_hits is None else self._validate_n_hits(n_hits)
        )
        effective_min_conf = (
            self.min_confidence if min_confidence is None else float(min_confidence)
        )

        query_list, extra_df = self._coerce_queries(queries, column)

        iterator = (
            tqdm(query_list, desc=f"Resolving {self.identifier_type.upper()}")
            if self.show_progress
            else query_list
        )

        # Each query yields a list of ranked hit dicts.  Track the source query
        # index so DataFrame/file extra columns can be broadcast across hits.
        rows: List[Dict[str, Any]] = []
        origin_index: List[int] = []
        for q_idx, q in enumerate(iterator):
            hits = self._resolve_single(q, effective_n_hits, effective_min_conf)
            for hit in hits:
                rows.append(hit)
                origin_index.append(q_idx)

        result_df = pd.DataFrame(rows)
        # Ensure all output columns are present (fill missing with None)
        for col in OUTPUT_COLUMNS:
            if col not in result_df.columns:
                result_df[col] = None
        ordered = list(OUTPUT_COLUMNS)
        if self.return_alternatives and "alternatives" in result_df.columns:
            ordered = ordered + ["alternatives"]
        result_df = result_df[ordered]

        # Broadcast extra columns from the original DataFrame across hit rows.
        if extra_df is not None and origin_index:
            extra_cols = [c for c in extra_df.columns if c not in result_df.columns]
            if extra_cols:
                broadcast = extra_df[extra_cols].iloc[origin_index].reset_index(drop=True)
                result_df = pd.concat(
                    [result_df.reset_index(drop=True), broadcast],
                    axis=1,
                )

        return result_df

    # ── Dataset enrichment ────────────────────────────────────────────────────

    def enrich(
        self,
        df: pd.DataFrame,
        column: str,
        *,
        prefix: str = "provesid_",
        n_hits: Optional[Union[int, str]] = None,
    ) -> pd.DataFrame:
        """Add resolved identifier columns to a DataFrame, searching each value once.

        Every *distinct* value in ``column`` is resolved once and the result is
        merged back onto every row that carries it. For measurement tables — where
        the same compound appears in many rows — this is far cheaper than
        resolving row by row, and it is the usual way to attach identifiers to an
        experimental dataset.

        Rows whose ``column`` value is empty, or which do not resolve, keep their
        original data and get empty identifier columns.

        Args:
            df: Input DataFrame. Returned unmodified; the result is a copy.
            column: Column holding the identifier to resolve. Its values are
                compared as stripped strings.
            prefix: Prepended to every added column, so the frame's own columns
                are never overwritten. Defaults to ``"provesid_"``.
            n_hits: Per-call override of the instance ``n_hits``. Leave at
                ``None`` (the default) unless you want more than one hit per
                query — with more than one, a query's rows are duplicated once
                per hit.

        Returns:
            A copy of ``df`` with the :data:`OUTPUT_COLUMNS` added under
            ``prefix``, in the original row order and with the original index.
            When ``n_hits`` yields more than one row per query the index is a
            fresh ``RangeIndex``, since rows no longer correspond one-to-one.

        Raises:
            KeyError: If ``column`` is not in ``df``.
            ValueError: If ``df`` already has columns starting with ``prefix``
                that would collide with the added ones.

        Example::

            import pandas as pd
            from provesid import Search

            #    8 rows, 3 distinct CAS numbers -> only 3 searches
            df = pd.DataFrame({
                "CAS": ["64-17-5", "64-17-5", "50-00-0", "50-78-2"],
                "boiling_point_C": [78.4, 78.2, -19.0, 140.0],
            })
            out = Search("cas").enrich(df, "CAS")
            out[["CAS", "boiling_point_C", "provesid_name", "provesid_InChIKey"]]
        """
        if column not in df.columns:
            raise KeyError(f"Column {column!r} is not in the DataFrame.")

        added = [f"{prefix}{c}" for c in OUTPUT_COLUMNS]
        collisions = [c for c in added if c in df.columns]
        if collisions:
            raise ValueError(
                f"DataFrame already has column(s) {collisions} that enrich() would "
                f"overwrite. Pass a different prefix."
            )

        # Normalise to stripped strings, with every missing form ("", None, NaN,
        # the literal "nan") collapsed to "" so it is never searched.
        key = df[column].map(lambda v: "" if _is_missing(v) else str(v).strip())
        queries = [q for q in key.unique().tolist() if q]

        if not queries:
            log.warning("Column %r has no non-empty values; nothing to resolve.", column)
            out = df.copy()
            for col in added:
                out[col] = None
            return out

        results = self.search(queries, n_hits=n_hits)

        lookup = results.add_prefix(prefix)
        lookup.insert(0, "_enrich_key", lookup[f"{prefix}query"].astype(str))
        if n_hits is None and self.n_hits == 1:
            # One row per query: guarantee a unique merge key so a left merge
            # cannot fan out the caller's rows.
            lookup = lookup.drop_duplicates(subset="_enrich_key", keep="first")

        out = df.copy()
        out["_enrich_key"] = key
        out = out.merge(lookup, on="_enrich_key", how="left").drop(columns="_enrich_key")

        # merge() returns a fresh RangeIndex; restore the caller's index unless
        # multi-hit results changed the row count.
        if len(out) == len(df):
            out.index = df.index
        return out

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

    def _resolve_single(
        self,
        query: str,
        n_hits: Union[int, str],
        min_confidence: float,
    ) -> List[Dict[str, Any]]:
        """Dispatch one query to the appropriate resolver and return ranked hits.

        Each resolver returns ``(base_template, pool, opsin_anchor)``; this
        method clusters the pool, ranks the clusters, and truncates to
        ``n_hits``.

        Args:
            query: A single identifier string.
            n_hits: Number of ranked hits to return (positive int or ``"all"``).
            min_confidence: Drop hits below this confidence before truncation.

        Returns:
            List of result dicts matching :data:`OUTPUT_COLUMNS` (length 1 when
            ``n_hits == 1``).
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
        base_template, pool, opsin_anchor = dispatch[self.identifier_type](query)
        return self._finalise_hits(base_template, pool, n_hits, min_confidence, opsin_anchor)

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
            "hit_rank": 0,
            "n_source_support": 0,
            "opsin_smiles": None,
        }

    # ── Pool construction helpers ─────────────────────────────────────────────

    @staticmethod
    def _tag_candidate(
        cand: Dict[str, Any],
        source_key: str,
        origin_rank: int,
        match_method: str,
        query_match_score: float,
    ) -> Dict[str, Any]:
        """Annotate a candidate record with pool/ranking metadata (in place).

        Args:
            cand: Candidate record from a ``_candidate_from_*`` helper.
            source_key: Originating source key (e.g. ``"chebi"``, ``"opsin"``).
            origin_rank: Rank position within the source's result list (0-based).
            match_method: How the candidate was found (key into
                :data:`_BASE_CONFIDENCE`).
            query_match_score: How well the candidate matches the query in
                [0, 1].

        Returns:
            The same candidate dict, mutated with the transient ``_`` keys.
        """
        cand["_source_key"] = source_key
        cand["_origin_rank"] = int(origin_rank)
        cand["_match_method"] = match_method
        cand["query_match_score"] = float(query_match_score)
        return cand

    def _pool_from_candidates_dict(
        self,
        candidates: Dict[str, Optional[Dict[str, Any]]],
        match_method: str,
        *,
        default_score: float = 1.0,
    ) -> List[Dict[str, Any]]:
        """Convert a per-source ``{key: candidate}`` dict into a tagged pool.

        Args:
            candidates: Mapping of source key → candidate (or ``None``).
            match_method: Match method to tag each candidate with.
            default_score: ``query_match_score`` assigned to every candidate
                (1.0 for exact-identifier matches; a similarity for fuzzy/
                Tanimoto matches).

        Returns:
            List of tagged candidate records.
        """
        pool: List[Dict[str, Any]] = []
        for key in self._SOURCE_KEYS:
            cand = candidates.get(key)
            if cand is None:
                continue
            self._tag_candidate(cand, key, 0, match_method, default_score)
            pool.append(cand)
        return pool

    def _name_score(self, query: str, cand: Dict[str, Any]) -> float:
        """Best similarity between the query name and a candidate's names.

        Compares the query against the candidate ``name``, ``IUPAC_name`` and
        each individual synonym using the configured fuzzy scorer (rapidfuzz)
        when available, falling back to :func:`_text_similarity`.

        Note:
            This is a ranking signal, not evidence of an exact match — use
            :func:`_matches_name_exactly` for that. The default scorer is
            ``ratio``; scorers with a partial-ratio term (``WRatio``,
            ``partial_ratio``) score a short name highly whenever it appears
            anywhere inside the query (``WRatio("caffiene", "ne") == 90``),
            which lets unrelated compounds past ``fuzzy_score_cutoff``.

        Args:
            query: Query name.
            cand: Candidate record.

        Returns:
            Best similarity in [0, 1].
        """
        names = _candidate_names(cand)
        if not names:
            return 0.0

        if RAPIDFUZZ_AVAILABLE and _fuzz is not None:
            scorer = getattr(_fuzz, self.fuzzy_scorer, _fuzz.ratio)
            try:
                return max(scorer(query, n) for n in names) / 100.0
            except Exception:
                pass
        return max(_text_similarity(query, n) for n in names)

    @staticmethod
    def _completeness_score(cand: Dict[str, Any]) -> float:
        """Fraction of key structural/identifier fields populated in [0, 1].

        Used as the query-agreement signal for formula matches (which have no
        name to compare against).

        Args:
            cand: Candidate record.

        Returns:
            Completeness fraction in [0, 1].
        """
        fields = ("SMILES", "InChIKey", "InChI", "molecular_mass", "name", "DTXSID")
        present = sum(1 for f in fields if not _is_missing(cand.get(f)))
        present += 1 if (cand.get("CAS_candidates") or []) else 0
        return present / (len(fields) + 1)

    def _inchikey_pool(
        self, inchikey: str, match_method: str, query_match_score: float
    ) -> List[Dict[str, Any]]:
        """Query every source by InChIKey and return a tagged candidate pool.

        Used by OPSIN anchoring to pull the structurally-correct compound from
        each source regardless of name spelling.

        Args:
            inchikey: Full InChIKey to look up.
            match_method: Match method to tag candidates with.
            query_match_score: Query-agreement score for the candidates.

        Returns:
            List of tagged candidate records (one per source that matched).
        """
        candidates: Dict[str, Optional[Dict[str, Any]]] = {k: None for k in self._SOURCE_KEYS}

        if self._chebi is not None:
            try:
                row = self._chebi.search_by_inchikey(inchikey)
                if row:
                    candidates["chebi"] = _candidate_from_chebi_row(row)
            except Exception as exc:
                log.warning("ChEBI OPSIN-InChIKey lookup failed: %s", exc)
        if self._comptox is not None:
            try:
                row = self._comptox.get_by_inchikey(inchikey)
                if row:
                    candidates["comptox"] = _candidate_from_comptox_row(row)
            except Exception as exc:
                log.warning("CompTox OPSIN-InChIKey lookup failed: %s", exc)
        if self._pubchem is not None:
            try:
                row = self._pubchem.get_by_inchikey(inchikey)
                if row:
                    candidates["pubchem"] = _candidate_from_pubchem_row(row)
            except Exception as exc:
                log.warning("PubChemID OPSIN-InChIKey lookup failed: %s", exc)
        if self._zeropm is not None:
            try:
                table = self._zeropm.get_id_table_from_inchikey(inchikey)
                candidates["zeropm"] = _candidate_from_zeropm_name_table(inchikey, table)
            except Exception as exc:
                log.warning("ZeroPM OPSIN-InChIKey lookup failed: %s", exc)
        if self._chembl is not None:
            try:
                row = self._chembl.search_by_inchikey(inchikey)
                if row:
                    candidates["chembl"] = _candidate_from_chembl_row(row, self._chembl)
            except Exception as exc:
                log.warning("ChEMBL OPSIN-InChIKey lookup failed: %s", exc)

        return self._pool_from_candidates_dict(
            candidates, match_method, default_score=query_match_score
        )

    # ── CAS resolver ─────────────────────────────────────────────────────────

    def _resolve_cas(self, cas: str) -> Tuple[Dict[str, Any], List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Resolve a CAS Registry Number into a unified identifier record.

        Queries ChEBI → CompTox → PubChemID → ZeroPM with waterfall priority,
        then enriches via ChEMBL.

        Args:
            cas: CAS Registry Number string.

        Returns:
            Tuple of (base result template, candidate pool, ``None``).
        """
        result = self._empty_result(cas, "CASRN")
        result["match_method"] = "exact_cas"
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

        pool = self._pool_from_candidates_dict(candidates, "exact_cas")
        return result, pool, None

    # ── Name resolver ─────────────────────────────────────────────────────────

    def _resolve_name(
        self, name: str
    ) -> Tuple[Dict[str, Any], List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Resolve a chemical name into a candidate pool for clustering.

        Pools the top-``k`` candidates from every source (exact name/synonym
        matches first; fuzzy-widened when ``self.fuzzy`` is set and exact
        matches are weak), then adds a PYOPSIN structure anchor and an
        InChIKey-driven structure lookup when ``self.use_opsin`` is enabled.

        Args:
            name: Chemical name string (common or IUPAC).

        Returns:
            Tuple of (base result template, candidate pool, OPSIN anchor dict).
        """
        result = self._empty_result(name, "name")
        result["match_method"] = "exact_name"
        pool = self._candidate_pool_from_name(name)

        # OPSIN structure anchoring (opt-in; needs Java).
        opsin_anchor: Optional[Dict[str, Any]] = None
        if self.use_opsin:
            opsin_anchor = self._opsin_anchor(name)
            if opsin_anchor is not None:
                anchor_cand = _make_candidate(
                    "OPSIN",
                    name=name,
                    smiles=opsin_anchor.get("smiles"),
                    inchikey=opsin_anchor.get("inchikey"),
                )
                self._tag_candidate(anchor_cand, "opsin", 0, "opsin", 1.0)
                pool.append(anchor_cand)
                # Pull the *correct* compound from each source by the OPSIN
                # InChIKey, even when the name spelling differs.
                if not _is_missing(opsin_anchor.get("inchikey")):
                    pool.extend(
                        self._inchikey_pool(str(opsin_anchor["inchikey"]), "opsin", 1.0)
                    )

        return result, pool, opsin_anchor

    def _candidate_pool_from_name(self, name: str) -> List[Dict[str, Any]]:
        """Build a flat, tagged candidate pool from a name query.

        Pulls up to ``self.top_k_per_source`` candidates from each source.
        When ``self.fuzzy`` is enabled and the exact pass yields no strong
        match, the search is widened with non-exact matching and ZeroPM's
        fuzzy ``query_similar_name``.

        Args:
            name: Chemical name to search.

        Returns:
            List of candidate records tagged with ``_source_key``,
            ``_origin_rank``, ``_match_method`` and ``query_match_score``.
        """
        pool: List[Dict[str, Any]] = []
        k = self.top_k_per_source

        def add(cand, source_key, rank, method):
            if cand is None:
                return
            self._tag_candidate(cand, source_key, rank, method, self._name_score(name, cand))
            pool.append(cand)

        # ── Exact pass ──────────────────────────────────────────────────────
        if self._chebi is not None:
            try:
                rows = self._chebi.search_by_name(name, exact=True) or []
                if not rows:
                    rows = self._chebi.search_by_synonym(name, exact=True) or []
                for rank, row in enumerate(rows[:k]):
                    add(_candidate_from_chebi_row(row), "chebi", rank, "exact_name")
            except Exception as exc:
                log.warning("ChEBI name lookup failed for %r: %s", name, exc)

        if self._comptox is not None:
            try:
                rows = self._comptox.search_by_name(name, exact=True, limit=k) or []
                if not rows:
                    row = self._comptox.get_by_name(name)
                    rows = [row] if row else []
                for rank, row in enumerate(rows[:k]):
                    add(_candidate_from_comptox_row(row), "comptox", rank, "exact_name")
            except Exception as exc:
                log.warning("CompTox name lookup failed for %r: %s", name, exc)

        if self._pubchem is not None:
            try:
                rows = self._pubchem.search_by_name(name, exact=True, limit=k) or []
                for rank, row in enumerate(rows[:k]):
                    add(_candidate_from_pubchem_row(row), "pubchem", rank, "exact_name")
            except Exception as exc:
                log.warning("PubChemID name lookup failed for %r: %s", name, exc)

        if self._zeropm is not None:
            try:
                table = self._zeropm.get_id_table_from_name(name)
                add(_candidate_from_zeropm_name_table(name, table), "zeropm", 0, "exact_name")
            except Exception as exc:
                log.warning("ZeroPM name lookup failed for %r: %s", name, exc)

        if self._chembl is not None:
            try:
                rows = self._chembl.search_by_name(name, limit=k, exact=True) or []
                for rank, row in enumerate(rows[:k]):
                    add(_candidate_from_chembl_row(row, self._chembl), "chembl", rank, "exact_name")
            except Exception as exc:
                log.warning("ChEMBL name lookup failed for %r: %s", name, exc)

        # ── Fuzzy widening ──────────────────────────────────────────────────
        # "Strong" means a candidate is genuinely *called* the query name, not
        # merely that it scored highly: WRatio gives a substring hit 85.7, so a
        # score-based test lets one spurious synonym match suppress the widening
        # that would find the right compound.
        cutoff = self.fuzzy_score_cutoff / 100.0
        strong = any(_matches_name_exactly(name, c) for c in pool)
        if self.fuzzy and not strong:
            norm_name = self._normalize_name(name)

            def add_fuzzy(cand, source_key, rank, score=None):
                """Add a fuzzy candidate, keeping only those at or above cutoff.

                ``score`` overrides the name-similarity estimate; pass it when
                the source already reported a true similarity, so it is not
                re-derived from a candidate whose recorded name is the query.
                """
                if cand is None:
                    return
                if score is None:
                    score = self._name_score(name, cand)
                if score < cutoff:
                    return
                self._tag_candidate(cand, source_key, rank, "fuzzy_name", score)
                pool.append(cand)

            if self._chebi is not None:
                try:
                    rows = self._chebi.search_by_name(norm_name, exact=False) or []
                    if not rows:
                        rows = self._chebi.search_by_synonym(norm_name, exact=False) or []
                    for rank, row in enumerate(rows[:k]):
                        add_fuzzy(_candidate_from_chebi_row(row), "chebi", rank)
                except Exception as exc:
                    log.warning("ChEBI fuzzy name lookup failed for %r: %s", name, exc)

            if self._comptox is not None:
                try:
                    rows = self._comptox.search_by_name(norm_name, exact=False, limit=k) or []
                    for rank, row in enumerate(rows[:k]):
                        add_fuzzy(_candidate_from_comptox_row(row), "comptox", rank)
                except Exception as exc:
                    log.warning("CompTox fuzzy name lookup failed for %r: %s", name, exc)

            if self._pubchem is not None:
                try:
                    rows = self._pubchem.search_by_name(norm_name, exact=False, limit=k) or []
                    for rank, row in enumerate(rows[:k]):
                        add_fuzzy(_candidate_from_pubchem_row(row), "pubchem", rank)
                except Exception as exc:
                    log.warning("PubChemID fuzzy name lookup failed for %r: %s", name, exc)

            # ZeroPM is the only source that does true fuzzy *retrieval* (the
            # others are substring-matched with exact=False), so it is the one
            # that can reach a typo like "asprin" -> "aspirin".
            if self._zeropm is not None:
                try:
                    table = self._zeropm.get_id_table_from_similar_name(
                        norm_name,
                        number_of_results=k,
                        score_cutoff=self.fuzzy_score_cutoff,
                    )
                    if table is not None and not table.empty:
                        # Label the candidate with what ZeroPM actually matched,
                        # not with the query, and use its reported similarity.
                        matched_name = str(table["matched_name"].iloc[0])
                        matched_score = float(table["match_score"].iloc[0]) / 100.0
                        add_fuzzy(
                            _candidate_from_zeropm_name_table(matched_name, table),
                            "zeropm",
                            0,
                            score=matched_score,
                        )
                except Exception as exc:
                    log.warning("ZeroPM fuzzy name lookup failed for %r: %s", name, exc)

            if self._chembl is not None:
                try:
                    rows = self._chembl.search_by_name(norm_name, limit=k, exact=False) or []
                    for rank, row in enumerate(rows[:k]):
                        add_fuzzy(_candidate_from_chembl_row(row, self._chembl), "chembl", rank)
                except Exception as exc:
                    log.warning("ChEMBL fuzzy name lookup failed for %r: %s", name, exc)

        return pool

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

    def _resolve_smiles(self, smiles: str) -> Tuple[Dict[str, Any], List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Resolve a SMILES string into a unified identifier record.

        Canonicalises the input, derives an InChIKey, and queries sources by
        InChIKey (ChEBI) or canonical SMILES.  Falls back to Tanimoto
        similarity search when ``self.similarity_threshold > 0`` and no exact
        match is found.

        Args:
            smiles: SMILES string.

        Returns:
            Tuple of (base result template, candidate pool, ``None``).
        """
        result = self._empty_result(smiles, "SMILES")
        result["match_method"] = "exact_smiles"
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
                score = tanimoto_score if tanimoto_score is not None else 0.0
                pool = self._pool_from_candidates_dict(
                    sim_candidates, "tanimoto", default_score=score
                )
                return result, pool, None

        pool = self._pool_from_candidates_dict(candidates, match_method)
        return result, pool, None

    # ── InChI resolver ────────────────────────────────────────────────────────

    def _resolve_inchi(self, inchi: str) -> Tuple[Dict[str, Any], List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Resolve an InChI string into a unified identifier record.

        Converts the InChI to InChIKey via RDKit and delegates to
        :meth:`_resolve_inchikey`.  Also queries sources that store InChI
        directly (ChEBI, CompTox, PubChemID).

        Args:
            inchi: InChI string (must start with ``"InChI="``).

        Returns:
            Tuple of (base result template, candidate pool, ``None``).
        """
        result = self._empty_result(inchi, "InChI")
        result["match_method"] = "inchi"
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

        pool = self._pool_from_candidates_dict(candidates, "inchi")
        return result, pool, None

    # ── InChIKey resolver ─────────────────────────────────────────────────────

    def _resolve_inchikey(self, inchikey: str) -> Tuple[Dict[str, Any], List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Resolve an InChIKey into a unified identifier record.

        Queries all offline sources by InChIKey.  Falls back to 14-character
        skeleton matching when ``self.inchikey_skeleton`` is True and no exact
        match is found.

        Args:
            inchikey: Full 27-character InChIKey
                (``XXXXXXXXXXXXXX-XXXXXXXXXX-X``).

        Returns:
            Tuple of (base result template, candidate pool, ``None``).
        """
        result = self._empty_result(inchikey, "InChIKey")
        result["match_method"] = "exact_inchikey"
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

        pool = self._pool_from_candidates_dict(candidates, match_method)
        return result, pool, None

    # ── DTXSID resolver ───────────────────────────────────────────────────────

    def _resolve_dtxsid(self, dtxsid: str) -> Tuple[Dict[str, Any], List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Resolve a CompTox DTXSID into a unified identifier record.

        Queries CompTox as the primary source, then cross-references other
        sources using the InChIKey derived from the CompTox result.

        Args:
            dtxsid: CompTox DTXSID string (e.g., ``"DTXSID7020182"``).

        Returns:
            Tuple of (base result template, candidate pool, ``None``).
        """
        result = self._empty_result(dtxsid, "DTXSID")
        result["match_method"] = "dtxsid"
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

        pool = self._pool_from_candidates_dict(candidates, "dtxsid")
        return result, pool, None

    # ── Formula resolver ──────────────────────────────────────────────────────

    def _resolve_formula(
        self, formula: str
    ) -> Tuple[Dict[str, Any], List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Resolve a molecular formula into a candidate pool.

        Formulas are not unique identifiers, so all sources may return many
        rows.  The top-``k`` rows per source are pooled and clustered; distinct
        compounds are returned ranked by completeness/consensus.  Confidence is
        capped low (base 0.30) because formula matches are ambiguous.

        Args:
            formula: Molecular formula string (e.g., ``"C9H8O4"``).

        Returns:
            Tuple of (base result template, candidate pool, ``None``).
        """
        result = self._empty_result(formula, "formula")
        result["match_method"] = "formula"
        result["molecular_formula"] = formula

        pool: List[Dict[str, Any]] = []
        k = self.top_k_per_source

        def add(cand, source_key, rank):
            if cand is None:
                return
            # Completeness drives the query_match_score for formula matches.
            score = self._completeness_score(cand)
            self._tag_candidate(cand, source_key, rank, "formula", score)
            pool.append(cand)

        if self._chebi is not None:
            try:
                rows = self._chebi.search_by_formula(formula) or []
                for rank, row in enumerate(_rank_rows_by_completeness(rows)[:k]):
                    add(_candidate_from_chebi_row(row), "chebi", rank)
            except Exception as exc:
                log.warning("ChEBI formula lookup failed for %r: %s", formula, exc)

        if self._comptox is not None:
            try:
                rows = self._comptox.search_by_formula(formula) or []
                for rank, row in enumerate(_rank_rows_by_completeness(rows)[:k]):
                    add(_candidate_from_comptox_row(row), "comptox", rank)
            except Exception as exc:
                log.warning("CompTox formula lookup failed for %r: %s", formula, exc)

        if self._pubchem is not None:
            try:
                rows = self._pubchem.search_by_formula(formula) or []
                for rank, row in enumerate(_rank_rows_by_completeness(rows)[:k]):
                    add(_candidate_from_pubchem_row(row), "pubchem", rank)
            except Exception as exc:
                log.warning("PubChemID formula lookup failed for %r: %s", formula, exc)

        return result, pool, None

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
        query_score: Optional[float] = None,
    ) -> float:
        """Compute the final confidence score for a result.

        The base confidence depends on the match method.  For fuzzy and
        Tanimoto methods, the raw similarity is used as the base.  The base is
        modulated by a query-agreement term (weighted by ``self.query_weight``)
        and by the cross-source consensus score, so that both a strong query
        match and multi-source agreement boost confidence.

        Formula::

            final = base
                  × (w_q × query_score + (1 − w_q))
                  × (0.5 + 0.5 × consensus_score)

        For exact-identifier methods ``query_score`` is 1.0, which collapses the
        middle term to 1.0 and reproduces the original
        ``base × (0.5 + 0.5 × consensus_score)`` behaviour.

        A ``consensus_score`` of exactly 0.0 short-circuits to 0.0 rather than
        following the formula. :func:`~provesid.tools._compute_consensus` only
        returns 0.0 when there were no candidates at all — one source scores 1.0,
        and even two fully disagreeing sources score 0.5 — so a zero consensus
        means nothing matched, and the formula's floor of ``0.5 × base`` would
        report a no-match row as half-confident.

        Args:
            match_method: One of the keys in :data:`_BASE_CONFIDENCE`.
            consensus_score: Cross-source consensus agreement in [0, 1].
            fuzzy_score: rapidfuzz similarity in [0, 1]; used when
                ``match_method == "fuzzy_name"``, scaled by the ``exact_name``
                base so a fuzzy match never outranks an exact one.
            tanimoto: Tanimoto similarity in [0, 1]; used when
                ``match_method == "tanimoto"``.
            query_score: Query-agreement signal in [0, 1] for name/formula
                methods.  Ignored (treated as 1.0) for fuzzy/Tanimoto where the
                similarity already lives in the base.

        Returns:
            Confidence value in [0, 1].
        """
        if consensus_score == 0.0:
            return 0.0

        base = _BASE_CONFIDENCE.get(match_method, 0.5)
        q = 1.0 if query_score is None else max(0.0, min(1.0, query_score))

        if match_method == "fuzzy_name":
            # Scaled by the exact-name base so an approximate name match can
            # never outrank an exact one: a perfect fuzzy score is worth exactly
            # as much as an exact name, and anything less is worth less.
            base = (
                fuzzy_score * _BASE_CONFIDENCE["exact_name"]
                if fuzzy_score is not None
                else 0.5
            )
            q = 1.0  # similarity already captured in base
        elif match_method == "tanimoto":
            base = (tanimoto * 0.85) if tanimoto is not None else 0.5
            q = 1.0

        w_q = max(0.0, min(1.0, self.query_weight))
        query_term = w_q * q + (1.0 - w_q)
        modulated = base * query_term * (0.5 + 0.5 * max(0.0, min(1.0, consensus_score)))
        return round(min(1.0, max(0.0, modulated)), 4)

    # ── Result finalisation ───────────────────────────────────────────────────

    def _finalise_hits(
        self,
        base_template: Dict[str, Any],
        pool: List[Dict[str, Any]],
        n_hits: Union[int, str],
        min_confidence: float,
        opsin_anchor: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Cluster a candidate pool, rank the clusters, and return ranked hits.

        Args:
            base_template: Empty result template (carries query/foundby and any
                pre-populated query fields).
            pool: Flat list of tagged candidate records.
            n_hits: Number of hits to return (positive int or ``"all"``).
            min_confidence: Drop hits below this confidence before truncation.
            opsin_anchor: Optional OPSIN structure anchor for this query.

        Returns:
            List of fully-populated result dicts ordered by descending
            confidence with ``hit_rank`` set.  Always contains at least one row
            (an empty/no-match row when nothing was found).
        """
        opsin_smiles = opsin_anchor.get("smiles") if opsin_anchor else None

        if not pool:
            # No source matched — still finalise structure/salt fields from any
            # pre-populated query fields (e.g. a SMILES/InChI query) via an empty
            # cluster, preserving the resolver's default match_method.
            empty = self._build_result_for_cluster(base_template, {"members": []}, opsin_smiles)
            empty["hit_rank"] = 0
            return [empty]

        clusters = _cluster_candidates(pool, by_skeleton=self.cluster_by_skeleton)
        hits = [
            self._build_result_for_cluster(base_template, cluster, opsin_smiles)
            for cluster in clusters
        ]

        # Rank: OPSIN match first, then confidence, support, query agreement,
        # and (lower) origin rank as a final tie-break.
        hits.sort(
            key=lambda h: (
                1 if h["_opsin_match"] else 0,
                h["confidence"],
                h["n_source_support"],
                h["_cluster_query_score"],
                -h["_min_origin_rank"],
            ),
            reverse=True,
        )

        filtered = [h for h in hits if h["confidence"] >= min_confidence]
        if not filtered:
            # Everything was below the floor — represent the query with a single
            # no-match row so it is not silently dropped.
            empty = self._build_result_for_cluster(base_template, {"members": []}, opsin_smiles)
            empty["hit_rank"] = 0
            return [empty]

        if n_hits != "all":
            filtered = filtered[: int(n_hits)]

        alternatives = None
        if self.return_alternatives and n_hits == 1 and len(hits) > 1:
            alternatives = [
                {
                    "name": h.get("name"),
                    "InChIKey": h.get("InChIKey"),
                    "confidence": h.get("confidence"),
                    "source": h.get("source"),
                }
                for h in hits[1:6]
            ]

        for rank, hit in enumerate(filtered):
            hit["hit_rank"] = rank
            if alternatives is not None and rank == 0:
                hit["alternatives"] = alternatives

        return filtered

    def _build_result_for_cluster(
        self,
        base_template: Dict[str, Any],
        cluster: Dict[str, Any],
        opsin_smiles: Optional[str],
    ) -> Dict[str, Any]:
        """Build one fully-populated result dict from a single structure cluster.

        All members of a cluster denote the same compound; this picks the best
        member per source, runs the existing consensus/merge machinery over
        them, normalises the structure, and computes confidence.

        Args:
            base_template: Empty result template to populate.
            cluster: A cluster dict with a ``members`` list of tagged candidates.
            opsin_smiles: OPSIN SMILES for the query (for the ``opsin_smiles``
                column), if any.

        Returns:
            A populated result dict, plus transient ``_``-prefixed ranking keys
            (dropped before output).
        """
        result = dict(base_template)
        members: List[Dict[str, Any]] = cluster["members"]

        opsin_match = any(m.get("_source_key") == "opsin" for m in members)

        # Best member per data source (lowest origin rank, then best query score).
        per_source: Dict[str, Dict[str, Any]] = {}
        for m in members:
            key = m.get("_source_key")
            if key in (None, "opsin"):
                continue
            cur = per_source.get(key)
            rank_tuple = (m.get("_origin_rank", 0), -m.get("query_match_score", 0.0))
            if cur is None or rank_tuple < (
                cur.get("_origin_rank", 0),
                -cur.get("query_match_score", 0.0),
            ):
                per_source[key] = m

        # Cluster match method = the strongest method among members; fall back
        # to the resolver's default (carried on the template) for empty clusters.
        cluster_method = max(
            (m.get("_match_method", "unknown") for m in members),
            key=lambda mm: _BASE_CONFIDENCE.get(mm, 0.5),
            default=base_template.get("match_method", "unknown"),
        )
        if opsin_match:
            cluster_method = "opsin"

        consensus_source, source_match_scores, match_score = _compute_consensus(per_source)
        consensus_candidate = per_source.get(consensus_source) if consensus_source else None

        for source_key in ["chebi", "comptox", "pubchem", "zeropm"]:
            candidate = per_source.get(source_key)
            if _candidate_compatible_with_consensus(
                candidate, consensus_candidate, self.consensus_compat_threshold
            ):
                _apply_candidate_to_result(result, candidate)

        chembl_cand = per_source.get("chembl")
        if _candidate_compatible_with_consensus(
            chembl_cand, consensus_candidate, self.consensus_compat_threshold
        ):
            _apply_candidate_to_result(result, chembl_cand)

        # OPSIN supplies a structure even when no source row carried one.
        if opsin_match and _is_missing(result.get("SMILES")) and not _is_missing(opsin_smiles):
            result["SMILES"] = opsin_smiles

        # Structure normalisation
        norm = normalize_structure(result.get("SMILES"))
        result["canonical_smiles"] = norm["canonical_smiles"]
        result["kekulized_smiles"] = norm["kekulized_smiles"]
        result["molecular_mass"] = _pick_first(result.get("molecular_mass"), norm["mol_weight"])

        rdkit_inchi = norm["inchi"]
        rdkit_ik = norm["inchikey"]
        if not _is_missing(result.get("InChIKey")) and not _is_missing(rdkit_ik):
            if result["InChIKey"] != rdkit_ik:
                log.debug(
                    "InChIKey mismatch for query %r: source=%r rdkit=%r",
                    result["query"], result["InChIKey"], rdkit_ik,
                )
        result["InChI"] = _pick_first(result.get("InChI"), rdkit_inchi)
        result["InChIKey"] = _pick_first(result.get("InChIKey"), rdkit_ik)

        result["name"] = _pick_first(result.get("name"), result.get("IUPAC_name"))
        result["IUPAC_name"] = _pick_first(result.get("IUPAC_name"), result.get("name"))
        result["source"] = _pick_first(
            result.get("source"),
            consensus_candidate.get("source") if consensus_candidate else None,
        )

        result["consensus_source"] = (
            consensus_candidate.get("source") if consensus_candidate else None
        )
        result["source_match_scores"] = {
            (per_source[src].get("source") if per_source.get(src) else src): round(score, 4)
            for src, score in source_match_scores.items()
        }
        result["match_score"] = round(match_score, 4)
        result["source_details"] = self._build_source_details(per_source)
        result["match_method"] = cluster_method

        cluster_query_score = max(
            (m.get("query_match_score", 0.0) for m in members), default=0.0
        )
        fuzzy_score = cluster_query_score if cluster_method == "fuzzy_name" else None
        tanimoto = cluster_query_score if cluster_method == "tanimoto" else None
        result["confidence"] = self._compute_confidence(
            cluster_method,
            match_score,
            fuzzy_score=fuzzy_score,
            tanimoto=tanimoto,
            query_score=cluster_query_score,
        )
        result["n_source_support"] = len(per_source)
        result["opsin_smiles"] = opsin_smiles

        # Salt stripping
        if self.strip_salts and not _is_missing(result.get("SMILES")):
            parent = strip_salts(result["SMILES"], self.salt_smarts or None)
            canonical = result.get("canonical_smiles")
            if not _is_missing(parent) and parent != canonical:
                result["parent_smiles"] = parent
                parent_norm = normalize_structure(parent)
                result["parent_inchikey"] = parent_norm["inchikey"]

        # Transient ranking metadata (dropped before DataFrame assembly).
        result["_opsin_match"] = opsin_match
        result["_cluster_query_score"] = cluster_query_score
        result["_min_origin_rank"] = min(
            (m.get("_origin_rank", 0) for m in members), default=0
        )
        return result


# ─────────────────────────────────────────────────────────────────────────────
# Cascade resolution for experimental datasets
# ─────────────────────────────────────────────────────────────────────────────

def mw_within(
    tolerance: float = 0.5,
    *,
    reference_column: str = "SMILES",
    name_column: Optional[str] = None,
):
    """Build an ``accept`` predicate that validates a hit by molecular weight.

    Most experimental datasets already carry *some* structure, which makes
    molecular weight a cheap and strict way to tell a correct identifier lookup
    from a plausible-looking wrong one: the same compound gives an exact match,
    so the default tolerance can be tight.

    The returned predicate accepts a hit only when the RDKit molecular weight of
    the hit's structure is within ``tolerance`` of the weight computed from the
    row's own ``reference_column``. It additionally *reports* — without requiring
    — agreement of the canonical SMILES and, when ``name_column`` is given, of
    the name, so :func:`resolve_cascade` can record how much evidence backed
    each row in its ``validated_by`` column.

    Args:
        tolerance: Maximum absolute difference in Da. Defaults to ``0.5``.
        reference_column: Column holding the row's own SMILES, used as the
            reference structure. Defaults to ``"SMILES"``.
        name_column: Optional column holding the row's own name. When given, a
            matching name is reported as an extra ``"name"`` check.

    Returns:
        A callable ``(hit, row) -> list[str]`` suitable for
        :func:`resolve_cascade`'s ``accept`` argument: the names of the checks
        that passed, or an empty list to reject the hit.

    Example::

        accept = mw_within(0.5, reference_column="canonical_SMILES", name_column="name")
        out = resolve_cascade(df, stages, accept=accept)
        out["provesid_validated_by"].value_counts()
        # mw+smiles+name    311
        # mw+smiles          64
        # mw                 12
    """
    def accept(hit: Dict[str, Any], row: Dict[str, Any]) -> List[str]:
        reference = normalize_structure(row.get(reference_column))
        hit_structure = normalize_structure(hit.get("SMILES"))

        reference_mass = reference["mol_weight"]
        hit_mass = hit_structure["mol_weight"]
        if reference_mass is None or hit_mass is None:
            return []
        if abs(hit_mass - reference_mass) > tolerance:
            return []

        passed = ["mw"]

        reference_smiles = reference["canonical_smiles"]
        hit_smiles = hit_structure["canonical_smiles"]
        if reference_smiles and hit_smiles and reference_smiles == hit_smiles:
            passed.append("smiles")

        if name_column is not None:
            wanted = row.get(name_column)
            if not _is_missing(wanted) and _matches_name_exactly(
                str(wanted), {"name": hit.get("name"), "IUPAC_name": hit.get("IUPAC_name")}
            ):
                passed.append("name")

        return passed

    return accept


def resolve_cascade(
    df: pd.DataFrame,
    stages: List[Tuple[str, "Search", str]],
    *,
    accept=None,
    fallback_column: Optional[str] = None,
    prefix: str = "provesid_",
) -> pd.DataFrame:
    """Resolve each row through a series of Search stages; the first hit wins.

    Experimental datasets are annotated unevenly — some rows have a CAS number,
    some only a name, some only a structure. This runs several
    :class:`Search` instances in order, passing to each stage only the rows that
    are still unresolved, so every row is resolved by the most reliable
    identifier it actually has.

    Each hit is checked with ``accept`` before it counts as resolved. A hit that
    fails leaves its row pending for the next stage, which is what stops a
    confident-but-wrong match from ending the cascade. Use :func:`mw_within` for
    the usual molecular-weight check.

    Args:
        df: Input DataFrame. Returned unmodified; the result is a copy.
        stages: Ordered list of ``(label, search, column)`` triples. ``label``
            names the stage in the output, ``search`` is a :class:`Search`
            instance, and ``column`` is the column it reads. Rows with an empty
            value in ``column`` skip that stage.
        accept: Optional ``(hit, row) -> bool | list[str]`` predicate, where
            ``hit`` is the Search result row and ``row`` the input row, both as
            dicts. Return ``True``, or the names of the checks that passed (they
            are joined into ``validated_by``); return ``False`` or an empty list
            to reject. When ``None``, any hit carrying an InChIKey is accepted.

            Both dicts come from DataFrame rows, so a missing field is ``NaN``
            rather than ``None`` — and ``bool(NaN)`` is ``True``. Test emptiness
            with :func:`pandas.isna` (or reuse :func:`mw_within`) rather than
            truthiness.
        fallback_column: Column holding a SMILES from which to derive identifiers
            for rows no stage resolved. Those rows get ``resolved_by="rdkit"``.
            When ``None``, unresolved rows are left empty.
        prefix: Prepended to every added column. Defaults to ``"provesid_"``.

    Returns:
        A copy of ``df`` with the :data:`OUTPUT_COLUMNS` added under ``prefix``,
        plus ``<prefix>resolved_by`` (the stage that resolved the row,
        ``"rdkit"``, or ``"none"``) and ``<prefix>validated_by``.

    Raises:
        KeyError: If a stage names a column that is not in ``df``.
        ValueError: If ``stages`` is empty.

    Example::

        from provesid import Search, resolve_cascade, mw_within

        out = resolve_cascade(
            df,
            stages=[
                ("cas",    Search("cas"),                  "CASRN"),
                ("name",   Search("name", use_opsin=True), "name"),
                ("smiles", Search("smiles"),               "SMILES"),
            ],
            accept=mw_within(0.5, reference_column="SMILES"),
            fallback_column="SMILES",
        )
        out["provesid_resolved_by"].value_counts()
        # cas       412
        # name       98
        # smiles     31
        # rdkit      14
        # none        2
    """
    if not stages:
        raise ValueError("stages must contain at least one (label, search, column).")
    for label, _, column in stages:
        if column not in df.columns:
            raise KeyError(f"Stage {label!r} reads column {column!r}, which is not in the DataFrame.")

    rows = df.reset_index(drop=True)
    pending = list(range(len(rows)))
    resolved: Dict[int, Dict[str, Any]] = {}

    def verdict(hit: Dict[str, Any], row: Dict[str, Any]) -> Optional[str]:
        """Run ``accept`` and return the validated_by text, or None to reject."""
        if accept is None:
            return "inchikey" if not _is_missing(hit.get("InChIKey")) else None
        outcome = accept(hit, row)
        if isinstance(outcome, bool):
            return "accept" if outcome else None
        checks = list(outcome or [])
        return "+".join(checks) if checks else None

    for label, searcher, column in stages:
        if not pending:
            break

        eligible = [i for i in pending if not _is_missing(rows.at[i, column])]
        if not eligible:
            continue

        queries = [str(rows.at[i, column]).strip() for i in eligible]
        hits = searcher.search(queries, n_hits=1).reset_index(drop=True)

        still_pending = []
        for position, i in enumerate(eligible):
            hit = hits.iloc[position].to_dict()
            validated_by = verdict(hit, rows.iloc[i].to_dict())
            if validated_by is None:
                still_pending.append(i)
                continue
            hit["resolved_by"] = label
            hit["validated_by"] = validated_by
            resolved[i] = hit

        eligible_set = set(eligible)
        pending = [i for i in pending if i not in eligible_set] + still_pending
        log.debug(
            "cascade stage %r: %d eligible, %d resolved, %d still pending",
            label, len(eligible), len(eligible) - len(still_pending), len(pending),
        )

    # Terminal RDKit fallback: derive what we can from the row's own structure.
    for i in list(pending):
        structure = (
            normalize_structure(rows.at[i, fallback_column])
            if fallback_column is not None
            else None
        )
        if structure is not None and structure["inchikey"] is not None:
            resolved[i] = {
                "SMILES": structure["canonical_smiles"],
                "canonical_smiles": structure["canonical_smiles"],
                "kekulized_smiles": structure["kekulized_smiles"],
                "InChI": structure["inchi"],
                "InChIKey": structure["inchikey"],
                "molecular_mass": structure["mol_weight"],
                "source": "RDKit",
                "resolved_by": "rdkit",
                "validated_by": "self (rdkit from the given structure)",
            }
        else:
            resolved[i] = {"resolved_by": "none", "validated_by": ""}

    enriched = pd.DataFrame(
        [resolved[i] for i in range(len(rows))],
        columns=OUTPUT_COLUMNS + ["resolved_by", "validated_by"],
    ).add_prefix(prefix)

    out = pd.concat([rows, enriched], axis=1)
    out.index = df.index
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Module-level private helpers
# ─────────────────────────────────────────────────────────────────────────────

def _candidate_names(cand: Dict[str, Any]) -> List[str]:
    """Every name a candidate is known by: name, IUPAC name and each synonym.

    Args:
        cand: Candidate record.

    Returns:
        List of non-empty name strings (may be empty).

    Example::

        _candidate_names({"name": "aspirin", "Synonyms": "ASA; 2-acetoxybenzoic acid"})
        # ["aspirin", "ASA", "2-acetoxybenzoic acid"]
    """
    names: List[str] = []
    for field in ("name", "IUPAC_name"):
        value = cand.get(field)
        if not _is_missing(value):
            names.append(str(value))
    synonyms = cand.get("Synonyms")
    if not _is_missing(synonyms):
        names.extend(s.strip() for s in str(synonyms).split(";") if s.strip())
    return names


def _matches_name_exactly(query: str, cand: Dict[str, Any]) -> bool:
    """Whether the query equals one of the candidate's names.

    Comparison ignores case, leading/trailing whitespace and repeated internal
    whitespace, but nothing else — this is an equality test, not a similarity
    test. It answers "is this actually what the compound is called?", which a
    fuzzy score cannot: ``WRatio("asprin", "Evasprin")`` is 85.7 even though
    the two are different compounds.

    Args:
        query: Query name.
        cand: Candidate record.

    Returns:
        True when one of the candidate's names equals the query.

    Example::

        _matches_name_exactly("Aspirin", {"name": "aspirin"})            # True
        _matches_name_exactly("asprin", {"Synonyms": "Evasprin"})        # False
    """
    def normalise(text: str) -> str:
        return re.sub(r"\s+", " ", str(text).strip().lower())

    target = normalise(query)
    if not target:
        return False
    return any(normalise(name) == target for name in _candidate_names(cand))


def _rank_rows_by_completeness(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Sort source rows by number of non-null fields (most complete first).

    Args:
        rows: List of raw source record dicts.

    Returns:
        New list ordered by descending completeness; stable for ties.
    """
    if not rows:
        return []
    return sorted(
        rows,
        key=lambda r: sum(1 for v in r.values() if not _is_missing(v)),
        reverse=True,
    )


def _candidate_cluster_keys(cand: Dict[str, Any], by_skeleton: bool) -> set:
    """Return the set of structure-identity keys a candidate belongs to.

    Two candidates are merged into the same cluster when their key sets
    intersect.  Priority: full InChIKey (plus 14-char skeleton when
    ``by_skeleton``) → canonical SMILES → normalised name.

    Args:
        cand: A tagged candidate record.
        by_skeleton: Whether to also emit the 14-char InChIKey skeleton key
            (merging stereo/charge/isotope variants).

    Returns:
        A set of hashable key tuples (empty when the candidate has no usable
        identity, signalling a unique singleton cluster).
    """
    keys: set = set()
    ik = cand.get("InChIKey")
    if not _is_missing(ik):
        ik_s = str(ik)
        keys.add(("ik", ik_s))
        if by_skeleton and len(ik_s) >= 14:
            keys.add(("skel", ik_s[:14]))
        return keys

    canon = cand.get("canonical_smiles") or cand.get("SMILES")
    if not _is_missing(canon):
        keys.add(("smi", str(canon)))
        return keys

    name = cand.get("name") or cand.get("IUPAC_name")
    if not _is_missing(name):
        keys.add(("name", str(name).strip().lower()))
    return keys


def _cluster_candidates(
    pool: List[Dict[str, Any]], by_skeleton: bool = True
) -> List[Dict[str, Any]]:
    """Group a candidate pool into structure-identity clusters.

    Uses a union-find over candidate identity keys so that, e.g., a candidate
    carrying a full InChIKey and one carrying only the matching skeleton end up
    in the same cluster.  Candidates with no usable identity become singletons.

    Args:
        pool: Flat list of tagged candidate records.
        by_skeleton: Merge stereo/charge/isotope variants via skeleton keys.

    Returns:
        List of cluster dicts, each ``{"members": [candidate, ...]}``.
    """
    parent: Dict[Any, Any] = {}

    def find(x: Any) -> Any:
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    def union(a: Any, b: Any) -> None:
        parent.setdefault(a, a)
        parent.setdefault(b, b)
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    cand_anchor: List[Any] = []
    for idx, cand in enumerate(pool):
        keys = _candidate_cluster_keys(cand, by_skeleton)
        if not keys:
            keys = {("uniq", idx)}  # unique singleton
        key_list = list(keys)
        for k in key_list:
            parent.setdefault(k, k)
        for k in key_list[1:]:
            union(key_list[0], k)
        cand_anchor.append(key_list[0])

    groups: Dict[Any, List[Dict[str, Any]]] = {}
    for cand, anchor in zip(pool, cand_anchor):
        groups.setdefault(find(anchor), []).append(cand)

    return [{"members": members} for members in groups.values()]


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
