"""PROVESID Search module — unified chemical identifier resolver.

Provides the [`Search`][provesid.search.Search] class for resolving chemical
identifiers across multiple offline databases (ChEBI, CompTox, PubChemID,
ChEMBL) with structure-aware matching, confidence scoring, fuzzy name search,
Tanimoto similarity search, InChIKey-skeleton matching, and salt/solvent
stripping.

The datasets these sources read are large --- ~21 GiB to download and ~6.7 GiB
installed, PubChem's being built from 14.3 GiB of FTP files that are deleted as
they are read, and ChEMBL being compacted from 27.7 GiB to 2.4 GiB as it
arrives --- and none of them is downloaded on the caller's behalf.
[`Search`][provesid.search.Search] queries whatever is installed and reports
what is missing; [`provesid.datasets`][provesid.datasets] installs them by
name.  Pass ``datasets="auto"`` to download what is missing, or
``datasets="required"`` to refuse to run on a partial set.

No socket is opened unless ``online_fallback=True``.  With it, a query that no
offline source answers is asked of PubChem's PUG-REST service and the NCI/CADD
resolver (CACTUS), and only such a query: offline first is a performance and
traffic decision, and it should not cost the answer.  Rows the network
supplied say so in ``source`` and ``source_details``, and
``df.attrs["online_fallbacks"]`` counts how many queries went online.

ZeroPM is **not** among the databases the resolver targets by default.  Its records
are harvested from regulatory inventories rather than curated compound-by-compound,
so its name→structure mappings are noisier than the other four sources and, being
counted as an independent vote, they used to push wrong structures up the
corroboration ranking.  The ZeroPM client itself is untouched and remains available
as [`ZeroPM`][provesid.zeropm.ZeroPM]; pass ``sources="all"`` (or a list
naming ``"zeropm"``) to let [`Search`][provesid.search.Search] query it again.

Supported identifier types:

- ``"cas"``     — CAS Registry Number
- ``"name"``    — Chemical name (common or IUPAC)
- ``"smiles"``  — SMILES string
- ``"inchi"``   — InChI string
- ``"inchikey"``— InChIKey
- ``"dtxsid"``  — CompTox DTXSID
- ``"formula"`` — Molecular formula

Examples:
    >>> from provesid import Search
    >>> df = Search("cas", show_progress=False).search(["50-00-0", "64-17-5"])
    >>> df[["query", "name", "canonical_smiles", "confidence"]]
         query          name canonical_smiles  confidence
    0  50-00-0  formaldehyde              C=O         0.9
    1  64-17-5       ethanol              CCO         0.9

    Typos: fuzzy names are retrieved mainly through ZeroPM, so the
    ``"recall"`` preset, which turns both on, is what rescues most of them.

    >>> df = Search("name", preset="recall", n_hits=1,
    ...             show_progress=False).search(["asprin", "caffiene"])
    >>> df["name"].tolist()
    ['Aspirin', 'caffeine']

    Salts, and InChIKeys that differ only in stereochemistry:

    >>> Search("smiles", strip_salts=True, show_progress=False).search(
    ...     "CC(=O)[O-].[Na+]")[["parent_smiles", "name"]].values.tolist()
    [['CC(=O)[O-]', 'sodium acetate']]
    >>> Search("inchikey", inchikey_skeleton=True, show_progress=False).search(
    ...     "BSYNRYMUTXBXSQ-UHFFFAOYSA-N")["CASRN"].tolist()
    ['50-78-2']
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from types import TracebackType
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Type, Union

import pandas as pd
from tqdm import tqdm

from .chebi_sdf import ChebiSDF
from .chembl import CheMBL
from .comptox import CompToxID
from .datasets import DATASETS, fetch_command, human_bytes, require
from .opsin import PYOPSIN
from .pubchem import PubChemAPI
from .pubchem_id import PubChemID
from .resolver import NCIChemicalIdentifierResolver
from .sources import LOOKUPS, ONLINE_SOURCE_KEYS, SOURCE_DISPLAY, SOURCE_KEYS, Query
from .sqlite_client import DatabaseClosedError
from .zeropm import ZeroPM
from .tools import (
    apply_candidate_to_result,
    candidate_compatible_with_consensus,
    candidate_from_chembl_row,
    candidate_from_pubchem_row,
    compute_consensus,
    formula_element_counts,
    identity_consistent,
    inchikey_from_smiles,
    is_missing,
    is_parsable_smiles,
    make_candidate,
    pick_casrn,
    pick_first,
    structure_formula,
    text_similarity,
)

# ── Optional RDKit ─────────────────────────────────────────────────────────────
try:
    from rdkit import Chem, DataStructs, rdBase
    from rdkit.Chem import AllChem, Descriptors
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from rdkit.Chem.SaltRemover import SaltRemover as _SaltRemover

    RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover
    Chem = None  # type: ignore[assignment]
    rdBase = None  # type: ignore[assignment]
    rdMolStandardize = None  # type: ignore[assignment]
    DataStructs = None  # type: ignore[assignment]
    AllChem = None  # type: ignore[assignment]
    Descriptors = None  # type: ignore[assignment]
    _SaltRemover = None  # type: ignore[assignment]
    RDKIT_AVAILABLE = False

# ── Optional rapidfuzz ─────────────────────────────────────────────────────────
try:
    from rapidfuzz import fuzz as _fuzz

    RAPIDFUZZ_AVAILABLE = True
except ImportError:  # pragma: no cover
    _fuzz = None  # type: ignore[assignment]
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

# Corroboration factor applied to the confidence, keyed by the number of
# *independent* databases that carry the same structure for the query.
#
# Cross-source agreement (``consensus_score``) cannot express this on its own: a
# lone source trivially agrees with itself and scores a perfect 1.0, so without
# this factor an uncorroborated hit outranks a structure that three databases
# agree on.  That is how a single ChEBI row used to beat a CompTox/PubChem/ZeroPM
# consensus and return the wrong compound for a CAS lookup.
#
# ``0`` sources means the cluster came from OPSIN alone (a deterministic
# name-to-structure parse rather than a database record), which is not a
# corroboration question and is left unpenalised.
_SUPPORT_FACTOR: Dict[int, float] = {
    0: 1.00,   # OPSIN-only cluster
    1: 0.85,   # single database, nothing corroborates it
    2: 0.95,   # two databases agree
}
_SUPPORT_FACTOR_MAX = 1.00  # three or more databases agree

# Confidence factor for a result that carries no usable structure.  A row with
# a name and a DTXSID but no structure is not as useful as one with a
# structure, and should not read as the same 0.765; it also makes a structured
# reading of the same query win a tie (CompTox lists eprinomectin as a
# structure-less mixture record, PubChem has a structure for it).
_NO_STRUCTURE_FACTOR = 0.80

# Match methods that claim the query *identifies* the compound.  Only clusters
# found this way can contradict one another; the clusters of a formula, a
# fuzzy name or a Tanimoto search are alternatives by design.
_EXACT_METHODS = frozenset(
    ["exact_inchikey", "exact_smiles", "inchi", "exact_cas", "dtxsid", "opsin", "exact_name"]
)

ON_CONFLICT_POLICIES = frozenset(["pick", "refuse"])
"""What [`Search`][provesid.search.Search] does when independent sources
return different compounds for one query: ``"pick"`` the best-supported one
and flag the row, or ``"refuse"`` to return any of them."""

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
    "structure_status",
    "identity_consistent",
    "source_inchikey",
    "source_formula",
    "n_clusters",
    "conflict",
    "runner_up_inchikey",
    "runner_up_sources",
    "runner_up_confidence",
]
"""The columns of the DataFrame ``Search.search`` returns, in this order.

``return_alternatives=True`` adds an ``alternatives`` column after them. A
column no candidate filled is still there, holding None.

The last four are quality control, and say whether the structure in a row can
be used as it is:

``structure_status``
    Where the structure came from: ``valid`` (the source's SMILES, as given),
    ``repaired:molfile`` / ``repaired:inchi`` / ``repaired:nh_repair`` (the
    source's SMILES did not parse and was rebuilt, with the source's InChIKey
    checked), ``invalid`` (a source gave a structure no one can read, and no
    other source gave a usable one) or ``none`` (no source gave one).  See
    [`STRUCTURE_STATUSES`][provesid.tools.STRUCTURE_STATUSES].
``identity_consistent``
    Whether the InChIKey and formula the sources state belong to the
    structure: True, False, or None when there is nothing to compare.  See
    [`identity_consistent`][provesid.tools.identity_consistent].
``source_inchikey``, ``source_formula``
    The InChIKey and formula as the sources state them.  When there is a
    structure, ``InChIKey``, ``InChI``, ``molecular_formula`` and
    ``molecular_mass`` are computed from it, so they always describe the
    ``SMILES`` of the row; these two keep what the sources said.
``n_clusters``
    How many distinct compounds the sources returned for the query.
``conflict``
    True when an independent source's best answer for the query is a
    different compound --- not a salt, solvate or stereoisomer of this one.
    The confidence of a contested row is multiplied by its share of the
    independent sources behind the two readings, so a one-against-one or
    two-against-two split scores half of what an unopposed answer would.
    ``preset="strict"`` (``on_conflict="refuse"``) returns no structure for
    a conflicted query.
``runner_up_inchikey``, ``runner_up_sources``, ``runner_up_confidence``
    The compound the row competed with: the best conflicting one when
    ``conflict`` is True, otherwise the best other one, if any.  For the top
    row that is the runner-up; for a lower row of ``n_hits="all"``, the
    winner.
"""

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

Hits = Dict[str, List[Dict[str, Any]]]
"""Source key -> that source's candidates, best first (``Search._collect``)."""


def _normalise_sources(sources: Union[str, Sequence[str]]) -> Tuple[str, ...]:
    """The offline sources a [`Search`][provesid.search.Search] is to query.

    Args:
        sources: ``"all"``, one key of
            [`SOURCE_KEYS`][provesid.sources.SOURCE_KEYS], or a sequence of
            them.  A string is one key, not a sequence of letters.

    Returns:
        The named sources, each once, in ``SOURCE_KEYS`` order.

    Raises:
        ValueError: If a name is not a key of ``SOURCE_KEYS``, or none is
            given.

    Examples:
        >>> _normalise_sources(["chembl", "chebi", "chebi"])
        ('chebi', 'chembl')
        >>> _normalise_sources("all")
        ('chebi', 'comptox', 'pubchem', 'zeropm', 'chembl')
    """
    if sources == "all":
        return tuple(SOURCE_KEYS)
    named = [sources] if isinstance(sources, str) else list(sources)
    unknown = [key for key in named if key not in SOURCE_KEYS]
    if unknown:
        raise ValueError(
            f"sources must be 'all' or names from {SOURCE_KEYS}; got {unknown}"
        )
    if not named:
        raise ValueError(f"sources must name at least one of {SOURCE_KEYS}")
    return tuple(key for key in SOURCE_KEYS if key in named)


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

    Examples:
        >>> rec = normalize_structure("c1ccccc1")
        >>> rec["canonical_smiles"], rec["kekulized_smiles"], rec["inchikey"]
        ('c1ccccc1', 'C1=CC=CC=C1', 'UHOVQNZJYSORNB-UHFFFAOYSA-N')
        >>> round(rec["mol_weight"], 3)
        78.114
        >>> normalize_structure("not a smiles")["inchikey"] is None
        True
    """
    empty: Dict[str, Any] = {
        "canonical_smiles": None,
        "kekulized_smiles": None,
        "inchi": None,
        "inchikey": None,
        "mol_weight": None,
        "mol": None,
    }
    if is_missing(smiles) or not RDKIT_AVAILABLE or Chem is None:
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

    Examples:
        >>> strip_salts("[Na+].[Cl-].CC(=O)O")
        'CC(=O)O'
        >>> strip_salts("CC(=O)[O-].[Na+]")    # the anion keeps its charge
        'CC(=O)[O-]'
    """
    if is_missing(smiles) or not RDKIT_AVAILABLE or Chem is None or _SaltRemover is None:
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
    DTXSID, or molecular formula — and queries ChEBI, CompTox, PubChemID and
    ChEMBL to build a harmonised result.  ZeroPM is excluded by default because
    its inventory-derived records are less reliable than the other four sources;
    ``sources="all"`` opts back in, and ``sources=[...]`` picks any subset.

    Features:

    - **Structure-aware matching**: canonicalisation and kekulisation via RDKit;
      InChIKey always derived and reported.
    - **Confidence scoring**: each result carries a ``confidence`` score in
      [0, 1] based on the match method, cross-source consensus, and how many
      independent databases corroborate the structure.
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
    - **No surprise downloads**: the offline datasets are ~21 GiB to fetch and
      ~6.7 GiB installed, and none of them is fetched on your behalf.  ``Search`` uses what is installed and
      reports the rest (``datasets="present"``, the default); install them
      deliberately with [`provesid.datasets.fetch`][provesid.datasets.fetch].
    - **Presets**: ``preset="balanced"`` (the default), ``"strict"`` or
      ``"recall"`` name a whole matching policy in one word; see
      [`PRESETS`][provesid.search.Search.PRESETS].
    - **Online fallback** (opt-in, ``online_fallback=True``): a query no
      offline source answers is retried against PubChem PUG-REST and CACTUS.

    Attributes:
        identifier_type (str): Input identifier type used for all queries.
        preset (str): The [`PRESETS`][provesid.search.Search.PRESETS] entry the
            instance started from.
        settings (dict): The matching and output settings in force, keyed as
            [`PRESETS`][provesid.search.Search.PRESETS]; explicit constructor
            arguments applied.
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
        min_source_support (int): Minimum number of databases that must carry a
            structure for it to be returned (0 disables the filter).
        on_conflict (str): ``"pick"`` or ``"refuse"``: what to return when
            independent sources name different compounds.
        use_opsin (bool): Enable PYOPSIN IUPAC→structure anchoring (needs Java).
        sources (tuple): The offline sources queried, in
            [`SOURCE_KEYS`][provesid.sources.SOURCE_KEYS] order (ZeroPM off by
            default).
        top_k_per_source (int): Candidates pulled per source before pooling.
        cluster_by_skeleton (bool): Merge stereo/charge variants when clustering.
        fuzzy_score_cutoff (float): Fuzzy score cut-off in [0, 100].
        fuzzy_scorer (str): rapidfuzz scorer name.
        consensus_compat_threshold (float): Min similarity to merge with anchor.
        query_weight (float): Weight of query agreement in the confidence score.
        return_alternatives (bool): Attach runner-up summaries when ``n_hits=1``.
        datasets (str): Dataset policy in force --- ``"present"`` (the default),
            ``"auto"`` or ``"required"``.  See the constructor.
        sources_available (list[str]): Source keys that initialised successfully,
            filled in on the first [`search`][provesid.search.Search.search]
            call.  Since corroboration drives confidence, a run missing a
            source scores lower than a full-source run; check this (or
            ``df.attrs["sources_available"]``) before comparing results across
            runs.
        sources_unavailable (list[str]): Source keys that failed to initialise.
        online_fallback (bool): Whether queries no offline source answers are
            retried online.  [`sources_available`][provesid.search.Search]
            lists offline sources only; the online ones are reported per row
            and in ``df.attrs``.

    Examples:
        >>> from provesid import Search
        >>> s = Search("cas", show_progress=False)
        >>> df = s.search(["50-00-0", "64-17-5"])
        >>> df[["CASRN", "name", "canonical_smiles", "confidence"]]
             CASRN          name canonical_smiles  confidence
        0  50-00-0  formaldehyde              C=O         0.9
        1  64-17-5       ethanol              CCO         0.9
        >>> df.attrs["sources_available"]
        ['chebi', 'comptox', 'pubchem', 'chembl']

        Named settings: ``"strict"`` returns only what two databases agree
        on, ``"recall"`` widens every way it can. Explicit arguments still
        win.

        >>> Search.PRESETS["strict"]["min_source_support"]
        2
        >>> Search("name", preset="strict", show_progress=False).search(
        ...     "atrazine")[["name", "CASRN"]].values.tolist()
        [['atrazine', '1912-24-9']]

        Every plausible reading of an ambiguous name:

        >>> df = Search("name", show_progress=False).search("xylene", n_hits="all")
        >>> df[["hit_rank", "name", "InChIKey"]]
           hit_rank      name                     InChIKey
        0         0  o-Xylene  CTQNGGLPUBDAKN-UHFFFAOYSA-N
        1         1  m-Xylene  IVSZLXZYQVIEFR-UHFFFAOYSA-N
        2         2  p-Xylene  URLKBWYHVLBVBO-UHFFFAOYSA-N

        Needing Java, or the network:

        >>> df = Search("name", use_opsin=True).search("2-(acetyloxy)benzoic acid")  # doctest: +SKIP
        >>> df = Search("cas", online_fallback=True).search(["50-78-2", "1912-24-9"])  # doctest: +SKIP
        >>> df.attrs["online_fallbacks"]    # queries that went online     # doctest: +SKIP
        0

        Hand the databases back when the run is over:

        >>> with Search("cas", show_progress=False) as s:
        ...     s.search("50-78-2")["CASRN"].tolist()
        ['50-78-2']
    """

    SUPPORTED_TYPES: frozenset = frozenset(
        ["cas", "name", "smiles", "inchi", "inchikey", "dtxsid", "formula"]
    )

    _SOURCE_DISPLAY: Dict[str, str] = SOURCE_DISPLAY

    # rapidfuzz scorer whitelist (name -> scorer callable resolved lazily).
    _FUZZY_SCORERS: frozenset = frozenset(
        ["WRatio", "ratio", "partial_ratio", "token_sort_ratio",
         "token_set_ratio", "QRatio"]
    )

    PRESETS: Dict[str, Dict[str, Any]] = {
        "balanced": {
            "fuzzy": False,
            "fuzzy_score_cutoff": 80.0,
            "fuzzy_scorer": "ratio",
            "inchikey_skeleton": False,
            "similarity_threshold": 0.0,
            # ZeroPM is a regulatory-inventory harvest rather than a curated
            # compound database, so its rows are kept out of the default
            # corroboration vote.
            "sources": ("chebi", "comptox", "pubchem", "chembl"),
            "top_k_per_source": 5,
            "cluster_by_skeleton": True,
            "consensus_compat_threshold": 0.35,
            "query_weight": 0.5,
            "n_hits": 1,
            "min_confidence": 0.0,
            "min_source_support": 0,
            "on_conflict": "pick",
        },
    }
    """Named settings for the arguments that decide what counts as a match
    and what is returned.  ``Search(..., preset=name)`` starts from one of
    these, and any of its keys passed explicitly overrides the preset's
    value.  ``"balanced"`` is the default and holds the constructor's
    defaults, so it is also the one place those defaults are written down.
    A preset is a name to cite: "resolved with ``Search('cas',
    preset='strict')``" says everything
    [`settings`][provesid.search.Search.settings] would.
    """
    # Precision first: only exact matches, only structures two independent
    # databases agree on, and nothing when independent databases disagree.
    PRESETS["strict"] = {
        **PRESETS["balanced"],
        "min_source_support": 2,
        "on_conflict": "refuse",
    }
    # Recall first: every widening the resolver has, every plausible
    # compound returned, and ZeroPM back in the pool because it is the only
    # source that retrieves by fuzzy name (see _candidate_pool_from_name).
    PRESETS["recall"] = {
        **PRESETS["balanced"],
        "fuzzy": True,
        "inchikey_skeleton": True,
        "similarity_threshold": 0.7,
        "sources": tuple(SOURCE_KEYS),
        "n_hits": "all",
    }

    DATASET_POLICIES: frozenset = frozenset(["present", "auto", "required"])
    """What to do about offline datasets that are not on disk.  ``"present"``
    is the default: a laptop should not spend ~32 GB on a first CAS lookup
    because a source client happens to default to ``auto_download=True``.
    """

    def __init__(
        self,
        identifier_type: str = "cas",
        *,
        preset: str = "balanced",
        strip_salts: bool = False,
        fuzzy: Optional[bool] = None,
        similarity_threshold: Optional[float] = None,
        inchikey_skeleton: Optional[bool] = None,
        show_progress: bool = True,
        salt_smarts: Optional[List[str]] = None,
        n_hits: Optional[Union[int, str]] = None,
        min_confidence: Optional[float] = None,
        min_source_support: Optional[int] = None,
        on_conflict: Optional[str] = None,
        use_opsin: bool = False,
        opsin_jar_fpath: str = "default",
        sources: Optional[Union[str, Sequence[str]]] = None,
        top_k_per_source: Optional[int] = None,
        cluster_by_skeleton: Optional[bool] = None,
        fuzzy_score_cutoff: Optional[float] = None,
        fuzzy_scorer: Optional[str] = None,
        consensus_compat_threshold: Optional[float] = None,
        query_weight: Optional[float] = None,
        return_alternatives: bool = False,
        validate: Optional[Dict[str, Any]] = None,
        online_fallback: bool = False,
        datasets: str = "present",
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
            preset: Named starting point for the matching and output
                settings, one of [`PRESETS`][provesid.search.Search.PRESETS]:

                ``"balanced"``
                    **The default.**  Exact matching only, uncorroborated hits
                    accepted, one row per query.
                ``"strict"``
                    As ``"balanced"``, but a structure is returned only when
                    at least two independent databases carry it
                    (``min_source_support=2``) and no independent database
                    names another compound (``on_conflict="refuse"``).
                    Fewer answers, fewer wrong ones.
                ``"recall"``
                    Fuzzy names, InChIKey-skeleton and Tanimoto (0.7)
                    widening, ZeroPM queried, and every plausible compound
                    returned (``n_hits="all"``).  Read ``confidence`` and
                    ``n_source_support`` before trusting a row.

                The arguments marked *preset* below default to ``None``, which
                takes the preset's value; passing one overrides the preset
                for that argument alone, so ``Search("name",
                preset="strict", n_hits=3)`` is strict with three hits.  The
                values in force are [`settings`][provesid.search.Search.settings].
            strip_salts: Strip salt/solvent fragments and populate
                ``parent_smiles`` / ``parent_inchikey`` columns.
            fuzzy: *Preset.*  Enable fuzzy name matching when an exact name
                match fails.  Requires rapidfuzz.  Balanced: ``False``.
            similarity_threshold: *Preset.*  Tanimoto similarity threshold in
                [0, 1].  When > 0 a Morgan-fingerprint similarity search is run
                as a fallback for SMILES queries with no exact match.  0.0
                disables the search entirely.  Balanced: ``0.0``.
            inchikey_skeleton: *Preset.*  When True, fall back to 14-character
                InChIKey prefix matching when an exact InChIKey match fails.
                Balanced: ``False``.
            show_progress: Display a tqdm progress bar during batch queries.
            salt_smarts: Additional SMARTS patterns passed to
                [`strip_salts`][provesid.search.strip_salts] when ``strip_salts=True``.
            n_hits: *Preset.*  Default number of ranked hits to return per
                query.  Either a positive integer or the literal ``"all"``.
                Balanced: ``1`` (one row per query).  Can be overridden
                per-call in [`search`][provesid.search.Search.search].
            min_confidence: *Preset.*  Drop hits whose confidence is below
                this value before truncating to ``n_hits``.  Balanced: ``0.0``.
            min_source_support: *Preset.*  Minimum number of independent
                databases that must carry a structure for it to be returned.
                ``0`` (balanced) accepts uncorroborated hits; ``2`` (strict)
                requires at least two databases to agree, trading recall for
                precision.  OPSIN-only clusters have no database support and
                are dropped by any value above ``0``.  A record found only
                through another source's answer (ChEMBL for a CAS number) is
                not counted.
            on_conflict: *Preset.*  What to do when independent sources
                return different compounds for the query (``conflict`` in
                the output; see [`OUTPUT_COLUMNS`][provesid.search.OUTPUT_COLUMNS]).
                ``"pick"`` (balanced, recall) returns the best-supported one,
                flagged and with its confidence lowered; ``"refuse"``
                (strict) returns no structure for the query, a no-match row
                with ``conflict=True``.
            use_opsin: Enable PYOPSIN IUPAC-name → structure anchoring for name
                queries.  Requires a Java runtime; falls back to plain name
                matching (with a one-time warning) when unavailable.  Defaults
                to ``False``.
            opsin_jar_fpath: ``jar_fpath`` passed to
                [`PYOPSIN`][provesid.opsin.PYOPSIN].
            sources: *Preset.*  The offline sources to query: any of
                [`SOURCE_KEYS`][provesid.sources.SOURCE_KEYS] (``"chebi"``,
                ``"comptox"``, ``"pubchem"``, ``"zeropm"``, ``"chembl"``), as
                a list, one key, or ``"all"``.  They are queried in
                ``SOURCE_KEYS`` order whatever order they are given in, so
                the answer does not depend on it.  Balanced and strict: all
                but ZeroPM; recall: ``"all"``.  ZeroPM aggregates regulatory
                inventories instead of curating compounds, so its
                name→structure rows are noisier than the other four's yet
                carry the same weight in the corroboration vote.  It is
                chiefly worth adding for fuzzy name queries, since it is the
                only source that does true fuzzy *retrieval* (see
                `_candidate_pool_from_name`).  A source left out is never
                opened, and a client passed for it is ignored with a warning.
                The online services are not listed here; see
                ``online_fallback``.
            top_k_per_source: *Preset.*  Number of candidate rows pulled from
                each source before pooling / clustering.  Balanced: ``5``.
            cluster_by_skeleton: *Preset.*  Merge stereo/charge/isotope
                variants when clustering candidates by structure (14-char
                InChIKey skeleton).  Balanced: ``True``.
            fuzzy_score_cutoff: *Preset.*  rapidfuzz / ZeroPM fuzzy score
                cut-off in [0, 100].  Balanced: ``80.0``.
            fuzzy_scorer: *Preset.*  rapidfuzz scorer name; one of ``WRatio``,
                ``ratio``, ``partial_ratio``, ``token_sort_ratio``,
                ``token_set_ratio``, ``QRatio``.  Balanced: ``"ratio"``.  Avoid ``WRatio`` and
                ``partial_ratio``: their partial-ratio term scores a short
                name highly whenever it appears anywhere inside the query, so
                ``fuzzy_score_cutoff`` stops discriminating (see
                `_name_score`).
            consensus_compat_threshold: *Preset.*  Minimum candidate
                similarity for a candidate to be merged with the consensus
                anchor.  Balanced: ``0.35``.
            query_weight: *Preset.*  Weight (in [0, 1]) of the query-agreement
                term versus the method base in the confidence formula.
                Balanced: ``0.5``.
            return_alternatives: When ``n_hits == 1``, attach compact runner-up
                summaries in an ``alternatives`` column.  Defaults to ``False``.
            validate: Check each hit against reference values in the input
                table, and pass over a hit that fails for the next-ranked
                compound --- what
                [`resolve_cascade`][provesid.search.resolve_cascade] does
                across stages, inside one search.  A dict naming columns of
                the DataFrame (or file) passed to
                [`search`][provesid.search.Search.search]:

                ``"mw"``
                    A reference molecular weight; the hit's must be within
                    ``"mw_tolerance"`` Da (default ``0.5``).
                ``"smiles"``
                    A reference structure; the hit's molecular weight must
                    be within the tolerance of its weight, and an identical
                    canonical SMILES is reported as ``"smiles"``.
                ``"formula"``
                    A reference formula; the hit's must have the same
                    element counts.
                ``"name"``
                    A reference name; reported as ``"name"`` when it is one
                    of the hit's names, never required.

                A row whose reference cell is empty skips that check.  A hit
                without a structure fails every required check it is put
                to.  The checks a hit passed are joined in a
                ``validated_by`` column (``"mw+smiles"``; ``""`` when the
                row had nothing to check against).  When every hit fails,
                the query gets a no-match row.  Defaults to ``None``: hits
                are returned unchecked.
            online_fallback: When True, a query that produced no candidate
                from any offline source --- and only such a query --- is asked
                of PubChem's PUG-REST service and of CACTUS, the NCI/CADD
                Chemical Identifier Resolver.  Their answers are pooled,
                clustered and scored exactly like offline ones, and each
                service is one more independent vote in ``n_source_support``.
                A row they supplied names them in ``source`` and
                ``source_details`` (``"PubChem (online)"``, ``"CACTUS"``), and
                ``df.attrs["online_fallbacks"]`` / ``["online_resolved"]``
                count the queries that went online and those it answered.

                Defaults to ``False``, so that a run opens no socket and a
                batch gives the same answer tomorrow as today.  Formula
                queries are never retried: a formula names thousands of
                PubChem compounds.  A query costs up to three PubChem
                requests (up to ``top_k_per_source + 2`` for a name) and two
                CACTUS requests, paced by the shared per-host limiter; results
                are cached as those clients cache them.  Each fallback is
                logged at DEBUG, and a service that fails is logged at WARNING
                and left out, as a failing database is.
            datasets: What to do about the offline datasets the sources read,
                when they are not on disk.  One of:

                ``"present"``
                    Use whatever is installed and say, once, which sources are
                    missing and what installing them would cost.  **The
                    default.**  Nothing is downloaded.
                ``"auto"``
                    Download whatever is missing, which on a clean machine is
                    ~21 GiB transferred and ~6.7 GiB installed for the four
                    default sources --- but up to ~37 GiB of free disk at the
                    worst moment, while ChEMBL's release is unpacked and
                    compacted.  This was the behaviour before the dataset
                    manager landed, and it happened without asking.
                ``"required"``
                    Raise
                    [`MissingDatasetError`][provesid.datasets.MissingDatasetError]
                    in the constructor, naming every missing dataset and the
                    exact ``provesid.datasets.fetch`` call that installs it.
                    Use this when a run on fewer sources would be worse than no
                    run at all --- confidence scores are not comparable across
                    different source sets.

                Install datasets deliberately with
                [`provesid.datasets.fetch`][provesid.datasets.fetch], and see
                what a download would cost with
                [`provesid.datasets.plan`][provesid.datasets.plan].
            data_dir: Optional shared data root used when lazily initialising
                source clients.
            redownload: If True, lazily initialised source clients force a
                fresh dataset download.  Requires ``datasets="auto"``, since
                the other two policies do not download at all.
            chebi: Pre-initialised [`ChebiSDF`][provesid.chebi_sdf.ChebiSDF]
                client.  Each queried source whose client is left ``None``
                is built on the first search, under the ``datasets`` policy,
                whether or not others were passed; a client that is passed is used
                as given and left open by
                [`close`][provesid.search.Search.close].  To leave a source
                out, leave it out of ``sources``.
            comptox: Pre-initialised [`CompToxID`][provesid.comptox.CompToxID] client.
            pubchem: Pre-initialised
                [`PubChemID`][provesid.pubchem_id.PubChemID] client.
            zeropm: Pre-initialised [`ZeroPM`][provesid.zeropm.ZeroPM] client.  Only
                used when ``sources`` includes ``"zeropm"``.
            chembl: Pre-initialised [`CheMBL`][provesid.chembl.CheMBL] client.

        Raises:
            ValueError: If ``identifier_type`` is not one of the supported
                values, ``preset`` is not a key of
                [`PRESETS`][provesid.search.Search.PRESETS], or ``datasets`` is
                not one of
                [`DATASET_POLICIES`][provesid.search.Search.DATASET_POLICIES],
                or ``redownload=True`` was combined with a policy that does not
                download, or ``sources`` names no source or an unknown one,
                or ``on_conflict`` is not one of
                [`ON_CONFLICT_POLICIES`][provesid.search.ON_CONFLICT_POLICIES].
            provesid.datasets.MissingDatasetError: If ``datasets="required"``
                and a dataset a queried source needs is not on disk.
        """
        if identifier_type not in self.SUPPORTED_TYPES:
            raise ValueError(
                f"identifier_type must be one of {sorted(self.SUPPORTED_TYPES)}, "
                f"got {identifier_type!r}"
            )

        if preset not in self.PRESETS:
            raise ValueError(
                f"preset must be one of {sorted(self.PRESETS)}, got {preset!r}"
            )
        # None means "not passed", so the preset supplies it; anything else
        # was asked for and wins.  No preset key legitimately takes None.
        explicit = {
            "fuzzy": fuzzy,
            "fuzzy_score_cutoff": fuzzy_score_cutoff,
            "fuzzy_scorer": fuzzy_scorer,
            "inchikey_skeleton": inchikey_skeleton,
            "similarity_threshold": similarity_threshold,
            "sources": sources,
            "top_k_per_source": top_k_per_source,
            "cluster_by_skeleton": cluster_by_skeleton,
            "consensus_compat_threshold": consensus_compat_threshold,
            "query_weight": query_weight,
            "n_hits": n_hits,
            "min_confidence": min_confidence,
            "min_source_support": min_source_support,
            "on_conflict": on_conflict,
        }
        chosen = dict(self.PRESETS[preset])
        chosen.update({k: v for k, v in explicit.items() if v is not None})

        self.identifier_type = identifier_type
        self.preset = preset
        self.strip_salts = strip_salts
        self.fuzzy = bool(chosen["fuzzy"])
        self.similarity_threshold = float(chosen["similarity_threshold"])
        self.inchikey_skeleton = bool(chosen["inchikey_skeleton"])
        self.show_progress = show_progress
        self.salt_smarts: List[str] = list(salt_smarts or [])

        # Multi-hit / tuning attributes
        self.n_hits = self._validate_n_hits(chosen["n_hits"])
        self.min_confidence = float(chosen["min_confidence"])
        self.min_source_support = max(0, int(chosen["min_source_support"]))
        if chosen["on_conflict"] not in ON_CONFLICT_POLICIES:
            raise ValueError(
                f"on_conflict must be one of {sorted(ON_CONFLICT_POLICIES)}, "
                f"got {chosen['on_conflict']!r}"
            )
        self.on_conflict: str = chosen["on_conflict"]
        self.use_opsin = bool(use_opsin)
        self.opsin_jar_fpath = opsin_jar_fpath
        self.sources: Tuple[str, ...] = _normalise_sources(chosen["sources"])
        self.top_k_per_source = max(1, int(chosen["top_k_per_source"]))
        self.cluster_by_skeleton = bool(chosen["cluster_by_skeleton"])
        self.fuzzy_score_cutoff = float(chosen["fuzzy_score_cutoff"])
        if chosen["fuzzy_scorer"] not in self._FUZZY_SCORERS:
            raise ValueError(
                f"fuzzy_scorer must be one of {sorted(self._FUZZY_SCORERS)}, "
                f"got {chosen['fuzzy_scorer']!r}"
            )
        self.fuzzy_scorer = chosen["fuzzy_scorer"]
        self.consensus_compat_threshold = float(chosen["consensus_compat_threshold"])
        self.query_weight = float(chosen["query_weight"])
        self.return_alternatives = bool(return_alternatives)
        self.validate: Optional[Dict[str, Any]] = _check_validate_spec(validate)
        self.online_fallback = bool(online_fallback)

        if datasets not in self.DATASET_POLICIES:
            raise ValueError(
                f"datasets must be one of {sorted(self.DATASET_POLICIES)}, "
                f"got {datasets!r}"
            )
        if redownload and datasets != "auto":
            # Silently ignoring it would be worse: the caller asked for a fresh
            # copy and would get a stale one with no indication.
            raise ValueError(
                f"redownload=True downloads, which datasets={datasets!r} does "
                "not permit. Pass datasets='auto' to re-download, or call "
                "provesid.datasets.fetch(..., force=True) yourself."
            )
        self.datasets = datasets

        self.data_dir = str(data_dir) if data_dir is not None else None
        self.redownload = redownload

        # OPSIN client — created lazily; disabled for the session on failure.
        self._opsin: Optional[PYOPSIN] = None
        self._opsin_available: bool = use_opsin

        self._SOURCE_KEYS: List[str] = list(self.sources)

        # A source left out of `sources` is not queried even when its client
        # is passed; otherwise which sources ran would depend on how the
        # caller happened to construct us.
        passed = {
            "chebi": chebi, "comptox": comptox, "pubchem": pubchem,
            "zeropm": zeropm, "chembl": chembl,
        }
        for key, client in passed.items():
            if client is not None and key not in self.sources:
                log.warning(
                    "A %s client was passed but %r is not in sources=%r; it will "
                    "not be queried.",
                    self._SOURCE_DISPLAY[key], key, list(self.sources),
                )
                passed[key] = None

        # Source key -> client, or None until _ensure_clients() builds it (or
        # for good, when it cannot be built).
        self._clients: Dict[str, Any] = {
            **passed,
            "pubchem_online": None,
            "cactus": None,
        }

        # Source keys whose client this instance constructed, and may
        # therefore close.  A client the caller passed in belongs to the
        # caller and outlives this Search; closing it would be closing
        # someone else's database.
        self._owned_clients: List[str] = []
        self._closed: bool = False

        # Whether _ensure_clients() has built the clients the caller did not
        # pass.  Passing some does not count: those are used as given, and the
        # rest are built on the first search as if none had been passed.
        self._clients_initialized: bool = False

        # The web services asked when every offline source missed.  Pooled
        # and reported after the offline sources, and not at all when the
        # fallback is off, so an offline run's source_details is unchanged.
        self._ONLINE_KEYS: List[str] = (
            list(ONLINE_SOURCE_KEYS) if self.online_fallback else []
        )
        self._online_clients_built: bool = False

        # Per search() call: queries retried online, and those it answered.
        self._online_fallbacks: int = 0
        self._online_resolved: int = 0

        # Sources that actually came up, filled in by _ensure_clients().
        self.sources_available: List[str] = []
        self.sources_unavailable: List[str] = []
        self._availability_logged: bool = False

        # "required" is checked here rather than on the first search, so the
        # run fails while the user is still looking at the line that started
        # it.  The check is a directory listing -- no client is constructed and
        # nothing is downloaded.
        if self.datasets == "required":
            require(self._datasets_needed(), self.data_dir)

    @property
    def settings(self) -> Dict[str, Any]:
        """The matching and output settings in force, keyed as
        [`PRESETS`][provesid.search.Search.PRESETS].

        The preset's values with any explicit constructor argument applied, as
        this instance will use them.  [`search`][provesid.search.Search.search]
        records the same dict, with its own per-call overrides applied, in
        ``df.attrs["settings"]``.

        Returns:
            A new dict with one entry per key of ``PRESETS["balanced"]``.

        Example::

            >>> s = Search("name", preset="strict", n_hits=3)
            >>> s.settings["min_source_support"], s.settings["n_hits"]
            (2, 3)
            >>> s.settings == Search.PRESETS["strict"]
            False
        """
        return {key: getattr(self, key) for key in self.PRESETS["balanced"]}

    def _provenance(self, **run_overrides: Any) -> Dict[str, Any]:
        """The ``df.attrs`` entries that say how a result frame was produced.

        Args:
            **run_overrides: Per-call values of
                [`settings`][provesid.search.Search.settings] keys, as
                [`search`][provesid.search.Search.search] resolved them.

        Returns:
            Dict of the preset, the settings in force for the call, the
            ``validate`` spec, the offline sources that backed it, the
            dataset files they read (`_dataset_versions`) and the
            online-fallback counters.
        """
        return {
            "preset": self.preset,
            "settings": {**self.settings, **run_overrides},
            "validate": dict(self.validate) if self.validate is not None else None,
            "sources_available": list(self.sources_available),
            "sources_unavailable": list(self.sources_unavailable),
            "dataset_versions": self._dataset_versions(),
            "online_fallbacks": self._online_fallbacks,
            "online_resolved": self._online_resolved,
        }

    def _dataset_versions(self) -> Dict[str, Any]:
        """Which build of each dataset, and of provesid, produced a frame.

        A result read months later can then be traced to the files behind
        it: a lookup made with the ChEBI index whose byte offsets were
        counted in text mode (fixed in 0.6.0) shows an index
        ``format_version`` of 1.  The files are only ``stat``-ed, never read.

        Returns:
            ``{"provesid": version, <source>: {"path", "bytes", "modified",
            "release"}, ...}`` for every offline source that backed the run;
            ChEBI adds the ``index`` metadata its cached index was built
            with.  ``modified`` is the file's UTC modification time, ISO 8601.
        """
        from datetime import datetime, timezone
        from importlib.metadata import PackageNotFoundError, version
        from .datasets import _release_of

        try:
            provesid_version = version("provesid")
        except PackageNotFoundError:  # pragma: no cover - running from a source tree
            provesid_version = None
        versions: Dict[str, Any] = {"provesid": provesid_version}
        for key in self.sources_available:
            client = self._clients.get(key)
            path = next(
                (value for value in (getattr(client, "sdf_path", None), getattr(client, "db_path", None))
                 if isinstance(value, (str, Path))),
                None,
            )
            if path is None or not Path(path).exists():
                versions[key] = {"path": str(path) if path is not None else None}
                continue
            stat = Path(path).stat()
            entry: Dict[str, Any] = {
                "path": str(path),
                "bytes": stat.st_size,
                "modified": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc)
                .isoformat(timespec="seconds"),
                "release": _release_of(key, str(path)) or None,
            }
            index = getattr(client, "index", None) if key == "chebi" else None
            index_meta = index.get("_meta") if isinstance(index, dict) else None
            if index_meta:
                entry["index"] = dict(index_meta)
            versions[key] = entry
        return versions

    # ── Client lifecycle ──────────────────────────────────────────────────────

    def _datasets_needed(self) -> List[str]:
        """Dataset names this instance would have to open on disk.

        The queried sources (`_SOURCE_KEYS`, from ``sources``) minus any
        whose client the caller constructed and
        passed in --- that client has already found its data, wherever it put
        it, so demanding a copy in the shared data directory would be wrong.

        Returns:
            Dataset names, in `_SOURCE_KEYS` order.
        """
        return [key for key in self._SOURCE_KEYS if self._clients[key] is None]

    def _ensure_clients(self) -> None:
        """Lazily initialise all offline source clients.

        Client construction is idempotent — it only runs once per Search
        instance.  Individual clients that fail to initialise are set to ``None``
        and a warning is logged; the search continues with the remaining sources
        and [`sources_available`][provesid.search.Search] /
        [`sources_unavailable`][provesid.search.Search] record which ones, so a
        three-source run stays distinguishable from a four-source one.

        Only the sources in `_SOURCE_KEYS` are constructed, so a source left
        out of ``sources`` (ZeroPM, by default) is never even opened.

        Whether a missing dataset is downloaded here is the ``datasets``
        policy's decision, and by default it is not: the clients are
        constructed with ``auto_download=False``, a missing one is reported
        with the size and the [`fetch`][provesid.datasets.fetch] call that would
        install it, and the search runs on the sources that are present.
        """
        if self._closed:
            raise DatabaseClosedError(
                "This Search was closed; the databases it opened are no longer "
                "available. Construct a new Search to query again."
            )

        if not self._clients_initialized:
            # Looked up here rather than held on the class, so a test that
            # patches ``provesid.search.PubChemID`` patches what is built.
            factories: Dict[str, Any] = {
                "chebi": ChebiSDF,
                "comptox": CompToxID,
                "pubchem": PubChemID,
                "zeropm": ZeroPM,
                "chembl": CheMBL,
            }
            auto = self.datasets == "auto"
            for key in self._SOURCE_KEYS:
                if self._clients[key] is None:
                    try:
                        self._clients[key] = factories[key](
                            data_dir=self.data_dir,
                            redownload=self.redownload,
                            auto_download=auto,
                        )
                        self._owned_clients.append(key)
                    except FileNotFoundError as exc:
                        if auto:
                            log.warning(
                                "Could not initialise offline source %s: %s", key, exc
                            )
                        else:
                            # Under datasets="present" an absent dataset is an
                            # ordinary state rather than a failure, so the line
                            # says what it would cost and how to install it
                            # instead of reading like an error.
                            dataset = DATASETS[key]
                            log.warning(
                                "%s is not installed, so the %s source is not "
                                "being queried (%s to download, %s on disk). "
                                "Install it with %s, or pass datasets='auto'.",
                                dataset.title, key,
                                human_bytes(dataset.download_bytes),
                                human_bytes(dataset.resident_bytes),
                                fetch_command(key),
                            )
                    except Exception as exc:
                        log.warning("Could not initialise offline source %s: %s", key, exc)

            self._clients_initialized = True

        self.sources_available = [
            key for key in self._SOURCE_KEYS if self._clients[key] is not None
        ]
        self.sources_unavailable = [
            key for key in self._SOURCE_KEYS if self._clients[key] is None
        ]

        # Corroboration drives confidence, so a missing source silently lowers
        # every score it would have voted on — say so once, loudly.
        if self.sources_unavailable and not self._availability_logged:
            log.warning(
                "Search is running with %d of %d sources; unavailable: %s. "
                "Confidence and min_source_support reflect the remaining sources "
                "only, so results are not comparable with a full-source run.",
                len(self.sources_available),
                len(self._SOURCE_KEYS),
                ", ".join(self._SOURCE_DISPLAY[k] for k in self.sources_unavailable),
            )
        self._availability_logged = True

    def _ensure_online_clients(self) -> None:
        """Build the web-service clients, on the first query that needs them.

        Not in `_ensure_clients`, because a run whose every query is
        answered offline should not construct them at all.  Nothing is
        contacted here; the clients only open a connection when asked.
        """
        if self._online_clients_built:
            return
        # Looked up at call time, like the offline factories, so a test that
        # patches ``provesid.search.PubChemAPI`` patches what is built.
        factories: Dict[str, Any] = {
            "pubchem_online": PubChemAPI,
            "cactus": NCIChemicalIdentifierResolver,
        }
        for key in self._ONLINE_KEYS:
            try:
                self._clients[key] = factories[key]()
                self._owned_clients.append(key)
            except Exception as exc:  # pragma: no cover - constructors do no I/O
                log.warning("Could not initialise online source %s: %s", key, exc)
        self._online_clients_built = True

    def close(self) -> None:
        """Close the source clients this instance constructed.

        A [`Search`][provesid.search.Search] may hold four SQLite databases
        open — CompTox, PubChemID, ChEMBL and, when ``sources`` names it, ZeroPM
        — totalling several gigabytes of mapped file.  Until this method
        existed there was no way to hand them back short of dropping the
        ``Search`` and waiting for the collector, which on Windows meant the
        files stayed locked.

        Only clients this instance built are closed.  One passed to the
        constructor belongs to the caller, who may still be using it, and
        closing it here would be closing someone else's database.

        Idempotent.  After it returns, [`search`][provesid.search.Search.search] raises
        [`DatabaseClosedError`][provesid.sqlite_client.DatabaseClosedError] rather than
        quietly running against whatever is left.

        Examples:
            >>> s = Search("cas", show_progress=False)
            >>> s.search("50-00-0")["name"].tolist()
            ['formaldehyde']
            >>> s.close()
            >>> s.search("50-00-0")
            Traceback (most recent call last):
            ...
            provesid.sqlite_client.DatabaseClosedError: ...
        """
        if self._closed:
            return
        self._closed = True

        for key in self._owned_clients:
            close = getattr(self._clients[key], "close", None)
            if close is not None:
                try:
                    close()
                except Exception as exc:  # pragma: no cover - close rarely fails
                    log.warning("Error closing the %s client: %s", key, exc)
            self._clients[key] = None

        self._owned_clients = []

    def __enter__(self) -> "Search":
        """Return the resolver, so ``with Search(...) as s`` binds it.

        Returns:
            (Search): ``self``.
        """
        return self

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc_value: Optional[BaseException],
        traceback: Optional[TracebackType],
    ) -> bool:
        """Close the clients this instance constructed, on the way out.

        Args:
            exc_type: Exception class, or None.
            exc_value: Exception instance, or None.
            traceback: Traceback, or None.

        Returns:
            (bool): False --- an exception raised in the block propagates.
        """
        self.close()
        return False

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
            A [`PYOPSIN`][provesid.opsin.PYOPSIN] instance, or ``None`` when OPSIN is
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
        if is_missing(smiles) or not str(smiles).strip():
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
        min_source_support: Optional[int] = None,
    ) -> pd.DataFrame:
        """Resolve one or more chemical identifiers and return a DataFrame.

        Args:
            queries: Input identifiers in any of the following forms:

                - A single string — returns a one-row DataFrame.
                - A list of strings — one row per query.
                - A `pandas.DataFrame` — the column given by ``column``
                  is used as the query list.  All other columns are preserved
                  in the output (broadcast across the hit rows of each query).
                - A file path (`pathlib.Path` or string ending in
                  ``.csv`` / ``.parquet``) — read into a DataFrame first;
                  ``column`` must be provided.

            column: Column name to read from a DataFrame or file input.
                Required when ``queries`` is a DataFrame or file path.
            n_hits: Per-call override of the instance ``n_hits`` (positive int
                or ``"all"``).  When ``None`` the instance default is used.
            min_confidence: Per-call override of the instance
                ``min_confidence``.  When ``None`` the instance default is used.
            min_source_support: Per-call override of the instance
                ``min_source_support``.  When ``None`` the instance default is
                used.

        Returns:
            DataFrame with columns defined in
            [`OUTPUT_COLUMNS`][provesid.search.OUTPUT_COLUMNS].  When ``n_hits
            == 1`` (the default) there is one row per query; otherwise up to
            ``n_hits`` ranked rows per query, ordered by descending confidence
            with a ``hit_rank`` column (0 = best).

            ``df.attrs["preset"]`` names the preset the instance was built
            from and ``df.attrs["settings"]`` holds the settings this call
            ran with ([`settings`][provesid.search.Search.settings] plus this
            call's ``n_hits``, ``min_confidence`` and ``min_source_support``),
            so a saved frame says how it was made.
            ``df.attrs["sources_available"]`` and
            ``df.attrs["sources_unavailable"]`` record which offline sources
            backed the run (see [`sources_available`][provesid.search.Search]).
             With ``online_fallback=True``, ``df.attrs["online_fallbacks"]``
            counts the queries no offline source answered, which were therefore
            asked online, and ``df.attrs["online_resolved"]`` those of them the
            online services answered.  Both are 0 when the fallback is off.

        Raises:
            ValueError: If a DataFrame/file input is given but ``column`` is
                not specified, or if ``n_hits`` is invalid.
            FileNotFoundError: If the given file path does not exist.

        Examples:
            >>> s = Search("cas", show_progress=False)
            >>> s.search(["50-00-0", "64-17-5"])["name"].tolist()
            ['formaldehyde', 'ethanol']
            >>> table = pd.DataFrame({"CAS": ["50-78-2"], "batch": ["A7"]})
            >>> s.search(table, column="CAS")[["CASRN", "batch"]].values.tolist()
            [['50-78-2', 'A7']]
            >>> df = s.search(Path("compounds.csv"), column="CAS")  # doctest: +SKIP
        """
        self._ensure_clients()
        self._online_fallbacks = 0
        self._online_resolved = 0

        effective_n_hits = (
            self.n_hits if n_hits is None else self._validate_n_hits(n_hits)
        )
        effective_min_conf = (
            self.min_confidence if min_confidence is None else float(min_confidence)
        )
        effective_min_support = (
            self.min_source_support
            if min_source_support is None
            else max(0, int(min_source_support))
        )

        query_list, extra_df = self._coerce_queries(queries, column)
        if self.validate is not None:
            columns = [col for key, col in self.validate.items() if key != "mw_tolerance"]
            if extra_df is None or any(col not in extra_df.columns for col in columns):
                raise ValueError(
                    f"validate= reads the columns {columns} of the input table; "
                    "pass a DataFrame or file holding them."
                )

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
            reference = (
                extra_df.iloc[q_idx].to_dict() if self.validate is not None else None
            )
            hits = self._resolve_single(
                q, effective_n_hits, effective_min_conf, effective_min_support,
                reference=reference,
            )
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
        if self.validate is not None:
            if "validated_by" not in result_df.columns:
                result_df["validated_by"] = None
            ordered = ordered + ["validated_by"]
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

        # Which sources backed this frame — a run degraded by a missing source
        # should not look like a full run afterwards.
        result_df.attrs.update(self._provenance(
            n_hits=effective_n_hits,
            min_confidence=effective_min_conf,
            min_source_support=effective_min_support,
        ))

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
            A copy of ``df`` with the
            [`OUTPUT_COLUMNS`][provesid.search.OUTPUT_COLUMNS] added under
            ``prefix``, in the original row order and with the original index.
            When ``n_hits`` yields more than one row per query the index is a
            fresh ``RangeIndex``, since rows no longer correspond one-to-one.
            ``df.attrs`` carries the same provenance
            [`search`][provesid.search.Search.search] records: the preset and
            settings, which offline sources backed the run and, with
            ``online_fallback=True``, how many queries went online.

        Raises:
            KeyError: If ``column`` is not in ``df``.
            ValueError: If ``df`` already has columns starting with ``prefix``
                that would collide with the added ones.

        Examples:
            >>> # 4 rows, 3 distinct CAS numbers -> only 3 searches
            >>> df = pd.DataFrame({
            ...     "CAS": ["64-17-5", "64-17-5", "50-00-0", "50-78-2"],
            ...     "boiling_point_C": [78.4, 78.2, -19.0, 140.0],
            ... })
            >>> out = Search("cas", show_progress=False).enrich(df, "CAS")
            >>> out[["CAS", "boiling_point_C", "provesid_name", "provesid_InChIKey"]]
                   CAS  boiling_point_C         provesid_name            provesid_InChIKey
            0  64-17-5             78.4               ethanol  LFQSCWFLJHTTHZ-UHFFFAOYSA-N
            1  64-17-5             78.2               ethanol  LFQSCWFLJHTTHZ-UHFFFAOYSA-N
            2  50-00-0            -19.0          formaldehyde  WSFSSNUMVMOOMR-UHFFFAOYSA-N
            3  50-78-2            140.0  acetylsalicylic acid  BSYNRYMUTXBXSQ-UHFFFAOYSA-N
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
        key = df[column].map(lambda v: "" if is_missing(v) else str(v).strip())
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

        # Carry the source provenance of the underlying search (merge drops attrs).
        # Read from the instance, not results.attrs, which a stubbed search()
        # need not set.
        run_overrides = {} if n_hits is None else {"n_hits": self._validate_n_hits(n_hits)}
        out.attrs.update(self._provenance(**run_overrides))
        return out

    # ── Input normalisation ───────────────────────────────────────────────────

    def _coerce_queries(
        self,
        queries: Union[str, List[str], pd.DataFrame, Path],
        column: Optional[str],
    ) -> Tuple[List[str], Optional[pd.DataFrame]]:
        """Convert the ``queries`` argument to a plain list of strings.

        Args:
            queries: Raw input from [`search`][provesid.search.Search.search].
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
        min_source_support: int = 0,
        reference: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Dispatch one query to the appropriate resolver and return ranked hits.

        Each resolver returns ``(base_template, pool, opsin_anchor)``; this
        method clusters the pool, ranks the clusters, and truncates to
        ``n_hits``.  An empty pool is where the online fallback happens, so
        that no resolver has to know about it (see `_online_pool`).

        Args:
            query: A single identifier string.
            n_hits: Number of ranked hits to return (positive int or ``"all"``).
            min_confidence: Drop hits below this confidence before truncation.
            min_source_support: Drop hits corroborated by fewer than this many
                databases before truncation.
            reference: The query's input row, for ``validate``.

        Returns:
            List of result dicts matching
            [`OUTPUT_COLUMNS`][provesid.search.OUTPUT_COLUMNS] (length 1 when
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
        if not pool and self._ONLINE_KEYS:
            pool = self._online_pool(query, base_template["match_method"])
        return self._finalise_hits(
            base_template,
            pool,
            n_hits,
            min_confidence,
            opsin_anchor,
            min_source_support=min_source_support,
            reference=reference,
        )

    # ── Empty result template ─────────────────────────────────────────────────

    def _empty_result(self, query: str, foundby: str) -> Dict[str, Any]:
        """Return a result dict with all fields initialised to None/defaults.

        Args:
            query: The original query string.
            foundby: The identifier type used for the search.

        Returns:
            Dict with all [`OUTPUT_COLUMNS`][provesid.search.OUTPUT_COLUMNS]
            keys present.
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
            "structure_status": "none",
            "identity_consistent": None,
            "source_inchikey": None,
            "source_formula": None,
            "n_clusters": 0,
            "conflict": False,
            "runner_up_inchikey": None,
            "runner_up_sources": None,
            "runner_up_confidence": None,
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
                `_BASE_CONFIDENCE`).
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

    def _collect(
        self,
        kind: str,
        value: str,
        *,
        label: Optional[str] = None,
        k: int = 1,
        sources: Optional[List[str]] = None,
    ) -> Hits:
        """Ask every available source one question from the lookup table.

        This is the only place a source is queried.  Each source that has a
        client and a row for ``kind`` in
        [`provesid.sources.LOOKUPS`][provesid.sources.LOOKUPS] is asked in
        turn.  One that raises is logged and left out, so a broken database
        costs its own vote rather than the query.

        Args:
            kind: Lookup kind, a key of [`LOOKUPS`][provesid.sources.LOOKUPS],
                such as ``"cas"`` or ``"fuzzy_name"``.
            value: The identifier to look up.
            label: Name for a ZeroPM candidate, when it should not be
                ``value`` (see [`Query`][provesid.sources.Query]).
            k: Candidates to take from each source.
            sources: Restrict the question to these source keys.  Defaults
                to every queried source.

        Returns:
            Source key -> that source's candidates, best first.  Sources that
            were asked and found nothing map to an empty list; sources that
            were not asked, or failed, are absent.
        """
        lookups = LOOKUPS[kind]
        query = Query(value, label=label, k=k, fuzzy_cutoff=self.fuzzy_score_cutoff)
        hits: Hits = {}
        for key in self._SOURCE_KEYS if sources is None else sources:
            client, lookup = self._clients.get(key), lookups.get(key)
            if client is None or lookup is None:
                continue
            try:
                hits[key] = lookup(client, query)
            except Exception as exc:
                # A held host was reported once, when the hold was recorded;
                # a batch should not repeat it for every query that follows.
                held = getattr(exc, "held_until", None) is not None
                log.log(
                    logging.DEBUG if held else logging.WARNING,
                    "%s %s lookup failed for %r: %s",
                    self._SOURCE_DISPLAY[key], kind, value, exc,
                )
        return hits

    def _pool(
        self,
        hits: Hits,
        match_method: str,
        score: Union[float, Callable[[Dict[str, Any]], float]] = 1.0,
    ) -> List[Dict[str, Any]]:
        """Flatten per-source hits into a tagged candidate pool.

        Candidates are pooled in `_SOURCE_KEYS` order, then the online
        services', and, within a source, in the order the source ranked them.

        Args:
            hits: Source key -> candidates, as `_collect` returns.
            match_method: Match method to tag each candidate with.
            score: The ``query_match_score`` of every candidate: a number
                (1.0 for exact-identifier matches), or a function of the
                candidate for matches whose quality varies, such as names.

        Returns:
            List of tagged candidate records.
        """
        pool: List[Dict[str, Any]] = []
        for key in self._SOURCE_KEYS + self._ONLINE_KEYS:
            for rank, cand in enumerate(hits.get(key) or []):
                cand_score = score(cand) if callable(score) else score
                pool.append(self._tag_candidate(cand, key, rank, match_method, cand_score))
        return pool

    def _name_score(self, query: str, cand: Dict[str, Any]) -> float:
        """Best similarity between the query name and a candidate's names.

        Compares the query against the candidate ``name``, ``IUPAC_name`` and
        each individual synonym using the configured fuzzy scorer (rapidfuzz)
        when available, falling back to
        [`text_similarity`][provesid.tools.text_similarity].

        Note:
            This is a ranking signal, not evidence of an exact match — use
            `_matches_name_exactly` for that. The default scorer is
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
        return max(text_similarity(query, n) for n in names)

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
        present = sum(1 for f in fields if not is_missing(cand.get(f)))
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
        return self._pool(self._collect("inchikey", inchikey), match_method, query_match_score)

    def _online_pool(self, query: str, match_method: str) -> List[Dict[str, Any]]:
        """Ask the online services a query no offline source answered.

        The question is the query itself, asked as its own identifier type:
        the identifier types and the lookup kinds share their names.  The
        cross-source routes the offline resolvers take (ChEMBL by the SMILES a
        CAS lookup found, and so on) have nothing to start from here, since
        nothing was found.  A kind with no online row --- ``formula`` ---
        asks nothing and is not counted as a fallback.

        Args:
            query: The query, as the user gave it.
            match_method: The resolver's match method, which the online
                candidates are tagged with: a CAS number PubChem knows is as
                exact a CAS match as one a database knows.

        Returns:
            The tagged candidate pool, empty when neither service answered.
        """
        kind = self.identifier_type
        if not any(key in LOOKUPS[kind] for key in self._ONLINE_KEYS):
            return []

        self._ensure_online_clients()
        self._online_fallbacks += 1
        log.debug(
            "No offline source answered %s %r; asking %s.", kind, query,
            ", ".join(self._SOURCE_DISPLAY[key] for key in self._ONLINE_KEYS),
        )

        k = self.top_k_per_source if kind == "name" else 1
        hits = self._collect(kind, query, k=k, sources=self._ONLINE_KEYS)
        score: Union[float, Callable[[Dict[str, Any]], float]] = (
            (lambda cand: self._name_score(query, cand)) if kind == "name" else 1.0
        )
        pool = self._pool(hits, match_method, score)

        if pool:
            self._online_resolved += 1
        log.debug(
            "Online fallback for %r: %s.", query,
            ", ".join(
                f"{self._SOURCE_DISPLAY[key]} {len(hits[key])}" for key in hits
            ) or "no service answered",
        )
        return pool

    # ── CAS resolver ─────────────────────────────────────────────────────────

    def _resolve_cas(self, cas: str) -> Tuple[Dict[str, Any], List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Resolve a CAS Registry Number into a unified identifier record.

        Queries ChEBI, CompTox and PubChemID (and ZeroPM when ``sources``
        names it) by CAS number.  ChEMBL records no CAS numbers, so it is
        asked for each structure the others found (see
        `_cross_references`); what it finds is an echo of their answer, not
        a vote.

        The template leaves ``CASRN`` empty, so a hit reports the compound's
        current number, as a name or structure search does, and not
        necessarily the one queried: atrazine found by the retired
        ``39400-72-1`` is reported as ``1912-24-9``. The number queried stays
        in ``query``.

        Args:
            cas: CAS Registry Number string.

        Returns:
            Tuple of (base result template, candidate pool, ``None``).
        """
        result = self._empty_result(cas, "CASRN")
        result["match_method"] = "exact_cas"

        hits = self._collect("cas", cas)
        hits.update(self._cross_references(hits, ["chembl"]))

        return result, self._pool(hits, "exact_cas"), None

    def _cross_references(
        self, hits: Hits, sources: List[str], *, label: Optional[str] = None
    ) -> Hits:
        """Look up in ``sources`` every structure the other sources found.

        For sources that do not hold the queried identifier.  Each distinct
        InChIKey among the candidates is looked up, so that every compound
        the sources proposed --- not only the first --- gets its record, and
        each record found is marked ``derived_from`` the source whose key
        found it.  Such a candidate fills fields the others left empty but
        does not vote: it agrees with the answer that found it by
        construction, and when it voted, a wrong ChEBI hit for a CAS number
        came with a free ChEMBL confirmation (see
        [`compute_consensus`][provesid.tools.compute_consensus]).

        Args:
            hits: The other sources' candidates, as `_collect` returns them.
            sources: Source keys to look the InChIKeys up in; those not
                queried by this search are skipped.
            label: Name for a ZeroPM candidate (see `_collect`).

        Returns:
            Source key -> derived candidates, for the sources that found
            something; ``{}`` when there was no InChIKey to look up.
        """
        sources = [
            key for key in sources
            if key in self._SOURCE_KEYS and self._clients.get(key) is not None
        ]
        origins: Dict[str, str] = {}
        for key in SOURCE_KEYS:
            for cand in hits.get(key) or []:
                inchikey = cand.get("InChIKey")
                if not is_missing(inchikey):
                    origins.setdefault(str(inchikey), self._SOURCE_DISPLAY[key])
        if not sources or not origins:
            return {}

        found: Hits = {}
        for inchikey, origin in origins.items():
            echoes = self._collect("inchikey", inchikey, label=label, sources=sources)
            for key, cands in echoes.items():
                for cand in cands:
                    cand["derived_from"] = origin
                    found.setdefault(key, []).append(cand)
        return found

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
                anchor_cand = make_candidate(
                    "OPSIN",
                    name=name,
                    smiles=opsin_anchor.get("smiles"),
                    inchikey=opsin_anchor.get("inchikey"),
                )
                self._tag_candidate(anchor_cand, "opsin", 0, "opsin", 1.0)
                pool.append(anchor_cand)
                # Pull the *correct* compound from each source by the OPSIN
                # InChIKey, even when the name spelling differs.
                if not is_missing(opsin_anchor.get("inchikey")):
                    pool.extend(
                        self._inchikey_pool(str(opsin_anchor["inchikey"]), "opsin", 1.0)
                    )

        return result, pool, opsin_anchor

    def _candidate_pool_from_name(self, name: str) -> List[Dict[str, Any]]:
        """Build a flat, tagged candidate pool from a name query.

        Pulls up to ``self.top_k_per_source`` candidates from each source.
        When ``self.fuzzy`` is enabled and the exact pass yields no strong
        match, the search is widened with non-exact matching and — only when
        ``sources`` names ZeroPM — its fuzzy ``get_id_table_from_similar_name``.

        Args:
            name: Chemical name to search.

        Returns:
            List of candidate records tagged with ``_source_key``,
            ``_origin_rank``, ``_match_method`` and ``query_match_score``.
        """
        k = self.top_k_per_source
        pool = self._pool(
            self._collect("name", name, k=k),
            "exact_name",
            lambda cand: self._name_score(name, cand),
        )

        # ── Fuzzy widening ──────────────────────────────────────────────────
        # "Strong" means a candidate is genuinely *called* the query name, not
        # merely that it scored highly: WRatio gives a substring hit 85.7, so a
        # score-based test lets one spurious synonym match suppress the widening
        # that would find the right compound.
        strong = any(_matches_name_exactly(name, c) for c in pool)
        if self.fuzzy and not strong:
            # ZeroPM is the only source that does true fuzzy *retrieval*, and
            # reports the similarity it matched on; that score is kept rather
            # than re-derived from the name ZeroPM's candidate was given.  It
            # is off unless sources names it, which is the cost of dropping it:
            # a typo that shares no substring with the real name stays
            # unresolved.
            def fuzzy_score(cand: Dict[str, Any]) -> float:
                reported = cand.get("query_match_score")
                return reported if reported is not None else self._name_score(name, cand)

            cutoff = self.fuzzy_score_cutoff / 100.0
            widened = self._pool(
                self._collect("fuzzy_name", self._normalize_name(name), k=k),
                "fuzzy_name",
                fuzzy_score,
            )
            pool.extend(c for c in widened if c["query_match_score"] >= cutoff)

        return pool

    # ── SMILES resolver ───────────────────────────────────────────────────────

    def _resolve_smiles(self, smiles: str) -> Tuple[Dict[str, Any], List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Resolve a SMILES string into a unified identifier record.

        Canonicalises the input, derives an InChIKey, and queries sources by
        SMILES (retrying CompTox and PubChemID with the canonical form) and
        ChEBI by the InChIKey.  Falls back to Tanimoto similarity search when
        ``self.similarity_threshold > 0`` and no exact match is found.

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
        inchikey = norm["inchikey"] or inchikey_from_smiles(smiles)

        hits = self._collect("smiles", smiles)
        if canonical != smiles:
            retry = [key for key in ("comptox", "pubchem") if hits.get(key) == []]
            hits.update(self._collect("smiles", canonical, sources=retry))
        if not is_missing(inchikey):
            hits.update(self._collect("inchikey", str(inchikey), sources=["chebi"]))

        # Tanimoto similarity fallback
        if not _any_candidate(hits) and self.similarity_threshold > 0:
            similar, tanimoto_score = self._tanimoto_candidates(smiles)
            if _any_candidate(similar):
                score = tanimoto_score if tanimoto_score is not None else 0.0
                return result, self._pool(similar, "tanimoto", score), None

        return result, self._pool(hits, "exact_smiles"), None

    # ── InChI resolver ────────────────────────────────────────────────────────

    def _resolve_inchi(self, inchi: str) -> Tuple[Dict[str, Any], List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Resolve an InChI string into a unified identifier record.

        Queries the sources that store InChI directly (ChEBI, PubChemID and,
        when enabled, ZeroPM), then CompTox by the InChIKey and ChEMBL by the
        SMILES that RDKit derives from the InChI.

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
        if not is_missing(inchikey):
            result["InChIKey"] = inchikey
        if not is_missing(smiles):
            result["SMILES"] = smiles

        hits = self._collect("inchi", inchi)
        if not is_missing(inchikey):
            hits.update(self._collect("inchikey", str(inchikey), sources=["comptox"]))
        if not is_missing(smiles):
            hits.update(self._collect("smiles", str(smiles), sources=["chembl"]))

        return result, self._pool(hits, "inchi"), None

    # ── InChIKey resolver ─────────────────────────────────────────────────────

    def _resolve_inchikey(self, inchikey: str) -> Tuple[Dict[str, Any], List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Resolve an InChIKey into a unified identifier record.

        Queries all offline sources by InChIKey.  Falls back to 14-character
        skeleton matching when ``self.inchikey_skeleton`` is True and no exact
        match is found.  The skeleton is the connectivity block, so it finds
        the compound regardless of stereochemistry, isotopes or charge.

        Args:
            inchikey: Full 27-character InChIKey
                (``XXXXXXXXXXXXXX-XXXXXXXXXX-X``).

        Returns:
            Tuple of (base result template, candidate pool, ``None``).
        """
        result = self._empty_result(inchikey, "InChIKey")
        result["match_method"] = "exact_inchikey"
        result["InChIKey"] = inchikey

        hits = self._collect("inchikey", inchikey)
        match_method = "exact_inchikey"

        if not _any_candidate(hits) and self.inchikey_skeleton:
            skeleton_hits = self._collect("inchikey_skeleton", inchikey)
            if _any_candidate(skeleton_hits):
                hits, match_method = skeleton_hits, "inchikey_skeleton"

        return result, self._pool(hits, match_method), None

    # ── DTXSID resolver ───────────────────────────────────────────────────────

    def _resolve_dtxsid(self, dtxsid: str) -> Tuple[Dict[str, Any], List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Resolve a CompTox DTXSID into a unified identifier record.

        CompTox is asked for the DTXSID; the other sources do not hold it.
        They are asked instead for the CAS number CompTox gives the
        substance, so that each one that holds that number answers from its
        own record: an independent vote on the structure, as for a CAS
        query.  The sources that do not hold the number, and ChEMBL, are
        then asked for each structure found, and what they find is an echo,
        not a vote (see `_cross_references`).

        Asking them for CompTox's InChIKey instead made every source agree
        with CompTox by construction: the wrong thiadiazole CompTox records
        for cefquinome sulfate (``DTXSID6048593``) came back with
        ``n_source_support = 2`` and no conflict, where its CAS number finds
        PubChem's aminothiazole.  A substance whose CAS number no other
        source holds is supported by CompTox alone, and ``strict`` refuses
        it.

        CompTox's candidate leads the pool, so a one-against-one split goes
        to the record the DTXSID names; the row is flagged ``conflict``.

        Args:
            dtxsid: CompTox DTXSID string (e.g., ``"DTXSID7020182"``).

        Returns:
            Tuple of (base result template, candidate pool, ``None``).
        """
        result = self._empty_result(dtxsid, "DTXSID")
        result["match_method"] = "dtxsid"
        result["DTXSID"] = dtxsid

        owner = self._collect("dtxsid", dtxsid, sources=["comptox"])
        comptox = owner.get("comptox") or []
        if not comptox:
            return result, [], None

        others = [key for key in self._SOURCE_KEYS if key not in ("comptox", "chembl")]
        # CompTox's own CASRN comes first among its numbers.
        casrn = next(iter(comptox[0].get("CAS_candidates") or []), None)
        votes: Hits = {}
        if not is_missing(casrn):
            votes = self._collect("cas", str(casrn), label=dtxsid, sources=others)
        silent = [key for key in others if not votes.get(key)] + ["chembl"]
        echoes = self._cross_references({**votes, **owner}, silent, label=dtxsid)

        pool = self._pool(owner, "dtxsid") + self._pool({**votes, **echoes}, "dtxsid")
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

        # Completeness drives the query_match_score for formula matches, which
        # have no name to compare against.
        hits = self._collect("formula", formula, k=self.top_k_per_source)
        return result, self._pool(hits, "formula", self._completeness_score), None

    # ── Fuzzy name search ─────────────────────────────────────────────────────

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
    ) -> Tuple[Hits, Optional[float]]:
        """Find structurally similar compounds using Tanimoto similarity.

        Computes a Morgan fingerprint for ``query_smiles`` and queries each
        source with its similarity search capabilities.  Returns candidates
        that meet ``self.similarity_threshold``.

        Args:
            query_smiles: Query SMILES string.

        Returns:
            Tuple of:

            - Source key -> candidates (the best match per source at or
              above threshold).
            - Best Tanimoto score observed, or ``None`` if RDKit is unavailable.

        Note:
            This is an initial implementation that uses per-source lookup; a
            future Parquet + vectorised fingerprint approach will be faster for
            large datasets.
        """
        candidates: Hits = {}

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
            if is_missing(smiles) or Chem is None:
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
        chembl = self._clients.get("chembl")
        if chembl is not None:
            try:
                row = chembl.search_by_smiles(query_smiles)
                if row:
                    t = _tanimoto_from_smiles(row.get("canonical_smiles"))
                    if t >= self.similarity_threshold:
                        candidates["chembl"] = [candidate_from_chembl_row(row, chembl)]
                        best_tanimoto = max(best_tanimoto, t)
            except Exception as exc:
                log.warning("ChEMBL Tanimoto search failed: %s", exc)

        # PubChemID — try canonical SMILES lookup as a proxy
        pubchem = self._clients.get("pubchem")
        if pubchem is not None:
            try:
                norm = normalize_structure(query_smiles)
                if not is_missing(norm["canonical_smiles"]):
                    row = pubchem.get_by_smiles(norm["canonical_smiles"])
                    if row:
                        t = _tanimoto_from_smiles(row.get("smiles") or row.get("canonical_smiles"))
                        if t >= self.similarity_threshold:
                            candidates["pubchem"] = [candidate_from_pubchem_row(row)]
                            best_tanimoto = max(best_tanimoto, t)
            except Exception as exc:
                log.warning("PubChemID Tanimoto search failed: %s", exc)

        return candidates, best_tanimoto if best_tanimoto > 0 else None

    # ── Source details ────────────────────────────────────────────────────────

    def _build_source_details(
        self, candidates: Dict[str, Optional[Dict[str, Any]]]
    ) -> Dict[str, Dict[str, Any]]:
        """Build a per-source traceability record from the candidates dict.

        For each source, records whether it was found and which output fields
        it has non-null values for.  The online services are listed only when
        ``online_fallback=True``.

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
        for key in self._SOURCE_KEYS + self._ONLINE_KEYS:
            display = self._SOURCE_DISPLAY[key]
            cand = candidates.get(key)
            if cand is None:
                details[display] = {"found": False, "fields": []}
            else:
                fields: List[str] = []
                for cand_field, out_field in _FIELD_MAP.items():
                    val = cand.get(cand_field)
                    if not is_missing(val):
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
        n_source_support: int = 0,
    ) -> float:
        """Compute the final confidence score for a result.

        The base confidence depends on the match method.  For fuzzy and
        Tanimoto methods, the raw similarity is used as the base.  The base is
        modulated by a query-agreement term (weighted by ``self.query_weight``),
        by the cross-source consensus score, and by how many independent
        databases carry the structure, so that a strong query match, agreement
        between sources, and corroboration all raise confidence.

        Formula::

            final = base
                  × (w_q × query_score + (1 − w_q))
                  × (0.5 + 0.5 × consensus_score)
                  × support_factor(n_source_support)

        For exact-identifier methods ``query_score`` is 1.0, which collapses the
        middle term to 1.0.

        The ``support_factor`` (`_SUPPORT_FACTOR`) is what keeps an
        uncorroborated hit from winning on provenance alone.  ``consensus_score``
        measures *how well* the sources that answered agree, not *how many*
        answered, and a lone source agrees with itself perfectly — so before this
        factor existed a single ChEBI row (0.90) outranked a structure that
        CompTox, PubChem and ZeroPM all agreed on (0.8777) and the resolver
        returned the wrong compound.

        A ``consensus_score`` of exactly 0.0 short-circuits to 0.0 rather than
        following the formula.
        [`compute_consensus`][provesid.tools.compute_consensus] only returns
        0.0 when there were no candidates at all — one source scores 1.0, and
        even two fully disagreeing sources score 0.5 — so a zero consensus
        means nothing matched, and the formula's floor of ``0.5 × base`` would
        report a no-match row as half-confident.

        Args:
            match_method: One of the keys in `_BASE_CONFIDENCE`.
            consensus_score: Cross-source consensus agreement in [0, 1].
            fuzzy_score: rapidfuzz similarity in [0, 1]; used when
                ``match_method == "fuzzy_name"``, scaled by the ``exact_name``
                base so a fuzzy match never outranks an exact one.
            tanimoto: Tanimoto similarity in [0, 1]; used when
                ``match_method == "tanimoto"``.
            query_score: Query-agreement signal in [0, 1] for name/formula
                methods.  Ignored (treated as 1.0) for fuzzy/Tanimoto where the
                similarity already lives in the base.
            n_source_support: Number of independent databases carrying this
                structure.  ``0`` means the cluster came from OPSIN alone and is
                left unpenalised.

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
        support_term = _SUPPORT_FACTOR.get(max(0, int(n_source_support)), _SUPPORT_FACTOR_MAX)
        modulated = (
            base
            * query_term
            * (0.5 + 0.5 * max(0.0, min(1.0, consensus_score)))
            * support_term
        )
        return round(min(1.0, max(0.0, modulated)), 4)

    # ── Result finalisation ───────────────────────────────────────────────────

    def _finalise_hits(
        self,
        base_template: Dict[str, Any],
        pool: List[Dict[str, Any]],
        n_hits: Union[int, str],
        min_confidence: float,
        opsin_anchor: Optional[Dict[str, Any]] = None,
        min_source_support: int = 0,
        reference: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Cluster a candidate pool, rank the clusters, and return ranked hits.

        Args:
            base_template: Empty result template (carries query/foundby and any
                pre-populated query fields).
            pool: Flat list of tagged candidate records.
            n_hits: Number of hits to return (positive int or ``"all"``).
            min_confidence: Drop hits below this confidence before truncation.
            opsin_anchor: Optional OPSIN structure anchor for this query.
            min_source_support: Drop hits carried by fewer than this many
                databases before truncation.  ``0`` disables the filter.
            reference: The query's input row, against which ``validate``
                checks each hit before truncation.

        Returns:
            List of fully-populated result dicts ordered by descending
            confidence with ``hit_rank`` set.  Always contains at least one row
            (an empty/no-match row when nothing was found).
        """
        opsin_smiles = opsin_anchor.get("smiles") if opsin_anchor else None

        # Drop group records (SMILES with an attachment point): they are never
        # the substance a query denotes.  A query that is itself a group SMILES
        # is exempt, since there the group *is* what was asked for.
        if not _has_attachment_point(base_template.get("query")):
            pool = [c for c in pool if not _has_attachment_point(c.get("SMILES"))]

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
        _rank_hits(hits)
        _mark_conflicts(hits)
        _rank_hits(hits)

        filtered = [
            h
            for h in hits
            if h["confidence"] >= min_confidence
            and h["n_source_support"] >= min_source_support
            and not (h["conflict"] and self.on_conflict == "refuse")
        ]
        if self.validate is not None:
            for h in filtered:
                checks = _validation_checks(h, reference or {}, self.validate)
                h["validated_by"] = None if checks is None else "+".join(checks)
            rejected = [h for h in filtered if h["validated_by"] is None]
            if rejected:
                log.debug(
                    "%r: validate= rejected %s.", base_template.get("query"),
                    ", ".join(f"{h.get('name')} ({h.get('molecular_formula')})" for h in rejected),
                )
            filtered = [h for h in filtered if h["validated_by"] is not None]
        if not filtered:
            # Everything was below the floor, or disputed under "refuse" ---
            # represent the query with a single no-match row so it is not
            # silently dropped, and say whether a dispute was the reason.
            empty = self._build_result_for_cluster(base_template, {"members": []}, opsin_smiles)
            empty["hit_rank"] = 0
            empty["n_clusters"] = len(hits)
            empty["conflict"] = any(h["conflict"] for h in hits)
            if empty["conflict"] and self.on_conflict == "refuse":
                log.info(
                    "%r: independent sources name different compounds (%s); "
                    "on_conflict='refuse' returns none of them.",
                    base_template.get("query"),
                    "; ".join(f"{h.get('name')} [{', '.join(h['_independent_sources'])}]"
                              for h in hits if h["conflict"]),
                )
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
            _warn_on_structure_quality(hit)

        return filtered

    @staticmethod
    def _set_structure_quality(
        result: Dict[str, Any],
        norm: Dict[str, Any],
        applied: List[Dict[str, Any]],
        query_structure_invalid: bool,
    ) -> None:
        """Fill a result's QC columns and make its structure authoritative.

        With a structure, ``InChIKey``, ``InChI``, ``molecular_formula`` and
        ``molecular_mass`` are computed from it, so that every identifier in
        the row describes the same molecule; what the sources stated is kept
        in ``source_inchikey`` and ``source_formula`` and checked against the
        structure (``identity_consistent``).  A mismatch is logged when the
        hit is returned (see `_finalise_hits`).  Without one, the sources'
        values stay where they are.

        Args:
            result: The result being built, modified in place.
            norm: [`normalize_structure`][provesid.search.normalize_structure]
                of ``result["SMILES"]``.
            applied: The candidates that filled the result, in the order
                they were applied.
            query_structure_invalid: Whether the query itself carried a
                structure that does not parse.
        """
        source_inchikey = result.get("InChIKey")
        source_formula = result.get("molecular_formula")
        result["source_inchikey"] = source_inchikey
        result["source_formula"] = source_formula

        if norm["mol"] is None:
            if not is_parsable_smiles(result.get("SMILES")):
                result["SMILES"] = None
            unusable = query_structure_invalid or any(
                cand.get("structure_status") == "invalid" for cand in applied
            )
            result["structure_status"] = "invalid" if unusable else "none"
            result["identity_consistent"] = None
            return

        provider = next(
            (cand for cand in applied if cand.get("SMILES") == result["SMILES"]), None
        )
        result["structure_status"] = (provider or {}).get("structure_status") or "valid"
        result["InChI"] = pick_first(norm["inchi"], result.get("InChI"))
        result["InChIKey"] = pick_first(norm["inchikey"], source_inchikey)
        result["molecular_formula"] = pick_first(
            structure_formula(result["SMILES"]), source_formula
        )
        result["molecular_mass"] = pick_first(norm["mol_weight"], result.get("molecular_mass"))
        result["identity_consistent"] = identity_consistent(
            result["SMILES"], source_inchikey, source_formula
        )

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

        # Candidates derived from another source's answer (ChEMBL for a CAS
        # query) are left out of the vote and the support count.
        consensus_source, source_match_scores, match_score = compute_consensus(per_source)
        consensus_candidate = per_source.get(consensus_source) if consensus_source else None
        independent = [key for key, cand in per_source.items() if not cand.get("derived_from")]

        # A query structure no one can read (a SMILES or InChI query) is not
        # kept: it would hold the field against the sources' structures.
        query_structure_invalid = (
            not is_missing(result.get("SMILES")) and not is_parsable_smiles(result.get("SMILES"))
        )
        if query_structure_invalid:
            result["SMILES"] = None

        # ChEMBL, then the online services, fill only what the others left.
        fill_order = [k for k in self._SOURCE_KEYS if k != "chembl"] + ["chembl"] + self._ONLINE_KEYS
        applied: List[Dict[str, Any]] = []
        for source_key in fill_order:
            candidate = per_source.get(source_key)
            if candidate_compatible_with_consensus(
                candidate, consensus_candidate, self.consensus_compat_threshold
            ):
                apply_candidate_to_result(result, candidate)
                applied.append(candidate)
        result["CASRN"] = pick_first(result.get("CASRN"), pick_casrn(applied))

        # OPSIN supplies a structure even when no source row carried one.
        if opsin_match and is_missing(result.get("SMILES")) and not is_missing(opsin_smiles):
            result["SMILES"] = opsin_smiles

        # Structure normalisation
        norm = normalize_structure(result.get("SMILES"))
        result["canonical_smiles"] = norm["canonical_smiles"]
        result["kekulized_smiles"] = norm["kekulized_smiles"]
        self._set_structure_quality(result, norm, applied, query_structure_invalid)

        result["name"] = pick_first(result.get("name"), result.get("IUPAC_name"))
        result["IUPAC_name"] = pick_first(result.get("IUPAC_name"), result.get("name"))
        result["source"] = pick_first(
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
        result["n_source_support"] = len(independent)
        result["confidence"] = self._compute_confidence(
            cluster_method,
            match_score,
            fuzzy_score=fuzzy_score,
            tanimoto=tanimoto,
            query_score=cluster_query_score,
            n_source_support=result["n_source_support"],
        )
        if result["structure_status"] in ("none", "invalid"):
            result["confidence"] = round(result["confidence"] * _NO_STRUCTURE_FACTOR, 4)
        result["opsin_smiles"] = opsin_smiles

        # Salt stripping
        if self.strip_salts and not is_missing(result.get("SMILES")):
            parent = strip_salts(result["SMILES"], self.salt_smarts or None)
            canonical = result.get("canonical_smiles")
            if not is_missing(parent) and parent != canonical:
                result["parent_smiles"] = parent
                parent_norm = normalize_structure(parent)
                result["parent_inchikey"] = parent_norm["inchikey"]

        # Transient ranking metadata (dropped before DataFrame assembly).
        result["_independent_sources"] = [
            per_source[key].get("source") for key in independent
        ] + (["OPSIN"] if opsin_match else [])
        # Whether some independent source's (or OPSIN's) *best* answer is this
        # compound, rather than one of its lower-ranked alternatives.
        result["_someones_first_answer"] = any(
            m.get("_origin_rank", 0) == 0 and not m.get("derived_from") for m in members
        )
        result["_opsin_match"] = opsin_match
        result["_cluster_query_score"] = cluster_query_score
        result["_min_origin_rank"] = min(
            (m.get("_origin_rank", 0) for m in members), default=0
        )
        return result


# ─────────────────────────────────────────────────────────────────────────────
# Cascade resolution for experimental datasets
# ─────────────────────────────────────────────────────────────────────────────

_VALIDATE_KEYS = frozenset(["mw", "smiles", "formula", "name", "mw_tolerance"])


def _check_validate_spec(validate: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """Check a ``Search(validate=...)`` spec and fill in its defaults.

    Args:
        validate: The spec, or None.

    Returns:
        A copy with ``mw_tolerance`` set, or None.

    Raises:
        ValueError: If a key is unknown, no check is named, or the tolerance
            is negative.
    """
    if validate is None:
        return None
    unknown = set(validate) - _VALIDATE_KEYS
    if unknown:
        raise ValueError(
            f"validate keys must be among {sorted(_VALIDATE_KEYS)}; got {sorted(unknown)}"
        )
    spec = {"mw_tolerance": 0.5, **validate}
    if not set(spec) - {"mw_tolerance"}:
        raise ValueError("validate must name at least one of 'mw', 'smiles', 'formula', 'name'")
    spec["mw_tolerance"] = float(spec["mw_tolerance"])
    if spec["mw_tolerance"] < 0:
        raise ValueError("validate['mw_tolerance'] must not be negative")
    return spec


def _validation_checks(
    hit: Dict[str, Any], row: Dict[str, Any], spec: Dict[str, Any]
) -> Optional[List[str]]:
    """Check a hit against the reference values of its input row.

    The semantics of [`mw_within`][provesid.search.mw_within], generalised
    to the references ``Search(validate=...)`` names: molecular weight
    (given, or computed from a reference structure) and formula must agree;
    an identical structure and a matching name are reported but never
    required.  A check whose reference cell is empty is skipped.

    Args:
        hit: The result to check.
        row: The query's input row.
        spec: The checked ``validate`` spec.

    Returns:
        The names of the checks that passed, in the order ``mw``,
        ``smiles``, ``formula``, ``name`` (empty when the row had nothing to
        check against), or None when a required check failed.
    """
    tolerance = spec["mw_tolerance"]
    hit_mass = hit.get("molecular_mass") if hit.get("structure_status") not in ("none", "invalid") else None
    passed: List[str] = []

    reference_masses: List[float] = []
    if "mw" in spec and not is_missing(row.get(spec["mw"])):
        try:
            reference_masses.append(float(row[spec["mw"]]))
        except (TypeError, ValueError):
            return None
    reference_structure = None
    if "smiles" in spec and not is_missing(row.get(spec["smiles"])):
        reference_structure = normalize_structure(row[spec["smiles"]])
        if reference_structure["mol_weight"] is None:
            return None
        reference_masses.append(reference_structure["mol_weight"])
    if reference_masses:
        if hit_mass is None or any(abs(float(hit_mass) - m) > tolerance for m in reference_masses):
            return None
        passed.append("mw")
    if reference_structure is not None and reference_structure["canonical_smiles"] == hit.get("canonical_smiles"):
        passed.append("smiles")

    if "formula" in spec and not is_missing(row.get(spec["formula"])):
        wanted = formula_element_counts(row[spec["formula"]])
        if wanted is None or wanted != formula_element_counts(hit.get("molecular_formula")):
            return None
        passed.append("formula")

    if "name" in spec and not is_missing(row.get(spec["name"])):
        if _matches_name_exactly(str(row[spec["name"]]), hit):
            passed.append("name")
    return passed


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
    the name, so [`resolve_cascade`][provesid.search.resolve_cascade] can
    record how much evidence backed each row in its ``validated_by`` column.

    Args:
        tolerance: Maximum absolute difference in Da. Defaults to ``0.5``.
        reference_column: Column holding the row's own SMILES, used as the
            reference structure. Defaults to ``"SMILES"``.
        name_column: Optional column holding the row's own name. When given, a
            matching name is reported as an extra ``"name"`` check.

    Returns:
        A callable ``(hit, row) -> list[str]`` suitable for
        [`resolve_cascade`][provesid.search.resolve_cascade]'s ``accept``
        argument: the names of the checks that passed, or an empty list to
        reject the hit.

    Examples:
        >>> accept = mw_within(0.5, reference_column="SMILES", name_column="name")
        >>> accept({"SMILES": "CCO", "name": "ethanol"}, {"SMILES": "OCC", "name": "Ethanol"})
        ['mw', 'smiles', 'name']
        >>> accept({"SMILES": "CCCO"}, {"SMILES": "CCO"})       # 60.1 vs 46.07 Da
        []
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
            if not is_missing(wanted) and _matches_name_exactly(
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
    [`Search`][provesid.search.Search] instances in order, passing to each
    stage only the rows that are still unresolved, so every row is resolved by
    the most reliable identifier it actually has.

    Each hit is checked with ``accept`` before it counts as resolved. A hit that
    fails leaves its row pending for the next stage, which is what stops a
    confident-but-wrong match from ending the cascade. Use
    [`mw_within`][provesid.search.mw_within] for the usual molecular-weight
    check.

    Args:
        df: Input DataFrame. Returned unmodified; the result is a copy.
        stages: Ordered list of ``(label, search, column)`` triples. ``label``
            names the stage in the output, ``search`` is a
            [`Search`][provesid.search.Search] instance, and ``column`` is the
            column it reads. Rows with an empty value in ``column`` skip that
            stage.
        accept: Optional ``(hit, row) -> bool | list[str]`` predicate, where
            ``hit`` is the Search result row and ``row`` the input row, both as
            dicts. Return ``True``, or the names of the checks that passed (they
            are joined into ``validated_by``); return ``False`` or an empty list
            to reject. When ``None``, any hit carrying an InChIKey is accepted.

            Both dicts come from DataFrame rows, so a missing field is ``NaN``
            rather than ``None`` — and ``bool(NaN)`` is ``True``. Test emptiness
            with `pandas.isna` (or reuse
            [`mw_within`][provesid.search.mw_within]) rather than truthiness.
        fallback_column: Column holding a SMILES from which to derive identifiers
            for rows no stage resolved. Those rows get ``resolved_by="rdkit"``.
            When ``None``, unresolved rows are left empty.
        prefix: Prepended to every added column. Defaults to ``"provesid_"``.

    Returns:
        A copy of ``df`` with the
        [`OUTPUT_COLUMNS`][provesid.search.OUTPUT_COLUMNS] added under
        ``prefix``, plus ``<prefix>resolved_by`` (the stage that resolved the
        row, ``"rdkit"``, or ``"none"``) and ``<prefix>validated_by``.

    Raises:
        KeyError: If a stage names a column that is not in ``df``.
        ValueError: If ``stages`` is empty.

    Examples:
        >>> data = pd.DataFrame({
        ...     "CASRN":  ["50-78-2", "", "0-00-0"],
        ...     "name":   ["", "caffeine", ""],
        ...     "SMILES": ["CC(=O)Oc1ccccc1C(=O)O", "Cn1c(=O)c2c(ncn2C)n(C)c1=O", "CC(C)O"],
        ... })
        >>> out = resolve_cascade(
        ...     data,
        ...     stages=[
        ...         ("cas",  Search("cas", show_progress=False),  "CASRN"),
        ...         ("name", Search("name", show_progress=False), "name"),
        ...     ],
        ...     accept=mw_within(0.5, reference_column="SMILES"),
        ...     fallback_column="SMILES",
        ... )
        >>> out[["provesid_name", "provesid_resolved_by", "provesid_validated_by"]].values.tolist()
        [['acetylsalicylic acid', 'cas', 'mw+smiles'], ['caffeine', 'name', 'mw+smiles'], [nan, 'rdkit', 'self (rdkit from the given structure)']]

        The third row's CAS number is not real, so no stage resolved it and
        its identifiers were derived from its own SMILES.
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
            return "inchikey" if not is_missing(hit.get("InChIKey")) else None
        outcome = accept(hit, row)
        if isinstance(outcome, bool):
            return "accept" if outcome else None
        checks = list(outcome or [])
        return "+".join(checks) if checks else None

    for label, searcher, column in stages:
        if not pending:
            break

        eligible = [i for i in pending if not is_missing(rows.at[i, column])]
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

def _rank_hits(hits: List[Dict[str, Any]]) -> None:
    """Sort cluster results best first, in place.

    OPSIN match first, then confidence, then a structure over none, support,
    query agreement, and (lower) origin rank as a final tie-break.
    Corroboration is folded into ``confidence`` itself (see
    `_SUPPORT_FACTOR`), so ``n_source_support`` here only breaks ties between
    equally confident clusters.  The sort is stable, so clusters still tied
    keep pool order, which is source reputation order.

    Args:
        hits: Results from ``Search._build_result_for_cluster``.
    """
    hits.sort(
        key=lambda h: (
            1 if h["_opsin_match"] else 0,
            h["confidence"],
            0 if h["structure_status"] in ("none", "invalid") else 1,
            h["n_source_support"],
            h["_cluster_query_score"],
            -h["_min_origin_rank"],
        ),
        reverse=True,
    )


def _parent_skeleton(canonical_smiles: Optional[str]) -> Optional[str]:
    """The InChIKey skeleton of a structure's neutral, canonical-tautomer parent.

    What two clusters are compared on to decide whether they are different
    compounds: a salt, a solvate, a charge state, a stereoisomer or a
    tautomer of one compound has the same parent skeleton.  The standard
    InChI already merges mobile ring hydrogens, but not keto/enol or
    azo/hydrazone forms, and CompTox and PubChem often store different ones:
    of 500 random CompTox CAS numbers, 6 of the 8 that looked contested were
    one compound drawn as two tautomers.  RDKit's tautomer canonicaliser is
    slow on large molecules, so this is computed only for queries that
    returned more than one compound.

    Args:
        canonical_smiles: The structure, or None.

    Returns:
        The 14-character skeleton, or None without a readable structure.
    """
    if is_missing(canonical_smiles):
        return None
    parent_smiles = strip_salts(canonical_smiles) or canonical_smiles
    mol = Chem.MolFromSmiles(parent_smiles) if RDKIT_AVAILABLE else None
    if mol is None:
        return None
    try:
        with rdBase.BlockLogs():
            mol = rdMolStandardize.Uncharger().uncharge(mol)
            mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
            inchikey = Chem.MolToInchiKey(mol)
    except Exception:
        inchikey = normalize_structure(parent_smiles)["inchikey"]
    return inchikey[:14] if inchikey else None


def _disputes(hit: Dict[str, Any], other: Dict[str, Any]) -> bool:
    """Whether ``other`` contradicts ``hit`` as the compound the query names.

    It does when both were found by an exact identifier match, ``other`` is
    the best answer of an independent source (not a lower-ranked
    alternative, and not a record found through another source's answer)
    and matches the query at least as well as ``hit`` does (a name search's
    top row that is not called by the name contests nothing), and the two
    are different compounds: their parent skeletons differ (see
    `_parent_skeleton`; the InChIKey skeleton stands in when there is no
    structure).

    Args:
        hit: One cluster's result.
        other: Another cluster's result for the same query.

    Returns:
        True when the two are rival answers.
    """
    def skeleton(h: Dict[str, Any]) -> Optional[str]:
        if "_parent_skeleton" not in h:
            h["_parent_skeleton"] = _parent_skeleton(h.get("canonical_smiles"))
        inchikey = h.get("InChIKey")
        return h["_parent_skeleton"] or (str(inchikey)[:14] if not is_missing(inchikey) else None)

    mine, theirs = skeleton(hit), skeleton(other)
    return (
        hit.get("match_method") in _EXACT_METHODS
        and other.get("match_method") in _EXACT_METHODS
        and other["_someones_first_answer"]
        and other["_cluster_query_score"] >= hit["_cluster_query_score"]
        and bool(other["_independent_sources"])
        and mine is not None
        and theirs is not None
        and mine != theirs
    )


def _mark_conflicts(hits: List[Dict[str, Any]]) -> None:
    """Fill the conflict columns of ranked cluster results, in place.

    Each result's rival is the best-ranked other result that disputes it
    (see `_disputes`).  A result with a rival is flagged ``conflict`` and
    its confidence multiplied by its share of the independent sources
    behind the two, ``own / (own + rival's)``: a one-against-one split
    halves it, three against one keeps three quarters.  The penalty is the
    point --- the support factor punishes a lone source, but before this a
    *contested* answer read the same as an unopposed one, and cefquinome's
    CAS resolved to a wrong CompTox record at 0.765 with no sign that
    PubChem said otherwise.

    ``runner_up_*`` describe the rival, or the best other result when there
    is none, with the confidence it has after its own penalty.

    Args:
        hits: Results ranked by `_rank_hits`; re-rank after this.
    """
    rivals = [
        next((other for other in hits if other is not hit and _disputes(hit, other)), None)
        for hit in hits
    ]
    for hit, rival in zip(hits, rivals):
        hit["n_clusters"] = len(hits)
        hit["conflict"] = rival is not None
        if rival is not None:
            own = max(1, len(hit["_independent_sources"]))
            theirs = max(1, len(rival["_independent_sources"]))
            hit["confidence"] = round(hit["confidence"] * own / (own + theirs), 4)
    for hit, rival in zip(hits, rivals):
        shown = rival or next((other for other in hits if other is not hit), None)
        if shown is not None:
            hit["runner_up_inchikey"] = shown.get("InChIKey")
            hit["runner_up_sources"] = list(shown["_independent_sources"])
            hit["runner_up_confidence"] = shown["confidence"]


def _warn_on_structure_quality(hit: Dict[str, Any]) -> None:
    """Log a returned hit whose structure cannot be used as it is.

    Logged at WARNING, because a row that looks resolved --- a name, an
    InChIKey, a confidence --- but has no readable structure, or a structure
    its stated InChIKey does not belong to, is exactly what a downstream
    pipeline takes at face value.

    Args:
        hit: A finished result dict.
    """
    if hit.get("structure_status") == "invalid":
        log.warning(
            "%r (%s): the structure %s gave cannot be read or repaired; "
            "the row has no structure.",
            hit.get("query"), hit.get("name"), hit.get("consensus_source"),
        )
    if hit.get("identity_consistent") is False:
        log.warning(
            "%r (%s): the stated InChIKey %s / formula %s do not belong to "
            "the structure %s (InChIKey %s).",
            hit.get("query"), hit.get("name"), hit.get("source_inchikey"),
            hit.get("source_formula"), hit.get("SMILES"), hit.get("InChIKey"),
        )


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
        if not is_missing(value):
            names.append(str(value))
    synonyms = cand.get("Synonyms")
    if not is_missing(synonyms):
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


def _has_attachment_point(smiles: Any) -> bool:
    """Whether a SMILES describes a *group* rather than a whole compound.

    An asterisk in SMILES is a dummy atom — an open valence where the group
    attaches to something else.  ChEBI indexes such groups (for example
    ``*C(=O)CCCC=CCC=CCCCCC``, "cis,cis-tetradeca-5,8-dienoyl group") alongside
    real compounds, and one of them being returned as the structure for a CAS
    lookup is always wrong: a registry number identifies a substance, never a
    substituent.  RDKit also cannot compute descriptors for these
    (``Unsupported in this mode element '*'``).

    Args:
        smiles: A SMILES string, or any missing-like value.

    Returns:
        True when the SMILES contains an attachment point.
    """
    return not is_missing(smiles) and "*" in str(smiles)


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
    if not is_missing(ik):
        ik_s = str(ik)
        keys.add(("ik", ik_s))
        if by_skeleton and len(ik_s) >= 14:
            keys.add(("skel", ik_s[:14]))
        return keys

    canon = cand.get("canonical_smiles") or cand.get("SMILES")
    if not is_missing(canon):
        keys.add(("smi", str(canon)))
        return keys

    name = cand.get("name") or cand.get("IUPAC_name")
    if not is_missing(name):
        keys.add(("name", str(name).strip().lower()))
    return keys


def _cluster_candidates(
    pool: List[Dict[str, Any]], by_skeleton: bool = True
) -> List[Dict[str, Any]]:
    """Group a candidate pool into structure-identity clusters.

    Uses a union-find over candidate identity keys so that, e.g., a candidate
    carrying a full InChIKey and one carrying only the matching skeleton end up
    in the same cluster.  Candidates with no usable identity become singletons.

    A record with no structure --- CompTox lists eprinomectin, gentamicin and
    other mixtures by name and DTXSID alone --- has no structure key to share,
    so it would form a cluster of its own and compete with the structured
    record of the same substance instead of supporting it.  Such a cluster is
    therefore merged into a structured one when exactly one structured
    cluster shares its DTXSID, else its name, else a CAS number (see
    `_structured_cluster_for`).  When several do, it is left alone.

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

    structured = [m for m in groups.values() if any(_has_structure_identity(c) for c in m)]
    clusters = list(structured)
    for members in groups.values():
        if any(_has_structure_identity(c) for c in members):
            continue
        target = _structured_cluster_for(members, structured)
        if target is None:
            clusters.append(members)
        else:
            target.extend(members)

    # Keep the clusters in the order their first member was pooled.
    first_seen = {id(cand): idx for idx, cand in enumerate(pool)}
    clusters.sort(key=lambda members: min(first_seen[id(c)] for c in members))
    return [{"members": members} for members in clusters]


def _has_structure_identity(cand: Dict[str, Any]) -> bool:
    """Whether a candidate carries a structure or an InChIKey to cluster on."""
    return not is_missing(cand.get("InChIKey")) or not is_missing(cand.get("canonical_smiles"))


def _structured_cluster_for(
    members: List[Dict[str, Any]], structured: List[List[Dict[str, Any]]]
) -> Optional[List[Dict[str, Any]]]:
    """The structured cluster a structure-less one describes, if exactly one.

    Tried in order, the first test that any structured cluster passes
    decides: the same DTXSID; a preferred name of one side among the names
    (synonyms included) of the other; a CAS number in common.  The decision
    is taken only when exactly one cluster passes it, so a structure-less
    record is never attached to one of two competing compounds.

    Args:
        members: The candidates of a cluster without structure.
        structured: The member lists of the clusters with one.

    Returns:
        The member list to merge into, or None.
    """
    def normalise(text: Any) -> str:
        return re.sub(r"\s+", " ", str(text).strip().lower())

    def preferred_names(cands: List[Dict[str, Any]]) -> set:
        return {
            normalise(c.get(field)) for c in cands for field in ("name", "IUPAC_name")
            if not is_missing(c.get(field))
        }

    def all_names(cands: List[Dict[str, Any]]) -> set:
        return {normalise(name) for c in cands for name in _candidate_names(c)}

    def values(cands: List[Dict[str, Any]], field: str) -> set:
        return {str(c.get(field)) for c in cands if not is_missing(c.get(field))}

    def cas_numbers(cands: List[Dict[str, Any]]) -> set:
        return {cas for c in cands for cas in (c.get("CAS_candidates") or [])}

    tests = (
        lambda other: bool(values(members, "DTXSID") & values(other, "DTXSID")),
        lambda other: bool(
            preferred_names(members) & all_names(other)
            or preferred_names(other) & all_names(members)
        ),
        lambda other: bool(cas_numbers(members) & cas_numbers(other)),
    )
    for passes in tests:
        matches = [other for other in structured if passes(other)]
        if len(matches) == 1:
            return matches[0]
        if matches:
            return None
    return None


def _any_candidate(hits: Hits) -> bool:
    """Return True if any source found at least one candidate.

    Args:
        hits: Source key -> candidates, as ``Search._collect`` returns.

    Returns:
        True when at least one source's list is non-empty.
    """
    return any(hits.values())
