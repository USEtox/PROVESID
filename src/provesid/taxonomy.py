"""Chemical taxonomy classification backends for PROVESID.

This module provides high-level "give me chemical-class labels for these
structures" capabilities on top of PROVESID's identifier resolution. The first
backend implemented here is :class:`ChebifierClassifier`, a wrapper around the
offline, AI-based ChEB-AI ``chebifier`` ensemble that assigns ChEBI ontology
classes to molecules.

``chebifier`` is a heavy, optional dependency (it pulls in PyTorch and, for the
graph models, the PyG stack). It is therefore **not** a core requirement of
PROVESID; install it with the ``chebifier`` extra plus the helper script::

    bash scripts/install_chebifier.sh

See ``docs/chebifier.md`` for the full installation story and known issues, and
``plans/2026-07-02-chemical-taxonomy-classyfire-chebifier.md`` (§10) for the
design rationale.

Key design points (mirrors the PROVESID conventions):

* **Optional + lazily imported.** Importing :mod:`provesid.taxonomy` never
  requires PyTorch. ``chebifier`` is imported only when a
  :class:`ChebifierClassifier` actually needs it, and a missing install raises a
  clear :class:`ChebifierMissingError` rather than a raw ``ModuleNotFoundError``.
* **Systemwide model storage.** Model weights are redirected to the shared
  per-user PROVESID dataset directory (:func:`provesid.utils.user_dataset_path`,
  overridable with ``PROVESID_DATA_DIR``) so a single copy is reused across all
  virtual environments on the machine, exactly like the other large PROVESID
  datasets.
* **InChIKey-keyed, resumable cache.** Each structure is classified once and
  cached on disk (keyed by InChIKey + chebifier version + configuration), so a
  re-run over the same chemicals hits the cache and never reloads the model.
* **Self-healing checkpoint compatibility.** chebifier 1.2.1's graph checkpoints
  require chebai-graph's property index vocabularies in a specific (older) state.
  :func:`ensure_v244_indices` restores them if a drifted version is installed.
"""

from __future__ import annotations

import importlib
import importlib.util
import logging
import os
from typing import Any, Dict, List, Optional, Sequence, Union

import pandas as pd

from .cache import _service_caches
from .search import normalize_structure
from .utils import user_dataset_path

logger = logging.getLogger(__name__)

#: chebifier release this backend is validated against (see docs/chebifier.md).
CHEBIFIER_PINNED_VERSION = "1.2.1"

#: Columns of the tidy taxonomy table returned by ``classify`` (shared with the
#: planned ClassyFire backend, hence the ClassyFire-only level columns).
TAXONOMY_COLUMNS = [
    "inchikey",
    "smiles",
    "kingdom",
    "superclass",
    "class",
    "subclass",
    "chebi_ids",
    "chebi_names",
    "source",
    "confidence",
]

# ---------------------------------------------------------------------------
# chebai-graph property-index compatibility (v244 checkpoints)
# ---------------------------------------------------------------------------
# chebifier 1.2.1 ships GNN checkpoints ("v244") whose edge/node feature widths
# are fixed by chebai-graph's one-hot property vocabularies. A chebai-graph
# commit (ea77f36, 2026-03-02, *after* chebifier 1.2.1 shipped) appended tokens
# to three of them, widening the graph feature vectors and breaking checkpoint
# loading with "mat1 and mat2 shapes cannot be multiplied". These are the
# pre-drift (commit 677d44b) contents that match the v244 checkpoints.
_V244_PROPERTY_INDICES: Dict[str, List[str]] = {
    "BondType": ["DATIVE", "SINGLE", "AROMATIC", "TRIPLE", "DOUBLE"],
    "AtomNumHs": ["0", "3", "2", "4", "1", "5", "6"],
    "NumAtomBonds": ["0", "1", "2", "4", "5", "3", "6", "8", "7", "10", "12"],
}


class ChebifierError(Exception):
    """Base class for errors raised by the chebifier taxonomy backend."""


class ChebifierMissingError(ChebifierError):
    """Raised when the optional ``chebifier`` dependency is not installed."""


def chebifier_available() -> bool:
    """Return whether the optional ``chebifier`` package is importable.

    This is a cheap check (it does not import PyTorch or load any model) suitable
    for feature-detection and for skipping tests when the extra is absent.

    Returns:
        ``True`` if ``chebifier`` can be imported, ``False`` otherwise.

    Example:
        >>> from provesid.taxonomy import chebifier_available
        >>> if chebifier_available():
        ...     ...  # safe to construct a ChebifierClassifier
    """
    return importlib.util.find_spec("chebifier") is not None


def _configure_chebifier_storage(data_dir: Optional[str] = None) -> str:
    """Redirect Hugging Face / torch model caches into the PROVESID dataset dir.

    chebifier downloads its model weights from the Hugging Face Hub on first use.
    By default those land in ``~/.cache/huggingface``; this points them at the
    shared per-user PROVESID dataset directory instead, so a single copy is
    reused across virtual environments. Environment variables already set by the
    user are respected (never overridden).

    Must be called **before** ``chebifier``/``torch``/``huggingface_hub`` are
    imported, since those libraries read these variables at import time.

    Args:
        data_dir: Explicit base directory for chebifier data. When ``None``,
            uses ``user_dataset_path("chebifier")`` (honors ``PROVESID_DATA_DIR``).

    Returns:
        The base directory used for chebifier storage.
    """
    base = data_dir or user_dataset_path("chebifier")
    os.environ.setdefault("HF_HOME", os.path.join(base, "huggingface"))
    os.environ.setdefault("HF_HUB_CACHE", os.path.join(base, "huggingface", "hub"))
    os.environ.setdefault("TORCH_HOME", os.path.join(base, "torch"))
    return base


def ensure_v244_indices() -> Dict[str, str]:
    """Ensure chebai-graph's property indices match chebifier's v244 checkpoints.

    Reverts the ``BondType``/``AtomNumHs``/``NumAtomBonds`` one-hot vocabularies
    inside the installed ``chebai_graph`` package to their pre-drift contents when
    a newer (drifted) version is present. Without this, chebifier 1.2.1's graph
    (GNN) models fail to load with a tensor-shape error. The operation is
    idempotent and a no-op when the indices already match or when ``chebai_graph``
    is not installed. See ``docs/chebifier.md`` for the full root-cause writeup.

    Returns:
        Mapping of property name to a status string: ``"patched"``, ``"ok"``
        (already matching), or ``"missing"`` (index file not found).
    """
    try:
        import chebai_graph  # noqa: F401 - only need its location
    except ImportError:
        logger.debug("chebai_graph not installed; skipping v244 index check")
        return {}

    bin_dir = os.path.join(
        os.path.dirname(chebai_graph.__file__), "preprocessing", "bin"
    )
    results: Dict[str, str] = {}
    for prop, tokens in _V244_PROPERTY_INDICES.items():
        path = os.path.join(bin_dir, prop, "indices_one_hot.txt")
        if not os.path.exists(path):
            logger.warning("chebai-graph index file not found: %s", path)
            results[prop] = "missing"
            continue
        with open(path, "r") as handle:
            current = [line.strip() for line in handle if line.strip()]
        if current == tokens:
            results[prop] = "ok"
            continue
        with open(path, "w") as handle:
            handle.write("\n".join(tokens) + "\n")
        logger.info(
            "Patched chebai-graph %s index to v244 state (%d -> %d tokens)",
            prop,
            len(current),
            len(tokens),
        )
        results[prop] = "patched"
    return results


def _load_chebifier():
    """Import and return the ``chebifier.BaseEnsemble`` class.

    Raises:
        ChebifierMissingError: If the optional ``chebifier`` extra is not
            installed, with instructions for installing it.
    """
    try:
        from chebifier import BaseEnsemble
    except ImportError as exc:  # pragma: no cover - exercised only without extra
        raise ChebifierMissingError(
            "The chebifier taxonomy backend requires the optional 'chebifier' "
            "extra and its model dependencies. Install with:\n"
            "    bash scripts/install_chebifier.sh\n"
            "or, for the transformer/rule-based models only:\n"
            "    pip install 'provesid[chebifier]'"
        ) from exc
    return BaseEnsemble


class ChebifierClassifier:
    """Classify molecules into ChEBI ontology classes with the chebifier ensemble.

    Wraps the ChEB-AI ``chebifier`` ensemble behind PROVESID conventions:
    systemwide model storage, an InChIKey-keyed resumable cache, and a tidy
    :class:`pandas.DataFrame` output. The (expensive) ensemble is constructed
    lazily on first use and reused for the lifetime of the instance.

    Args:
        data_dir: Base directory for chebifier model storage. Defaults to
            ``user_dataset_path("chebifier")`` (honors ``PROVESID_DATA_DIR``).
        use_cache: When ``True`` (default), classified structures are cached on
            disk by InChIKey and reused on subsequent calls.
        resolve_names: When ``True``, resolve predicted ChEBI IDs to names via the
            online :class:`provesid.ChEBI` client (cached). Defaults to ``False``
            to keep classification fast and dependency-light.
        model_configs: Optional custom chebifier ensemble configuration (path or
            dict) passed straight to ``BaseEnsemble``. When ``None``, chebifier's
            default ensemble is used.
        patch_indices: When ``True`` (default), call :func:`ensure_v244_indices`
            before loading the ensemble so the graph models load correctly.

    Raises:
        ChebifierMissingError: If ``chebifier`` is not installed (raised when the
            ensemble is first constructed).

    Example:
        >>> from provesid.taxonomy import ChebifierClassifier
        >>> clf = ChebifierClassifier()
        >>> df = clf.classify(["c1ccccc1", "OCC1OC(O)C(O)C(O)C1O"])
        >>> df[["smiles", "chebi_ids"]]  # doctest: +SKIP
    """

    _CACHE_FUNC_NAME = "provesid.taxonomy.ChebifierClassifier"

    def __init__(
        self,
        data_dir: Optional[str] = None,
        use_cache: bool = True,
        resolve_names: bool = False,
        model_configs: Optional[Union[str, Dict[str, Any]]] = None,
        patch_indices: bool = True,
    ) -> None:
        self.data_dir = _configure_chebifier_storage(data_dir)
        self.use_cache = use_cache
        self.resolve_names = resolve_names
        self.model_configs = model_configs
        self.patch_indices = patch_indices
        self._ensemble = None
        self._chebi = None
        self._cache = _service_caches["chebifier"]

    # -- lazy resources ----------------------------------------------------
    @property
    def chebifier_version(self) -> str:
        """Installed chebifier version (falls back to the pinned version)."""
        try:
            return importlib.import_module("chebifier").__version__
        except Exception:
            return CHEBIFIER_PINNED_VERSION

    @property
    def ensemble(self):
        """The lazily-constructed, reused ``BaseEnsemble`` instance."""
        if self._ensemble is None:
            base_ensemble_cls = _load_chebifier()
            if self.patch_indices:
                ensure_v244_indices()
            logger.info(
                "Loading chebifier ensemble (weights cache: %s). First run "
                "downloads model weights.",
                os.environ.get("HF_HOME", "<default>"),
            )
            if self.model_configs is not None:
                self._ensemble = base_ensemble_cls(model_configs=self.model_configs)
            else:
                self._ensemble = base_ensemble_cls()
        return self._ensemble

    def _get_chebi(self):
        """Lazily construct an online ChEBI client for name resolution."""
        if self._chebi is None:
            from .chebi import ChEBI

            self._chebi = ChEBI()
        return self._chebi

    # -- caching helpers ---------------------------------------------------
    def _cache_key(self, inchikey: str) -> str:
        """Build the cache key for a structure (InChIKey + version + config)."""
        config_sig = "default" if self.model_configs is None else "custom"
        return f"{inchikey}::{self.chebifier_version}::{config_sig}"

    @staticmethod
    def _normalize_prediction(prediction: Any) -> Dict[str, Optional[float]]:
        """Normalise one chebifier prediction into ``{chebi_id: confidence}``.

        chebifier's ``predict_smiles_list`` returns, per molecule, either a list
        of predicted ChEBI id strings or a mapping of id to score (depending on
        version/config); ``None`` for structures it could not classify. This
        collapses all shapes into an ordered ``{id: confidence}`` dict (with
        ``None`` confidence when the backend does not provide one).
        """
        if prediction is None:
            return {}
        if isinstance(prediction, dict):
            return {str(k): (float(v) if v is not None else None)
                    for k, v in prediction.items()}
        # list/tuple/set of ids
        return {str(cid): None for cid in prediction}

    # -- core API ----------------------------------------------------------
    def classify(
        self,
        smiles: Union[str, Sequence[str]],
        inchikeys: Optional[Sequence[Optional[str]]] = None,
    ) -> pd.DataFrame:
        """Classify one or more structures into ChEBI ontology classes.

        Args:
            smiles: A SMILES string or a sequence of SMILES strings.
            inchikeys: Optional pre-computed InChIKeys aligned with ``smiles``.
                When omitted, InChIKeys are derived from the SMILES with RDKit.
                The InChIKey is the cache key, so supplying canonical values makes
                caching consistent across equivalent SMILES.

        Returns:
            A :class:`pandas.DataFrame` with one row per input structure and the
            columns in :data:`TAXONOMY_COLUMNS`. For this backend the ClassyFire
            level columns (``kingdom``/``superclass``/``class``/``subclass``) are
            ``None``; ``chebi_ids`` holds the ``|``-joined predicted ChEBI ids,
            ``chebi_names`` the ``|``-joined names when ``resolve_names`` is set,
            ``source`` is ``"chebifier"``, and ``confidence`` the ``|``-joined
            per-label scores when available.

        Raises:
            ChebifierMissingError: If ``chebifier`` is not installed.
        """
        if isinstance(smiles, str):
            smiles_list: List[str] = [smiles]
        else:
            smiles_list = list(smiles)

        if inchikeys is not None and len(inchikeys) != len(smiles_list):
            raise ValueError(
                "inchikeys must be the same length as smiles "
                f"({len(inchikeys)} != {len(smiles_list)})"
            )

        # Resolve InChIKeys (cache keys) for every input.
        resolved_keys: List[Optional[str]] = []
        for idx, smi in enumerate(smiles_list):
            key = inchikeys[idx] if inchikeys is not None else None
            if not key:
                key = normalize_structure(smi).get("inchikey")
            resolved_keys.append(key)

        # Look up cache; collect the structures that still need the model.
        predictions: List[Optional[Dict[str, Optional[float]]]] = [None] * len(smiles_list)
        to_run_idx: List[int] = []
        for idx, key in enumerate(resolved_keys):
            if self.use_cache and key:
                found, value = self._cache.get(
                    self._CACHE_FUNC_NAME, (self._cache_key(key),), {}
                )
                if found:
                    predictions[idx] = value
                    continue
            to_run_idx.append(idx)

        # Run the ensemble on the cache misses (single batched call).
        if to_run_idx:
            batch = [smiles_list[i] for i in to_run_idx]
            raw = self.ensemble.predict_smiles_list(batch)
            for pos, idx in enumerate(to_run_idx):
                norm = self._normalize_prediction(raw[pos])
                predictions[idx] = norm
                key = resolved_keys[idx]
                if self.use_cache and key:
                    self._cache.set(
                        self._CACHE_FUNC_NAME, (self._cache_key(key),), {}, norm
                    )

        rows = [
            self._build_row(smiles_list[i], resolved_keys[i], predictions[i] or {})
            for i in range(len(smiles_list))
        ]
        return pd.DataFrame(rows, columns=TAXONOMY_COLUMNS)

    def _build_row(
        self,
        smiles: str,
        inchikey: Optional[str],
        prediction: Dict[str, Optional[float]],
    ) -> Dict[str, Any]:
        """Assemble one tidy-schema row from a normalised prediction."""
        chebi_ids = list(prediction.keys())
        confidences = [prediction[cid] for cid in chebi_ids]
        names = self._resolve_names(chebi_ids) if self.resolve_names else None
        has_conf = any(c is not None for c in confidences)
        return {
            "inchikey": inchikey,
            "smiles": smiles,
            "kingdom": None,
            "superclass": None,
            "class": None,
            "subclass": None,
            "chebi_ids": "|".join(chebi_ids) if chebi_ids else None,
            "chebi_names": "|".join(names) if names else None,
            "source": "chebifier",
            "confidence": (
                "|".join("" if c is None else f"{c:g}" for c in confidences)
                if has_conf
                else None
            ),
        }

    def _resolve_names(self, chebi_ids: Sequence[str]) -> List[str]:
        """Resolve ChEBI ids to names via the ChEBI client (best-effort)."""
        names: List[str] = []
        for cid in chebi_ids:
            names.append(self._resolve_single_name(cid))
        return names

    def _resolve_single_name(self, chebi_id: str) -> str:
        """Resolve one ChEBI id to a name, caching and degrading gracefully."""
        cid = chebi_id if str(chebi_id).upper().startswith("CHEBI:") else f"CHEBI:{chebi_id}"
        found, value = self._cache.get("provesid.taxonomy.chebi_name", (cid,), {})
        if found:
            return value
        name = cid
        try:
            entity = self._get_chebi().get_compound(cid)
            if isinstance(entity, dict):
                name = entity.get("chebiAsciiName") or entity.get("name") or cid
        except Exception as exc:  # network / lookup failures must not break classify
            logger.debug("ChEBI name lookup failed for %s: %s", cid, exc)
        self._cache.set("provesid.taxonomy.chebi_name", (cid,), {}, name)
        return name

    @staticmethod
    def to_labels(df: pd.DataFrame, level: str = "chebi_ids") -> Dict[str, Optional[str]]:
        """Collapse a taxonomy table to a ``{inchikey: label}`` mapping.

        Args:
            df: A taxonomy table as returned by :meth:`classify`.
            level: The column to use as the label (e.g. ``"chebi_ids"`` or
                ``"chebi_names"``).

        Returns:
            Mapping from InChIKey to the chosen label (``None`` when absent).

        Raises:
            KeyError: If ``level`` is not a column of ``df``.
        """
        if level not in df.columns:
            raise KeyError(f"Unknown level {level!r}; available: {list(df.columns)}")
        return dict(zip(df["inchikey"], df[level]))


def classify_chebifier(
    smiles: Union[str, Sequence[str]],
    inchikeys: Optional[Sequence[Optional[str]]] = None,
    *,
    data_dir: Optional[str] = None,
    use_cache: bool = True,
    resolve_names: bool = False,
) -> pd.DataFrame:
    """Classify structures with the chebifier ensemble (convenience wrapper).

    Thin functional wrapper over :class:`ChebifierClassifier` for one-off calls.
    For repeated calls, construct a :class:`ChebifierClassifier` once and reuse it
    so the model is loaded a single time.

    Args:
        smiles: A SMILES string or sequence of SMILES strings.
        inchikeys: Optional pre-computed InChIKeys aligned with ``smiles``.
        data_dir: Base directory for chebifier model storage.
        use_cache: Whether to use the on-disk InChIKey cache.
        resolve_names: Whether to resolve ChEBI ids to names.

    Returns:
        A tidy taxonomy :class:`pandas.DataFrame` (see :meth:`ChebifierClassifier.classify`).

    Raises:
        ChebifierMissingError: If ``chebifier`` is not installed.
    """
    classifier = ChebifierClassifier(
        data_dir=data_dir, use_cache=use_cache, resolve_names=resolve_names
    )
    return classifier.classify(smiles, inchikeys=inchikeys)
