#!/usr/bin/env bash
#
# install_chebifier.sh — install the optional chebifier backend for PROVESID,
# including the graph (GNN) models, on Linux/CPU.
#
# See: plans/2026-07-02-chemical-taxonomy-classyfire-chebifier.md  (§10)
#
# This installs the ChEB-AI `chebifier` ensemble (chebifier 1.2.1) so that ALL
# of its models run on CPU: the Electra transformer, the graph models
# (gat-aug, resgated-aug), and the rule-based ChemLog + C3P models.
#
# Two non-obvious things this script handles, both verified by actually running
# the ensemble in a throwaway env:
#
#  1. PIN TORCH TO 2.11. The graph models need chebai-graph, which requires the
#     PyG compiled extensions (torch_scatter/sparse/cluster/pyg_lib). torch_cluster
#     has NO prebuilt CPU wheel for torch 2.12 — the newest with wheels is 2.11.
#     Also, the wheels must come from the PyG index matching the exact torch
#     version + platform, and be installed BEFORE chebai-graph.
#
#  2. PATCH THE PROPERTY INDEX FILES. chebifier 1.2.1 ships GNN checkpoints
#     ("v244") trained against chebai-graph's one-hot property vocabularies as
#     they were at commit 677d44b. A later commit (ea77f36, 2026-03-02, AFTER
#     chebifier 1.2.1 shipped) APPENDED tokens to BondType / AtomNumHs /
#     NumAtomBonds, which widens the graph models' edge/node feature vectors by
#     one and makes them fail to load ("mat1 and mat2 shapes cannot be
#     multiplied"). We revert those three index files inside the installed
#     chebai-graph to their pre-drift (v244-matching) contents. This is why the
#     deployed chebifier web app works (frozen env predating the drift) but a
#     fresh unpinned install does not.
#
# Overridable via environment variables:
#   TORCH_VERSION      PyTorch version (default 2.11.0; must be 2.9-2.11 for torch_cluster)
#   CHEBIFIER_VERSION  chebifier release to pin      (default 1.2.1)
#   CHEBAI_VERSION     chebai release to pin         (default 1.3.0)
#   CHEBAI_GRAPH_VERSION  chebai-graph release       (default 1.1.0)
#   PIP                installer command             (default: uv pip, else python -m pip)
#
# Usage:
#   bash scripts/install_chebifier.sh
#
set -euo pipefail

TORCH_VERSION="${TORCH_VERSION:-2.11.0}"
CHEBIFIER_VERSION="${CHEBIFIER_VERSION:-1.2.1}"
CHEBAI_VERSION="${CHEBAI_VERSION:-1.3.0}"
CHEBAI_GRAPH_VERSION="${CHEBAI_GRAPH_VERSION:-1.1.0}"
TORCH_CPU_INDEX="https://download.pytorch.org/whl/cpu"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --- pick an installer (uv preferred, per project convention) ----------------
if [[ -n "${PIP:-}" ]]; then
    INSTALL=(${PIP})
elif command -v uv >/dev/null 2>&1; then
    INSTALL=(uv pip)
elif command -v python3 >/dev/null 2>&1; then
    INSTALL=(python3 -m pip)
else
    echo "ERROR: no installer found (need 'uv' or 'python3'). Set \$PIP." >&2
    exit 1
fi

PY="$(command -v python3 || command -v python || true)"
if [[ -z "${PY}" ]]; then
    echo "ERROR: no python interpreter on PATH." >&2
    exit 1
fi

log() { printf '\n==> %s\n' "$*"; }

# --- 0. platform sanity check ------------------------------------------------
if [[ "$(uname -s)" != "Linux" ]]; then
    echo "WARNING: this script targets Linux; on other platforms the PyG wheel" >&2
    echo "         tags (and torch_cluster availability) differ." >&2
fi

# --- 1. CPU PyTorch (pinned; must stay <= 2.11 for torch_cluster wheels) ------
if EXISTING_TORCH="$("${PY}" -c 'import torch; print(torch.__version__)' 2>/dev/null)"; then
    EXISTING_VER="${EXISTING_TORCH%%+*}"
    EX_MINOR="$(printf '%s' "${EXISTING_VER}" | cut -d. -f1-2)"
    case "${EX_MINOR}" in
        2.9|2.10|2.11)
            TORCH_VERSION="${EXISTING_VER}"
            log "Found existing PyTorch ${EXISTING_TORCH}; matching PyG wheels to ${TORCH_VERSION}"
            ;;
        *)
            echo "ERROR: existing torch ${EXISTING_TORCH} is outside 2.9-2.11, which the" >&2
            echo "       graph models' PyG wheels (esp. torch_cluster) require. Install" >&2
            echo "       torch 2.11 in a clean env, or set TORCH_VERSION and remove torch." >&2
            exit 1
            ;;
    esac
else
    log "Installing CPU PyTorch ${TORCH_VERSION}"
    "${INSTALL[@]}" install "torch==${TORCH_VERSION}" --index-url "${TORCH_CPU_INDEX}"
fi

PYG_WHL="https://data.pyg.org/whl/torch-${TORCH_VERSION}+cpu.html"

# --- 2. PyG stack (CPU wheels for this torch) --------------------------------
# torch_geometric (pure python) + the compiled extensions pyg_lib / torch_scatter
# / torch_sparse / torch_cluster, all from the matching wheel index.
log "Checking PyG CPU wheel index: ${PYG_WHL}"
if ! curl -sfIL "${PYG_WHL}" >/dev/null 2>&1; then
    echo "ERROR: no PyG wheel index for torch ${TORCH_VERSION} (cpu)." >&2
    echo "       Use a torch version in 2.9-2.11 with wheels at https://data.pyg.org/whl/" >&2
    exit 1
fi
log "Installing PyG stack (torch_geometric pyg_lib torch_scatter torch_sparse torch_cluster)"
"${INSTALL[@]}" install torch_geometric pyg_lib torch_scatter torch_sparse torch_cluster -f "${PYG_WHL}"

# --- 3. chebifier + chebai + graph backend + rule-based backends -------------
# chebi-utils is an UNDECLARED import of chebai-graph, so we install it explicitly.
log "Installing chebifier==${CHEBIFIER_VERSION}, chebai==${CHEBAI_VERSION}, chebai-graph==${CHEBAI_GRAPH_VERSION}, chebi-utils"
"${INSTALL[@]}" install \
    "chebifier==${CHEBIFIER_VERSION}" \
    "chebai==${CHEBAI_VERSION}" \
    "chebai-graph==${CHEBAI_GRAPH_VERSION}" \
    chebi-utils

log "Installing chemlog-extra (chemlog_element, chemlog_organox models)"
"${INSTALL[@]}" install "git+https://github.com/ChEB-AI/chemlog-extra.git"

log "Installing c3p (c3p model; cross-platform fork)"
"${INSTALL[@]}" install "git+https://github.com/sfluegel05/c3p.git"

# --- 4. the provesid[chebifier] extra ----------------------------------------
if [[ -f "${REPO_ROOT}/pyproject.toml" ]]; then
    log "Installing provesid[chebifier] (editable, from ${REPO_ROOT})"
    "${INSTALL[@]}" install -e "${REPO_ROOT}[chebifier]"
else
    log "Installing provesid[chebifier]"
    "${INSTALL[@]}" install "provesid[chebifier]"
fi

# --- 5. patch chebai-graph property indices to the v244-matching state --------
# Revert BondType / AtomNumHs / NumAtomBonds one-hot vocabularies to their
# pre-drift contents (chebai-graph commit 677d44b) so the chebifier 1.2.1 GNN
# checkpoints load. Contents are embedded (tiny + immutable given the pin), so
# no network is needed for the patch.
log "Patching chebai-graph property index files (v244 checkpoint compatibility)"
"${PY}" - <<'PYPATCH'
import os, chebai_graph
bin_dir = os.path.join(os.path.dirname(chebai_graph.__file__), "preprocessing", "bin")

# Pre-drift (commit 677d44b) one-hot vocabularies matching the v244 checkpoints.
PRE_DRIFT = {
    "BondType":     ["DATIVE", "SINGLE", "AROMATIC", "TRIPLE", "DOUBLE"],
    "AtomNumHs":    ["0", "3", "2", "4", "1", "5", "6"],
    "NumAtomBonds": ["0", "1", "2", "4", "5", "3", "6", "8", "7", "10", "12"],
}

for prop, tokens in PRE_DRIFT.items():
    path = os.path.join(bin_dir, prop, "indices_one_hot.txt")
    if not os.path.exists(path):
        print(f"  WARNING: {path} not found (chebai-graph layout changed?) — skipping")
        continue
    with open(path) as fh:
        current = [ln.strip() for ln in fh if ln.strip()]
    if current == tokens:
        print(f"  {prop}: already v244-matching ({len(tokens)} tokens)")
        continue
    with open(path, "w") as fh:
        fh.write("\n".join(tokens) + "\n")
    print(f"  {prop}: patched {len(current)} -> {len(tokens)} tokens")
PYPATCH

# --- 6. verify the full stack imports (no model download) --------------------
log "Verifying imports"
"${PY}" - <<'PYCODE'
import importlib
mods = ["torch", "torch_scatter", "torch_sparse", "torch_cluster",
        "torch_geometric", "chebai", "chebai_graph", "chebi_utils", "chebifier"]
for m in mods:
    mod = importlib.import_module(m)
    print(f"  {m:16s} {getattr(mod, '__version__', 'ok')}")
print("\nchebifier stack ready (all models incl. GNN). Model weights download on")
print("first BaseEnsemble() use into the shared PROVESID dataset dir")
print("(set PROVESID_DATA_DIR to override).")
PYCODE

log "Done."
