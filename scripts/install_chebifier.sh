#!/usr/bin/env bash
#
# install_chebifier.sh — install the optional chebifier backend for PROVESID,
# including the graph (GNN) models, on Linux/CPU.
#
# See: docs/chebifier.md  and
#      plans/2026-07-02-chemical-taxonomy-classyfire-chebifier.md  (§10)
#
# Since chebifier 1.2.2 this is essentially two pip commands, because upstream
# now ships a `models` extra that pins the whole model stack:
#
#   uv pip install "chebifier[models]"
#   uv pip install torch==2.12.0 torch_scatter torch_geometric \
#       -f https://data.pyg.org/whl/torch-2.12.0+cpu.html
#
# This script does exactly that (plus `provesid[chebifier]` and a verify step),
# and by default pulls torch from PyTorch's CPU index: measured 1.6 GB of
# site-packages instead of 5.4 GB, since plain PyPI torch adds 2.7 GB of CUDA
# wheels plus triton. Set TORCH_INDEX_URL="" to use plain PyPI torch.
#
# What the `models` extra buys us (verified by running the ensemble):
#
#  * chebai-graph is pinned to 1.0.0, whose one-hot property vocabularies still
#    match chebifier's "v244" GNN checkpoints. The index-file drift that used to
#    break the graph models ("mat1 and mat2 shapes cannot be multiplied") came
#    from chebai-graph 1.1.0; with the pin there is nothing left to patch.
#    `provesid.taxonomy.ensure_v244_indices()` still runs as a safety net and is
#    a no-op on a clean install.
#  * chemlog-extra and c3p now come from PyPI, so the old git installs are gone.
#  * Only torch_scatter (plus pure-python torch_geometric) is needed — not
#    torch_sparse/torch_cluster/pyg_lib. torch_cluster was what previously
#    forced torch <= 2.11, so torch 2.12 is now fine.
#
# Overridable via environment variables:
#   TORCH_VERSION      PyTorch version              (default 2.12.0)
#   CHEBIFIER_VERSION  chebifier release to pin     (default 1.2.2)
#   TORCH_INDEX_URL    torch wheel index            (default PyTorch CPU index;
#                                                    set to "" for plain PyPI)
#   PIP                installer command            (default: uv pip, else python -m pip)
#
# Usage:
#   bash scripts/install_chebifier.sh
#
set -euo pipefail

TORCH_VERSION="${TORCH_VERSION:-2.12.0}"
CHEBIFIER_VERSION="${CHEBIFIER_VERSION:-1.2.2}"
TORCH_INDEX_URL="${TORCH_INDEX_URL-https://download.pytorch.org/whl/cpu}"

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

# Verify with the interpreter we install *into*: when VIRTUAL_ENV is exported but
# the env is not activated, a bare `python3` is the system interpreter and the
# verification below would fail even though the install succeeded.
if [[ -n "${VIRTUAL_ENV:-}" && -x "${VIRTUAL_ENV}/bin/python" ]]; then
    PY="${VIRTUAL_ENV}/bin/python"
else
    PY="$(command -v python3 || command -v python || true)"
fi
if [[ -z "${PY}" ]]; then
    echo "ERROR: no python interpreter on PATH." >&2
    exit 1
fi

log() { printf '\n==> %s\n' "$*"; }

# --- 0. platform sanity check ------------------------------------------------
if [[ "$(uname -s)" != "Linux" ]]; then
    echo "WARNING: this script targets Linux; on other platforms the PyG wheel" >&2
    echo "         tags differ (see https://data.pyg.org/whl/)." >&2
fi

PYG_WHL="https://data.pyg.org/whl/torch-${TORCH_VERSION}+cpu.html"

if ! curl -sfIL "${PYG_WHL}" >/dev/null 2>&1; then
    echo "ERROR: no PyG wheel index for torch ${TORCH_VERSION} (cpu)." >&2
    echo "       Pick a TORCH_VERSION listed at https://data.pyg.org/whl/" >&2
    exit 1
fi

# --- 1. torch, pinned to the version the PyG wheels are built for ------------
# Installed first (and by default from the CPU index) so that chebifier's own
# resolution does not first pull a multi-GB CUDA torch that we then replace.
log "Installing torch ${TORCH_VERSION}${TORCH_INDEX_URL:+ from ${TORCH_INDEX_URL}}"
if [[ -n "${TORCH_INDEX_URL}" ]]; then
    "${INSTALL[@]}" install "torch==${TORCH_VERSION}" --index-url "${TORCH_INDEX_URL}"
else
    "${INSTALL[@]}" install "torch==${TORCH_VERSION}"
fi

# --- 2. chebifier + the whole model stack (`models` extra) -------------------
# chebifier[models] pins chebai, chebai-graph, chemlog-extra and c3p, all PyPI.
log "Installing chebifier[models]==${CHEBIFIER_VERSION} (chebai, chebai-graph, chemlog-extra, c3p)"
"${INSTALL[@]}" install "chebifier[models]==${CHEBIFIER_VERSION}"

# --- 3. the PyG bits the graph models need -----------------------------------
# torch_scatter is a compiled extension with no source install: it must come
# from the PyG wheel index matching the exact torch version and platform.
log "Installing torch_scatter + torch_geometric for torch ${TORCH_VERSION}"
"${INSTALL[@]}" install torch_scatter torch_geometric -f "${PYG_WHL}"

# --- 4. the provesid[chebifier] extra ----------------------------------------
if [[ -f "${REPO_ROOT}/pyproject.toml" ]]; then
    log "Installing provesid[chebifier] (editable, from ${REPO_ROOT})"
    "${INSTALL[@]}" install -e "${REPO_ROOT}[chebifier]"
else
    log "Installing provesid[chebifier]"
    "${INSTALL[@]}" install "provesid[chebifier]"
fi

# --- 5. verify the full stack imports (no model download) --------------------
log "Verifying imports"
"${PY}" - <<'PYCODE'
import importlib
mods = ["torch", "torch_scatter", "torch_geometric", "chebai", "chebai_graph",
        "chemlog", "chemlog_extra", "c3p", "chebifier"]
for m in mods:
    mod = importlib.import_module(m)
    print(f"  {m:16s} {getattr(mod, '__version__', 'ok')}")

# Safety net only: with chebai-graph pinned by chebifier[models] this reports
# every property index as "ok" and changes nothing.
from provesid.taxonomy import ensure_v244_indices
print(f"  chebai-graph property indices: {ensure_v244_indices()}")

print("\nchebifier stack ready (all models incl. GNN). Model weights download on")
print("first BaseEnsemble() use into the shared PROVESID dataset dir")
print("(set PROVESID_DATA_DIR to override).")
PYCODE

log "Done."
