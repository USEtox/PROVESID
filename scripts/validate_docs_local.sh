#!/usr/bin/env bash
set -euo pipefail

# Local documentation validation helper:
# 1) strict MkDocs build (non-executed)
# 2) MyST tutorial round-trip conversion to .ipynb
# Optional: execute round-tripped notebooks with nbclient (--execute)

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

EXECUTE_NOTEBOOKS=0
if [[ "${1:-}" == "--execute" ]]; then
  EXECUTE_NOTEBOOKS=1
fi

TUTORIAL_FILES=(
  "docs/examples/CCC/CAS_Common_Chemistry_tutorial.md"
  "docs/examples/ChEBI/ChEBI_tutorial.md"
  "docs/examples/ChEBI/chebi_sdf_tutorial.md"
  "docs/examples/ClassyFire/classyfire_tutorial.md"
  "docs/examples/OPSIN/opsin_tutorial.md"
  "docs/examples/pubchem/pubchem_tutorial.md"
  "docs/examples/pubchemview/pubchem_view_tutorial.md"
  "docs/examples/resolver/chem_id_resolver_tutorial.md"
  "docs/examples/chembl/chembl_tutorial.md"
  "docs/examples/zeropm/zeropm-example.md"
)

echo "[1/3] Running strict docs build"
uv run --extra docs mkdocs build --strict

mkdir -p /tmp/provesid-docs-validate

echo "[2/3] Validating MyST -> notebook round-trip"
for tutorial in "${TUTORIAL_FILES[@]}"; do
  out_nb="/tmp/provesid-docs-validate/$(basename "${tutorial%.md}").ipynb"
  echo "  - $tutorial -> $out_nb"
  uv run --with jupytext jupytext --to notebook "$tutorial" --output "$out_nb"
done

if [[ "$EXECUTE_NOTEBOOKS" -eq 1 ]]; then
  echo "[3/3] Executing round-tripped notebooks with nbclient"
  uv run --with nbclient --with nbformat python - <<'PY'
from pathlib import Path
import nbformat
from nbclient import NotebookClient

tmp_dir = Path("/tmp/provesid-docs-validate")
for nb_path in sorted(tmp_dir.glob("*.ipynb")):
    print(f"  - executing {nb_path}")
    nb = nbformat.read(nb_path, as_version=4)
    client = NotebookClient(nb, timeout=600, kernel_name="python3", allow_errors=False)
    client.execute()
PY
else
  echo "[3/3] Skipping execution. Re-run with --execute to run code cells locally."
fi

echo "Docs validation complete."