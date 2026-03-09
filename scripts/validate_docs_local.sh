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
  "docs/quickstart.md"
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
import os
import nbformat
from nbclient import NotebookClient

tmp_dir = Path("/tmp/provesid-docs-validate")
cas_key_available = bool(os.getenv("CCC_API_KEY") or os.getenv("CAS_API_KEY"))
failures = []

for nb_path in sorted(tmp_dir.glob("*.ipynb")):
  if nb_path.name == "CAS_Common_Chemistry_tutorial.ipynb" and not cas_key_available:
    print(f"  - skipping {nb_path} (CCC_API_KEY/CAS_API_KEY not set)")
    continue

  print(f"  - executing {nb_path}")
  nb = nbformat.read(nb_path, as_version=4)
  client = NotebookClient(nb, timeout=600, kernel_name="python3", allow_errors=False)
  try:
    client.execute()
  except Exception as exc:
    failures.append((str(nb_path), str(exc)))

if failures:
  print("Execution failures detected:")
  for path, err in failures:
    print(f"  - {path}: {err}")
  raise SystemExit(1)
PY
else
  echo "[3/3] Skipping execution. Re-run with --execute to run code cells locally."
fi

echo "Docs validation complete."