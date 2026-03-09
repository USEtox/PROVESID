#!/usr/bin/env bash
set -euo pipefail

# Local documentation validation helper:
# 1) strict MkDocs build (non-executed)
# 2) MyST tutorial round-trip conversion to .ipynb
# Optional: strict MkDocs build with notebook execution enabled (--execute)

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

echo "[1/3] Running strict docs build (non-executed)"
PROVESID_DOCS_EXECUTE=false uv run --extra docs mkdocs build --strict

mkdir -p /tmp/provesid-docs-validate

echo "[2/3] Validating MyST -> notebook round-trip"
for tutorial in "${TUTORIAL_FILES[@]}"; do
  out_nb="/tmp/provesid-docs-validate/$(basename "${tutorial%.md}").ipynb"
  echo "  - $tutorial -> $out_nb"
  uv run --with jupytext jupytext --to notebook "$tutorial" --output "$out_nb"
done

if [[ "$EXECUTE_NOTEBOOKS" -eq 1 ]]; then
  echo "[3/3] Running strict docs build with notebook execution enabled"
  PROVESID_DOCS_EXECUTE=true uv run --extra docs mkdocs build --strict

  echo "      Verifying rendered pages contain notebook output blocks"
  for page in \
    "site/quickstart/index.html" \
    "site/examples/pubchem/pubchem_tutorial/index.html" \
    "site/examples/chembl/chembl_tutorial/index.html"; do
    if [[ ! -f "$page" ]]; then
      echo "Missing built page: $page"
      exit 1
    fi

    if ! grep -q "jp-OutputArea" "$page"; then
      echo "Expected executed output block not found in: $page"
      exit 1
    fi
  done
else
  echo "[3/3] Skipping execution. Re-run with --execute to run code cells locally."
fi

echo "Docs validation complete."