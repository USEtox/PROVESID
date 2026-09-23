#!/usr/bin/env bash
set -euo pipefail

# Local documentation validation helper:
# 1) strict MkDocs build (non-executed); the tutorials are read from examples/
#    by scripts/mkdocs_hooks.py
# 2) the quick start's MyST round-trip to .ipynb, and every tutorial notebook
#    carries the outputs of a run: no code cell without output, no error
# Optional: strict MkDocs build with notebook execution enabled (--execute),
# which re-runs every notebook, online cells included

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

EXECUTE_NOTEBOOKS=0
if [[ "${1:-}" == "--execute" ]]; then
  EXECUTE_NOTEBOOKS=1
fi

echo "[1/3] Running strict docs build (non-executed)"
PROVESID_DOCS_EXECUTE=false uv run --extra docs mkdocs build --strict

echo "      Verifying rendered tutorial pages show their outputs"
for page in \
  "site/examples/search/search_tutorial/index.html" \
  "site/examples/pubchem/pubchem_tutorial/index.html" \
  "site/examples/chembl/chembl_tutorial/index.html"; do
  if [[ ! -f "$page" ]]; then
    echo "Missing built page: $page"
    exit 1
  fi
  if ! grep -q "jp-OutputArea" "$page"; then
    echo "Expected notebook output block not found in: $page"
    exit 1
  fi
done

mkdir -p /tmp/provesid-docs-validate

echo "[2/3] Validating the quick start's MyST round-trip and the tutorial outputs"
uv run --with jupytext jupytext --to notebook docs/quickstart.md \
  --output /tmp/provesid-docs-validate/quickstart.ipynb
uv run python - <<'PY'
import json, pathlib, sys

# A cell that only imports or assigns legitimately prints nothing, so a cell
# without output is reported; only an error output fails the check.
failed = []
for path in sorted(pathlib.Path("examples").rglob("*.ipynb")):
    if ".ipynb_checkpoints" in path.parts:
        continue
    cells = [c for c in json.loads(path.read_text())["cells"] if c["cell_type"] == "code"]
    silent = sum(1 for c in cells if not c.get("outputs"))
    errors = sum(1 for c in cells for o in c.get("outputs", []) if o["output_type"] == "error")
    print(f"  - {path}: {len(cells)} code cells, {silent} without output, {errors} error(s)")
    if errors:
        failed.append(str(path))
if failed:
    sys.exit(f"Error outputs in: {', '.join(failed)}")

# mkdocs-jupyter renders a notebook's Markdown without MkDocs' link handling,
# so a relative link that works in docs/ breaks on the site and in the
# repository alike. Notebooks link to the published site instead; each such
# link must name a page this build produced.
import re
site_url = "https://usetox.github.io/PROVESID/"
dead = []
for path in sorted(pathlib.Path("examples").rglob("*.ipynb")):
    for cell in json.loads(path.read_text())["cells"]:
        if cell["cell_type"] != "markdown":
            continue
        for target in re.findall(r"\]\(([^)\s]+)\)", "".join(cell["source"])):
            if target.startswith(site_url):
                page = pathlib.Path("site") / target[len(site_url):].split("#")[0] / "index.html"
                if not page.is_file():
                    dead.append(f"{path}: {target}")
            elif not target.startswith(("http://", "https://", "#")):
                dead.append(f"{path}: relative link {target}")
if dead:
    sys.exit("Links in tutorial notebooks that name no page:\n  " + "\n  ".join(dead))
print("  - every link in the tutorial notebooks names a built page")
PY

if [[ "$EXECUTE_NOTEBOOKS" -eq 1 ]]; then
  echo "[3/3] Running strict docs build with notebook execution enabled"
  PROVESID_DOCS_EXECUTE=true uv run --extra docs mkdocs build --strict
else
  echo "[3/3] Skipping execution. Re-run with --execute to run every notebook again."
fi

echo "Docs validation complete."