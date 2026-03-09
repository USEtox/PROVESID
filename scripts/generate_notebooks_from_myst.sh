#!/usr/bin/env bash
set -euo pipefail

# Generate .ipynb notebooks from MyST Markdown tutorial sources under examples/.
# Requires: jupytext CLI available in active environment.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

if ! command -v jupytext >/dev/null 2>&1; then
  echo "jupytext is required. Install with: pip install jupytext"
  exit 1
fi

mapfile -t MYST_FILES < <(find examples -type f -name "*.md" | sort)

if [[ ${#MYST_FILES[@]} -eq 0 ]]; then
  echo "No Markdown tutorial files found under examples/."
  exit 0
fi

for md in "${MYST_FILES[@]}"; do
  nb="${md%.md}.ipynb"
  echo "Generating $nb from $md"
  jupytext --to notebook "$md" --output "$nb"
done

echo "Done. Notebook files generated from Markdown sources."
