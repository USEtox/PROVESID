#!/usr/bin/env bash
set -euo pipefail

# Convert tutorial notebooks under examples/ to MyST Markdown sources.
# Requires: jupytext CLI available in active environment.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

if ! command -v jupytext >/dev/null 2>&1; then
  echo "jupytext is required. Install with: pip install jupytext"
  exit 1
fi

mapfile -t NOTEBOOKS < <(find examples -type f -name "*.ipynb" | sort)

if [[ ${#NOTEBOOKS[@]} -eq 0 ]]; then
  echo "No notebooks found under examples/."
  exit 0
fi

for nb in "${NOTEBOOKS[@]}"; do
  md="${nb%.ipynb}.md"
  echo "Converting $nb -> $md"
  jupytext --to myst "$nb" --output "$md"
done

echo "Done. Review generated .md files before removing .ipynb sources."
