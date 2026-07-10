#!/usr/bin/env bash
# Assemble docs/manuscript/*.md into a single preprint PDF.
#   markdown --(preprocess)--> pandoc --citeproc (Nature CSL) --> tectonic --> PDF
# Re-runnable; touches nothing under docs/manuscript/ except this build/ dir.
set -euo pipefail
export PATH="/opt/homebrew/bin:$PATH"

BUILD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$BUILD/../../.." && pwd)"        # worktree / repo root
OUT="$BUILD/manuscript.pdf"

echo "[1/2] preprocess -> _body.md"
python3 "$BUILD/preprocess.py"

echo "[2/2] pandoc + tectonic -> manuscript.pdf"
cd "$REPO"                                    # so docs/figures/*.png resolve
pandoc "$BUILD/_body.md" \
  --metadata-file="$BUILD/meta.yaml" \
  --citeproc \
  --bibliography="$BUILD/refs.bib" \
  --csl="$BUILD/nature.csl" \
  --pdf-engine=tectonic \
  --from=markdown+tex_math_dollars+raw_tex \
  -V papersize=a4 \
  -o "$OUT"

echo "done: $OUT ($(du -h "$OUT" | cut -f1))"
