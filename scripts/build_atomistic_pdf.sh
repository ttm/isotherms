#!/usr/bin/env bash
# Rebuild the molecular forward-theory manuscript from its source.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT/article"
pandoc atomistic_surface_potential.md \
  -o atomistic_surface_potential.pdf \
  --pdf-engine=pdflatex \
  -V documentclass=article \
  -V geometry:margin=0.8in \
  -V fontsize=10pt \
  -V linestretch=1.13 \
  -V colorlinks=true \
  --resource-path=.
echo "Built article/atomistic_surface_potential.pdf"
