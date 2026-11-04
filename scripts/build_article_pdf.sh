#!/usr/bin/env bash
set -euo pipefail

# Rebuild the article PDF from the canonical Markdown source.
# Requires pandoc, xelatex, and the Noto Serif font.
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root/article"
pandoc molecular_surface_potential.md \
  --from=markdown-implicit_figures \
  --to=pdf \
  --pdf-engine=xelatex \
  -V geometry:margin=2.1cm \
  -V 'mainfont=Noto Serif' \
  -V fontsize=10pt \
  -V linkcolor=blue \
  -V urlcolor=blue \
  -o molecular_surface_potential.pdf
