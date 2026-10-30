#!/usr/bin/env bash
# Regenerate the whole review (comparison tables, figures and REVIEW.pdf) from the
# current pipeline outputs in ../../cr2sub/. Run after main_pipeline.py.
# REVIEW.pdf is built from REVIEW.md (style: review_pdf_header.tex).
# Needs: .venv with pandas, matplotlib, geopandas; R with terra/sf; pandoc + xelatex;
# .env with CR2SUB_HEAVY_DATA (CLSoilMaps rasters). Memory capped at 8 GB.
set -euo pipefail
ulimit -v 8000000
HERE="$(cd "$(dirname "$0")" && pwd)"
REVIEW="$(dirname "$HERE")"
ROOT="$(dirname "$REVIEW")"
PY="$ROOT/.venv/bin/python"

# 1. Snapshot of the current outputs (v1.1 snapshot in versions/v1.1 is kept as is)
cp "$ROOT"/cr2sub/cr2sub_v1.2_*.csv "$REVIEW/versions/v1.2/"

# 2. Series and metadata comparison, then CLSoilMaps re-extraction (repo root)
(cd "$HERE" && python3 compare_outputs.py)
(cd "$ROOT" && Rscript review_dga_xlsx_migration/scripts/validate_clsoilmaps_extraction.R)

# 3. Attribute / CLSoilMaps comparison and figures
(cd "$HERE" && "$PY" compare_attributes.py && "$PY" compare_plots.py && "$PY" compare_vs_dga.py && "$PY" plot_series_extension.py && "$PY" plot_coverage_by_region.py)

# 4. PDF (style in review_pdf_header.tex)
cd "$REVIEW"
pandoc REVIEW.md -o REVIEW.pdf --pdf-engine=xelatex --resource-path=. \
  -H scripts/review_pdf_header.tex \
  -V geometry:margin=2cm -V fontsize=10pt -V mainfont="Latin Modern Sans" \
  -V monofont="DejaVu Sans Mono" -V monofontoptions="Scale=0.82" \
  -V colorlinks=true -V linkcolor=cr2blue -V urlcolor=cr2blue
