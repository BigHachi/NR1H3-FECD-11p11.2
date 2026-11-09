#!/bin/bash
# ============================================================
# Step 139b (10 Nov 2026, review): Figures 1-7 as separate files for the journal, written by the plotting
#   steps themselves with FIG_DIR set (fig_common.R save_fig): figures/FigureN.pdf (vector, fonts embedded;
#   printed from the SVG by headless Google Chrome/Chromium, skipped with a message if neither is installed),
#   figures/FigureN.svg and figures/FigureN.tif (600 dpi, LZW). Re-plots from existing outputs; no analysis is rerun.
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
cd "$FECD_BASE"
export FIG_DIR="$FECD_BASE/figures"
mkdir -p "$FIG_DIR"
Rscript "$SCRIPTS_DIR/69b_fig1_fig2_v2.R"                     # Figures 1-2
Rscript "$SCRIPTS_DIR/81c_fecd_coloc_heatmap_fullcis.R"       # Figure 3 (rewrites S6 output identically)
Rscript "$SCRIPTS_DIR/87_nr1h3_fecd_mr_forest.R" EUR          # Figure 4
Rscript "$SCRIPTS_DIR/76e_fig5_glaucoma_label.R"              # Figure 5
Rscript "$SCRIPTS_DIR/89b_nr1h3_cismr_standardised_perSNP_N.R" # Figure 6
Rscript "$SCRIPTS_DIR/70b_flowchart_NR1H3_v2.R"               # Figure 7
for i in 1 2 3 4 5 6 7; do
  for e in pdf svg tif; do [ -s "$FIG_DIR/Figure$i.$e" ] && printf "Figure%s.%s  " "$i" "$e" || printf "Figure%s.%s MISSING  " "$i" "$e"; done
  echo
done
