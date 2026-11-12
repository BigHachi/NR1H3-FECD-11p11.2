#!/bin/bash
# ============================================================
# Step 76h (29 Oct 2026): genome-wide FUSION whole-blood TWAS + conditional/joint analysis for
#   glaucoma (FinnGen R12 H7_GLAUCOMA), with the input as reported and after the allele-matching
#   correction of step 55c applied to every chromosome. Recomputes the genome-wide ranking of the
#   conditionally independent signals that step 76f reads from the June 2025 workstation run.
#   Weights: GTEx v8 EUR Whole_Blood FUSION (step 39). LD: 1000G EUR (FUSION LDREF).
#   12 Nov 2026: steps 1-3 run through scripts/fusion_genomewide.sh, the procedure used for every trait (FECD step
#   82, IOP and refractive error step 76i): per-chromosome inputs restricted to the LD-reference SNPs,
#   FUSION.assoc_test, models in the MHC kept aside (they were in the conditional analysis up to release v14.5),
#   genome-wide Bonferroni threshold, FUSION.post_process on the significant models (100-kb locus window).
#   Then the ranking of the jointly included genes and the comparison with the June 2025 run
#   (76h_glaucoma_wb_genomewide_rank.py). Step 76d takes the chromosome-11 results for S8 from the allele-corrected run.
#   10 Nov 2026: allele_corrected/ is the input reported from release v14 (step 60b's allele-matched file);
#   as_reported/ keeps its name but is the first-record input reported up to release v13 (step 60).
#   Outputs: results/twas_wb_genomewide/{as_reported,allele_corrected}/
#            analysis/glaucoma_wb_genomewide_joint_ranking.tsv
#            analysis/glaucoma_wb_genomewide_rerun_summary.tsv
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
cd "$FECD_BASE"
OUT=results/twas_wb_genomewide
WPOS=data/weights_eur/GTExv8.EUR.Whole_Blood.pos
mkdir -p "$OUT" logs

bash "$SCRIPTS_DIR/fusion_genomewide.sh" glaucoma_wb data/gwas/glaucoma_gwas_fusion.tsv "$WPOS" "$OUT/allele_corrected"
bash "$SCRIPTS_DIR/fusion_genomewide.sh" glaucoma_wb data/gwas/glaucoma_gwas_fusion_firstrecord.tsv.gz "$WPOS" "$OUT/as_reported"
"${PY:-python3}" "$SCRIPTS_DIR/76h_glaucoma_wb_genomewide_rank.py"
