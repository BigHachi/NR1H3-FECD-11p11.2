#!/bin/bash
# ============================================================
# Step 76i (12 Nov 2026): genome-wide FUSION whole-blood TWAS + conditional/joint analysis for the other two
#   cross-trait phenotypes, IOP (Bonnemaijer et al. 2019, GCST009413) and refractive error (Hysi et al. 2020,
#   GCST010003 public subset), as step 76h does for glaucoma. Step 76d runs all three traits on chromosome 11
#   (S8, Figure 5); this step adds the genome-wide context for IOP and refractive error, as the June 2025
#   workstation run had for all three traits.
#   Weights: GTEx v8 EUR Whole_Blood FUSION (step 39). LD: 1000G EUR (FUSION LDREF).
#   Runs scripts/fusion_genomewide.sh, the procedure used for every trait (FECD step 82, glaucoma step 76h), on
#   data/gwas/{iop,myopia}_gwas_fusion.tsv (step 60): per-chromosome inputs restricted to the LD-reference SNPs,
#   FUSION.assoc_test, models in the MHC kept aside, genome-wide Bonferroni threshold, FUSION.post_process on the
#   significant models (100-kb locus window). Then the NR1H3 status and the jointly included genes, compared with
#   the June 2025 run when the non-public workbook is present (76i_crosstrait_wb_genomewide_summary.py).
#   Step 76d takes the chromosome-11 results for S8 from these runs.
#   Outputs: results/twas_wb_genomewide_crosstrait/{iop,myopia}/
#            analysis/crosstrait_wb_genomewide_joint_ranking.tsv
#            analysis/crosstrait_wb_genomewide_summary.tsv
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
cd "$FECD_BASE"
OUT=results/twas_wb_genomewide_crosstrait
WPOS=data/weights_eur/GTExv8.EUR.Whole_Blood.pos
mkdir -p "$OUT" logs

for T in iop myopia; do
  bash "$SCRIPTS_DIR/fusion_genomewide.sh" "${T}_wb" "data/gwas/${T}_gwas_fusion.tsv" "$WPOS" "$OUT/$T"
done
"${PY:-python3}" "$SCRIPTS_DIR/76i_crosstrait_wb_genomewide_summary.py"
