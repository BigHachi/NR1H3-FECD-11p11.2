#!/bin/bash
# ============================================================
# Step 77: Glaucoma S-PrediXcan + S-MultiXcan — whole-blood panel
#   Panel: Cells_Cultured_fibroblasts + Brain_Cortex + Whole_Blood
#   (skin replaced by whole blood). Fibroblast + brain S-PrediXcan
#   already exist; runs whole blood, then S-MultiXcan across the 3.
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -euo pipefail
cd "$FECD_BASE"
METAXCAN="$METAXCAN_DIR"
GWAS=data/gwas/glaucoma_gwas_sprediXcan.txt.gz
PREDICTDB=data/predictdb
SPX=results/sprediXcan
STAGED=results/sprediXcan_glaucoma_wb_staged
OUT_DIR=results/multixcan
LOG=logs/77_glaucoma_wholeblood_panel.log
mkdir -p "$SPX" "$STAGED" "$OUT_DIR" logs
exec > >(tee "$LOG") 2>&1

# ---- 1. S-PrediXcan: Whole_Blood ----
echo "=== S-PrediXcan: Whole_Blood ==="
WB_OUT=${SPX}/glaucoma_whole_blood.csv
if [ ! -f "$WB_OUT" ]; then
  $METAXCAN_PY ${METAXCAN}/SPrediXcan.py \
    --model_db_path            ${PREDICTDB}/mashr_Whole_Blood.db \
    --covariance               ${PREDICTDB}/mashr_Whole_Blood.txt.gz \
    --gwas_file                ${GWAS} \
    --snp_column               panel_variant_id \
    --effect_allele_column     effect_allele \
    --non_effect_allele_column non_effect_allele \
    --beta_column              effect_size \
    --se_column                standard_error \
    --pvalue_column            pvalue \
    --freq_column              frequency \
    --keep_non_rsid --model_db_snp_key varID --remove_ens_version \
    --output_file              "$WB_OUT"
fi
echo "  Whole_Blood genes: $(tail -n +2 "$WB_OUT" | wc -l | tr -d ' ')"

# ---- 2. Stage panel ----
echo ""
echo "=== Staging fibroblast + brain + whole blood ==="
rm -f ${STAGED}/*.csv
cp ${SPX}/glaucoma_cells_cultured_fibroblasts.csv ${STAGED}/mashr_Cells_Cultured_fibroblasts.csv
cp ${SPX}/glaucoma_brain_cortex.csv               ${STAGED}/mashr_Brain_Cortex.csv
cp ${SPX}/glaucoma_whole_blood.csv                ${STAGED}/mashr_Whole_Blood.csv
echo "  staged 3 tissues"

# ---- 3. S-MultiXcan ----
echo ""
echo "=== S-MultiXcan: fibroblast + brain + whole blood ==="
COVAR=${PREDICTDB}/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz
OUT_FILE=${OUT_DIR}/glaucoma_smultixcan_wb.txt
$METAXCAN_PY ${METAXCAN}/SMulTiXcan.py \
  --models_folder                    "$PREDICTDB" \
  --models_name_filter               "mashr_Cells_Cultured_fibroblasts.db|mashr_Brain_Cortex.db|mashr_Whole_Blood.db" \
  --models_name_pattern              "(.*).db" \
  --snp_covariance                   "$COVAR" \
  --metaxcan_folder                  "$STAGED" \
  --metaxcan_filter                  "mashr_.*csv" \
  --metaxcan_file_name_parse_pattern "(mashr_.*)\\.csv" \
  --gwas_file                        "$GWAS" \
  --snp_column                       panel_variant_id \
  --effect_allele_column             effect_allele \
  --non_effect_allele_column         non_effect_allele \
  --beta_column                      effect_size \
  --se_column                        standard_error \
  --pvalue_column                    pvalue \
  --freq_column                      frequency \
  --keep_non_rsid --trimmed_ensemble_id --model_db_snp_key varID \
  --cutoff_condition_number 30 \
  --output                           "$OUT_FILE"
echo ""
echo "=== Done ==="
echo "S-MultiXcan (whole-blood panel) genes: $(tail -n +2 "$OUT_FILE" | wc -l | tr -d ' ')"
