#!/bin/bash
# ============================================================
# Step 75: Glaucoma S-PrediXcan + S-MultiXcan — brain panel
#   New panel (drops esophagus, adds brain cortex):
#     Cells_Cultured_fibroblasts  (TM/outflow proxy)
#     Brain_Cortex                (RGC/CNS neuronal proxy)
#     Skin_Not_Sun_Exposed_Suprapubic (connective-tissue proxy)
#   Reuses existing fibroblast + skin S-PrediXcan; runs brain;
#   then S-MultiXcan across the 3-tissue panel.
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -euo pipefail
cd "$FECD_BASE"

METAXCAN="$METAXCAN_DIR"
GWAS=data/gwas/glaucoma_gwas_sprediXcan.txt.gz
PREDICTDB=data/predictdb
SPX=results/sprediXcan
STAGED=results/sprediXcan_glaucoma_brain_staged
OUT_DIR=results/multixcan
LOG=logs/75_glaucoma_brain_panel.log
mkdir -p "$SPX" "$STAGED" "$OUT_DIR" logs
exec > >(tee "$LOG") 2>&1

PANEL=(Cells_Cultured_fibroblasts Brain_Cortex Skin_Not_Sun_Exposed_Suprapubic)

# ---- 1. S-PrediXcan for Brain_Cortex (fibroblast/skin already done) ----
echo "=== S-PrediXcan: Brain_Cortex ==="
BC_OUT=${SPX}/glaucoma_brain_cortex.csv
if [ ! -f "$BC_OUT" ]; then
  $METAXCAN_PY ${METAXCAN}/SPrediXcan.py \
    --model_db_path            ${PREDICTDB}/mashr_Brain_Cortex.db \
    --covariance               ${PREDICTDB}/mashr_Brain_Cortex.txt.gz \
    --gwas_file                ${GWAS} \
    --snp_column               panel_variant_id \
    --effect_allele_column     effect_allele \
    --non_effect_allele_column non_effect_allele \
    --beta_column              effect_size \
    --se_column                standard_error \
    --pvalue_column            pvalue \
    --freq_column              frequency \
    --keep_non_rsid --model_db_snp_key varID --remove_ens_version \
    --output_file              "$BC_OUT"
fi
echo "  Brain_Cortex genes: $(tail -n +2 "$BC_OUT" | wc -l | tr -d ' ')"

# ---- 2. Stage the 3 panel tissues for S-MultiXcan ----
echo ""
echo "=== Staging panel for S-MultiXcan ==="
rm -f ${STAGED}/*.csv
for T in "${PANEL[@]}"; do
  case "$T" in
    Cells_Cultured_fibroblasts)      SRC=${SPX}/glaucoma_cells_cultured_fibroblasts.csv ;;
    Brain_Cortex)                    SRC=${SPX}/glaucoma_brain_cortex.csv ;;
    Skin_Not_Sun_Exposed_Suprapubic) SRC=${SPX}/glaucoma_skin_not_sun_exposed_suprapubic.csv ;;
  esac
  cp "$SRC" "${STAGED}/mashr_${T}.csv"
  echo "  staged mashr_${T}.csv"
done

# ---- 3. S-MultiXcan across the brain panel ----
echo ""
echo "=== S-MultiXcan: fibroblast + brain + skin ==="
COVAR=${PREDICTDB}/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz
OUT_FILE=${OUT_DIR}/glaucoma_smultixcan_brain.txt
$METAXCAN_PY ${METAXCAN}/SMulTiXcan.py \
  --models_folder                    "$PREDICTDB" \
  --models_name_filter               "mashr_Cells_Cultured_fibroblasts.db|mashr_Brain_Cortex.db|mashr_Skin_Not_Sun_Exposed_Suprapubic.db" \
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
echo "S-MultiXcan (brain panel) genes: $(tail -n +2 "$OUT_FILE" | wc -l | tr -d ' ')"
