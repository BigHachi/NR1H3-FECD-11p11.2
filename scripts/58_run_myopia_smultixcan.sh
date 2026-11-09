#!/bin/bash
# ============================================================
# Myopia Analysis - Step 48: Run S-MultiXcan (cross-tissue)
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -e
cd "$FECD_BASE"

METAXCAN="$METAXCAN_DIR"
GWAS=data/gwas/myopia_gwas_sprediXcan.txt.gz
PREDICTDB=data/predictdb
SPREDIXCAN_DIR=results/sprediXcan
STAGED_DIR=results/sprediXcan_myopia_staged
OUT_DIR=results/multixcan
OUT_FILE="${OUT_DIR}/myopia_smultixcan.txt"

mkdir -p "$OUT_DIR" "$STAGED_DIR" logs

COVAR="${PREDICTDB}/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz"

echo "=== S-MultiXcan: Myopia GWAS (cross-tissue) ==="
echo "Tissues: Esophagus_Mucosa, Skin_Not_Sun_Exposed_Suprapubic, Cells_Cultured_fibroblasts"
echo "GWAS:    $GWAS"
echo "Output:  $OUT_FILE"
echo ""

# Stage myopia S-PrediXcan outputs with mashr_ prefix for SMulTiXcan
for TISSUE in Esophagus_Mucosa Skin_Not_Sun_Exposed_Suprapubic Cells_Cultured_fibroblasts; do
  case "$TISSUE" in
    Esophagus_Mucosa)                LOWER="esophagus_mucosa" ;;
    Skin_Not_Sun_Exposed_Suprapubic) LOWER="skin_not_sun_exposed_suprapubic" ;;
    Cells_Cultured_fibroblasts)      LOWER="cells_cultured_fibroblasts" ;;
  esac
  SRC="${SPREDIXCAN_DIR}/myopia_${LOWER}.csv"
  DST="${STAGED_DIR}/mashr_${TISSUE}.csv"
  if [ ! -f "$SRC" ]; then
    echo "ERROR: S-PrediXcan output missing: $SRC"
    echo "Run 47_run_myopia_sprediXcan.sh first."
    exit 1
  fi
  cp "$SRC" "$DST"
  echo "  Staged: $DST"
done

# Check prerequisites
for f in \
  "${PREDICTDB}/mashr_Esophagus_Mucosa.db" \
  "${PREDICTDB}/mashr_Skin_Not_Sun_Exposed_Suprapubic.db" \
  "${PREDICTDB}/mashr_Cells_Cultured_fibroblasts.db" \
  "$COVAR"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: Required file missing: $f"; exit 1
  fi
done

echo ""
$METAXCAN_PY ${METAXCAN}/SMulTiXcan.py \
  --models_folder                    "$PREDICTDB" \
  --models_name_filter               "mashr_Esophagus_Mucosa.db|mashr_Skin_Not_Sun_Exposed_Suprapubic.db|mashr_Cells_Cultured_fibroblasts.db" \
  --models_name_pattern              "(.*).db" \
  --snp_covariance                   "$COVAR" \
  --metaxcan_folder                  "$STAGED_DIR" \
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
  --keep_non_rsid \
  --trimmed_ensemble_id \
  --model_db_snp_key                 varID \
  --cutoff_condition_number          30 \
  --output                           "$OUT_FILE" \
  2>&1 | tee logs/myopia_smultixcan.log

echo ""
echo "=== S-MultiXcan complete ==="
if [ -f "$OUT_FILE" ]; then
  NGENES=$(tail -n +2 "$OUT_FILE" | wc -l | tr -d ' ')
  echo "Output: $OUT_FILE"
  echo "Genes tested: $NGENES"
else
  echo "WARNING: No output — check logs/myopia_smultixcan.log"
  exit 1
fi
