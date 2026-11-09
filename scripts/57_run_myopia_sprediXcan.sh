#!/bin/bash
# ============================================================
# Myopia Analysis - Step 47: Run S-PrediXcan (per tissue)
# GTEx v8 MASHR: Esophagus_Mucosa, Skin_Not_Sun_Exposed_Suprapubic,
#                Cells_Cultured_fibroblasts
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -e
cd "$FECD_BASE"

METAXCAN="$METAXCAN_DIR"
GWAS=data/gwas/myopia_gwas_sprediXcan.txt.gz
PREDICTDB=data/predictdb
OUT_DIR=results/sprediXcan

mkdir -p "$OUT_DIR" logs

if [ ! -f "$GWAS" ]; then
  echo "ERROR: $GWAS not found — run 46_prep_myopia_sprediXcan.sh first"
  exit 1
fi

echo "=== S-PrediXcan: Myopia GWAS x GTEx v8 MASHR ==="
echo "GWAS: $GWAS"
echo ""

for TISSUE in Esophagus_Mucosa Skin_Not_Sun_Exposed_Suprapubic Cells_Cultured_fibroblasts; do
  case "$TISSUE" in
    Esophagus_Mucosa)                LOWER="esophagus_mucosa" ;;
    Skin_Not_Sun_Exposed_Suprapubic) LOWER="skin_not_sun_exposed_suprapubic" ;;
    Cells_Cultured_fibroblasts)      LOWER="cells_cultured_fibroblasts" ;;
  esac

  DB_FILE="${PREDICTDB}/mashr_${TISSUE}.db"
  COV_FILE="${PREDICTDB}/mashr_${TISSUE}.txt.gz"
  OUT_FILE="${OUT_DIR}/myopia_${LOWER}.csv"

  echo "--- Tissue: ${TISSUE} ---"

  if [ ! -f "$DB_FILE" ]; then
    echo "  ERROR: Model not found: $DB_FILE"; exit 1
  fi
  if [ ! -f "$COV_FILE" ]; then
    echo "  ERROR: Covariance not found: $COV_FILE"; exit 1
  fi

  $METAXCAN_PY ${METAXCAN}/SPrediXcan.py \
    --model_db_path            "$DB_FILE" \
    --covariance               "$COV_FILE" \
    --gwas_file                "$GWAS" \
    --snp_column               panel_variant_id \
    --effect_allele_column     effect_allele \
    --non_effect_allele_column non_effect_allele \
    --beta_column              effect_size \
    --se_column                standard_error \
    --pvalue_column            pvalue \
    --freq_column              frequency \
    --keep_non_rsid \
    --model_db_snp_key         varID \
    --remove_ens_version \
    --output_file              "$OUT_FILE" \
    2>&1 | tee "logs/myopia_sprediXcan_${LOWER}.log"

  if [ -f "$OUT_FILE" ]; then
    NGENES=$(tail -n +2 "$OUT_FILE" | wc -l | tr -d ' ')
    echo "  Done: ${NGENES} genes"
  else
    echo "  WARNING: No output — check logs/myopia_sprediXcan_${LOWER}.log"
  fi
  echo ""
done

echo "=== S-PrediXcan complete ==="
ls -lh ${OUT_DIR}/myopia_*.csv 2>/dev/null || echo "  (no outputs found)"
echo ""
echo "Gene counts:"
for f in ${OUT_DIR}/myopia_*.csv; do
  [ -f "$f" ] && echo "  $(basename $f): $(tail -n +2 $f | wc -l) genes"
done
