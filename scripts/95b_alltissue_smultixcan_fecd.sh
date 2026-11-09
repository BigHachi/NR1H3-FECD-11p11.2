#!/bin/bash
# ============================================================
# Step 95b: All-tissue S-MultiXcan — FECD x 49 GTEx v8 MASHR tissues
#   Review point 7: show the NR1H3 signal is not a consequence of
#   selecting 3 proxy tissues. Same settings as the published 3-tissue
#   run (step 19): condition-number cutoff 30, MASHR SNP covariance.
#   Input: per-tissue S-PrediXcan from step 95 run on the remapped GWAS
#          (step 95a), e.g.
#     GWAS_FILE=data/gwas/fecd_alltissue_sprediXcan.txt.gz \
#     OUT_DIR=results/sprediXcan_alltissue_remap bash scripts/95_alltissue_sprediXcan_fecd.sh
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -eu
cd "$FECD_BASE"

METAXCAN="$METAXCAN_DIR"
M=data/predictdb_all/eqtl/mashr
GWAS=${GWAS_FILE:-data/gwas/fecd_alltissue_sprediXcan.txt.gz}
SPX_DIR=${SPX_DIR:-results/sprediXcan_alltissue_remap}
STAGED=${SPX_DIR}_staged
OUT=${OUT_FILE:-results/multixcan/fecd_smultixcan_alltissue49.txt}
COVAR=data/predictdb/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz

mkdir -p "$STAGED" "$(dirname "$OUT")" logs
n=0
for f in "$SPX_DIR"/fecd_*.csv; do
  T=$(basename "$f" .csv); T=${T#fecd_}
  cp "$f" "$STAGED/mashr_${T}.csv"; n=$((n+1))
done
[ "$n" -eq 49 ] || { echo "ERROR: expected 49 S-PrediXcan outputs, found $n"; exit 1; }
rm -f "$OUT"

$METAXCAN_PY ${METAXCAN}/SMulTiXcan.py \
  --models_folder                    "$M" \
  --models_name_filter               "mashr_.*.db" \
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
  --keep_non_rsid \
  --trimmed_ensemble_id \
  --model_db_snp_key                 varID \
  --cutoff_condition_number          30 \
  --output                           "$OUT" \
  > logs/fecd_smultixcan_alltissue49.log 2>&1
echo "=== done: $OUT ($(($(wc -l < "$OUT") - 1)) genes) ==="
