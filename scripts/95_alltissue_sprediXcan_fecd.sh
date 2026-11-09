#!/bin/bash
# ============================================================
# Step 95: All-tissue S-PrediXcan — FECD x GTEx v8 MASHR (49 tissues)
#   Extends the 3 proxy tissues to the full GTEx v8 panel so NR1H3's
#   tissue profile can be read against a genome-wide background.
#   Models: PredictDB GTEx v8 MASHR (Zenodo 3518299).
#   Logs the % of each model's SNPs matched in the GWAS - a low value
#   means a varID/allele-orientation mismatch, not a biological result.
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -u
cd "$FECD_BASE"

METAXCAN="$METAXCAN_DIR"
M=data/predictdb_all/eqtl/mashr
GWAS=${GWAS_FILE:-data/gwas/fecd_gwas_sprediXcan.txt.gz}
OUT_DIR=${OUT_DIR:-results/sprediXcan_alltissue}
TAG=${TAG:-fecd}

mkdir -p "$OUT_DIR" logs
SNPLOG="$OUT_DIR/_snp_usage.tsv"
printf "tissue\tpct_snps_used\tn_genes\n" > "$SNPLOG"

n=0
for DB in "$M"/mashr_*.db; do
  T=$(basename "$DB" .db); T=${T#mashr_}
  n=$((n+1))
  OUT="$OUT_DIR/${TAG}_${T}.csv"
  LOG="logs/${TAG}_alltissue_${T}.log"
  printf "[%2d/49] %-45s " "$n" "$T"
  rm -f "$OUT"          # S-PrediXcan silently skips if the output already exists
  $METAXCAN_PY ${METAXCAN}/SPrediXcan.py \
    --model_db_path            "$DB" \
    --covariance               "$M/mashr_${T}.txt.gz" \
    --gwas_file                "$GWAS" \
    --snp_column               panel_variant_id \
    --effect_allele_column     effect_allele \
    --non_effect_allele_column non_effect_allele \
    --beta_column              effect_size \
    --se_column                standard_error \
    --pvalue_column            pvalue \
    --freq_column              frequency \
    --keep_non_rsid --model_db_snp_key varID --remove_ens_version \
    --output_file              "$OUT" > "$LOG" 2>&1
  PCT=$(grep -o "[0-9]\+ % of model's snps used" "$LOG" | grep -o "^[0-9]\+" | tail -1)
  NG=$([ -f "$OUT" ] && echo $(( $(wc -l < "$OUT") - 1 )) || echo 0)
  printf "snps_used=%3s%%  genes=%s\n" "${PCT:-NA}" "$NG"
  printf "%s\t%s\t%s\n" "$T" "${PCT:-NA}" "$NG" >> "$SNPLOG"
done
echo "=== done: $OUT_DIR ==="
