#!/bin/bash
# ============================================================
# Step 76d: FUSION whole-blood TWAS + conditional/joint analysis
#   for the cross-trait ocular phenotypes (Figure 5 source).
#   Exposure weights: GTEx v8 EUR Whole_Blood FUSION (downloaded from
#     the same S3 bucket as the FECD proxy-tissue weights, step 39).
#   For each trait: FUSION.assoc_test (chr11) -> FUSION.post_process
#     (conditional/joint), then record NR1H3 marginal TWAS.Z and its
#     conditional status (conditionally independent vs conditioned away).
#   Glaucoma input: allele-matched (step 60b, 10 Nov 2026); the first-record input used up to
#     release v13 is rerun into results/twas_wb_conditional/previous_first_record/ for comparison.
#   Outputs: results/twas_wb_conditional/{trait}_wb_chr11.dat + .joint.*
#            analysis/nr1h3_crosstrait_fusion_conditional.tsv
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -euo pipefail
cd "$FECD_BASE"
BASE="$FECD_BASE"
FU="$FUSION_DIR"              # FUSION.assoc_test.R + FUSION.post_process.R (+ utils/plink_utils.R)
WGT_DIR=data/weights_eur
LD=data/ld_ref/LDREF/1000G.EUR.                 # per-chromosome 1000G EUR PLINK ref
OUT=results/twas_wb_conditional
mkdir -p "$OUT" logs

# ---- 0. whole-blood EUR FUSION weights (once) ----
if [ ! -f "$WGT_DIR/GTExv8.EUR.Whole_Blood.pos" ]; then
  curl -sSL -o "$WGT_DIR/GTExv8.EUR.Whole_Blood.tar.gz" \
    "https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Whole_Blood.tar.gz"
  ( cd "$WGT_DIR" && tar -xzf GTExv8.EUR.Whole_Blood.tar.gz )
fi
# chr11 subset of the WB weight positions (NR1H3 = 11p11.2)
head -1 "$WGT_DIR/GTExv8.EUR.Whole_Blood.pos" > /tmp/wb_chr11.pos
awk -F'\t' '$4==11' "$WGT_DIR/GTExv8.EUR.Whole_Blood.pos" >> /tmp/wb_chr11.pos

# ---- 1. FUSION-format GWAS restricted to chr11 SNPs (speed) ----
#   (data/gwas/{trait}_gwas_fusion.tsv built by step 60: SNP A1 A2 Z N)
awk '{print $2}' data/ld_ref/LDREF/1000G.EUR.11.bim | sort -u > /tmp/chr11_snps.txt
for T in glaucoma iop myopia; do
  awk 'NR==FNR{k[$1]=1;next} FNR==1{print;next} ($1 in k)' \
      /tmp/chr11_snps.txt data/gwas/${T}_gwas_fusion.tsv > data/gwas/${T}_gwas_fusion_chr11.tsv
done
# glaucoma as analysed up to release v13 (first FinnGen record per rsID, step 60), for the
# comparison with the allele-matched input of step 60b (additional output, step 139)
gzip -cd data/gwas/glaucoma_gwas_fusion_firstrecord.tsv.gz | awk 'NR==FNR{k[$1]=1;next} FNR==1{print;next} ($1 in k)' \
    /tmp/chr11_snps.txt - > data/gwas/glaucoma_gwas_fusion_firstrecord_chr11.tsv

# ---- 2. assoc_test + conditional/joint per trait ----
for T in glaucoma iop myopia; do
  Rscript "$FU/FUSION.assoc_test.R" \
    --sumstats data/gwas/${T}_gwas_fusion_chr11.tsv \
    --weights /tmp/wb_chr11.pos --weights_dir "$WGT_DIR" \
    --ref_ld_chr "$LD" --chr 11 \
    --out ${OUT}/${T}_wb_chr11.dat > logs/fus_${T}.log 2>&1
  Rscript "$FU/FUSION.post_process.R" \
    --sumstats data/gwas/${T}_gwas_fusion_chr11.tsv \
    --input ${OUT}/${T}_wb_chr11.dat \
    --out ${OUT}/${T}_wb_chr11.joint \
    --ref_ld_chr "$LD" --chr 11 --locus_win 100000 >> logs/fus_${T}.log 2>&1
done

# glaucoma, first-record input (releases up to v13): same commands, separate folder
PREV=${OUT}/previous_first_record; mkdir -p "$PREV"
Rscript "$FU/FUSION.assoc_test.R" \
  --sumstats data/gwas/glaucoma_gwas_fusion_firstrecord_chr11.tsv \
  --weights /tmp/wb_chr11.pos --weights_dir "$WGT_DIR" \
  --ref_ld_chr "$LD" --chr 11 \
  --out ${PREV}/glaucoma_wb_chr11.dat > logs/fus_glaucoma_firstrecord.log 2>&1
Rscript "$FU/FUSION.post_process.R" \
  --sumstats data/gwas/glaucoma_gwas_fusion_firstrecord_chr11.tsv \
  --input ${PREV}/glaucoma_wb_chr11.dat \
  --out ${PREV}/glaucoma_wb_chr11.joint \
  --ref_ld_chr "$LD" --chr 11 --locus_win 100000 >> logs/fus_glaucoma_firstrecord.log 2>&1

# ---- 3. consolidate NR1H3 marginal + conditional status ----
SUM=analysis/nr1h3_crosstrait_fusion_conditional.tsv
printf "trait\ttissue\tmethod\tTWAS_Z\tTWAS_P\tJOINT_Z\tJOINT_P\tconditional_status\n" > "$SUM"
for T in glaucoma iop myopia; do
  inc=${OUT}/${T}_wb_chr11.joint.joint_included.dat
  drp=${OUT}/${T}_wb_chr11.joint.joint_dropped.dat
  if grep -q ENSG00000025434 "$inc"; then r=$(grep ENSG00000025434 "$inc"); st=conditionally_independent
  else r=$(grep ENSG00000025434 "$drp"); st=conditioned_away; fi
  # cols: FILE ID TWAS.Z TWAS.P BETA BETA.SE (JOINT|COND).Z (JOINT|COND).P
  printf "%s\twhole_blood\tFUSION\t%s\t%s\t%s\t%s\t%s\n" "$T" \
    "$(echo "$r"|awk '{print $3}')" "$(echo "$r"|awk '{print $4}')" \
    "$(echo "$r"|awk '{print $7}')" "$(echo "$r"|awk '{print $8}')" "$st" >> "$SUM"
done
# FECD disease-proxy tissues: FUSION EUR marginal (from step 40/82 output)
printf "FECD_fibroblast\tcultured_fibroblasts\tFUSION_EUR\t-3.74\t1.85e-4\tNA\tNA\tconditionally_independent\n" >> "$SUM"
printf "FECD_skin\tskin_not_sun_exposed\tFUSION_EUR\t-3.53\t4.21e-4\tNA\tNA\tNA\n" >> "$SUM"
echo "=== NR1H3 cross-trait FUSION conditional summary ==="; column -t -s$'\t' "$SUM"
