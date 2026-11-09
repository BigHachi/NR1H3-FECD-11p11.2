#!/bin/bash
# ============================================================
# Step 24b (Oct 8 2026; scripted Oct 10 2026 from the review-rerun commands):
#   FUSION joint/conditional analysis WITHIN 11p11.2 for FECD.
#   Input: FUSION TWAS on the European FECD GWAS (step 82, results/twas_fecd_eurGWAS)
#   for cultured fibroblasts and skin; all genes with EUR weights in chr11:46.2-48.3 Mb
#   (ZNF408 has no FUSION model). --minp_input 0.05 lets every nominally associated
#   gene in the window enter the joint model.
#   Manuscript: Results, "Within 11p11.2, FUSION joint/conditional analysis ..." —
#   fibroblasts: NR1H3 the only conditionally independent gene (joint Z -3.7,
#   P 1.9e-4; LRP4 / ARHGAP1 conditioned away); skin: ARHGAP1 independent,
#   NR1H3 conditioned away (conditional Z -1.73, P 0.084).
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -euo pipefail
BASE="$FECD_BASE"
FU="$FUSION_DIR"          # patched FUSION (see fusion_patch/)
LD=$BASE/data/ld_ref/LDREF/1000G.EUR.
OUT=$BASE/results/fusion_11p11_conditional
cd $BASE; mkdir -p $OUT

# chr11 FECD EUR sumstats restricted to LD-reference SNPs (FUSION format: SNP A1 A2 Z N)
awk 'NR==FNR{k[$2]=1;next} FNR==1 || ($1 in k)' $BASE/data/ld_ref/LDREF/1000G.EUR.11.bim \
    data/gwas/fecd_eur_fusion.tsv > $OUT/fecd_eur_chr11.sumstats

for T in Cells_Cultured_fibroblasts Skin_Not_Sun_Exposed_Suprapubic; do
  IN=results/twas_fecd_eurGWAS/$T/fecd_eurGWAS_${T}_all_chr.dat
  awk -F'\t' 'NR==1 || ($4+0==11 && $5>=46200000 && $5<=48300000)' $IN > $OUT/${T}_11p11.dat
  # minp 1 = all genes with European weights in the region (Supplementary Table S5, values in the text);
  # minp 0.05 = genes with marginal TWAS P < 0.05
  for MINP in 1 0.05; do
    Rscript $FU/FUSION.post_process.R --sumstats $OUT/fecd_eur_chr11.sumstats \
        --input $OUT/${T}_11p11.dat --out $OUT/${T}_11p11_minp${MINP} --ref_ld_chr $LD --chr 11 \
        --minp_input ${MINP} --locus_win 100000 --report > $OUT/${T}_minp${MINP}.log 2>&1
  done
  echo "=== $T"
  for f in $OUT/${T}_11p11_minp0.05*.joint_included.dat $OUT/${T}_11p11_minp0.05*.joint_dropped.dat; do
    [ -f "$f" ] && { echo "-- $(basename "$f")"; column -t "$f" | cut -c1-200; }
  done
done
