#!/bin/bash
# ============================================================
# Step 83a-TSPAN10 (13 Nov 2026: copy of the NR1H3 script with only the gene changed to TSPAN10, ENSG00000182612, 17q25.3,
#   and the LD reference to chromosome 17; run by run_full_pipeline.sh from release v14.8; S7 TSPAN10 rows): eQTLGen whole-blood cis-eQTLs for TSPAN10
#   eQTLGen phase I cis-eQTLs (FDR < 0.05, probe level; N up to 31,684), filtered to
#   ENSG00000182612. Output data/eqtl_eqtlgen/tspan10_eqtlgen.txt is the exposure file
#   for 91_mr_tspan10_eur.R and 83_mr_tspan10_eqtlgen_fecd.R.
#   Also writes data/eqtl_eqtlgen/tspan10_eqtlgen_af.tsv (eQTLGen allele frequencies, 10 Nov 2026).
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -euo pipefail
cd "$FECD_BASE"; mkdir -p data/eqtl_eqtlgen
F=data/eqtl_eqtlgen/eqtlgen_cis_FDR05.txt.gz
[ -s $F ] || curl -L --fail -o $F \
  https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
gzip -cd $F | awk -F'\t' 'NR==1 || $8=="ENSG00000182612"' > data/eqtl_eqtlgen/tspan10_eqtlgen.txt
echo "TSPAN10 cis-eQTL rows: $(($(wc -l < data/eqtl_eqtlgen/tspan10_eqtlgen.txt)-1))"

# eQTLGen allele frequencies (10 Nov 2026, review): AlleleB frequency over the eQTLGen cohorts with
# genotype counts. Every cis-MR step standardises the exposure Z with the frequency of the assessed allele
# in this file (exposure source); up to release v13 the outcome GWAS frequency was used instead.
AF=data/eqtl_eqtlgen/eqtlgen_SNP_AF_AlleleB.txt.gz
[ -s $AF ] || curl -L --fail -o $AF \
  https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz
gzip -cd $AF | awk -F'\t' 'NR==FNR{if(FNR>1)k[$2]=1;next} FNR==1||($1 in k)' data/eqtl_eqtlgen/tspan10_eqtlgen.txt - \
  > data/eqtl_eqtlgen/tspan10_eqtlgen_af.tsv
echo "eQTLGen allele frequencies for TSPAN10 cis-eQTL variants: $(($(wc -l < data/eqtl_eqtlgen/tspan10_eqtlgen_af.tsv)-1))"

# FECD GWAS rows for those eQTLGen variants (outcome files for cis-MR: 83, 88 multi-ancestry; 91 European).
# European BETA/SE from the OR and 95% CI (SE = (ln CIU - ln CIL) / (2 x qnorm(0.975))).
E=data/eqtl_eqtlgen/tspan10_eqtlgen.txt
gzip -cd data/gwas/fecd_gwas_full.tsv.gz | awk -F'\t' -v OFS='\t' 'NR==FNR{if(FNR>1)k[$2]=1;next}
  FNR==1{print "SNP","A1","A2","BETA","SE","P","FREQ";next} ($1 in k){print $1,$4,$5,$6,$7,$10,$9}' $E - \
  > data/gwas/fecd_tspan10_eqtlgen_region.tsv
gzip -cd data/gwas/fecd_eur_GCST90475865.tsv.gz | awk -F'\t' -v OFS='\t' 'NR==FNR{if(FNR>1)k[$2]=1;next}
  FNR==1{print "SNP","A1","A2","BETA","SE","P","FREQ";next}
  ($9 in k) && $5+0>0 && $10+0>0 && $11+0>0 {b=log($5+0); se=(log($10+0)-log($11+0))/(2*1.959963984540054); if(se>0) print $9,$3,$4,b,se,$8,$7}' $E - \
  > data/gwas/fecd_eur_tspan10_eqtlgen_region.tsv
echo "FECD rows at TSPAN10 eQTLGen variants: multi-ancestry $(($(wc -l < data/gwas/fecd_tspan10_eqtlgen_region.tsv)-1)), European $(($(wc -l < data/gwas/fecd_eur_tspan10_eqtlgen_region.tsv)-1))"
