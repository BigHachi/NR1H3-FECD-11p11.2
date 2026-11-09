#!/bin/bash
# ============================================================
# Step 83a (scripted Oct 10 2026): eQTLGen whole-blood cis-eQTLs for NR1H3
#   eQTLGen phase I cis-eQTLs (FDR < 0.05, probe level; N up to 31,684), filtered to
#   ENSG00000025434. Output data/eqtl_eqtlgen/nr1h3_eqtlgen.txt is the exposure file
#   for every cis-MR step (83, 88, 89b, 89c, 91) and SMR (90-93).
#   Verified Oct 10 2026: byte-identical to the file used for the manuscript.
#   Also writes data/eqtl_eqtlgen/nr1h3_eqtlgen_af.tsv (eQTLGen allele frequencies, 10 Nov 2026).
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -euo pipefail
cd "$FECD_BASE"; mkdir -p data/eqtl_eqtlgen
F=data/eqtl_eqtlgen/eqtlgen_cis_FDR05.txt.gz
[ -s $F ] || curl -L --fail -o $F \
  https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
gzip -cd $F | awk -F'\t' 'NR==1 || $8=="ENSG00000025434"' > data/eqtl_eqtlgen/nr1h3_eqtlgen.txt
echo "NR1H3 cis-eQTL rows: $(($(wc -l < data/eqtl_eqtlgen/nr1h3_eqtlgen.txt)-1))"

# eQTLGen allele frequencies (10 Nov 2026, review): AlleleB frequency over the eQTLGen cohorts with
# genotype counts. Every cis-MR step standardises the exposure Z with the frequency of the assessed allele
# in this file (exposure source); up to release v13 the outcome GWAS frequency was used instead.
AF=data/eqtl_eqtlgen/eqtlgen_SNP_AF_AlleleB.txt.gz
[ -s $AF ] || curl -L --fail -o $AF \
  https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz
gzip -cd $AF | awk -F'\t' 'NR==FNR{if(FNR>1)k[$2]=1;next} FNR==1||($1 in k)' data/eqtl_eqtlgen/nr1h3_eqtlgen.txt - \
  > data/eqtl_eqtlgen/nr1h3_eqtlgen_af.tsv
echo "eQTLGen allele frequencies for NR1H3 cis-eQTL variants: $(($(wc -l < data/eqtl_eqtlgen/nr1h3_eqtlgen_af.tsv)-1))"

# FECD GWAS rows for those eQTLGen variants (outcome files for cis-MR: 83, 88 multi-ancestry; 91 European).
# European BETA/SE from the OR and 95% CI (SE = (ln CIU - ln CIL) / (2 x qnorm(0.975))).
E=data/eqtl_eqtlgen/nr1h3_eqtlgen.txt
gzip -cd data/gwas/fecd_gwas_full.tsv.gz | awk -F'\t' -v OFS='\t' 'NR==FNR{if(FNR>1)k[$2]=1;next}
  FNR==1{print "SNP","A1","A2","BETA","SE","P","FREQ";next} ($1 in k){print $1,$4,$5,$6,$7,$10,$9}' $E - \
  > data/gwas/fecd_nr1h3_eqtlgen_region.tsv
gzip -cd data/gwas/fecd_eur_GCST90475865.tsv.gz | awk -F'\t' -v OFS='\t' 'NR==FNR{if(FNR>1)k[$2]=1;next}
  FNR==1{print "SNP","A1","A2","BETA","SE","P","FREQ";next}
  ($9 in k) && $5+0>0 && $10+0>0 && $11+0>0 {b=log($5+0); se=(log($10+0)-log($11+0))/(2*1.959963984540054); if(se>0) print $9,$3,$4,b,se,$8,$7}' $E - \
  > data/gwas/fecd_eur_nr1h3_eqtlgen_region.tsv
echo "FECD rows at NR1H3 eQTLGen variants: multi-ancestry $(($(wc -l < data/gwas/fecd_nr1h3_eqtlgen_region.tsv)-1)), European $(($(wc -l < data/gwas/fecd_eur_nr1h3_eqtlgen_region.tsv)-1))"
