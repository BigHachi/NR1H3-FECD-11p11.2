#!/bin/bash
# ============================================================
# Step 112a: European LD panel for fine-mapping 11p11.2 (steps 112, 113).
#   1000 Genomes phase 3, GRCh38 biallelic SNV/indel release (20190312), chr11:46.2-48.3 Mb,
#   503 European samples, MAF >= 1%, variant IDs chr_pos_allele1_allele2 (alleles in ASCII order).
#   The HapMap3-only FUSION LD reference is too sparse for variant-level fine-mapping.
#   Output: data/ld_1000g_grch38/eur_chr11_46.2-48.3Mb.{bed,bim,fam}
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
cd "$FECD_BASE"
D=data/ld_1000g_grch38; mkdir -p $D; cd $D
U=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr11.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
[ -s integrated_call_samples_v3.20130502.ALL.panel ] || \
  curl -s --fail -o integrated_call_samples_v3.20130502.ALL.panel http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
awk '$3=="EUR"{print $1, $1}' integrated_call_samples_v3.20130502.ALL.panel > eur_keep.txt
[ -s chr11_46.2-48.3Mb.vcf.gz ] || { "$TABIX" -h $U 11:46200000-48300000 | gzip > chr11_46.2-48.3Mb.vcf.gz; rm -f ALL.chr11*.tbi; }
"$PLINK" --vcf chr11_46.2-48.3Mb.vcf.gz --double-id --keep eur_keep.txt --maf 0.01 --biallelic-only strict \
  --set-missing-var-ids '@_#_$1_$2' --make-bed --out eur_chr11_46.2-48.3Mb --silent
echo "European samples: $(wc -l < eur_keep.txt); variants (MAF >= 1%): $(wc -l < eur_chr11_46.2-48.3Mb.bim)"

# European GWAS for the same region (step 112, 113)
cd "$FECD_BASE"
mkdir -p analysis
gzip -cd data/gwas/fecd_eur_GCST90475865.tsv.gz | awk -F'\t' 'NR==1 || ($1==11 && $2>=46200000 && $2<=48300000)' \
  > analysis/fecd_eur_11p11_region.tsv
echo "European GWAS variants in region: $(( $(wc -l < analysis/fecd_eur_11p11_region.tsv) - 1 ))"
