#!/bin/bash
# ============================================================
# Step 100: fetch FULL cis summary statistics (all variant-gene pairs,
#           significant AND non-significant) for coloc.abf / coloc.susie.
#
# Source: eQTL Catalogue, release 6 = GTEx **v8** (r8 is GTEx v10 - do
#         NOT mix v10 eQTLs with v8 MASHR/FUSION weights).
#         Study QTS000015. Files are tabix-indexed, so a region is
#         fetched remotely: no 3.16 GB download per tissue.
#
# Usage:  bash 100_fetch_eqtl_catalogue_region.sh <ENSG> <chr> <start> <end> [outdir]
#   e.g.  bash 100_fetch_eqtl_catalogue_region.sh ENSG00000025434 11 46248300 48248400
#
# Requires tabix (htslib); set TABIX in config.local.sh if it is not on PATH
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -eu
ENSG=${1:?ENSG id required}; CHR=${2:?chr required}
START=${3:?start required}; END=${4:?end required}
OUT=${5:-data/eqtl_catalogue}
cd "$FECD_BASE"; mkdir -p "$OUT"

BASE=https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015
HDR="molecular_trait_id\tchromosome\tposition\tref\talt\tvariant\tma_samples\tmaf\tpvalue\tbeta\tse\ttype\tac\tan\tr2\tmolecular_trait_object_id\tgene_id\tmedian_tpm\trsid"

# dataset_id:label:N   (GTEx v8 eQTL sample sizes, from the release metadata;
# the 'an' column in the data equals 2N, so each file self-documents its N)
DATASETS="QTD000216:fibroblast:483 QTD000311:skin_not_sun_exposed:517 QTD000241:esophagus_mucosa:493 QTD000356:blood:670"

for D in $DATASETS; do
  ID=${D%%:*}; REST=${D#*:}; LAB=${REST%%:*}; N=${REST##*:}
  F="$OUT/${ENSG}_${LAB}_v8_allpairs.tsv"
  printf "%-24s " "$LAB"
  printf "$HDR\n" > "$F"
  if $TABIX "$BASE/$ID/$ID.all.tsv.gz" "${CHR}:${START}-${END}" 2>/dev/null \
       | awk -F'\t' -v g="$ENSG" '$1==g' >> "$F"; then
    NV=$(( $(wc -l < "$F") - 1 ))
    BEST=$(awk -F'\t' 'NR>1{if(b==""||$9<b){b=$9;r=$19}} END{print b" ("r")"}' "$F")
    AN=$(awk -F'\t' 'NR>1{print $14/2; exit}' "$F")
    printf "variants=%-6s N_from_data=%-5s best_P=%s\n" "$NV" "${AN:-NA}" "$BEST"
    [ "${AN:-0}" != "$N" ] && echo "    NOTE: N in data ($AN) != metadata N ($N)"
  else
    echo "FETCH FAILED"
  fi
done

cat <<'NOTE'

For coloc.susie with matched LD you do NOT need to supply LD yourself -
eQTL Catalogue publishes per-dataset SuSiE log-Bayes-factor files:
  https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000015/<QTD>/<QTD>.lbf_variable.txt.gz
Pair those with coloc::coloc.susie() (or coloc.bf_bf) on the eQTL side.
NOTE
