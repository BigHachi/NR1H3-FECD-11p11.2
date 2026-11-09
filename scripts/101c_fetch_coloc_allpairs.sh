#!/bin/bash
# ============================================================
# Step 101c: as 101, gene list from $GENES (default genes.tsv). Used with genes_coloc.tsv
#   = Table 2 genes + every gene with a PP4 in Supp Table S6 (review point 2).
#           for the 21 FECD EUR S-MultiXcan genes ranked 1..NR1H3 (rank 21),
#           in the 3 TWAS tissues, from eQTL Catalogue r6 (= GTEx v8).
#   Region = TSS +/- 1.05 Mb, filtered to the gene's molecular_trait_id,
#   so the variant set is exactly the cis window eQTL Catalogue tested.
#   Gene list/coords: data/tophits21/genes.tsv (Ensembl REST GRCh38;
#   PANO1 ENSG00000274897 is retired in current Ensembl -> approx TSS; its eQTL Catalogue ID is mapped below).
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, TABIX
set -eu
cd "$FECD_BASE"
OUT=data/tophits21/allpairs; mkdir -p "$OUT"
BASE=https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015
HDR="molecular_trait_id\tchromosome\tposition\tref\talt\tvariant\tma_samples\tmaf\tpvalue\tbeta\tse\ttype\tac\tan\tr2\tmolecular_trait_object_id\tgene_id\tmedian_tpm\trsid"

fetch() {  # ensg gene chr tss qtd label
  set -o pipefail   # a failed or interrupted tabix download must be retried, never saved as "no data"
  local F="$OUT/$2_$6.tsv"
  [ -s "$F" ] && [ "$(wc -l < "$F")" -gt 1 ] && return 0   # resumable; header-only files are re-checked
  local s=$(( $4 - 1050000 )); [ $s -lt 1 ] && s=1
  local e=$(( $4 + 1050000 ))
  local tmp="$F.tmp"
  for try in 1 2 3 4 5; do
    printf "$HDR\n" > "$tmp"
    if (cd "$OUT" && $TABIX "$BASE/$5/$5.all.tsv.gz" "$3:$s-$e") 2>/dev/null \
         | awk -F'\t' -v g="$1" '$1==g' >> "$tmp"; then
      mv "$tmp" "$F"
      echo "$2 $6 variants=$(( $(wc -l < "$F") - 1 ))"   # 0 = gene not tested in that tissue
      return 0
    fi
    sleep $(( try * 15 ))
  done
  rm -f "$tmp"; echo "FAILED: $2 $6 (tabix failed 5 times)" >&2; return 1
}
export -f fetch; export OUT BASE TABIX HDR
for q in QTD000216 QTD000311 QTD000241; do
  for try in 1 2 3; do (cd "$OUT" && $TABIX -l "$BASE/$q/$q.all.tsv.gz" > /dev/null 2>&1) && break; sleep 15; done
done

tail -n +2 ${GENES:-data/tophits21/genes.tsv} | while IFS=$'\t' read rank ensg gene p chr st en strand tss bt; do
  # PANO1 is ENSG00000274897 in GENCODE v26 (GTEx models, S-MultiXcan) but ENSG00000288675 in eQTL Catalogue r6
  # (formerly fetched separately by step 101b)
  [ "$ensg" = ENSG00000274897 ] && ensg=ENSG00000288675
  for T in QTD000216:fibroblast QTD000311:skin_not_sun_exposed QTD000241:esophagus_mucosa; do
    echo "$ensg $gene $chr $tss ${T%%:*} ${T##*:}"
  done
done | xargs -P 8 -n 6 bash -c 'fetch "$@"' _
