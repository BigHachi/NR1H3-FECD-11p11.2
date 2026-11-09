#!/bin/bash
# ============================================================
# FECD Analysis - Step 17: Prepare FECD GWAS for S-PrediXcan
# Reformats step 11's output, keeping its panel_variant_id (chr_pos_A2_A1_b38).
# That ID matches a model varID only when the effect allele is ALT, so step 95a
# re-matches every variant to each model by position and allele pair; this is the
# file that produced the reported results. (An rsID->varID remap tried here earlier
# was superseded by 95a and dropped one indel, so it has been removed.)
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -e
cd "$FECD_BASE"

IN=data/gwas/fecd_gwas_full.tsv.gz
OUT=data/gwas/fecd_gwas_sprediXcan.txt.gz

[ -s "$IN" ] || { echo "ERROR: $IN not found (run step 11 first)"; exit 1; }
echo "=== Preparing FECD GWAS for S-PrediXcan ==="
echo "Input:  $IN"

# Build the S-PrediXcan input
echo "Converting FECD GWAS to S-PrediXcan format..."
$METAXCAN_PY - <<'PYEOF'
import gzip

in_file  = "data/gwas/fecd_gwas_full.tsv.gz"
out_file = "data/gwas/fecd_gwas_sprediXcan.txt.gz"

# Columns in fecd_gwas_full.tsv.gz:
# SNP(0), CHR(1), BP(2), A1(3), A2(4), BETA(5), SE(6), OR(7), FREQ(8), P(9),
# ZSCORE(10), panel_variant_id(11)

header = "panel_variant_id\trsid\tchromosome\tposition\teffect_allele\tnon_effect_allele\tzscore\tpvalue\tsample_size\teffect_size\tstandard_error\tfrequency\n"

rows_written = 0
seen_pvids = set()

with gzip.open(in_file, "rt") as fin, gzip.open(out_file, "wt") as fout:
    fout.write(header)
    next(fin)  # skip input header
    for line in fin:
        parts = line.rstrip("\n").split("\t")
        snp  = parts[0]   # rsid
        chrom = parts[1]
        bp    = parts[2]
        a1    = parts[3]   # effect allele
        a2    = parts[4]   # non-effect allele
        beta  = parts[5]
        se    = parts[6]
        freq  = parts[8]
        pval  = parts[9]
        zscore = parts[10]
        pvid  = parts[11]  # original panel_variant_id

        # Filter: SE > 0 and ZSCORE is finite
        try:
            se_f = float(se)
            zs_f = float(zscore)
        except (ValueError, TypeError):
            continue
        if se_f <= 0 or not (-1e38 < zs_f < 1e38):
            continue

        # Deduplicate on panel_variant_id
        if pvid in seen_pvids:
            continue
        seen_pvids.add(pvid)

        fout.write(f"{pvid}\t{snp}\t{chrom}\t{bp}\t{a1}\t{a2}\t{zscore}\t{pval}\t570321\t{beta}\t{se}\t{freq}\n")
        rows_written += 1

print(f"Rows written: {rows_written}")
PYEOF

echo ""
echo "=== FECD S-PrediXcan prep complete ==="
echo "Output: $(du -sh $OUT | cut -f1)  ($OUT)"
[ -s "$OUT" ] || { echo "ERROR: $OUT was not written"; exit 1; }
