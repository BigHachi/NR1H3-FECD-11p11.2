#!/bin/bash
# ============================================================
# Glaucoma Analysis - Step 46: Prepare GWAS for S-PrediXcan
# Uses rsid->varID remapping (same approach as FECD, script 17)
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -e
cd "$FECD_BASE"

IN=data/gwas/glaucoma_gwas_full.tsv
OUT=data/gwas/glaucoma_gwas_sprediXcan.txt.gz
PREDICTDB=data/predictdb

echo "=== Preparing Glaucoma GWAS for S-PrediXcan ==="
echo "Input:  $IN"

if [ ! -f "$IN" ]; then
  echo "ERROR: $IN not found — run 45_download_glaucoma_gwas.sh first"
  exit 1
fi

# Step 1: Build rsid->varID mapping from prediction model DBs
echo "Building rsid->varID mapping from prediction models..."
cat > /tmp/glaucoma_build_map.py << 'PYEOF'
import sqlite3, json, os

predictdb = "data/predictdb"
rsid_to_varid = {}

for f in os.listdir(predictdb):
    if not f.endswith(".db"):
        continue
    path = os.path.join(predictdb, f)
    con = sqlite3.connect(path)
    rows = con.execute("SELECT rsid, varID FROM weights").fetchall()
    con.close()
    for rsid, varid in rows:
        if rsid not in rsid_to_varid:
            rsid_to_varid[rsid] = varid

print(f"rsid->varID entries: {len(rsid_to_varid)}")
with open("/tmp/glaucoma_rsid_varid_map.json", "w") as fh:
    json.dump(rsid_to_varid, fh)
print("Map saved to /tmp/glaucoma_rsid_varid_map.json")
PYEOF
$METAXCAN_PY /tmp/glaucoma_build_map.py

echo ""

# Step 2: Convert to S-PrediXcan format with varID remap
echo "Converting to S-PrediXcan format..."
cat > /tmp/glaucoma_convert.py << 'PYEOF'
import gzip, json, sys

with open("/tmp/glaucoma_rsid_varid_map.json") as fh:
    rsid_to_varid = json.load(fh)
print(f"Loaded {len(rsid_to_varid)} rsid->varID mappings")

in_file  = "data/gwas/glaucoma_gwas_full.tsv"
out_file = "data/gwas/glaucoma_gwas_sprediXcan.txt.gz"

# Input columns: SNP CHR BP A1 A2 BETA SE P FREQ
# A1=effect/alt, A2=non-effect/ref (FinnGen convention)
SAMPLE_SIZE = 500348  # total N: 26,591 cases + 473,757 controls (FinnGen R12)

header = "panel_variant_id\trsid\tchromosome\tposition\teffect_allele\tnon_effect_allele\tzscore\tpvalue\tsample_size\teffect_size\tstandard_error\tfrequency\n"

rows_written  = 0
rows_remapped = 0
seen_pvids    = set()

with open(in_file, "r") as fin, gzip.open(out_file, "wt") as fout:
    fout.write(header)
    next(fin)  # skip header
    for line in fin:
        parts = line.rstrip("\n").split("\t")
        snp   = parts[0]
        chrom = parts[1]
        bp    = parts[2]
        a1    = parts[3]
        a2    = parts[4]
        beta  = parts[5]
        se    = parts[6]
        pval  = parts[7]
        freq  = parts[8]

        try:
            se_f   = float(se)
            beta_f = float(beta)
            zscore = beta_f / se_f
        except (ValueError, ZeroDivisionError):
            continue
        if se_f <= 0 or not (-1e38 < zscore < 1e38):
            continue

        pvid = f"chr{chrom}_{bp}_{a2}_{a1}_b38"
        if snp in rsid_to_varid:
            pvid = rsid_to_varid[snp]
            rows_remapped += 1

        if pvid in seen_pvids:
            continue
        seen_pvids.add(pvid)

        fout.write(f"{pvid}\t{snp}\t{chrom}\t{bp}\t{a1}\t{a2}\t{zscore:.6f}\t{pval}\t{SAMPLE_SIZE}\t{beta}\t{se}\t{freq}\n")
        rows_written += 1

print(f"Rows written:          {rows_written}")
print(f"Rows with varID remap: {rows_remapped}")
PYEOF
$METAXCAN_PY /tmp/glaucoma_convert.py

echo ""
echo "=== Glaucoma S-PrediXcan prep complete ==="
echo "Output: $(du -sh $OUT | cut -f1)  ($OUT)"
