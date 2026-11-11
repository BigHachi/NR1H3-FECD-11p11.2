#!/bin/bash
# ============================================================
# IOP Analysis - Step 50: Convert IOP GWAS + liftover GRCh37->GRCh38
# Source  : International Glaucoma Genetics Consortium IOP meta-analysis released with
#           Bonnemaijer et al., Commun Biol 2019 (GCST009413, PMID 31798171):
#           N = 31,269 European-ancestry individuals (GWAS Catalog; = maximum per-variant Ntot)
# Download: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009413/Bonnemaijer_Intraocular_pressure_pmid.txt.gz
# Uses pyliftover (streaming Python) for minimal RAM usage
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -e
cd "$FECD_BASE"

GWAS_RAW=data/gwas/iop_gwas_raw.txt.gz
GWAS_FULL=data/gwas/iop_gwas_full.tsv
CHAIN=data/liftover/hg19ToHg38.over.chain.gz

mkdir -p data/liftover

if [ ! -f "$GWAS_RAW" ]; then
  echo "=== Downloading IOP GWAS (Bonnemaijer 2019, GCST009413) ==="
  curl -L --fail --retry 3 -o "$GWAS_RAW" \
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009413/Bonnemaijer_Intraocular_pressure_pmid.txt.gz"
  echo "  Downloaded: $(du -sh $GWAS_RAW | cut -f1)"
else
  echo "  IOP GWAS already downloaded: $(du -sh $GWAS_RAW | cut -f1)"
fi

if [ ! -f "$CHAIN" ]; then
  echo "=== Downloading hg19->hg38 chain file ==="
  curl -L --fail --retry 3 -o "$CHAIN" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
  echo "  Downloaded: $(du -sh $CHAIN | cut -f1)"
else
  echo "  Chain file already present"
fi

echo ""
echo "=== Converting IOP GWAS + liftover GRCh37->GRCh38 (streaming Python) ==="

cat > /tmp/iop_liftover_convert.py << 'PYEOF'
import gzip, sys
from pyliftover import LiftOver

print("Loading hg19->hg38 chain (one-time cost ~200MB)...")
lo = LiftOver("data/liftover/hg19ToHg38.over.chain.gz")
print("Chain loaded.")

# Input columns: MarkerName CHR POS SNP EffectAllele NonEffectAllele
#                MAF MAFSE MinMAF MaxMAF Effect StdErr P-value ...
in_file  = "data/gwas/iop_gwas_raw.txt.gz"
out_file = "data/gwas/iop_gwas_full.tsv"

written = unmapped = skipped = 0
seen_snps = set()

with gzip.open(in_file, "rt") as fin, open(out_file, "w") as fout:
    fout.write("SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\tFREQ\n")
    header = next(fin).strip().split("\t")
    col = {c: i for i, c in enumerate(header)}

    for i, line in enumerate(fin):
        if i % 1000000 == 0 and i > 0:
            print(f"  {i:,} rows processed | written={written:,} unmapped={unmapped:,} skipped={skipped:,}")

        parts = line.rstrip("\n").split("\t")
        try:
            chrom_raw = str(int(parts[col["CHR"]]))
            pos_37    = int(parts[col["POS"]])
            snp       = parts[col["SNP"]]
            a1        = parts[col["EffectAllele"]].upper()
            a2        = parts[col["NonEffectAllele"]].upper()
            beta      = float(parts[col["Effect"]])
            se        = float(parts[col["StdErr"]])
            pval      = parts[col["P-value"]]
            freq      = parts[col["MAF"]]
        except (ValueError, IndexError, KeyError):
            skipped += 1
            continue

        if se <= 0 or not (-1e38 < beta/se < 1e38):
            skipped += 1
            continue
        if int(chrom_raw) not in range(1, 23):
            skipped += 1
            continue

        # Liftover position
        result = lo.convert_coordinate(f"chr{chrom_raw}", pos_37 - 1)  # pyliftover is 0-based
        if not result:
            unmapped += 1
            continue

        chrom_38 = result[0][0].replace("chr", "")
        pos_38   = result[0][1] + 1  # back to 1-based

        # Deduplicate on SNP id
        if snp in seen_snps:
            skipped += 1
            continue
        seen_snps.add(snp)

        fout.write(f"{snp}\t{chrom_38}\t{pos_38}\t{a1}\t{a2}\t{beta}\t{se}\t{pval}\t{freq}\n")
        written += 1

print(f"\nDone.")
print(f"  Written:  {written:,}")
print(f"  Unmapped: {unmapped:,}")
print(f"  Skipped:  {skipped:,}")
PYEOF

$METAXCAN_PY /tmp/iop_liftover_convert.py

echo ""
echo "=== IOP conversion complete ==="
echo "Output: $GWAS_FULL  ($(wc -l < $GWAS_FULL) lines)"
