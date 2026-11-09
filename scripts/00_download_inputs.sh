#!/bin/bash
# ============================================================
# Step 00: download every public input and the external tools into $FECD_BASE.
#   Re-running skips anything already present. Files with a published or
#   recorded checksum are verified, so a fresh download is byte-identical to
#   the copy analysed for the manuscript.
#   Inputs fetched by later steps are not repeated here: glaucoma (45), IOP (50),
#   refractive error (55), GTEx v8 EUR FUSION weights (39), eQTLGen cis-eQTLs (83a),
#   FinnGen corneal endpoints (84b), eQTL Catalogue (100), GWAS Catalog (109a).
#   Disk: ~3 GB for this step; ~12 GB for all downloads plus the manual BESD files.
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -euo pipefail
mkdir -p "$FECD_BASE" && cd "$FECD_BASE"
mkdir -p data/gwas/raw data/predictdb data/predictdb_all data/ld_ref data/weights \
         data/endothelium_rnaseq data/magma data/besd software logs analysis results

md5of () { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
# fetch URL DEST [MD5]  — download once (resumable), then verify the checksum if one is given
fetch () {
  local url=$1 dest=$2 sum=${3:-}
  if [ ! -s "$dest" ]; then
    echo "  downloading $(basename "$dest")"
    curl -L --fail --retry 3 -C - -o "$dest.part" "$url" && mv "$dest.part" "$dest"
  else
    echo "  present     $(basename "$dest")"
  fi
  if [ -n "$sum" ] && [ "$(md5of "$dest")" != "$sum" ]; then
    echo "ERROR: checksum mismatch for $dest (expected $sum)"; exit 1
  fi
}

echo "=== 1. FECD GWAS (Verma et al. 2024, MVP; GWAS Catalog) ==="
FTP=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
fetch $FTP/GCST90480001-GCST90481000/GCST90480062/GCST90480062.tsv.gz data/gwas/raw/GCST90480062.tsv.gz e167ea1a27ea4e1bd313641b293a49eb
fetch $FTP/GCST90475001-GCST90476000/GCST90475865/GCST90475865.tsv.gz data/gwas/fecd_eur_GCST90475865.tsv.gz 17f5b6b9bcc46c6414075dddd5589877

echo "=== 2. GTEx v8 MASHR expression-prediction models (PredictDB, Zenodo 10.5281/zenodo.3518299) ==="
Z=https://zenodo.org/api/records/3518299/files
fetch $Z/mashr_eqtl.tar/content data/predictdb_all/mashr_eqtl.tar 87f3470bf2676043c748b684fb35fa7d
fetch $Z/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz/content \
      data/predictdb/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz dda0eedeb842cfc272e76ad432753d73
[ -d data/predictdb_all/eqtl/mashr ] || tar -xf data/predictdb_all/mashr_eqtl.tar -C data/predictdb_all
# the five tissues used by the 3-tissue FECD and glaucoma analyses
for t in Cells_Cultured_fibroblasts Skin_Not_Sun_Exposed_Suprapubic Esophagus_Mucosa Brain_Cortex Whole_Blood; do
  for ext in db txt.gz; do
    [ -e data/predictdb/mashr_$t.$ext ] || ln -s ../predictdb_all/eqtl/mashr/mashr_$t.$ext data/predictdb/mashr_$t.$ext
  done
done

echo "=== 3. FUSION LD reference (1000 Genomes Phase 3 EUR, HapMap3 SNPs) ==="
fetch https://alkesgroup.broadinstitute.org/FUSION/LDREF.tar.bz2 data/ld_ref/LDREF.tar.bz2
[ -s data/ld_ref/LDREF/1000G.EUR.22.bim ] || tar -xjf data/ld_ref/LDREF.tar.bz2 -C data/ld_ref

echo "=== 4. FUSION GTEx v8 weights, all ancestries, cultured fibroblasts (gene positions for coloc) ==="
fetch https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/ALL/GTExv8.ALL.Cells_Cultured_fibroblasts.tar.gz \
      data/weights/GTExv8.ALL.Cells_Cultured_fibroblasts.tar.gz 57270d1a3af4093d5699587731836ebe
[ -e data/weights/GTExv8.ALL.Cells_Cultured_fibroblasts.pos ] || tar -xzf data/weights/GTExv8.ALL.Cells_Cultured_fibroblasts.tar.gz -C data/weights

echo "=== 5. Corneal endothelium RNA-seq (GEO GSE112201, GSE101872) ==="
GEO=https://ftp.ncbi.nlm.nih.gov/geo/series
fetch $GEO/GSE112nnn/GSE112201/suppl/GSE112201_GeneCounts.tsv.gz data/endothelium_rnaseq/GSE112201_GeneCounts.tsv.gz 5046c0e88702da239dea86e4770f64a8
fetch $GEO/GSE101nnn/GSE101872/suppl/GSE101872_GeneCounts.tsv.gz data/endothelium_rnaseq/GSE101872_GeneCounts.tsv.gz 30dd16148029ebcd5ddaa2c9324778ee
fetch $GEO/GSE112nnn/GSE112201/matrix/GSE112201-GPL11154_series_matrix.txt.gz data/endothelium_rnaseq/GSE112201-GPL11154_series_matrix.txt.gz
fetch $GEO/GSE112nnn/GSE112201/matrix/GSE112201-GPL20301_series_matrix.txt.gz data/endothelium_rnaseq/GSE112201-GPL20301_series_matrix.txt.gz
fetch $GEO/GSE101nnn/GSE101872/matrix/GSE101872_series_matrix.txt.gz data/endothelium_rnaseq/GSE101872_series_matrix.txt.gz

echo "=== 6. External tools (pinned versions) ==="
command -v git >/dev/null || { echo "ERROR: git is required"; exit 1; }
if [ ! -s "$METAXCAN_DIR/SPrediXcan.py" ]; then        # MetaXcan v0.8.1-1 (commit 9b9f76b), unmodified
  git clone -q https://github.com/hakyimlab/MetaXcan software/MetaXcan
  git -C software/MetaXcan checkout -q 9b9f76b
fi
if [ ! -s "$FUSION_DIR/FUSION.assoc_test.R" ]; then    # FUSION e1ba5f7 + post-processing from 9346c12 + local patch
  git clone -q https://github.com/gusevlab/fusion_twas software/fusion_twas
  git -C software/fusion_twas checkout -q e1ba5f7
  for f in FUSION.post_process.R utils/plink_utils.R; do
    git -C software/fusion_twas show 9346c12:$f > software/fusion_twas/$f
  done
  patch -s software/fusion_twas/FUSION.assoc_test.R < "$SCRIPTS_DIR/../fusion_patch/FUSION.assoc_test.coloc.patch"
fi
# FUSION.post_process.R sources here("utils","plink_utils.R"), i.e. relative to the project folder
[ -s utils/plink_utils.R ] || { mkdir -p utils && cp "$FUSION_DIR/utils/plink_utils.R" utils/; }
if [ ! -x "$PLINK" ]; then                             # PLINK v1.90b7.2 (11 Dec 2023)
  case "$(uname -s)-$(uname -m)" in
    Darwin-*)       Z1=plink_mac_20231211.zip ;;
    Linux-x86_64)   Z1=plink_linux_x86_64_20231211.zip ;;
    *) echo "ERROR: download PLINK 1.9 for this platform to $PLINK"; exit 1 ;;
  esac
  mkdir -p "$(dirname "$PLINK")"
  fetch https://s3.amazonaws.com/plink1-assets/$Z1 "$(dirname "$PLINK")/$Z1"
  unzip -oq "$(dirname "$PLINK")/$Z1" plink -d "$(dirname "$PLINK")"
fi


echo ""
echo "All public inputs and tools are in place."
