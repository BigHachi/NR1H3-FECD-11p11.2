#!/bin/bash
# ============================================================
# Step 39: Download GTEx v8 EUR-only FUSION weights
#
# Motivation: the original FECD/BK TWAS used ALL-ancestry weights
# (GTExv8.ALL.*) with a European-only 1000G LD reference and a
# European-dominant GWAS.  Using ALL-ancestry eQTL models with an
# EUR LD panel can introduce bias because LD between instruments
# differs across ancestries.  This script downloads EUR-only weights
# for the three analysis tissues so that script 40 can rerun the
# FUSION TWAS with an ancestry-matched model + LD reference.
#
# Tissues:
#   Cells_Cultured_fibroblasts        (primary — used in original TWAS)
#   Esophagus_Mucosa                  (used in S-PrediXcan/MultiXcan)
#   Skin_Not_Sun_Exposed_Suprapubic   (used in S-PrediXcan/MultiXcan)
#
# EUR weights URL pattern:
#   https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.{TISSUE}.tar.gz
#
# NOTE: Not all tissues have EUR-specific weights.  This script
# checks each download and reports which tissues are available.
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -euo pipefail
cd "$FECD_BASE"
mkdir -p data/weights_eur logs

LOG=logs/39_download_eur_weights.log
exec > >(tee -a "$LOG") 2>&1

echo "=== Download EUR FUSION weights: $(date) ==="
echo ""

BASE_URL="https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR"

TISSUES=(
  "Cells_Cultured_fibroblasts"
  "Esophagus_Mucosa"
  "Skin_Not_Sun_Exposed_Suprapubic"
)

SUCCESS=()
FAILED=()

for TISSUE in "${TISSUES[@]}"; do
  TAR="GTExv8.EUR.${TISSUE}.tar.gz"
  URL="${BASE_URL}/${TAR}"
  OUT="data/weights_eur/${TAR}"
  DEST="data/weights_eur"

  echo "--- ${TISSUE} ---"

  # Skip if already extracted (pos file present)
  POS="data/weights_eur/GTExv8.EUR.${TISSUE}.pos"
  if [[ -f "$POS" ]]; then
    echo "  Already extracted: ${POS} — skipping download"
    SUCCESS+=("${TISSUE}")
    continue
  fi

  # Try download (HEAD first to check existence without downloading full file)
  HTTP_CODE=$(curl -s -o /dev/null -w "%{http_code}" --head "${URL}")
  if [[ "$HTTP_CODE" != "200" ]]; then
    echo "  WARNING: EUR weights not available for ${TISSUE} (HTTP ${HTTP_CODE})"
    echo "  URL tried: ${URL}"
    FAILED+=("${TISSUE}")
    continue
  fi

  echo "  Downloading ${TAR}..."
  curl -L --fail --retry 3 -o "${OUT}" "${URL}"

  echo "  Extracting..."
  tar -xzf "${OUT}" -C "${DEST}/"
  rm "${OUT}"   # remove tar after extraction to save disk

  if [[ -f "$POS" ]]; then
    N_GENES=$(wc -l < "$POS")
    echo "  OK — ${N_GENES} gene models: ${POS}"
    SUCCESS+=("${TISSUE}")
  else
    echo "  WARNING: pos file not found after extraction — check tar structure"
    FAILED+=("${TISSUE}")
  fi
  echo ""
done

echo "======================================"
echo "Download summary"
echo "======================================"
echo ""
echo "Available EUR tissues (${#SUCCESS[@]}):"
for T in "${SUCCESS[@]}"; do
  POS="data/weights_eur/GTExv8.EUR.${T}.pos"
  N=$( [[ -f "$POS" ]] && wc -l < "$POS" || echo "?" )
  echo "  ${T}: ${N} gene models"
done

if [[ ${#FAILED[@]} -gt 0 ]]; then
  echo ""
  echo "NOT available as EUR-only (${#FAILED[@]}):"
  for T in "${FAILED[@]}"; do
    echo "  ${T}  — will fall back to ALL weights for this tissue"
  done
  echo ""
  echo "  Fallback strategy (script 40):"
  echo "  For tissues without EUR weights, the ALL-ancestry weights are"
  echo "  used with a note that this is a limitation.  Results for those"
  echo "  tissues should be interpreted cautiously."
fi

echo ""
echo "EUR weights complete: $(date)"
echo "Next step: bash scripts/40_run_eur_twas.sh"
