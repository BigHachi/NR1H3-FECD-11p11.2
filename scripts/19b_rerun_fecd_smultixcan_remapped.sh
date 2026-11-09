#!/bin/bash
# ============================================================
# Step 19b: Re-run 3-tissue FECD S-MultiXcan with the REMAPPED GWAS
#   Step 19 passed the original GWAS (data/gwas/fecd_gwas_sprediXcan.txt.gz)
#   as --gwas_file. Its panel_variant_id orientation matches only ~33% of
#   fibroblast model varIDs, so S-MultiXcan dropped/mis-specified genes
#   (15,854 of 16,032 rows with a P; 1,046 fibroblast genes missing).
#   z-scores AND --gwas_file must come from the SAME GWAS file: mixing the old
#   per-tissue z (results/sprediXcan/) with the remapped GWAS inflated S-MultiXcan
#   (genes with max single-tissue |z|<1.7 at P~1e-10; kept as
#   fecd_smultixcan_remapped.INVALID_mixed_inputs.txt). So the per-tissue inputs
#   are the remapped step-95 outputs for the 3 proxy tissues.
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -eu
cd "$FECD_BASE"
METAXCAN="$METAXCAN_DIR"
PREDICTDB=data/predictdb
GWAS=data/gwas/fecd_alltissue_sprediXcan.txt.gz
STAGED=results/sprediXcan_fecd_remap_staged
OUT=results/multixcan/fecd_smultixcan_remapped.txt
COVAR=${PREDICTDB}/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz
mkdir -p "$STAGED"
for T in Esophagus_Mucosa Skin_Not_Sun_Exposed_Suprapubic Cells_Cultured_fibroblasts; do
  cp results/sprediXcan_alltissue_remap/fecd_${T}.csv "$STAGED/mashr_${T}.csv"
done
rm -f "$OUT"
$METAXCAN_PY ${METAXCAN}/SMulTiXcan.py \
  --models_folder "$PREDICTDB" \
  --models_name_filter "mashr_Esophagus_Mucosa.db|mashr_Skin_Not_Sun_Exposed_Suprapubic.db|mashr_Cells_Cultured_fibroblasts.db" \
  --models_name_pattern "(.*).db" --snp_covariance "$COVAR" \
  --metaxcan_folder "$STAGED" --metaxcan_filter "mashr_.*csv" --metaxcan_file_name_parse_pattern "(mashr_.*)\\.csv" \
  --gwas_file "$GWAS" --snp_column panel_variant_id --effect_allele_column effect_allele \
  --non_effect_allele_column non_effect_allele --beta_column effect_size --se_column standard_error \
  --pvalue_column pvalue --freq_column frequency --keep_non_rsid --trimmed_ensemble_id \
  --model_db_snp_key varID --cutoff_condition_number 30 --output "$OUT" > logs/fecd_smultixcan_remapped.log 2>&1
echo "done: $OUT"
