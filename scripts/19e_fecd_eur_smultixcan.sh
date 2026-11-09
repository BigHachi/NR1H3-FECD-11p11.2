#!/bin/bash
# ============================================================
# Step 19e (30 Oct 2026): European-ancestry FECD S-PrediXcan in the three proxy tissues + 3-tissue
#   S-MultiXcan, with the same commands as the multi-ancestry analysis (steps 95 and 19b). Source of the
#   S_MultiXcan_rank column of Supplementary Table S6 (step 19f).
#   Input : data/gwas/fecd_eur_alltissue_sprediXcan.txt.gz (step 19d, then step 95a: model variants matched by
#           chromosome, GRCh38 position and allele pair)
#   Output: results/sprediXcan_fecd_eur_remap_staged/mashr_{tissue}.csv,
#           results/multixcan/fecd_eur_smultixcan_remapped.txt
#   (The 4 Jul 2026 run matched model variants by rsID; its outputs, results/sprediXcan_fecd_eur_staged/ and
#   results/multixcan/fecd_eur_smultixcan.txt, are kept for reference.)
#   S-PrediXcan does not overwrite an existing output file (it skips the tissue), so outputs are removed first.
#   EUR_STAGED / EUR_SMX override the output paths (e.g. to rerun without replacing the reported outputs).
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, METAXCAN_DIR, METAXCAN_PY
set -eu
cd "$FECD_BASE"
GWAS=data/gwas/fecd_eur_alltissue_sprediXcan.txt.gz
P=data/predictdb
ST=${EUR_STAGED:-results/sprediXcan_fecd_eur_remap_staged}
OUT=${EUR_SMX:-results/multixcan/fecd_eur_smultixcan_remapped.txt}
mkdir -p "$ST" "$(dirname "$OUT")"
for T in Cells_Cultured_fibroblasts Skin_Not_Sun_Exposed_Suprapubic Esophagus_Mucosa; do
  rm -f "$ST/mashr_$T.csv"
  $METAXCAN_PY "$METAXCAN_DIR/SPrediXcan.py" --model_db_path $P/mashr_$T.db --covariance $P/mashr_$T.txt.gz \
    --gwas_file $GWAS --snp_column panel_variant_id --effect_allele_column effect_allele \
    --non_effect_allele_column non_effect_allele --beta_column effect_size --se_column standard_error \
    --pvalue_column pvalue --freq_column frequency --keep_non_rsid --model_db_snp_key varID --throw \
    --remove_ens_version --output_file "$ST/mashr_$T.csv"
done
rm -f "$OUT"
$METAXCAN_PY "$METAXCAN_DIR/SMulTiXcan.py" --models_folder $P \
  --models_name_filter "mashr_Esophagus_Mucosa.db|mashr_Skin_Not_Sun_Exposed_Suprapubic.db|mashr_Cells_Cultured_fibroblasts.db" \
  --models_name_pattern "(.*).db" --snp_covariance $P/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
  --metaxcan_folder "$ST" --metaxcan_filter "mashr_.*csv" --metaxcan_file_name_parse_pattern "(mashr_.*)\\.csv" \
  --gwas_file $GWAS --snp_column panel_variant_id --effect_allele_column effect_allele \
  --non_effect_allele_column non_effect_allele --beta_column effect_size --se_column standard_error \
  --pvalue_column pvalue --freq_column frequency --keep_non_rsid --trimmed_ensemble_id \
  --model_db_snp_key varID --cutoff_condition_number 30 --output "$OUT"
echo "done: $OUT"
