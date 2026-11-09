#!/bin/bash
# Step 75b (review point 3): glaucoma S-MultiXcan WITHOUT brain cortex.
#   NR1H3 is not a GTEx v8 Brain_Cortex eGene, so its cortex prediction is
#   borrowed through mashr sharing. Panel = fibroblast + skin (the step-75
#   panel minus Brain_Cortex). Same inputs/settings as step 75 otherwise.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -euo pipefail
cd "$FECD_BASE"
METAXCAN="$METAXCAN_DIR"
PREDICTDB=data/predictdb
STAGED=results/sprediXcan_glaucoma_nobrain_staged; mkdir -p $STAGED; rm -f $STAGED/*.csv
cp results/sprediXcan_glaucoma_brain_staged/mashr_Cells_Cultured_fibroblasts.csv $STAGED/
cp results/sprediXcan_glaucoma_brain_staged/mashr_Skin_Not_Sun_Exposed_Suprapubic.csv $STAGED/
$METAXCAN_PY ${METAXCAN}/SMulTiXcan.py \
  --models_folder "$PREDICTDB" \
  --models_name_filter "mashr_Cells_Cultured_fibroblasts.db|mashr_Skin_Not_Sun_Exposed_Suprapubic.db" \
  --models_name_pattern "(.*).db" \
  --snp_covariance ${PREDICTDB}/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
  --metaxcan_folder $STAGED --metaxcan_filter "mashr_.*csv" --metaxcan_file_name_parse_pattern "(mashr_.*)\\.csv" \
  --gwas_file data/gwas/glaucoma_gwas_sprediXcan.txt.gz --snp_column panel_variant_id \
  --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele \
  --beta_column effect_size --se_column standard_error --pvalue_column pvalue --freq_column frequency \
  --keep_non_rsid --trimmed_ensemble_id --model_db_snp_key varID --cutoff_condition_number 30 \
  --output results/multixcan/glaucoma_smultixcan_nobrain.txt > logs/75b_glaucoma_nobrain.log 2>&1
echo done
