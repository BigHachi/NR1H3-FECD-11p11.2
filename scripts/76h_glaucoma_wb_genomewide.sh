#!/bin/bash
# ============================================================
# Step 76h (29 Oct 2026): genome-wide FUSION whole-blood TWAS + conditional/joint analysis for
#   glaucoma (FinnGen R12 H7_GLAUCOMA), with the input as reported and after the allele-matching
#   correction of step 55c applied to every chromosome. Recomputes the genome-wide ranking of the
#   conditionally independent signals that step 76f reads from the June 2025 workstation run.
#   Weights: GTEx v8 EUR Whole_Blood FUSION (step 76d). LD: 1000G EUR (FUSION LDREF).
#   1. per-chromosome inputs (76h_glaucoma_wb_genomewide_inputs.py)
#   2. FUSION.assoc_test per chromosome, both inputs
#   3. FUSION.post_process per chromosome (MHC models included) on the models with TWAS P below
#      the genome-wide Bonferroni threshold (0.05 / models with a TWAS P), 100-kb locus window
#   4. ranking of the jointly included genes and comparison with the June 2025 run
#      (76h_glaucoma_wb_genomewide_rank.py)
#   10 Nov 2026: allele_corrected/ is the input reported from release v14 (step 60b's allele-matched file);
#   as_reported/ keeps its name but is the first-record input reported up to release v13 (step 60).
#   Outputs: results/twas_wb_genomewide/{as_reported,allele_corrected}/
#            analysis/glaucoma_wb_genomewide_joint_ranking.tsv
#            analysis/glaucoma_wb_genomewide_rerun_summary.tsv
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
cd "$FECD_BASE"
FU="$FUSION_DIR"
WGT_DIR=data/weights_eur
LD=data/ld_ref/LDREF/1000G.EUR.
OUT=results/twas_wb_genomewide
JOBS=${JOBS:-4}
mkdir -p "$OUT" logs

# ---- 1. inputs
if [ ! -f "$OUT/inputs/replaced_records.tsv" ]; then
  "${PY:-python3}" "$SCRIPTS_DIR/76h_glaucoma_wb_genomewide_inputs.py" "$OUT/inputs" > logs/76h_inputs.log 2>&1
fi

# ---- 2. assoc_test (whole-blood weights of one chromosome at a time)
assoc() {   # $1 = input set, $2 = chromosome
  local d="$OUT/$1"
  mkdir -p "$d"
  head -1 "$WGT_DIR/GTExv8.EUR.Whole_Blood.pos" > "$d/wb_chr$2.pos"
  awk -F'\t' -v c="$2" '$4==c' "$WGT_DIR/GTExv8.EUR.Whole_Blood.pos" >> "$d/wb_chr$2.pos"
  Rscript "$FU/FUSION.assoc_test.R" --sumstats "$OUT/inputs/$1/glaucoma_fusion_chr$2.tsv" \
    --weights "$d/wb_chr$2.pos" --weights_dir "$WGT_DIR" --ref_ld_chr "$LD" --chr "$2" \
    --out "$d/glaucoma_wb_chr$2.dat" > "$d/assoc_chr$2.log" 2>&1
}
export -f assoc
export OUT WGT_DIR FU LD
for S in as_reported allele_corrected; do for C in $(seq 1 22); do
  [ -s "$OUT/$S/glaucoma_wb_chr$C.dat" ] || echo "$S $C"
done; done | xargs -P "$JOBS" -n 2 bash -c 'assoc "$0" "$1"'
[ "${ASSOC_ONLY:-0}" = 1 ] && exit 0

# ---- 3. conditional/joint per chromosome on Bonferroni-significant models
for S in as_reported allele_corrected; do
  d="$OUT/$S"
  # FUSION writes MHC models (chr6) to a separate .dat.MHC file; they are analysed with the rest, as in June 2025
  for C in $(seq 1 22); do
    cat "$d/glaucoma_wb_chr$C.dat" > "$d/glaucoma_wb_chr$C.with_mhc.dat"
    if [ -f "$d/glaucoma_wb_chr$C.dat.MHC" ]; then
      tail -n +2 "$d/glaucoma_wb_chr$C.dat.MHC" >> "$d/glaucoma_wb_chr$C.with_mhc.dat"
    fi
  done
  awk 'FNR==1 && NR>1{next} {print}' "$d"/glaucoma_wb_chr{1..22}.with_mhc.dat > "$d/glaucoma_wb_all_chr.dat"
  NTEST=$(awk -F'\t' 'NR>1 && $NF !~ /NA/' "$d/glaucoma_wb_all_chr.dat" | wc -l | tr -d ' ')
  BONF=$(awk -v n="$NTEST" 'BEGIN{printf "%.6e", 0.05/n}')
  echo -e "models_with_TWAS_P\t$NTEST\nbonferroni_threshold\t$BONF" > "$d/bonferroni.txt"
  for C in $(seq 1 22); do
    awk -F'\t' -v b="$BONF" 'NR==1 || ($NF !~ /NA/ && $NF+0 < b+0)' "$d/glaucoma_wb_chr$C.with_mhc.dat" > "$d/glaucoma_wb_chr$C.sig.dat"
    if [ "$(wc -l < "$d/glaucoma_wb_chr$C.sig.dat")" -le 1 ]; then continue; fi
    Rscript "$FU/FUSION.post_process.R" --sumstats "$OUT/inputs/$S/glaucoma_fusion_chr$C.tsv" \
      --input "$d/glaucoma_wb_chr$C.sig.dat" --out "$d/glaucoma_wb_chr$C.joint" \
      --ref_ld_chr "$LD" --chr "$C" --locus_win 100000 > "$d/post_chr$C.log" 2>&1
  done
done

# ---- 4. ranking and comparison
"${PY:-python3}" "$SCRIPTS_DIR/76h_glaucoma_wb_genomewide_rank.py"
