#!/bin/bash
# ============================================================
# fusion_genomewide.sh (12 Nov 2026): the one genome-wide FUSION TWAS + conditional/joint procedure used for every
#   trait: FECD (step 82: cultured fibroblasts, skin, oesophagus mucosa), glaucoma (step 76h) and IOP and
#   refractive error (step 76i), the three cross-trait phenotypes in whole blood.
#   usage: fusion_genomewide.sh NAME SUMSTATS WEIGHTS_POS OUTDIR
#     NAME         output prefix, e.g. glaucoma_wb
#     SUMSTATS     FUSION-format summary statistics (SNP A1 A2 Z N; plain text or .gz)
#     WEIGHTS_POS  FUSION .pos file of the weight set; its folder is the weights directory
#     OUTDIR       output folder
#   1. per-chromosome inputs restricted to the 1000G EUR LD-reference SNPs of that chromosome (FUSION uses no
#      other SNP): OUTDIR/inputs/NAME_chr{1..22}.tsv
#   2. FUSION.assoc_test per chromosome, JOBS at a time: OUTDIR/NAME_chr{C}.dat
#   3. OUTDIR/NAME_all_chr.dat: every model outside the MHC. FUSION writes models in the MHC (chr6) to a separate
#      file because of the region's extended LD; they are collected in OUTDIR/NAME_MHC.dat and not used further.
#   4. genome-wide Bonferroni threshold, 0.05 / models with a TWAS P in NAME_all_chr.dat (OUTDIR/bonferroni.txt)
#   5. FUSION.post_process per chromosome on the models below that threshold (100-kb locus window):
#      OUTDIR/NAME_chr{C}.joint.joint_included.dat / .joint_dropped.dat. FUSION keeps a model in the joint model
#      at P < 0.05 / (models entered on that chromosome).
#   Reruns: inputs are rebuilt when SUMSTATS is newer than them; a chromosome is re-tested when its input is newer
#   than its result; steps 4-5 always rerun.
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
cd "$FECD_BASE"
[ $# -eq 4 ] || { echo "usage: fusion_genomewide.sh NAME SUMSTATS WEIGHTS_POS OUTDIR" >&2; exit 2; }
NAME=$1; SUMSTATS=$2; WPOS=$3; OUT=$4
WDIR="$(dirname "$WPOS")/"   # trailing slash as step 82 always passed it: FUSION writes FILE as <dir>/<model>, and S3 reports FILE
LD=data/ld_ref/LDREF/1000G.EUR.
FU="$FUSION_DIR"
JOBS=${JOBS:-4}
for f in "$SUMSTATS" "$WPOS" "${LD}22.bim"; do [ -s "$f" ] || { echo "ERROR: $f not found" >&2; exit 1; }; done
mkdir -p "$OUT/inputs"

# ---- 1. per-chromosome inputs (LD-reference SNPs only)
if [ ! -s "$OUT/inputs/${NAME}_chr22.tsv" ] || [ "$SUMSTATS" -nt "$OUT/inputs/${NAME}_chr22.tsv" ]; then
  "${PY:-python3}" - "$SUMSTATS" "$OUT/inputs" "$NAME" <<'PYEOF'
import gzip, os, sys
src, out, name = sys.argv[1:4]
chrom = {}
for c in range(1, 23):
    for line in open(f"data/ld_ref/LDREF/1000G.EUR.{c}.bim"):
        chrom[line.split("\t", 2)[1]] = c
tmp = {c: f"{out}/{name}_chr{c}.tsv.tmp" for c in range(1, 23)}
fh = {c: open(p, "w") for c, p in tmp.items()}
n = dict.fromkeys(range(1, 23), 0)
with (gzip.open(src, "rt") if src.endswith(".gz") else open(src)) as f:
    hdr = f.readline()
    for h in fh.values():
        h.write(hdr)
    for line in f:
        c = chrom.get(line.split("\t", 1)[0])
        if c:
            fh[c].write(line); n[c] += 1
for c, h in fh.items():
    h.close(); os.replace(tmp[c], tmp[c][:-4])
print(f"{name}: LD-reference SNPs in the input: {sum(n.values())} ({' '.join(f'{c}:{n[c]}' for c in range(1, 23))})")
PYEOF
fi

# ---- 2. FUSION.assoc_test per chromosome
assoc() {   # $1 = chromosome
  local c=$1 pos="$OUT/inputs/${NAME}_weights_chr$1.pos"
  head -1 "$WPOS" > "$pos"
  awk -F'\t' -v c="$c" '$4==c' "$WPOS" >> "$pos"
  rm -f "$OUT/${NAME}_chr$c.dat" "$OUT/${NAME}_chr$c.dat.MHC"
  Rscript "$FU/FUSION.assoc_test.R" --sumstats "$OUT/inputs/${NAME}_chr$c.tsv" \
    --weights "$pos" --weights_dir "$WDIR" --ref_ld_chr "$LD" --chr "$c" \
    --out "$OUT/${NAME}_chr$c.dat" > "$OUT/${NAME}_assoc_chr$c.log" 2>&1
}
export -f assoc
export OUT NAME WPOS WDIR LD FU
for C in $(seq 1 22); do
  if [ ! -s "$OUT/${NAME}_chr$C.dat" ] || [ "$OUT/inputs/${NAME}_chr$C.tsv" -nt "$OUT/${NAME}_chr$C.dat" ]; then echo "$C"; fi
done | xargs -P "$JOBS" -n 1 bash -c 'assoc "$0"'
for C in $(seq 1 22); do
  [ -s "$OUT/${NAME}_chr$C.dat" ] || { echo "ERROR: FUSION.assoc_test failed for $NAME chr$C (see $OUT/${NAME}_assoc_chr$C.log)" >&2; exit 1; }
done

# ---- 3. all models outside the MHC; MHC models kept aside
awk 'FNR==1 && NR>1{next} {print}' "$OUT"/${NAME}_chr{1..22}.dat > "$OUT/${NAME}_all_chr.dat"
MHC=$(ls "$OUT"/${NAME}_chr*.dat.MHC 2>/dev/null || true)
if [ -n "$MHC" ]; then awk 'FNR==1 && NR>1{next} {print}' $MHC > "$OUT/${NAME}_MHC.dat"; else rm -f "$OUT/${NAME}_MHC.dat"; fi

# ---- 4. genome-wide Bonferroni threshold
NTEST=$(awk -F'\t' 'NR>1 && $NF !~ /NA/' "$OUT/${NAME}_all_chr.dat" | wc -l | tr -d ' ')
BONF=$(awk -v n="$NTEST" 'BEGIN{printf "%.6e", 0.05/n}')
echo -e "models_with_TWAS_P\t$NTEST\nbonferroni_threshold\t$BONF" > "$OUT/bonferroni.txt"

# ---- 5. conditional/joint analysis per chromosome on the Bonferroni-significant models
NSIG=0
for C in $(seq 1 22); do
  rm -f "$OUT/${NAME}_chr$C.joint".* "$OUT/${NAME}_chr$C.sig.dat"
  awk -F'\t' -v b="$BONF" 'NR==1 || ($NF !~ /NA/ && $NF+0 < b+0)' "$OUT/${NAME}_chr$C.dat" > "$OUT/${NAME}_chr$C.sig.dat"
  k=$(($(wc -l < "$OUT/${NAME}_chr$C.sig.dat") - 1))
  if [ "$k" -eq 0 ]; then rm -f "$OUT/${NAME}_chr$C.sig.dat"; continue; fi
  NSIG=$((NSIG + k))
  Rscript "$FU/FUSION.post_process.R" --sumstats "$OUT/inputs/${NAME}_chr$C.tsv" \
    --input "$OUT/${NAME}_chr$C.sig.dat" --out "$OUT/${NAME}_chr$C.joint" \
    --ref_ld_chr "$LD" --chr "$C" --locus_win 100000 > "$OUT/${NAME}_post_chr$C.log" 2>&1
done
echo -e "bonferroni_significant_models\t$NSIG" >> "$OUT/bonferroni.txt"
echo "$NAME: $NTEST models with a TWAS P outside the MHC; Bonferroni $BONF; $NSIG significant; $(cat "$OUT"/${NAME}_chr*.joint.joint_included.dat 2>/dev/null | grep -vc '^FILE' || true) jointly included"
