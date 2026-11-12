#!/bin/bash
# ============================================================
# Step 76d: chromosome-11 whole-blood FUSION results for the cross-trait phenotypes (S8, Figure 5 source).
#   12 Nov 2026 (release v14.6): taken from the genome-wide runs of the procedure shared by every trait
#   (scripts/fusion_genomewide.sh): glaucoma from step 76h (allele-matched input of step 60b; the first-record input
#   of releases up to v13 into previous_first_record/ for comparison), IOP and refractive error from step 76i.
#   The conditional/joint status is therefore the genome-wide one: chromosome-11 models below the trait's genome-wide
#   Bonferroni threshold entered FUSION.post_process (100-kb locus window); the others are not in the joint model.
#   Up to release v14.5 this step ran FUSION on chromosome 11 alone with every chr11 model entering the joint model
#   (joint threshold 0.05/463), which conditioned NR1H3 away for refractive error; with the genome-wide threshold
#   NR1H3 stays in the joint model for refractive error, alongside correlated neighbours (see README).
#   Outputs: results/twas_wb_conditional/{trait}_wb_chr11.dat (all chr11 models), .sig.dat (models entered in the
#              joint model), .joint.joint_included.dat / .joint_dropped.dat, .bonferroni.txt
#            data/gwas/{trait}_gwas_fusion_chr11.tsv (the chr11 input; audited by step 55c)
#            analysis/nr1h3_crosstrait_fusion_conditional.tsv
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
cd "$FECD_BASE"
OUT=results/twas_wb_conditional
GW=results/twas_wb_genomewide
CT=results/twas_wb_genomewide_crosstrait
mkdir -p "$OUT/previous_first_record"
JHDR_INC=$'FILE\tID\tTWAS.Z\tTWAS.P\tJOINT.BETA\tJOINT.BETA.SE\tJOINT.Z\tJOINT.P'
JHDR_DRP=$'FILE\tID\tTWAS.Z\tTWAS.P\tCOND.BETA\tCOND.BETA.SE\tCOND.Z\tCOND.P'

take() {   # $1 = source prefix (genome-wide run, NAME), $2 = source folder, $3 = destination prefix
  local src="$2/$1" dst=$3
  for f in "$src"_chr11.dat "$2/bonferroni.txt" "$2/inputs/$1_chr11.tsv"; do
    [ -s "$f" ] || { echo "ERROR: $f not found (run steps 76h and 76i first)" >&2; exit 1; }
  done
  cp "$src"_chr11.dat "$dst"_chr11.dat
  cp "$2/bonferroni.txt" "$dst"_chr11.bonferroni.txt
  if [ -s "$src"_chr11.sig.dat ]; then cp "$src"_chr11.sig.dat "$dst"_chr11.sig.dat; else head -1 "$src"_chr11.dat > "$dst"_chr11.sig.dat; fi
  if [ -s "$src"_chr11.joint.joint_included.dat ]; then cp "$src"_chr11.joint.joint_included.dat "$dst"_chr11.joint.joint_included.dat
  else echo "$JHDR_INC" > "$dst"_chr11.joint.joint_included.dat; fi
  if [ -s "$src"_chr11.joint.joint_dropped.dat ]; then cp "$src"_chr11.joint.joint_dropped.dat "$dst"_chr11.joint.joint_dropped.dat
  else echo "$JHDR_DRP" > "$dst"_chr11.joint.joint_dropped.dat; fi
}
take glaucoma_wb "$GW/allele_corrected" "$OUT/glaucoma_wb"
take glaucoma_wb "$GW/as_reported"      "$OUT/previous_first_record/glaucoma_wb"
take iop_wb      "$CT/iop"              "$OUT/iop_wb"
take myopia_wb   "$CT/myopia"           "$OUT/myopia_wb"
cp "$GW/allele_corrected/inputs/glaucoma_wb_chr11.tsv" data/gwas/glaucoma_gwas_fusion_chr11.tsv
cp "$GW/as_reported/inputs/glaucoma_wb_chr11.tsv"      data/gwas/glaucoma_gwas_fusion_firstrecord_chr11.tsv
cp "$CT/iop/inputs/iop_wb_chr11.tsv"                   data/gwas/iop_gwas_fusion_chr11.tsv
cp "$CT/myopia/inputs/myopia_wb_chr11.tsv"             data/gwas/myopia_gwas_fusion_chr11.tsv

# ---- NR1H3 marginal + conditional status ----
# marginal TWAS Z from the full-precision results file (the joint-model files round it); P = 2*pnorm(-|Z|), FUSION's own
# formula (the results file prints P to 3 significant figures); JOINT/COND Z and P from the joint-model files
SUM=analysis/nr1h3_crosstrait_fusion_conditional.tsv
"${PY:-python3}" - "$OUT" "$SUM" <<'PYEOF'
import csv, math, sys
out, summ = sys.argv[1:3]
def rows(p):
    return [{k.strip(): (v or "").strip() for k, v in r.items()} for r in csv.DictReader(open(p), delimiter="\t")]
with open(summ, "w") as fo:
    fo.write("trait\ttissue\tmethod\tTWAS_Z\tTWAS_P\tJOINT_Z\tJOINT_P\tconditional_status\n")
    for t in ("glaucoma", "iop", "myopia"):
        m = next(r for r in rows(f"{out}/{t}_wb_chr11.dat") if r["ID"].startswith("ENSG00000025434"))
        z = float(m["TWAS.Z"])
        inc = [r for r in rows(f"{out}/{t}_wb_chr11.joint.joint_included.dat") if r["ID"].startswith("ENSG00000025434")]
        drp = [r for r in rows(f"{out}/{t}_wb_chr11.joint.joint_dropped.dat") if r["ID"].startswith("ENSG00000025434")]
        if inc:
            jz, jp, st = inc[0]["JOINT.Z"], inc[0]["JOINT.P"], "conditionally_independent"
        elif drp:
            jz, jp, st = drp[0]["COND.Z"], drp[0]["COND.P"], "conditioned_away"
        else:   # below the genome-wide threshold: not in the conditional analysis
            jz, jp, st = "NA", "NA", "not_in_joint_model"
        fo.write(f"{t}\twhole_blood\tFUSION\t{z}\t{math.erfc(abs(z) / math.sqrt(2)):.4g}\t{jz}\t{jp}\t{st}\n")
PYEOF
# FECD disease-proxy tissues: FUSION EUR marginal (step 82 output); fibroblast status from the 11p11.2 conditional (step 24b)
printf "FECD_fibroblast\tcultured_fibroblasts\tFUSION_EUR\t-3.74\t1.85e-4\tNA\tNA\tconditionally_independent\n" >> "$SUM"
printf "FECD_skin\tskin_not_sun_exposed\tFUSION_EUR\t-3.53\t4.21e-4\tNA\tNA\tNA\n" >> "$SUM"
echo "=== NR1H3 cross-trait FUSION conditional summary ==="; column -t -s$'\t' "$SUM"
