#!/bin/bash
# ============================================================
# Step 113: gene-level fine-mapping of the 11p11.2 TWAS signal with FOCUS (Mancuso et al. 2019, v0.803).
#   FOCUS models the TWAS z-scores of all genes in the region jointly, using the correlation of their
#   predicted expression (weights x LD), and returns a posterior inclusion probability (PIP) per gene
#   model plus a 90% credible gene set; a null model ("no causal gene among those tested") is included.
#   Weights: GTEx v8 MASHR models (the S-PrediXcan models; FUSION has no ZNF408 model) for cultured
#            fibroblasts, skin (not sun exposed) and oesophagus mucosa, restricted to genes with weights
#            in chr11:45.0-49.5 Mb.
#   GWAS   : European FECD GWAS (GCST90475865; SE from the 95% CI), chr11:46.2-48.3 Mb.
#   LD     : 1000 Genomes phase 3 Europeans (GRCh38), variant IDs set to rsIDs (step 112a panel).
#   Region : chr11:46,200,000-48,300,000; --p-threshold 1e-4 (the regional minimum P is 3.7e-5, so the
#            default genome-wide threshold would skip the locus).
#   Output : results/focus_11p11/
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
cd "$FECD_BASE"
FOCUS=${FOCUS:-$FECD_BASE/software/focus_venv/bin/focus}
OUT=results/focus_11p11; mkdir -p $OUT
LDP=data/ld_1000g_grch38/eur_chr11_46.2-48.3Mb

# 1. region-restricted copies of the MASHR databases, imported into one FOCUS database
rm -f $OUT/mashr_11p11.db
for T in Cells_Cultured_fibroblasts Skin_Not_Sun_Exposed_Suprapubic Esophagus_Mucosa; do
  python3 - "$T" "$OUT" <<'PY'
import sqlite3, sys, shutil, os
T, OUT = sys.argv[1], sys.argv[2]
src, dst = f"data/predictdb/mashr_{T}.db", f"{OUT}/mashr_{T}_11p11.db"
if os.path.exists(dst): os.remove(dst)
shutil.copy(src, dst); c = sqlite3.connect(dst)
genes = {g for g, v in c.execute("select gene, varID from weights")
         if v.split("_")[0] == "chr11" and 45_000_000 <= int(v.split("_")[1]) <= 49_500_000}
c.execute("create temp table keep(gene text)"); c.executemany("insert into keep values (?)", [(g,) for g in genes])
c.execute("delete from weights where gene not in (select gene from keep)")
c.execute("delete from extra where gene not in (select gene from keep)")
c.commit(); c.execute("vacuum"); c.close(); print(f"{T}: {len(genes)} gene models in chr11:45.0-49.5 Mb")
PY
  "$FOCUS" import "$OUT/mashr_${T}_11p11.db" predixcan --tissue "$T" --name GTExv8 --assay rnaseq \
    --use-ens-id --from-gencode --predixcan-method mashr -o $OUT/mashr_11p11 -q
done

# 2. LD panel with rsIDs (from the European GWAS and the MASHR model variants, matched on position + alleles)
python3 - "$LDP" "$OUT" <<'PY'
import sys, sqlite3, glob
LDP, OUT = sys.argv[1], sys.argv[2]
key = lambda p, a, b: (int(p), *sorted((a, b)))
rs = {}
for db in glob.glob(f"{OUT}/mashr_*_11p11.db"):
    for r, v in sqlite3.connect(db).execute("select rsid, varID from weights"):
        c, p, a, b, _ = v.split("_"); rs[key(p, a, b)] = r
with open("analysis/fecd_eur_11p11_region.tsv") as f:
    h = f.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
    for l in f:
        x = l.rstrip("\n").split("\t")
        if x[ix["rsid"]].startswith("rs"): rs.setdefault(key(x[ix["base_pair_location"]], x[ix["effect_allele"]], x[ix["other_allele"]]), x[ix["rsid"]])
n = 0; seen = set()
with open(f"{LDP}.bim") as f, open(f"{OUT}/ld_rsid.bim", "w") as o:
    for l in f:
        c, i, cm, p, a1, a2 = l.split()
        r = rs.get(key(p, a1, a2), i)
        if r in seen: r = i
        seen.add(r); n += r.startswith("rs")
        o.write(f"{c}\t{r}\t{cm}\t{p}\t{a1}\t{a2}\n")
print(f"LD panel variants labelled with rsIDs: {n}")
PY
cp $LDP.bed $OUT/ld_rsid.bed; cp $LDP.fam $OUT/ld_rsid.fam

# 3. European GWAS in FOCUS format (SNP CHR BP A1 A2 Z N; A1 = effect allele)
python3 - "$OUT" <<'PY'
import sys, math, gzip
OUT = sys.argv[1]; N = 448712
with open("analysis/fecd_eur_11p11_region.tsv") as f, gzip.open(f"{OUT}/fecd_eur_11p11.sumstats.gz", "wt") as o:
    h = f.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
    o.write("SNP\tCHR\tBP\tA1\tA2\tZ\tN\n"); n = 0
    for l in f:
        x = l.rstrip("\n").split("\t")
        try:
            b = math.log(float(x[ix["odds_ratio"]]))
            se = float(x[ix["standard_error"]]) if x[ix["standard_error"]] not in ("#NA", "NA", "") else \
                 (math.log(float(x[ix["ci_upper"]])) - math.log(float(x[ix["ci_lower"]]))) / (2 * 1.959963984540054)
        except (ValueError, ZeroDivisionError): continue
        if se <= 0 or not x[ix["rsid"]].startswith("rs"): continue
        o.write(f"{x[ix['rsid']]}\t11\t{x[ix['base_pair_location']]}\t{x[ix['effect_allele']]}\t{x[ix['other_allele']]}\t{b/se:.6f}\t{N}\n"); n += 1
print(f"GWAS variants written: {n}")
PY

# 4. fine-map. Primary: one region spanning both NR1H3 and ZNF408 (FOCUS's default GRCh38 European LD
#    blocks split at 46,984,586, between the two genes). Sensitivity: the default blocks.
#    FOCUS 0.802 needs pandas 1.5 / numpy 1.26 / xarray 2023.1 and two one-line patches (np.warnings ->
#    warnings; DataFrame.pivot keyword arguments); see bundle README.
printf "chrom\tstart\tstop\n11\t46200000\t48300000\n" > $OUT/region_11p11.bed
"$FOCUS" finemap $OUT/fecd_eur_11p11.sumstats.gz $OUT/ld_rsid $OUT/mashr_11p11.db \
  --locations $OUT/region_11p11.bed --p-threshold 1e-4 --trait FECD_EUR -o $OUT/fecd_eur_11p11 -q
BLOCKS=$("$(dirname "$FOCUS")/python" -c "import pyfocus, os; print(os.path.join(os.path.dirname(pyfocus.__file__), 'data/ld_blocks/grch38.eur.loci.bed'))")
"$FOCUS" finemap $OUT/fecd_eur_11p11.sumstats.gz $OUT/ld_rsid $OUT/mashr_11p11.db \
  --locations "$BLOCKS" --chr 11 --p-threshold 1e-4 --trait FECD_EUR -o $OUT/fecd_eur_11p11_defaultblocks -q

# 5. summary: per-model PIPs, and per-gene summed PIP, for both runs
python3 - "$OUT" <<'PY'
import csv, sys, collections
OUT = sys.argv[1]
with open(f"{OUT}/focus_11p11_summary.tsv", "w") as fo:
    fo.write("analysis\tregion\tgene\ttissue\ttwas_z\tPIP\tin_90pct_credible_set\tgene_summed_PIP\n")
    for lab, f in [("primary (chr11:46.2-48.3 Mb)", "fecd_eur_11p11"), ("default LD blocks", "fecd_eur_11p11_defaultblocks")]:
        R = list(csv.DictReader(open(f"{OUT}/{f}.focus.tsv"), delimiter="\t"))
        tot = collections.defaultdict(float)
        for r in R: tot[(r["block"], r["mol_name"])] += float(r["pips_pop1"])
        for r in sorted(R, key=lambda r: (r["block"], -float(r["pips_pop1"]))):
            fo.write(f"{lab}\t{r['block']}\t{r['mol_name']}\t{r['tissue']}\t{r['twas_z_pop1']}\t{float(r['pips_pop1']):.4f}\t"
                     f"{r['in_cred_set_pop1']}\t{tot[(r['block'], r['mol_name'])]:.4f}\n")
        top = sorted({(b, g): v for (b, g), v in tot.items()}.items(), key=lambda x: -x[1])[:4]
        print(lab, "| top genes (summed PIP):", ", ".join(f"{g} {v:.2f}" for (b, g), v in top))
PY
