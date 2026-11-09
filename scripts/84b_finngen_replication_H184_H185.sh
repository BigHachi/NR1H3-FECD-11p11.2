#!/bin/bash
# ============================================================
# Step 84b (Oct 8 2026; scripted Oct 10 2026 from the review-rerun commands):
#   NR1H3 / ZNF408 S-PrediXcan look-up in two FinnGen R12 corneal endpoints.
#     H7_CORNEALDYSTROPHY  491 cases  = ICD-10 H18.4 corneal DEGENERATION (despite its name)
#     H7_KERATOCONUS     1,076 cases  = ICD-10 H18.5 hereditary corneal DYSTROPHIES (includes FECD)
#   Controls 473,095 for both (FinnGen Risteys). Supersedes step 84 (H18.4 only).
#   Only the 11p11.2 region (chr11:45-49.5 Mb, GRCh38) is needed; H18.5 is fetched
#   by tabix from the public FinnGen bucket, H18.4 from the full download (step 84).
#   GWAS variants are matched to the GTEx v8 MASHR model varIDs by position + allele pair.
#   Output: analysis/nr1h3_finngen_replication_H184_vs_H185.tsv  (Supplementary Table S4)
#   Manuscript: H18.4 fib Z -2.10 P 0.035, skin Z -2.39 P 0.017;
#               H18.5 fib Z -0.32 P 0.75, skin Z -0.41 P 0.68.  Power: step 84c.
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -euo pipefail
BASE="$FECD_BASE"; cd $BASE
TABIX=${TABIX:-tabix}
PY=${PY:-$METAXCAN_PY}                  # env with MetaXcan dependencies
M=data/predictdb_all/eqtl/mashr                       # GTEx v8 MASHR models (49 tissues)
R=data/gwas/finngen_h185_region; mkdir -p $R results/sprediXcan_finngen_repl logs
U=https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release

# ---- 1. 11p11.2 region of each endpoint ----
[ -s data/gwas/finngen_R12_H7_CORNEALDYSTROPHY.gz ] || curl -L --fail -o data/gwas/finngen_R12_H7_CORNEALDYSTROPHY.gz $U/finngen_R12_H7_CORNEALDYSTROPHY.gz
# Each region file is written to .part and renamed only when it has data rows, so a failed
# fetch is retried next time instead of leaving a header-only file that would be skipped.
K=H7_KERATOCONUS_chr11_45-49.5Mb.tsv
[ -s $R/$K ] || {
  ( cd $R
    { { gzip -cd ../finngen_R12_H7_CORNEALDYSTROPHY.gz || true; } | head -1   # header; gzip stops early when head closes the pipe
      $TABIX $U/finngen_R12_H7_KERATOCONUS.gz 11:45000000-49500000; } > $K.part
    [ "$(wc -l < $K.part)" -gt 1 ] && mv $K.part $K ); }
C=H7_CORNEALDYSTROPHY_chr11_45-49.5Mb.tsv
[ -s $R/$C ] || {
  gzip -cd data/gwas/finngen_R12_H7_CORNEALDYSTROPHY.gz | \
    awk -F'\t' 'NR==1 || ($1==11 && $2>=45000000 && $2<=49500000)' > $R/$C.part
  [ "$(wc -l < $R/$C.part)" -gt 1 ] && mv $R/$C.part $R/$C; }

# ---- 2. convert to S-PrediXcan input keyed on MASHR varIDs ----
python3 - <<'E'
import sqlite3, glob, os, gzip
B = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
want = {}
for db in glob.glob(f'{B}/data/predictdb_all/eqtl/mashr/mashr_*.db'):
    con = sqlite3.connect(db)
    for (v,) in con.execute("SELECT DISTINCT varID FROM weights WHERE varID LIKE 'chr11_%'"):
        c, p, a1, a2, _ = v.split('_'); want.setdefault((p, frozenset((a1, a2))), set()).add(v)
for ep, ncase in [('H7_KERATOCONUS', 1076), ('H7_CORNEALDYSTROPHY', 491)]:
    N = ncase + 473095
    fi = open(f'{B}/data/gwas/finngen_h185_region/{ep}_chr11_45-49.5Mb.tsv')
    h = fi.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(h)}
    with gzip.open(f'{B}/data/gwas/finngen_{ep}_chr11_sprediXcan.txt.gz', 'wt') as fo:
        fo.write('panel_variant_id\trsid\tchromosome\tposition\teffect_allele\tnon_effect_allele\tzscore\tpvalue\tsample_size\teffect_size\tstandard_error\tfrequency\n')
        n = 0
        for line in fi:
            x = line.rstrip('\n').split('\t')
            try: b = float(x[c['beta']]); se = float(x[c['sebeta']])
            except ValueError: continue
            if se <= 0: continue
            for v in want.get((x[c['pos']], frozenset((x[c['ref']], x[c['alt']]))), ()):
                fo.write(f"{v}\t{x[c['rsids']]}\t11\t{x[c['pos']]}\t{x[c['alt']]}\t{x[c['ref']]}\t{b/se:.6f}\t{x[c['pval']]}\t{N}\t{b}\t{se}\t{x[c['af_alt']]}\n"); n += 1
    print(ep, 'model variants matched in region:', n)
E

# ---- 3. S-PrediXcan, 3 proxy tissues ----
for ep in H7_KERATOCONUS H7_CORNEALDYSTROPHY; do
  for T in Cells_Cultured_fibroblasts Skin_Not_Sun_Exposed_Suprapubic Esophagus_Mucosa; do
    O=results/sprediXcan_finngen_repl/${ep}_${T}.csv; rm -f $O          # S-PrediXcan skips existing output
    $PY software/MetaXcan/software/SPrediXcan.py --model_db_path $M/mashr_$T.db --covariance $M/mashr_$T.txt.gz \
      --gwas_file data/gwas/finngen_${ep}_chr11_sprediXcan.txt.gz --snp_column panel_variant_id \
      --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele \
      --beta_column effect_size --se_column standard_error --pvalue_column pvalue --freq_column frequency \
      --keep_non_rsid --model_db_snp_key varID --remove_ens_version --output_file $O > logs/repl_${ep}_$T.log 2>&1
  done
done

# ---- 4. NR1H3 / ZNF408 table ----
python3 - <<'E'
import csv, os
B = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS")); rows = []
for ep, lab, ca in [('H7_CORNEALDYSTROPHY', 'H18.4 corneal degeneration', 491), ('H7_KERATOCONUS', 'H18.5 hereditary corneal dystrophies', 1076)]:
    for T in ['Cells_Cultured_fibroblasts', 'Skin_Not_Sun_Exposed_Suprapubic', 'Esophagus_Mucosa']:
        for x in csv.DictReader(open(f'{B}/results/sprediXcan_finngen_repl/{ep}_{T}.csv')):
            if x['gene_name'] in ('NR1H3', 'ZNF408'):
                rows.append([ep, lab, ca, 473095, T, x['gene_name'], x['zscore'], x['pvalue'], x['n_snps_used'], x['n_snps_in_model']])
with open(f'{B}/analysis/nr1h3_finngen_replication_H184_vs_H185.tsv', 'w') as f:
    f.write('endpoint\tICD10\tcases\tcontrols\ttissue\tgene\tzscore\tpvalue\tn_snps_used\tn_snps_in_model\n')
    for r in rows: f.write('\t'.join(map(str, r)) + '\n')
for r in rows:
    if r[5] == 'NR1H3': print(f"{r[0]:20} {r[4]:32} z={float(r[6]):6.2f} P={float(r[7]):.3g}")
E
