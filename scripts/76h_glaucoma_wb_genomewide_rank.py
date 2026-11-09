#!/usr/bin/env python3
# ============================================================
# Step 76h-b (29 Oct 2026): genome-wide ranking of the conditionally independent glaucoma
#   whole-blood FUSION signals, from the step 76h reruns (input as reported; input after the
#   allele-matching correction), compared with the June 2025 workstation run read by step 76f
#   (data/derived/glaucoma_wb_conditional_joint_included.tsv).
#   Genes are ranked by JOINT.P, as in step 76f.
#   Outputs: analysis/glaucoma_wb_genomewide_joint_ranking.tsv (one row per gene jointly included
#            in any of the three runs), analysis/glaucoma_wb_genomewide_rerun_summary.tsv
# ============================================================
import csv, glob, os, sqlite3

B = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
os.chdir(B)
GW = "results/twas_wb_genomewide"
RUNS = ("June2025", "as_reported", "allele_corrected")


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# gene symbols: June workbook names, then step 76g's chr11 symbols, then GENCODE v26 names from the MASHR models
SYM = {}
for db in sorted(glob.glob("data/predictdb_all/eqtl/mashr/mashr_*.db")) + sorted(glob.glob("data/predictdb/mashr_*.db")):
    con = sqlite3.connect(db)
    for g, n in con.execute("SELECT gene, genename FROM extra"):
        SYM.setdefault(g.split(".")[0], n)
    con.close()
for r in tsv("data/derived/s8_chr11_wb_gene_symbols.tsv"):
    SYM[r["ENSG"].split(".")[0]] = r["symbol"]
JUNE = tsv("data/derived/glaucoma_wb_conditional_joint_included.tsv")
for r in JUNE:
    SYM[r["ID"].split(".")[0]] = r["gene"]
if os.path.exists("data/gene_symbol_display.tsv"):
    for r in tsv("data/gene_symbol_display.tsv"):
        SYM[r["ENSG"]] = r["display_symbol"]


def sym(i):
    return SYM.get(i.split(".")[0], i)


def ranked(rows):
    rows = sorted(rows, key=lambda r: (num(r["JOINT.P"]) if num(r["JOINT.P"]) is not None else 2,
                                       -abs(num(r["JOINT.Z"]) or 0)))
    return {r["ID"]: (i, r) for i, r in enumerate(rows, 1)}


joint, marg, chrom = {"June2025": ranked(JUNE)}, {}, {}
for s in RUNS[1:]:
    rows = []
    for f in sorted(glob.glob(f"{GW}/{s}/glaucoma_wb_chr*.joint.joint_included.dat")):
        rows += [{k: v.strip() for k, v in r.items()} for r in tsv(f)]
    joint[s] = ranked(rows)
    marg[s] = {}
    for r in tsv(f"{GW}/{s}/glaucoma_wb_all_chr.dat"):
        marg[s][r["ID"]] = r
        chrom[r["ID"]] = (int(r["CHR"]), int(r["P0"]))

ids = sorted(set().union(*joint.values()), key=lambda i: chrom.get(i, (99, 0)))
out = []
for i in ids:
    row = [sym(i), i, *chrom.get(i, ("", ""))]
    for s in RUNS:
        rk, r = joint[s].get(i, ("", None))
        if s == "June2025":
            row += [rk, r["TWAS.Z"] if r else "", r["TWAS.P"] if r else ""]
        else:
            m = marg[s].get(i, {})
            row += [rk, m.get("TWAS.Z", "").strip(), m.get("TWAS.P", "").strip()]
        row += [r["JOINT.Z"] if r else "", r["JOINT.P"] if r else ""]
    out.append(row)
out.sort(key=lambda r: (r[4] if r[4] != "" else 99, r[9] if r[9] != "" else 99))
os.makedirs("analysis", exist_ok=True)
with open("analysis/glaucoma_wb_genomewide_joint_ranking.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["gene", "ID", "CHR", "P0"] + [f"{s}_{k}" for s in RUNS for k in
                                              ("joint_rank", "TWAS.Z", "TWAS.P", "JOINT.Z", "JOINT.P")])
    w.writerows(out)

# ---- summary
NR = next(i for i in ids if i.startswith("ENSG00000025434"))
summ = []
for s in RUNS:
    top = sorted(joint[s].values(), key=lambda t: t[0])
    b = dict(l.rstrip("\n").split("\t") for l in open(f"{GW}/{s}/bonferroni.txt")) if s != "June2025" else {}
    summ.append([s, b.get("models_with_TWAS_P", ""), b.get("bonferroni_threshold", ""), len(top),
                 joint[s][NR][0] if NR in joint[s] else "not jointly included",
                 ", ".join(f"{sym(r['ID'])} ({r['JOINT.Z']})" for _, r in top[:5])])
dz = [(abs(num(marg["allele_corrected"][i]["TWAS.Z"]) - num(marg["as_reported"][i]["TWAS.Z"])), i)
      for i in marg["as_reported"] if num(marg["as_reported"][i]["TWAS.Z"]) is not None
      and num(marg["allele_corrected"][i]["TWAS.Z"]) is not None]
changed = [t for t in dz if t[0] > 0]
big = max(dz)
with open("analysis/glaucoma_wb_genomewide_rerun_summary.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["run", "models_with_TWAS_P", "bonferroni_threshold", "jointly_included_genes", "NR1H3_rank",
                "top5_by_JOINT.P (JOINT.Z)"])
    w.writerows(summ)
    w.writerow([])
    w.writerow(["models_with_changed_TWAS.Z_after_correction", len(changed)])
    w.writerow(["largest_abs_change_TWAS.Z", f"{big[0]:.4f}", sym(big[1])])
    same = set(joint["as_reported"]) == set(joint["allele_corrected"])
    w.writerow(["same_jointly_included_genes_reported_vs_corrected", same])
    w.writerow(["same_jointly_included_genes_June2025_vs_rerun_as_reported",
                set(joint["June2025"]) == set(joint["as_reported"])])
for l in open("analysis/glaucoma_wb_genomewide_rerun_summary.tsv"):
    print(l.rstrip())
