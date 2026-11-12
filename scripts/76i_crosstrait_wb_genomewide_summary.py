#!/usr/bin/env python3
# ============================================================
# Step 76i-b (12 Nov 2026): summary of the genome-wide whole-blood FUSION runs for IOP and refractive error
#   (step 76i). For each trait: the genome-wide Bonferroni threshold, the jointly included genes ranked by
#   JOINT.P (as steps 76f and 76h rank glaucoma), the NR1H3 status (jointly included, conditioned away, or below
#   the genome-wide threshold and so not in the conditional analysis) and, when the non-public June 2025 workbook
#   is present (data/private/), a comparison with the June 2025 workstation run.
#   Outputs: analysis/crosstrait_wb_genomewide_joint_ranking.tsv, analysis/crosstrait_wb_genomewide_summary.tsv
# ============================================================
import csv, glob, os, sqlite3

B = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
os.chdir(B)
GW = "results/twas_wb_genomewide_crosstrait"
TRAITS = {"iop": ("IOP", "IOP"), "myopia": ("Refractive error", "myopia")}   # trait -> (label, June 2025 sheet key)
JUNE_XLSX = os.environ.get("GLAUCOMA_COND_XLSX", "data/private/Eye Conditional results 12_6_25 FinnGEN GTEV8 whole blood.xlsx")   # config.sh, as step 76f
NR = "ENSG00000025434"


def tsv(path):
    return [{k.strip(): (v or "").strip() for k, v in r.items()} for r in csv.DictReader(open(path), delimiter="\t")]


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


SYM = {}
for db in sorted(glob.glob("data/predictdb_all/eqtl/mashr/mashr_*.db")) + sorted(glob.glob("data/predictdb/mashr_*.db")):
    con = sqlite3.connect(db)
    for g, n in con.execute("SELECT gene, genename FROM extra"):
        SYM.setdefault(g.split(".")[0], n)
    con.close()
for r in tsv("data/derived/s8_chr11_wb_gene_symbols.tsv"):
    SYM[r["ENSG"].split(".")[0]] = r["symbol"]
if os.path.exists("data/gene_symbol_display.tsv"):
    for r in tsv("data/gene_symbol_display.tsv"):
        SYM[r["ENSG"]] = r["display_symbol"]


def sym(i):
    return SYM.get(i.split(".")[0], i)


def ranked(rows):
    rows = sorted(rows, key=lambda r: (num(r["JOINT.P"]) if num(r["JOINT.P"]) is not None else 2,
                                       -abs(num(r["JOINT.Z"]) or 0)))
    return {r["ID"]: (i, r) for i, r in enumerate(rows, 1)}


def june(key):
    """jointly included genes of the June 2025 sheet for this trait, or None when the workbook is absent"""
    if not os.path.exists(JUNE_XLSX):
        return None
    import openpyxl
    wb = openpyxl.load_workbook(JUNE_XLSX, read_only=True, data_only=True)
    name = next(s for s in wb.sheetnames if key.lower() in s.lower())
    it = wb[name].iter_rows(values_only=True)
    hdr = [str(h).strip() if h is not None else "" for h in next(it)]
    rows = []
    for v in it:
        r = dict(zip(hdr, v))
        if r.get("ID") and str(r["ID"]).startswith("ENSG"):
            rows.append({k: ("" if r.get(k) is None else str(r[k])) for k in ("ID", "TWAS.Z", "TWAS.P", "JOINT.Z", "JOINT.P", "Genename")})
    for r in rows:
        if r["Genename"]:
            SYM.setdefault(r["ID"].split(".")[0], r["Genename"])
    return ranked(rows)


rank_rows, summ = [], []
for t, (label, key) in TRAITS.items():
    d = f"{GW}/{t}"
    b = dict(l.rstrip("\n").split("\t") for l in open(f"{d}/bonferroni.txt"))
    marg = {r["ID"]: r for r in tsv(f"{d}/{t}_wb_all_chr.dat")}             # models outside the MHC
    inc, drop = [], []
    for f in sorted(glob.glob(f"{d}/{t}_wb_chr*.joint.joint_included.dat")):
        inc += tsv(f)
    for f in sorted(glob.glob(f"{d}/{t}_wb_chr*.joint.joint_dropped.dat")):
        drop += tsv(f)
    joint = ranked(inc)
    dropped = {r["ID"]: r for r in drop}
    jn = june(key)
    nr = next(i for i in marg if i.startswith(NR))
    m = marg[nr]
    if nr in joint:
        k, r = joint[nr]
        status = f"conditionally independent (rank {k} of {len(joint)} by JOINT.P; JOINT.Z {r['JOINT.Z']}, JOINT.P {r['JOINT.P']})"
    elif nr in dropped:
        r = dropped[nr]
        status = f"conditioned away (COND.Z {r.get('COND.Z', '')}, COND.P {r.get('COND.P', '')})"
    else:
        status = "below the genome-wide threshold, so not in the conditional analysis"
    summ += [[label, "models_with_TWAS_P", b["models_with_TWAS_P"]],
             [label, "bonferroni_threshold", b["bonferroni_threshold"]],
             [label, "jointly_included_genes", len(joint)],
             [label, "top5_by_JOINT.P", ", ".join(f"{sym(i)} ({r['JOINT.Z']})" for i, (k, r) in sorted(joint.items(), key=lambda x: x[1][0])[:5])],
             [label, "NR1H3_TWAS.Z_TWAS.P", f"{m['TWAS.Z']}, {m['TWAS.P']}"],
             [label, "NR1H3_genome_wide_status", status]]
    if jn is not None:
        nj = jn.get(nr)
        summ += [[label, "June2025_jointly_included_genes", len(jn)],
                 [label, "June2025_NR1H3_status", f"conditionally independent (rank {nj[0]}; JOINT.Z {nj[1]['JOINT.Z']}, JOINT.P {nj[1]['JOINT.P']})" if nj else "not jointly included"],
                 [label, "jointly_included_in_both_runs", len(set(jn) & set(joint))],
                 [label, "only_in_rerun", ", ".join(sorted(sym(i) for i in set(joint) - set(jn)))],
                 [label, "only_in_June2025", ", ".join(sorted(sym(i) for i in set(jn) - set(joint)))]]
    for i in sorted(set(joint) | set(jn or {}), key=lambda i: (joint.get(i, (9999,))[0], (jn or {}).get(i, (9999,))[0])):
        k, r = joint.get(i, ("", {}))
        mm = marg.get(i, {})
        kj, rj = (jn or {}).get(i, ("", {}))
        rank_rows.append([label, sym(i), i, mm.get("CHR", ""), mm.get("P0", ""), k, mm.get("TWAS.Z", ""), mm.get("TWAS.P", ""),
                          r.get("JOINT.Z", ""), r.get("JOINT.P", ""), kj, rj.get("TWAS.Z", ""), rj.get("JOINT.Z", ""), rj.get("JOINT.P", "")])

os.makedirs("analysis", exist_ok=True)
with open("analysis/crosstrait_wb_genomewide_joint_ranking.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["trait", "gene", "ID", "CHR", "P0", "rerun_joint_rank", "TWAS.Z", "TWAS.P", "JOINT.Z", "JOINT.P",
                "June2025_joint_rank", "June2025_TWAS.Z", "June2025_JOINT.Z", "June2025_JOINT.P"])
    w.writerows(rank_rows)
with open("analysis/crosstrait_wb_genomewide_summary.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["trait", "item", "value"])
    w.writerows(summ)
for r in summ:
    print("\t".join(map(str, r)))
