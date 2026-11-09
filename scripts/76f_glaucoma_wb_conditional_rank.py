#!/usr/bin/env python3
# ============================================================
# Step 76f (Oct 9 2026; scripted Oct 10 2026): genome-wide rank of NR1H3 among the
#   conditionally independent FUSION whole-blood TWAS signals for glaucoma.
#   Source: genome-wide FUSION GTEx v8 EUR whole-blood TWAS + FUSION.post_process
#   (conditional/joint) on FinnGen H7_GLAUCOMA, run 12 June 2025 on the University
#   workstation; results kept in "Eye Conditional results 12_6_25 FinnGEN GTEV8 whole
#   blood.xlsx" (sheet 2 = jointly included genes). The chr11 re-run in step 76d
#   reproduces the NR1H3 values (TWAS Z -7.7, JOINT P 1.2e-14).
#   Manuscript: NR1H3 is the third-strongest conditionally independent glaucoma signal,
#   behind LOXL1-AS1 and AFAP1.
# ============================================================
import csv, os, sys

B = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
X = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "GLAUCOMA_COND_XLSX", f"{B}/data/private/Eye Conditional results 12_6_25 FinnGEN GTEV8 whole blood.xlsx")
DERIVED = f"{B}/data/derived/glaucoma_wb_conditional_joint_included.tsv"   # bundled fallback
COLS = ["FILE", "ID", "TWAS.Z", "TWAS.P", "JOINT.BETA", "JOINT.BETA.SE", "JOINT.Z", "JOINT.P", "gene"]
if os.path.exists(X):          # original workbook (analysis Mac): read sheet 2 and refresh the bundled table
    import openpyxl
    ws = openpyxl.load_workbook(X, read_only=True, data_only=True).worksheets[1]
    rows = [r for r in list(ws.iter_rows(values_only=True))[1:] if r[1] not in (None, "None")]
    os.makedirs(os.path.dirname(DERIVED), exist_ok=True)
    with open(DERIVED, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t"); w.writerow(COLS[1:])
        for r in rows: w.writerow([("" if v is None else v) for v in r[1:9]])
else:                          # anyone else: the jointly included genes, as exported from the workbook
    rows = [[None] + [r[c] for c in COLS[1:]] for r in csv.DictReader(open(DERIVED), delimiter="\t")]
def num(v):
    try: return float(v)
    except (TypeError, ValueError): return 1.0
rows.sort(key=lambda r: num(r[7]))                     # JOINT.P
print(f"jointly included genes: {len(rows)}")
for i, r in enumerate(rows[:5], 1):
    print(f"  {i}\t{r[8]}\tTWAS.Z={r[2]}\tTWAS.P={r[3]}\tJOINT.Z={r[6]}\tJOINT.P={r[7]}")
print("NR1H3 rank:", next(i for i, r in enumerate(rows, 1) if r[8] == "NR1H3"))
