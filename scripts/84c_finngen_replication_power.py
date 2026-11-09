#!/usr/bin/env python3
# ============================================================
# Step 84c: power of the FinnGen look-up (step 84b), per gene and tissue
#   Expected look-up |Z| = discovery |Z| x sqrt(Neff_look-up / Neff_discovery),
#   Neff = 4 / (1/cases + 1/controls); discovery = MVP multi-ancestry (3,928 / 566,393).
#   Discovery Z = S-PrediXcan Z in the same tissue (step 95 outputs, results/sprediXcan_alltissue_remap).
#   Threshold: two-sided P < 0.025 (0.05 / 2 tissues for the prioritised gene).
#   Oct 21 2026: power is computed for each gene x tissue from its own discovery Z (earlier versions used the
#   NR1H3 Z for every row, so the ZNF408 rows of S4 showed NR1H3's power) and written to a table.
#   Manuscript (NR1H3): H18.4 17% (fibroblasts) / 21% (skin); H18.5 36% / 45%.
#   Output: analysis/nr1h3_finngen_lookup_power.tsv
# ============================================================
import csv, os
from math import sqrt
from statistics import NormalDist

B = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
N = NormalDist()
def neff(ca, co): return 4 / (1 / ca + 1 / co)
disc = neff(3928, 566393)
crit = N.inv_cdf(1 - 0.025 / 2)
GENES = {"ENSG00000025434": "NR1H3", "ENSG00000175213": "ZNF408"}
TISSUES = ["Cells_Cultured_fibroblasts", "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Mucosa"]
ENDPOINTS = [("H7_CORNEALDYSTROPHY", "H18.4", 491, 473095), ("H7_KERATOCONUS", "H18.5", 1076, 473095)]

z = {}
for t in TISSUES:
    for r in csv.DictReader(open(f"{B}/results/sprediXcan_alltissue_remap/fecd_{t}.csv")):
        g = GENES.get(r["gene"].split(".")[0])
        if g and r["zscore"] not in ("", "NA"): z[(g, t)] = float(r["zscore"])
print(f"discovery Neff={disc:,.0f}; two-sided critical |Z| at alpha=0.025: {crit:.3f}")
with open(f"{B}/analysis/nr1h3_finngen_lookup_power.tsv", "w") as fo:
    fo.write("gene\tendpoint\tICD10\ttissue\tdiscovery_z\texpected_abs_z\tpower\n")
    for ep, icd, ca, co in ENDPOINTS:
        ratio = sqrt(neff(ca, co) / disc)
        for g in ("ZNF408", "NR1H3"):
            for t in TISSUES:
                if (g, t) not in z: continue
                ez = abs(z[(g, t)]) * ratio
                pw = N.cdf(ez - crit) + N.cdf(-ez - crit)
                fo.write(f"{g}\t{ep}\t{icd}\t{t}\t{z[(g, t)]:.4f}\t{ez:.4f}\t{pw:.4f}\n")
                print(f"{ep:20} {g:7} {t:32} Neff={neff(ca, co):6,.0f} expected |Z|={ez:.2f} power={pw:.2f}")
