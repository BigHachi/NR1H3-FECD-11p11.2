#!/usr/bin/env python3
# ============================================================
# Step 19f (30 Oct 2026): ranking of the European FECD S-MultiXcan (step 19e) and the S_MultiXcan_rank
#   column of Supplementary Table S6.
#   - ranks every gene with an S-MultiXcan P-value by ascending P -> analysis/fecd_eur_smultixcan_ranking.tsv
#   - S6 genes (data/tophits21/genes_coloc.tsv): every gene with European S-MultiXcan P < 1e-4 (23 genes; rule
#     stated this way from 13 Nov 2026, the gene list is unchanged), plus genes kept from earlier versions of S6.
#     Checks that each gene with P < 1e-4 is in the list, then writes each gene's European rank and P into the
#     list's rank and smx_p columns (NA where the gene has no S-MultiXcan P-value); steps 102c and 81c carry
#     the rank into S6.
#     -> data/tophits21/genes_coloc.tsv (rank, smx_p), analysis/fecd_eur_smultixcan_S6_genes.tsv
#   Table 1 and Supplementary Table S2 use the multi-ancestry S-MultiXcan (step 19b), not this ranking.
# Usage: python3 19f_fecd_eur_smultixcan_rank.py [S-MultiXcan output]
# ============================================================
import csv, os, sys

B = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
os.chdir(B)
SMX = sys.argv[1] if len(sys.argv) > 1 else "results/multixcan/fecd_eur_smultixcan_remapped.txt"
GENES = "data/tophits21/genes_coloc.tsv"
P_SEL = 1e-4   # genes with a European S-MultiXcan P below this are colocalised


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


smx = tsv(SMX)
rows = [r for r in smx if num(r["pvalue"]) is not None]
rows.sort(key=lambda r: num(r["pvalue"]))
rank = {r["gene"].split(".")[0]: (i, r) for i, r in enumerate(rows, 1)}
os.makedirs("analysis", exist_ok=True)
with open("analysis/fecd_eur_smultixcan_ranking.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["rank", "gene", "gene_name", "pvalue", "n", "n_indep"])
    for i, r in enumerate(rows, 1):
        w.writerow([i, r["gene"], r["gene_name"], r["pvalue"], r["n"], r["n_indep"]])
print(f"genes with an S-MultiXcan P: {len(rows)} (of {len(smx)})")

genes = tsv(GENES)
cols = list(genes[0].keys())
nr = sum(num(r["pvalue"]) < P_SEL for r in rows)
listed = {g["ensg"] for g in genes}
missing = [(i, r["gene_name"]) for e, (i, r) in rank.items() if i <= nr and e not in listed]
print(f"genes with European S-MultiXcan P < {P_SEL:g}: {nr}; missing from the S6 list: {missing or 'none'}")
assert not missing, missing
for g in genes:
    i, r = rank.get(g["ensg"], ("NA", None))
    g["rank"], g["smx_p"] = i, (r["pvalue"] if r else "NA")
genes.sort(key=lambda g: (g["rank"] == "NA", g["rank"] if g["rank"] != "NA" else 0))
with open(GENES, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(genes)
with open("analysis/fecd_eur_smultixcan_S6_genes.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["gene", "ensg", "European_rank", "European_P", "selection", "n_genes_ranked"])
    for g in genes:
        sel = "P < 1e-4" if g["rank"] != "NA" and g["rank"] <= nr else "kept from earlier S6"
        w.writerow([g["gene"], g["ensg"], g["rank"], g["smx_p"], sel, len(rows)])
        if sel != "P < 1e-4":
            print(f"  {g['gene']:16s} European rank {g['rank']} ({sel})")
