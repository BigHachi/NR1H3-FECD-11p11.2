#!/usr/bin/env python3
# ============================================================
# Step 19c (Oct 10 2026): Tables 1-2 and the TWAS numbers quoted in the Results,
#   from the canonical remapped runs (step 95 S-PrediXcan via 95a; step 19b S-MultiXcan).
#   Table 2 (from manuscript v3.9; file table1_*): NR1H3 per-tissue S-PrediXcan (Z, beta, SE = |beta/Z|, P, n SNPs) + BH FDR
#            within tissue (manuscript: FDR 0.36 fib, 0.053 skin, 0.95 oesophagus;
#            ZNF408 0.029 fib, 0.022 skin).
#   Table 1 (from manuscript v3.9; file table2_*): S-MultiXcan genes ranked 1..NR1H3, with "P value of locus" = smallest
#            multi-ancestry FECD GWAS P within +/-1 Mb of the gene TSS (TSS from step 109).
#            Cytoband and FECD locus label are curated annotations (below).
#   Gene labels follow data/gene_symbol_display.tsv (same as Figures 1-3).
# ============================================================
import csv, gzip, os

B = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
TISSUES = ["Cells_Cultured_fibroblasts", "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Mucosa"]
disp = {r["ENSG"]: r["display_symbol"] for r in csv.DictReader(open(f"{B}/data/gene_symbol_display.tsv"), delimiter="\t")}
ens = lambda g: g.split(".")[0]

def bh(ps):
    n = len(ps); order = sorted(range(n), key=lambda i: ps[i]); q = [0.0] * n; m = 1.0
    for rank, i in reversed(list(enumerate(order, 1))):
        m = min(m, ps[i] * n / rank); q[i] = m
    return q

# ---- manuscript Table 2 (per-tissue NR1H3) + FDR ----
t1 = []
for T in TISSUES:
    rows = [r for r in csv.DictReader(open(f"{B}/results/sprediXcan_alltissue_remap/fecd_{T}.csv")) if r["pvalue"] not in ("", "NA")]
    q = bh([float(r["pvalue"]) for r in rows])
    for r, qq in zip(rows, q):
        if r["gene_name"] in ("NR1H3", "ZNF408"):
            z, b = float(r["zscore"]), float(r["effect_size"])
            t1.append([T, r["gene_name"], r["n_snps_used"], f"{z:.2f}", f"{b:.3f}", f"{abs(b / z):.3f}", f"{float(r['pvalue']):.2g}", f"{qq:.2g}", len(rows)])
with open(f"{B}/analysis/table1_nr1h3_znf408_sprediXcan.tsv", "w") as f:
    f.write("tissue\tgene\tn_snps\tZ\tbeta\tSE\tP\tFDR_within_tissue\tgenes_tested\n")
    for r in t1: f.write("\t".join(map(str, r)) + "\n")
print("Table 2 (per-tissue NR1H3) / FDR"); [print("  ", *r) for r in t1]

# ---- manuscript Table 1 (S-MultiXcan ranking) ----
CURATED = {  # ENSG -> (cytoband, FECD locus label)
    "ENSG00000267325": ("18q21.2", "TCF4"), "ENSG00000196628": ("18q21.2", "TCF4"), "ENSG00000166510": ("18q21.2", "TCF4"),
    "ENSG00000135862": ("1q25.3", "LAMC1"), "ENSG00000058085": ("1q25.3", "LAMC1"),
    "ENSG00000255142": ("11p15.5", "PIDD1"), "ENSG00000177666": ("11p15.5", "PIDD1"), "ENSG00000177595": ("11p15.5", "PIDD1"),
    "ENSG00000255284": ("11p15.5", "PIDD1"), "ENSG00000274897": ("11p15.5", "PIDD1"), "ENSG00000184524": ("11p15.5", "PIDD1"),
    "ENSG00000177156": ("11p15.5", "PIDD1"), "ENSG00000177225": ("11p15.5", "PIDD1"), "ENSG00000177685": ("11p15.5", "PIDD1"),
    "ENSG00000177542": ("11p15.5", "PIDD1"),
    "ENSG00000079819": ("6q23.2", "AKAP7/EPB41L2"), "ENSG00000118507": ("6q23.2", "AKAP7/EPB41L2"),
    "ENSG00000175213": ("11p11.2", "NR1H3"), "ENSG00000025434": ("11p11.2", "NR1H3"),
    "ENSG00000237515": ("16p13.12", "SHISA9"), "ENSG00000182612": ("17q25.3", "TSPAN10"),
}
smx = [r for r in csv.DictReader(open(f"{B}/results/multixcan/fecd_smultixcan_remapped.txt"), delimiter="\t") if r["pvalue"] not in ("", "NA")]
smx.sort(key=lambda r: float(r["pvalue"]))
top = []
for r in smx:
    top.append(r)
    if r["gene_name"] == "NR1H3": break
tss = {r["ENSG"]: (r["chr"], int(r["tss"])) for r in csv.DictReader(open(f"{B}/analysis/fecd_twas_gene_locus_class.tsv"), delimiter="\t") if r["tss"] not in ("", "NA", "None")}
win = {ens(r["gene"]): tss[ens(r["gene"])] for r in top if ens(r["gene"]) in tss}
best = {e: None for e in win}
with gzip.open(f"{B}/data/gwas/fecd_gwas_full.tsv.gz", "rt") as f:
    h = f.readline().rstrip("\n").split("\t"); iC, iB, iP = h.index("CHR"), h.index("BP"), h.index("P")
    for line in f:
        x = line.split("\t")
        for e, (c, t) in win.items():
            if x[iC] == c and abs(int(x[iB]) - t) <= 1_000_000:
                try: p = float(x[iP])
                except ValueError: continue
                if best[e] is None or p < best[e]: best[e] = p
with open(f"{B}/analysis/table2_smultixcan_ranking.tsv", "w") as f:
    f.write("rank\tsymbol\tENSG\tsmultixcan_P\tcytoband\tFECD_locus\tlocus_min_GWAS_P_1Mb\n")
    print("Table 1 (S-MultiXcan ranking)")
    for i, r in enumerate(top, 1):
        e = ens(r["gene"]); cb, loc = CURATED.get(e, ("NA", "NA")); lp = best.get(e)
        # GWAS P-values below 1e-300 sit at the floor of double precision (stored as ~1.1e-308 or 0), so they are
        # reported as "<1e-300" rather than as a spurious number (Oct 21 2026)
        lps = "NA" if lp is None else ("<1e-300" if lp < 1e-300 else f"{lp:.2g}")
        row = [i, disp.get(e, r["gene_name"]), e, f"{float(r['pvalue']):.3g}", cb, loc, lps]
        f.write("\t".join(map(str, row)) + "\n"); print("  ", *row)
print(f"S-MultiXcan genes with a result: {len(smx)}; Bonferroni 0.05/{len(smx)} = {0.05/len(smx):.2g}")
