#!/usr/bin/env python3
# ============================================================
# Step 109b (Oct 9 2026; scripted Oct 10 2026): 11p11.2 in the Insawang et al. 2026
#   European FECD meta-analysis (MVP + FinnGen + UKB/AoU; summary statistics shared by
#   B. Insawang, not public: H185_MVP_FinnGen_UKB+AoU_SNPsinatleast2studies.gz).
#   (a) the MVP 11p11.2 lead variant rs4752929; (b) the smallest P among variants
#   within 1 Mb of the NR1H3 TSS (positions from the MVP multi-ancestry GWAS, as the
#   meta-analysis file is rsID-only).
#   Manuscript (Discussion): rs4752929 P = 2.2e-3, same direction in all three
#   cohorts (+++); smallest P within 1 Mb of NR1H3 = 1.5e-5.
#   Not an independent replication (the meta-analysis includes MVP).
# ============================================================
import gzip, os

B = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
H185 = os.environ.get("H185_FILE", f"{B}/data/private/H185_MVP_FinnGen_UKB+AoU_SNPsinatleast2studies.gz")   # not public
OUT  = f"{B}/data/derived/insawang2026_11p11_lookup.tsv"   # bundled fallback (the values reported in the Discussion)
NR1H3_TSS, WIN = 47_248_421, 1_000_000

pos = {}
with gzip.open(f"{B}/data/gwas/fecd_gwas_full.tsv.gz", "rt") as f:
    h = f.readline().rstrip("\n").split("\t"); iS, iC, iB = h.index("SNP"), h.index("CHR"), h.index("BP")
    for line in f:
        x = line.split("\t")
        if x[iC] == "11" and abs(int(x[iB]) - NR1H3_TSS) <= WIN: pos[x[iS]] = int(x[iB])

if os.path.exists(H185):
    lead, hits = None, []
    with gzip.open(H185, "rt") as f:
        h = f.readline().rstrip("\n").split("\t"); c = {k: i for i, k in enumerate(h)}
        for line in f:
            x = line.rstrip("\n").split("\t"); s = x[c["MarkerName"]]
            if s == "rs4752929": lead = x
            if s in pos: hits.append((float(x[c["Pvalue"]]), s, pos[s], x[c["Direction"]]))
    hits.sort()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fo:          # refresh the bundled table from the shared summary statistics
        fo.write("row\trsid\tpos_GRCh38\tAllele1\tAllele2\tEffect\tStdErr\tPvalue\tDirection\tn_variants_within_1Mb\n")
        fo.write("\t".join(["MVP_lead", "rs4752929", "", *(lead[c[k]] for k in ("Allele1", "Allele2", "Effect", "StdErr", "Pvalue", "Direction")), str(len(hits))]) + "\n")
        for k, (p, s, b, d) in enumerate(hits[:3], 1):
            fo.write("\t".join([f"min_within_1Mb_{k}", s, str(b), "", "", "", "", f"{p:.3g}", d, str(len(hits))]) + "\n")
elif os.path.exists(OUT):
    print(f"Insawang summary statistics not available; reading {OUT}")
else:
    # 10 Nov 2026: the derived look-up table is not distributed with the public release
    print("Insawang summary statistics and the derived look-up table are not available: the 11p11.2 look-up "
          "reported in the Discussion needs the authors' file (H185_FILE); step skipped")
    raise SystemExit(0)
import csv
R = list(csv.DictReader(open(OUT), delimiter="\t"))
L = R[0]
print("rs4752929:", {k: L[k] for k in ("Allele1", "Allele2", "Effect", "StdErr", "Pvalue", "Direction")})
print(f"variants within 1 Mb of NR1H3 in both files: {L['n_variants_within_1Mb']}")
for r in R[1:]: print(f"  {r['rsid']}\tchr11:{r['pos_GRCh38']}\tP={r['Pvalue']}\t{r['Direction']}")
