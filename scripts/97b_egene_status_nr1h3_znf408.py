#!/usr/bin/env python3
# ============================================================
# Step 97b (Oct 9 2026; scripted Oct 10 2026): GTEx v8 eGene status of NR1H3 and ZNF408
#   Significant variant-gene pairs from the GTEx Portal API (singleTissueEqtl,
#   gtex_v8). A gene is an eGene in a tissue if any significant pair is returned.
#   Also flags genome-wide-significant (P < 5e-8) eQTLs and whether each is present
#   in the 1000 Genomes EUR LD reference used for MR clumping (HapMap3 SNPs only).
#   Manuscript: NR1H3 is an eGene in fibroblasts (lead P 2.4e-6), skin (1.2e-8),
#   oesophagus mucosa and whole blood (3.2e-7), not brain cortex; ZNF408 only in
#   oesophagus mucosa (and brain cortex). The skin GW-significant variants
#   (rs7119517, rs2203712) are absent from the LD reference.
#   Output: analysis/egene_status_nr1h3_znf408.tsv
# ============================================================
import json, os, urllib.request, collections

B = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
GENES = {"NR1H3": "ENSG00000025434.18", "ZNF408": "ENSG00000175213.2"}
TISSUES = ["Cells_Cultured_fibroblasts", "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Mucosa", "Whole_Blood", "Brain_Cortex"]
ldref = {l.split()[1] for l in open(f"{B}/data/ld_ref/LDREF/1000G.EUR.11.bim")}

out = []
for g, gid in GENES.items():
    url = f"https://gtexportal.org/api/v2/association/singleTissueEqtl?gencodeId={gid}&datasetId=gtex_v8&itemsPerPage=100000"
    rows = json.load(urllib.request.urlopen(url, timeout=120))["data"]
    by = collections.defaultdict(list)
    for r in rows: by[r["tissueSiteDetailId"]].append(r)
    for t in TISSUES:
        rr = sorted(by.get(t, []), key=lambda r: r["pValue"])
        gw = [r["snpId"] for r in rr if r["pValue"] < 5e-8]
        out.append([g, t, len(rr), "yes" if rr else "no",
                    rr[0]["snpId"] if rr else "", f"{rr[0]['pValue']:.2g}" if rr else "",
                    ",".join(gw), ",".join(s for s in gw if s not in ldref)])
with open(f"{B}/analysis/egene_status_nr1h3_znf408.tsv", "w") as f:
    f.write("gene\ttissue\tn_sig_pairs\teGene\tlead_snp\tlead_P\tGW_sig_snps\tGW_sig_snps_not_in_LDref\n")
    for r in out: f.write("\t".join(map(str, r)) + "\n")
for r in out: print(*r, sep="\t")
