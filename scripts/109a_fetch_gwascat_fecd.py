#!/usr/bin/env python3
# ============================================================
# Step 109a (Oct 8 2026; scripted Oct 10 2026): FECD associations from the GWAS Catalog
#   REST API v2 for every published FECD GWAS analysis, written to
#   data/gwascat/gwascat_fecd.json (input to step 109, Supplementary Table S12).
#   Note: the Catalog is updated continuously; re-fetching later may add associations.
# ============================================================
import json, os, urllib.request

B = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS")); os.makedirs(f"{B}/data/gwascat", exist_ok=True)
studies = {"GCST000790": "Baratz 2010", "GCST004217": "Afshari 2017",
           "GCST90428877": "Gorman 2024 (EUR)", "GCST90428878": "Gorman 2024 (multi)",
           "GCST90428879": "Gorman 2024 (AFR)", "GCST90428880": "Gorman 2024 (HIS)",
           "GCST90475865": "Verma 2024 (EUR)", "GCST90476531": "Verma 2024 (AFR)", "GCST90480062": "Verma 2024 (multi)"}
rows = []
for acc, lab in studies.items():
    page = 0
    while True:
        u = f"https://www.ebi.ac.uk/gwas/rest/api/v2/associations?accession_id={acc}&size=200&page={page}"
        d = json.load(urllib.request.urlopen(u, timeout=90))
        for a in d.get("_embedded", {}).get("associations", []):
            rows.append([lab, acc, ";".join(a.get("snp_effect_allele") or []), a.get("p_value"),
                         ";".join(a.get("mapped_genes") or []), ";".join(str(x) for x in (a.get("locations") or []))])
        page += 1
        if page >= d.get("page", {}).get("totalPages", 1): break
    print(lab, acc, sum(1 for r in rows if r[1] == acc))
json.dump(rows, open(f"{B}/data/gwascat/gwascat_fecd.json", "w"))
