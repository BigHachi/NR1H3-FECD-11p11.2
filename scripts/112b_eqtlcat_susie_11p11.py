#!/usr/bin/env python3
# ============================================================
# Step 112b: fine-mapped cis-eQTLs at 11p11.2 in eQTL Catalogue (r6, GTEx v8 re-analysis).
#   eQTL Catalogue fine-maps every gene with SuSiE using in-sample LD and publishes the 95% credible sets.
#   For cultured fibroblasts, skin (not sun exposed), oesophagus mucosa and whole blood, list every gene
#   with a credible set containing a variant in chr11:46.2-48.3 Mb, and record whether NR1H3 or ZNF408 has one.
#   Output: analysis/eqtlcat_susie_11p11.tsv
# ============================================================
import csv, gzip, json, os, socket, sys, time, urllib.request
socket.setdefaulttimeout(120)   # a stalled FTP connection must fail, not hang

def download(url, f):
    """Download to f.part and rename, so an interrupted download is never mistaken for a complete file."""
    for i in range(5):
        try:
            urllib.request.urlretrieve(url, f + ".part"); os.replace(f + ".part", f); return
        except OSError as e:
            print(f"download failed ({e}); retry {i + 1}/5", flush=True); time.sleep(15 * (i + 1))
    sys.exit(f"could not download {url}")
BASE = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
D = f"{BASE}/data/eqtl_catalogue/susie"; os.makedirs(D, exist_ok=True)
DATASETS = {"QTD000216": "Fibroblasts", "QTD000311": "Skin (not sun exposed)",
            "QTD000241": "Oesophagus mucosa", "QTD000356": "Whole blood"}
URL = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000015/{q}/{q}.credible_sets.tsv.gz"
TARGETS = {"ENSG00000025434": "NR1H3", "ENSG00000175213": "ZNF408"}
LO, HI = 46_200_000, 48_300_000

rows, names = [], {}
for q, tissue in DATASETS.items():
    f = f"{D}/{q}.credible_sets.tsv.gz"
    if not os.path.exists(f): download(URL.format(q=q), f)
    cs = {}
    with gzip.open(f, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c, p = r["variant"].split("_")[0], int(r["variant"].split("_")[1])
            if c == "chr11" and LO <= p <= HI:
                k = (r["gene_id"], r["cs_id"])
                cs.setdefault(k, {"n": int(r["cs_size"]), "top": None, "pip": -1.0, "minp": 1.0})
                if float(r["pip"]) > cs[k]["pip"]:
                    cs[k].update(top=r["rsid"] or r["variant"], pip=float(r["pip"]))
                cs[k]["minp"] = min(cs[k]["minp"], float(r["pvalue"]))
    for (g, c), v in cs.items(): rows.append([tissue, g, c, v["n"], v["top"], v["pip"], v["minp"]])
    for g, s in TARGETS.items():
        if not any(r[0] == tissue and r[1] == g for r in rows):
            rows.append([tissue, g, "none", 0, "", "", ""])

genes = sorted({r[1] for r in rows})
for i in range(0, len(genes), 50):          # gene symbols from Ensembl (POST lookup)
    body = json.dumps({"ids": genes[i:i + 50]}).encode()
    req = urllib.request.Request("https://rest.ensembl.org/lookup/id", data=body,
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
    for _ in range(5):
        try: names.update({k: (v or {}).get("display_name", k) for k, v in json.load(urllib.request.urlopen(req, timeout=60)).items()}); break
        except Exception: time.sleep(5)

with open(f"{BASE}/analysis/eqtlcat_susie_11p11.tsv", "w") as fo:
    fo.write("tissue\tgene\tENSG\tcredible_set\tcs_size\ttop_variant\ttop_PIP\tmin_eQTL_P\n")
    for r in sorted(rows, key=lambda r: (r[0], names.get(r[1], r[1]))):
        fo.write(f"{r[0]}\t{names.get(r[1], r[1])}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\n")
for t in DATASETS.values():
    have = sorted({names.get(r[1], r[1]) for r in rows if r[0] == t and r[2] != "none"})
    miss = [s for g, s in TARGETS.items() if not any(r[0] == t and r[1] == g and r[2] != "none" for r in rows)]
    print(f"{t}: genes with credible sets: {', '.join(have) or 'none'} | no credible set: {', '.join(miss)}")
