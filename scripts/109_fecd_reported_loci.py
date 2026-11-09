#!/usr/bin/env python3
# ============================================================
# Step 109 (review point 1): every reported FECD GWAS locus + classification
#   of each TWAS gene against it.
#   Loci sources (GRCh38):
#     - GWAS Catalog associations (P<5e-8) for Baratz 2010, Afshari 2017,
#       Gorman 2024 (4 analyses), Verma 2024 (3 analyses)   [scratch json]
#     - Genome-wide loci re-derived from summary statistics, because the
#       Catalog curates only lead SNPs:  Verma/MVP multi-ancestry + EUR, and
#       the Insawang et al. 2026 meta-analysis (H185 file from B. Insawang;
#       rsID-only -> positions from the MVP file).
#   Locus = greedy clumping by distance: best SNP P<5e-8, drop +/-1 Mb, repeat.
#   A TWAS gene is "locus-recapturing" if its TSS lies within 1 Mb of any
#   reported lead SNP; otherwise "not previously reported".
# ============================================================
import gzip, json, csv, os, sys, urllib.error, urllib.request
BASE = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
SCR  = sys.argv[1] if len(sys.argv) > 1 else f"{BASE}/data/gwascat"   # gwascat_fecd.json from step 109a
H185 = os.environ.get("H185_FILE", f"{BASE}/data/private/H185_MVP_FinnGen_UKB+AoU_SNPsinatleast2studies.gz")   # not public
H185_LEADS = os.environ.get("H185_LEADS_FILE", f"{BASE}/data/derived/insawang2026_clumped_lead_variants.tsv")   # derived from the shared file (not in the public release)
WIN, GWS = 1_000_000, 5e-8

def clump(snps):  # snps: list of (chr, pos, p, rsid)
    out = []
    for s in sorted(snps, key=lambda x: x[2]):
        if all(not (s[0] == o[0] and abs(s[1] - o[1]) <= WIN) for o in out):
            out.append(s)
    return out

loci = []  # (source, chr, pos, lead_rsid, P, mapped_gene)
# ---- GWAS Catalog ----
for lab, acc, snp, p, genes, locs in json.load(open(f"{SCR}/gwascat_fecd.json")):
    if p is None or float(p) >= GWS or not locs: continue
    c, pos = locs.split(";")[0].split(":")
    loci.append((f"{lab} [GWAS Catalog {acc}]", c, int(pos), snp.split("-")[0], float(p), genes))

# ---- MVP multi-ancestry (CHR BP P SNP) ----
gw = []
with gzip.open(f"{BASE}/data/gwas/fecd_gwas_full.tsv.gz", "rt") as f:
    h = f.readline().rstrip("\n").split("\t"); iC, iB, iP, iS = (h.index(k) for k in ("CHR", "BP", "P", "SNP"))
    for line in f:
        x = line.split("\t")
        try: p = float(x[iP])
        except ValueError: continue
        if p < GWS: gw.append((x[iC], int(x[iB]), p, x[iS]))
for s in clump(gw): loci.append(("Verma 2024 multi-ancestry [summary statistics]", s[0], s[1], s[3], s[2], ""))

# ---- MVP EUR ----
gw = []
with gzip.open(f"{BASE}/data/gwas/fecd_eur_GCST90475865.tsv.gz", "rt") as f:
    h = f.readline().rstrip("\n").split("\t")
    iC, iB, iP, iS = (h.index(k) for k in ("chromosome", "base_pair_location", "p_value", "rsid"))
    for line in f:
        x = line.rstrip("\n").split("\t")
        try: p = float(x[iP])
        except ValueError: continue
        if p < GWS: gw.append((x[iC], int(x[iB]), p, x[iS]))
for s in clump(gw): loci.append(("Verma 2024 EUR [summary statistics]", s[0], s[1], s[3], s[2], ""))

# ---- Insawang 2026 (H185 meta; rsID only) ----
# With the shared summary statistics, clump them and refresh the bundled lead-variant table;
# without them (anyone else), read that table, which gives the same loci.
if os.path.exists(H185):
    h185 = {}
    with gzip.open(H185, "rt") as f:
        h = f.readline().rstrip("\n").split("\t"); iS, iP = h.index("MarkerName"), h.index("Pvalue" if "Pvalue" in h else "P-value")
        for line in f:
            x = line.rstrip("\n").split("\t")
            try: p = float(x[iP])
            except ValueError: continue
            if p < GWS: h185[x[iS]] = p
    pos = {}
    with gzip.open(f"{BASE}/data/gwas/fecd_gwas_full.tsv.gz", "rt") as f:
        h = f.readline().rstrip("\n").split("\t"); iC, iB, iS = (h.index(k) for k in ("CHR", "BP", "SNP"))
        for line in f:
            x = line.split("\t")
            for r in x[iS].split(","):
                if r in h185 and r not in pos: pos[r] = (x[iC], int(x[iB]))
    gw = [(pos[r][0], pos[r][1], p, r) for r, p in h185.items() if r in pos]
    print(f"Insawang/H185: {len(h185)} GW SNPs, {len(gw)} positioned via MVP file", file=sys.stderr)
    leads = clump(gw)
    os.makedirs(os.path.dirname(H185_LEADS), exist_ok=True)
    with open(H185_LEADS, "w") as fo:
        fo.write("chr\tpos_GRCh38\tlead_rsid\tP\n")
        for c_, b_, p_, r_ in leads: fo.write(f"{c_}\t{b_}\t{r_}\t{p_:.3g}\n")
elif os.path.exists(H185_LEADS):
    leads = [(r["chr"], int(r["pos_GRCh38"]), float(r["P"]), r["lead_rsid"])
             for r in csv.DictReader(open(H185_LEADS), delimiter="\t")]
    print(f"Insawang/H185: summary statistics not available; using {len(leads)} lead variants from {H185_LEADS}", file=sys.stderr)
else:
    # 10 Nov 2026: the lead-variant table derived from the authors' file is not distributed with the public release.
    # Without either file the loci re-derived from the Insawang et al. summary statistics are left out; their
    # published loci still enter through the GWAS Catalog (step 109a).
    leads = []
    print("Insawang/H185: neither the summary statistics nor the derived lead-variant table is available; "
          "loci re-derived from them are left out (Supplementary Table S12 lists them)", file=sys.stderr)
for s in leads: loci.append(("Insawang 2026 [H185 summary statistics]", s[0], s[1], s[3], s[2], ""))

# ---- consolidate into locus regions (union across sources) ----
regions = []
for L in sorted(loci, key=lambda x: (x[1].zfill(2), x[2])):
    for R in regions:
        if R["chr"] == L[1] and R["start"] - WIN <= L[2] <= R["end"] + WIN:
            R["start"] = min(R["start"], L[2]); R["end"] = max(R["end"], L[2]); R["hits"].append(L); break
    else:
        regions.append({"chr": L[1], "start": L[2], "end": L[2], "hits": [L]})
print(f"{len(loci)} lead associations -> {len(regions)} reported FECD loci", file=sys.stderr)

# ---- gene coordinates for TWAS genes (S-MultiXcan Bonferroni + Table 2 + per-tissue Bonferroni) ----
genes = {}
smx = [x for x in csv.DictReader(open(f"{BASE}/results/multixcan/fecd_smultixcan_remapped.txt"), delimiter="\t") if x["pvalue"] not in ("", "NA")]
smx.sort(key=lambda x: float(x["pvalue"])); bonf = 0.05 / len(smx)
for k, x in enumerate(smx, 1):
    if float(x["pvalue"]) < bonf or k <= 21: genes[x["gene"]] = {"symbol": x["gene_name"], "smx_rank": k, "smx_p": float(x["pvalue"]), "spx": []}
for T in ("Cells_Cultured_fibroblasts", "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Mucosa"):
    R = [x for x in csv.DictReader(open(f"{BASE}/results/sprediXcan_alltissue_remap/fecd_{T}.csv")) if x["pvalue"] not in ("", "NA")]
    b = 0.05 / len(R)
    for x in R:
        if float(x["pvalue"]) < b:
            g = genes.setdefault(x["gene"], {"symbol": x["gene_name"], "smx_rank": None, "smx_p": None, "spx": []})
            g["spx"].append(T.split("_")[0])
req = urllib.request.Request("https://rest.ensembl.org/lookup/id", data=json.dumps({"ids": list(genes)}).encode(),
                             headers={"Content-Type": "application/json", "Accept": "application/json"})
def ensembl_lookup(req, tries=int(os.environ.get("ENSEMBL_TRIES", 8))):
    """Ensembl REST is occasionally unavailable (HTTP 429/5xx, timeouts): retry with back-off (10 Nov 2026)"""
    import time
    for k in range(tries):
        try:
            return json.load(urllib.request.urlopen(req, timeout=120))
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as e:
            if k == tries - 1 or (isinstance(e, urllib.error.HTTPError) and e.code not in (429, 500, 502, 503, 504)):
                raise
            wait = min(300, 15 * 2 ** k)
            print(f"Ensembl REST: {e}; retrying in {wait} s", file=sys.stderr)
            time.sleep(wait)
# Fallback (10 Nov 2026): the gene positions this release used (Ensembl REST lookup of 10 Nov 2026, GRCh38), read
# when Ensembl cannot be reached; a live lookup can also differ after an Ensembl release changes a gene model
TSS_FILE = f"{BASE}/data/derived/ensembl_gene_tss.tsv"
try:
    coord = ensembl_lookup(req)
except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as e:
    print(f"Ensembl REST unavailable ({e}); using the gene positions in {TSS_FILE}", file=sys.stderr)
    coord = {r["ENSG"]: {"seq_region_name": r["chr"], "start": int(r["tss_GRCh38"]), "end": int(r["tss_GRCh38"]), "strand": 1}
             for r in csv.DictReader(open(TSS_FILE), delimiter="\t")}
GTEX_FALLBACK = {"ENSG00000274897": ("11", 818000), "ENSG00000279672": ("11", 779617)}  # retired IDs

with open(f"{BASE}/analysis/fecd_reported_loci.tsv", "w") as fo:
    fo.write("locus_id\tchr\tstart_lead_bp\tend_lead_bp\tn_lead_associations\tsources\tlead_rsids\tbest_P\tmapped_genes\n")
    for i, R in enumerate(regions, 1):
        R["id"] = f"L{i:02d}"
        src = sorted({h[0].split(" [")[0] for h in R["hits"]})
        fo.write(f"{R['id']}\t{R['chr']}\t{R['start']}\t{R['end']}\t{len(R['hits'])}\t{'; '.join(src)}\t"
                 f"{';'.join(sorted({h[3] for h in R['hits']}))}\t{min(h[4] for h in R['hits']):.2g}\t"
                 f"{';'.join(sorted({g for h in R['hits'] for g in h[5].split(';') if g}))}\n")
with open(f"{BASE}/analysis/fecd_twas_gene_locus_class.tsv", "w") as fo:
    fo.write("ENSG\tsymbol\tchr\ttss\tsmx_rank\tsmx_p\tsprediXcan_bonferroni_tissues\tlocus_id\tdistance_to_nearest_lead_kb\tnearest_lead_sources\tclassification\n")
    for e, g in sorted(genes.items(), key=lambda kv: (kv[1]["smx_rank"] or 1e9)):
        c = coord.get(e)
        if c: ch, tss = c["seq_region_name"], (c["start"] if c["strand"] == 1 else c["end"])
        else: ch, tss = GTEX_FALLBACK.get(e, ("NA", None))
        best = None
        for R in regions:
            if R["chr"] != ch or tss is None: continue
            for h in R["hits"]:
                d = abs(h[2] - tss)
                if best is None or d < best[0]: best = (d, R, h)
        inloc = best is not None and best[0] <= WIN
        fo.write(f"{e}\t{g['symbol']}\t{ch}\t{tss}\t{g['smx_rank'] or ''}\t{g['smx_p'] or ''}\t{','.join(g['spx'])}\t"
                 f"{best[1]['id'] if inloc else ''}\t{round(best[0]/1000) if best else ''}\t"
                 f"{'; '.join(sorted({h[0].split(' [')[0] for h in best[1]['hits']})) if inloc else ''}\t"
                 f"{'locus-recapturing' if inloc else 'not previously reported'}\n")
print("wrote analysis/fecd_reported_loci.tsv and analysis/fecd_twas_gene_locus_class.tsv", file=sys.stderr)
