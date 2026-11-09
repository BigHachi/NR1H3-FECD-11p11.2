#!/usr/bin/env python3
# ============================================================
# Step 139 (30 Oct 2026): additional outputs requested in peer review (29 Oct 2026), collected in one folder with
#   the tables derived for the request, so that they ship with the release (analysis/additional_outputs/):
#   glaucoma_allele_matching/         (10 Nov 2026) the glaucoma FUSION input is allele-matched from release v14
#                                     (step 60b); this folder compares it with the first-record input reported up
#                                     to v13. allele_matching_records.tsv = every record replaced or dropped;
#     chr11/                          whole-blood FUSION + conditional analysis on chr11 (step 76d; S8): reported
#                                     (allele-matched) and first-record runs; chr11_first_record_vs_reported.tsv
#                                     puts all 463 models side by side
#     genomewide/                     step 76h, all 22 chromosomes, both inputs: marginal results, jointly
#                                     included and dropped genes, Bonferroni threshold, logs
#   fecd_european_smultixcan/         the European S-MultiXcan behind S6 (steps 19d-19f) with its per-tissue
#                                     S-PrediXcan results; the 4 Jul 2026 rsID-matched run (if present, for
#                                     reference); S6_rank_lookup.tsv = each S6 gene's rank in both runs and in the
#                                     multi-ancestry S-MultiXcan (Table 1)
#   rnaseq_batch_metadata/            basis for assigning the GSE101872 controls to batch 1 (step 110): GEO sample
#                                     annotations, gene-count identity of the five shared FECD samples, and the SRA
#                                     run information (fetched from NCBI E-utilities on first run, then cached)
#   Every file here is copied from, or derived from, outputs of earlier steps or the downloaded inputs; nothing
#   is re-estimated.
# Usage: python3 139_additional_outputs.py
# ============================================================
import csv, glob, gzip, os, shutil, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET

B = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
OUT = os.environ.get("ADDITIONAL_OUT", f"{B}/analysis/additional_outputs")
TISSUES = ("Cells_Cultured_fibroblasts", "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Mucosa")
if os.path.exists(OUT):
    shutil.rmtree(OUT)


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def write(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    print(f"{os.path.relpath(path, OUT)}: {len(rows)} rows")


def copy(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)


def fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------- glaucoma chr11, first record (<= v13) vs reported
A = f"{OUT}/glaucoma_allele_matching"
D = f"{A}/chr11"
copy(f"{B}/analysis/glaucoma_fusion_allele_matching.tsv", f"{A}/allele_matching_records.tsv")
REP, COR = (f"{B}/results/twas_wb_conditional/previous_first_record/glaucoma_wb_chr11",
            f"{B}/results/twas_wb_conditional/glaucoma_wb_chr11")
for s in (".dat", ".joint.joint_included.dat", ".joint.joint_dropped.dat"):
    copy(REP + s, f"{D}/first_record_v13/glaucoma_wb_chr11{s}")
    copy(COR + s, f"{D}/reported_allele_matched/glaucoma_wb_chr11{s}")
sym = {r["ENSG"]: r["symbol"] for r in tsv(f"{B}/data/derived/s8_chr11_wb_gene_symbols.tsv")}


def run(prefix):
    dat = {r["ID"]: r for r in tsv(prefix + ".dat")}
    st = {}
    for r in tsv(prefix + ".joint.joint_included.dat"):
        st[r["ID"]] = ("jointly_included", r["JOINT.Z"].strip(), r["JOINT.P"].strip())
    for r in tsv(prefix + ".joint.joint_dropped.dat"):
        st[r["ID"]] = ("conditioned_or_not_tested", r["COND.Z"].strip(), r["COND.P"].strip())
    return dat, st


rep, rep_st = run(REP)
cor, cor_st = run(COR)
assert set(rep) == set(cor)
rows = []
for gid in rep:
    a, b = rep[gid], cor[gid]
    za, zb = fl(a["TWAS.Z"]), fl(b["TWAS.Z"])
    rows.append([gid, sym.get(gid, ""), a["CHR"], a["P0"], a["P1"], a["NSNP"], b["NSNP"], a["MODEL"],
                 a["TWAS.Z"].strip(), b["TWAS.Z"].strip(), "" if za is None or zb is None else f"{zb - za:.5f}",
                 a["TWAS.P"].strip(), b["TWAS.P"].strip(), *rep_st.get(gid, ("", "", "")), *cor_st.get(gid, ("", "", ""))])
rows.sort(key=lambda r: fl(r[11]) if fl(r[11]) is not None else 2)
write(f"{D}/chr11_first_record_vs_reported.tsv",
      ["ID", "symbol", "CHR", "P0", "P1", "NSNP_first_record", "NSNP_reported", "MODEL", "TWAS.Z_first_record",
       "TWAS.Z_reported", "delta_Z", "TWAS.P_first_record", "TWAS.P_reported", "status_first_record",
       "joint_or_cond_Z_first_record", "joint_or_cond_P_first_record", "status_reported", "joint_or_cond_Z_reported",
       "joint_or_cond_P_reported"], rows)

# ---------------------------------------------------------------- glaucoma genome-wide (step 76h), both inputs
G = f"{A}/genomewide"


def cat(files, out):
    with open(out, "w") as fo:
        for i, f in enumerate(files):
            with open(f) as fi:
                h = fi.readline()
                if i == 0:
                    fo.write(h)
                fo.writelines(fi)


for s, lab in (("as_reported", "first_record_v13"), ("allele_corrected", "reported_allele_matched")):
    d, o = f"{B}/results/twas_wb_genomewide/{s}", f"{G}/{lab}"
    chrs = lambda suf: [f"{d}/glaucoma_wb_chr{c}{suf}" for c in range(1, 23) if os.path.exists(f"{d}/glaucoma_wb_chr{c}{suf}")]
    copy(f"{d}/glaucoma_wb_all_chr.dat", f"{o}/glaucoma_wb_all_chr.dat")
    copy(f"{d}/bonferroni.txt", f"{o}/bonferroni.txt")
    cat(chrs(".joint.joint_included.dat"), f"{o}/glaucoma_wb_all_chr.joint_included.dat")
    cat(chrs(".joint.joint_dropped.dat"), f"{o}/glaucoma_wb_all_chr.joint_dropped.dat")
    with open(f"{o}/post_process_logs.txt", "w") as fo:
        for f in sorted(glob.glob(f"{d}/post_chr*.log"), key=lambda x: int(x.rsplit("chr", 1)[1][:-4])):
            fo.write(f"==== {os.path.basename(f)}\n" + open(f).read())
    print(f"{os.path.relpath(o, OUT)}/: {len(os.listdir(o))} files")

# ---------------------------------------------------------------- European S-MultiXcan (steps 19d-19f)
F = f"{OUT}/fecd_european_smultixcan"
copy(f"{B}/results/multixcan/fecd_eur_smultixcan_remapped.txt", f"{F}/fecd_eur_smultixcan_remapped.txt")
for t in TISSUES:
    copy(f"{B}/results/sprediXcan_fecd_eur_remap_staged/mashr_{t}.csv",
         f"{F}/per_tissue_S-PrediXcan/fecd_eur_sprediXcan_remapped_{t}.csv")
JUL = f"{B}/results/multixcan/fecd_eur_smultixcan.txt"      # 4 Jul 2026 run (rsID-matched input); not rerun by this release
if os.path.exists(JUL):
    P = f"{F}/previous_run_4Jul2026_rsID_matched"
    copy(JUL, f"{P}/fecd_eur_smultixcan.txt")
    for t in TISSUES:
        copy(f"{B}/results/sprediXcan_fecd_eur_staged/mashr_{t}.csv", f"{P}/per_tissue_S-PrediXcan/fecd_eur_sprediXcan_{t}.csv")
else:
    print("4 Jul 2026 European run not present: previous-run folder and its rank columns left empty")


def ranks(path):
    if not os.path.exists(path):
        return {}
    rr = [r for r in tsv(path) if r["pvalue"] not in ("", "NA", "nan")]
    rr.sort(key=lambda r: float(r["pvalue"]))
    return {r["gene"].split(".")[0]: (i, r["pvalue"]) for i, r in enumerate(rr, 1)}


jul = ranks(JUL)
eur = ranks(f"{B}/results/multixcan/fecd_eur_smultixcan_remapped.txt")
mul = ranks(f"{B}/results/multixcan/fecd_smultixcan_remapped.txt")
s6 = {}
for r in tsv(f"{B}/analysis/fecd_coloc_S6_fullcis.tsv"):            # step 81c (S6)
    s6.setdefault(r["ENSG"], (r["SYMBOL"], r["S_MultiXcan_rank"]))
rows = []
for e, (sym6, rk) in s6.items():
    rk = "" if rk in ("", "NA") else int(rk)
    assert rk == eur.get(e, ("",))[0], (e, rk)
    jr = jul.get(e, ("no P",))[0] if jul else ""
    old = (jr if isinstance(jr, int) and jr <= 21 else 0) if jul else ""     # S6 rank in v2.15/v2.16: ranks 1-21, else 0
    rows.append([sym6, e, old, jr, rk, eur.get(e, ("no P", ""))[0], eur.get(e, ("", ""))[1],
                 mul.get(e, ("no P", ""))[0], mul.get(e, ("", ""))[1]])
rows.sort(key=lambda r: (r[4] == "", r[4] if r[4] != "" else 0))
write(f"{F}/S6_rank_lookup.tsv",
      ["S6_SYMBOL", "ENSG", "S6_rank_v2.15_v2.16", "European_rank_4Jul2026_rsID_matched", "S6_rank_v2.17",
       f"European_rank_position_allele_matched_of_{len(eur)}", "European_P", f"multiancestry_rank_of_{len(mul)}",
       "multiancestry_P"], rows)

# ---------------------------------------------------------------- RNA-seq batch basis (step 110)
R = f"{OUT}/rnaseq_batch_metadata"
os.makedirs(R)
E = f"{B}/data/endothelium_rnaseq"
meta = []
for f in ("GSE101872_series_matrix.txt.gz", "GSE112201-GPL11154_series_matrix.txt.gz",
          "GSE112201-GPL20301_series_matrix.txt.gz"):
    fields = {}
    for line in gzip.open(f"{E}/{f}", "rt", errors="ignore"):
        if not line.startswith("!Sample_"):
            continue
        p = [x.strip('"') for x in line.rstrip("\n").split("\t")]
        fields.setdefault(p[0], []).append(p[1:])
    for i in range(len(fields["!Sample_geo_accession"][0])):
        ch = [c[i] for c in fields.get("!Sample_characteristics_ch1", [])]
        get = lambda k: next((c.split(": ", 1)[1] for c in ch if c.startswith(k + ":")), "")
        sra = next((x[i].split("term=")[-1] for x in fields.get("!Sample_relation", []) if "SRA" in x[i]), "")
        meta.append(["GSE101872" if f.startswith("GSE101872") else "GSE112201", fields["!Sample_geo_accession"][0][i],
                     fields["!Sample_title"][0][i], fields["!Sample_platform_id"][0][i],
                     fields["!Sample_instrument_model"][0][i], get("batch") or "not given",
                     get("group") or get("subject status"), get("re"), fields["!Sample_submission_date"][0][i], sra])
write(f"{R}/geo_sample_metadata.tsv",
      ["GEO_series", "GSM", "GEO_title", "platform", "instrument", "GEO_batch_label", "GEO_group", "GEO_re",
       "GEO_submission_date", "SRA_experiment"], meta)


def counts(f):
    r = csv.reader(gzip.open(f"{E}/{f}", "rt"), delimiter="\t")
    h = next(r)
    return h, {tuple(x[:3]): x for x in r}


h1, c1 = counts("GSE101872_GeneCounts.tsv.gz")
h2, c2 = counts("GSE112201_GeneCounts.tsv.gz")
keys = sorted(set(c1) & set(c2))
rows = []
for a, b in [("2011-020", "RNA10"), ("2011-024", "RNA16"), ("2011-038", "RNA20"), ("2011-041", "RNA15"), ("6004", "RNA19")]:
    i, j = h1.index(a), h2.index(b)
    rows.append([a, b, len(keys), sum(c1[k][i] == c2[k][j] for k in keys),
                 sum(int(c1[k][i]) for k in keys), sum(int(c2[k][j]) for k in keys)])
write(f"{R}/shared_donor_count_identity.tsv",
      ["GSE101872_sample", "GSE112201_batch1_sample", "genes_compared", "genes_with_identical_counts",
       "total_counts_GSE101872", "total_counts_GSE112201"], rows)

# SRA run information for every sample's SRA experiment (NCBI E-utilities, cached in data/endothelium_rnaseq/)
X = f"{E}/sra_runinfo.xml"
if not os.path.exists(X):
    U = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    term = " OR ".join(f"{m[9]}[accn]" for m in meta if m[9])
    q = urllib.request.urlopen(f"{U}/esearch.fcgi?db=sra&retmax=200&term={urllib.parse.quote(term)}", timeout=60).read()
    ids = [i.text for i in ET.fromstring(q).iter("Id")]
    time.sleep(1)
    xml = urllib.request.urlopen(f"{U}/efetch.fcgi?db=sra&rettype=runinfo&id={','.join(ids)}", timeout=120).read()
    open(X, "wb").write(xml)
geo = {m[9]: m for m in meta if m[9]}
FLD = ["Experiment", "Run", "Model", "LibraryLayout", "avgLength", "spots", "LoadDate", "SampleName", "BioSample"]
rows = []
for r in ET.parse(X).getroot().findall("Row"):
    d = {c.tag: (c.text or "") for c in r}
    g = geo.get(d.get("Experiment"))
    if g:
        rows.append([g[0], g[2], g[5]] + [d.get(f, "") for f in FLD])
rows.sort(key=lambda x: x[3])
assert len({x[3] for x in rows}) == len(geo), "SRA run information missing for some samples"
write(f"{R}/sra_runinfo_extract.tsv", ["GEO_series", "GEO_title", "GEO_batch_label"] + FLD, rows)
print("additional outputs ->", OUT)
