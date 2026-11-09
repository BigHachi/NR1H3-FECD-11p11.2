#!/usr/bin/env python3
# ============================================================
# Step 55c: Input audit of every GWAS used by the TWAS, conditional TWAS and cis-MR analyses
#   FECD multi-ancestry (GCST90480062): S-PrediXcan / S-MultiXcan input (95a), cis-MR sensitivity (83)
#   FECD European (GCST90475865): FUSION (82), 11p11.2 conditional FUSION (24b), primary cis-MR (91)
#   Glaucoma, IOP, refractive error: S-PrediXcan (46/51/56), whole-blood FUSION + chr11 conditional (60, 76d),
#     cross-trait cis-MR (89b/89c). Refractive error starts from the identity-merged file of step 55
#     (checked against the original METAL records by step 55b).
# Every analysis input written by the pipeline is linked back to the source GWAS records and to the
# reference variant set it is matched against:
#   FUSION       1000G EUR LD reference SNPs (rsID; FUSION aligns alleles to the reference)
#   S-PrediXcan  GTEx v8 MASHR model variants (varID chr_pos_ref_alt_b38)
#   cis-MR       eQTLGen NR1H3 cis-eQTL variants (rsID; alleles aligned by TwoSampleMR)
# and every reference variant is counted as matched, missing, ambiguous (>1 source record for the
# pipeline's key; the pipeline keeps the first) or allele-inconsistent. Each input record is also checked
# against the source record it came from (alleles and Z).
# The quantities behind the reported results are then re-derived from the audited inputs:
#   - prediction-model variants for NR1H3 and ZNF408 (S-PrediXcan, every tissue analysed) and the
#     weights of every FUSION model in the 11p11.2 conditional region (found, missing = imputed by
#     FUSION, allele-inconsistent), with the region rebuilt as FUSION.post_process.R builds it;
#   - the cis-MR instrument funnel (eQTLGen variants in the GWAS -> eQTLGen MAF > 0.01 -> |Z| > 5.45 -> PLINK
#     clump -> harmonised) and the instruments' source records.
# Outputs (analysis/input_audit/):
#   input_audit_summary.tsv      counts per trait x analysis x category
#   input_audit_records.tsv      every ambiguous-conflicting, allele- or position-inconsistent or unverified
#                                record used by a reported analysis (11p11.2 FUSION models, NR1H3/ZNF408
#                                S-PrediXcan models, eQTLGen NR1H3 cis-variants)
#   twas_model_coverage.tsv      per model (S-PrediXcan NR1H3/ZNF408; FUSION 11p11.2 region + NR1H3)
#   fusion_conditional_regions.tsv  the rebuilt 11p11.2 conditional regions
#   cis_instrument_coverage.tsv  instrument funnel per outcome (variants; "_rows" = rows after the rsID merge,
#                                which the cis-MR steps print in their logs)
#   mr_instruments_verified.tsv  each instrument: source record(s) vs the cis-MR input record
#   fusion_allele_sensitivity.tsv  only when a whole-blood conditional region holds an allele-inconsistent weight SNP:
#                                chr11 FUSION + conditional rerun (step 76d's commands) with every record whose
#                                alleles disagree with the LD reference replaced by an allele-consistent source record
#                                for the same rsID, compared with the reported run (FUSION outputs in
#                                results/input_audit_sensitivity/)
# Needs: numpy, scipy; Rscript (reads the FUSION .wgt.RDat files); PLINK 1.9.
# ============================================================
import glob, gzip, math, os, sqlite3, subprocess, tempfile, time
from array import array
from collections import Counter, defaultdict, namedtuple

import numpy as np
from scipy.stats import norm

BASE = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
os.chdir(BASE)
PLINK = os.environ.get("PLINK", f"{BASE}/software/plink/plink")
FUSION_DIR = os.environ.get("FUSION_DIR", f"{BASE}/software/fusion_twas")
OUT = "analysis/input_audit"
os.makedirs(OUT, exist_ok=True)
LDREF = "data/ld_ref/LDREF/1000G.EUR"
EQTLGEN = "data/eqtl_eqtlgen/nr1h3_eqtlgen.txt"
MASHR49 = sorted(glob.glob("data/predictdb_all/eqtl/mashr/mashr_*.db"))
MASHR5 = "data/predictdb"
NR1H3, ZNF408 = "ENSG00000025434", "ENSG00000175213"
AUTO = {str(c) for c in range(1, 23)}
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.0f}s]", *a, flush=True)


# ---------------------------------------------------------------- allele helpers
_COMP = str.maketrans("ACGT", "TGCA")
ACGT = set("ACGT")


def comp(a):
    return a.translate(_COMP)


def fusion_drops(a1, a2):
    """FUSION allele.qc keep = FALSE: strand-ambiguous pair or not a single A/C/G/T"""
    return not (a1 in ACGT and a2 in ACGT) or comp(a1) == a2


def orient(a1, a2, r1, r2):
    """+1 same (or complement), -1 swapped (or complement-swapped), 0 allele-inconsistent"""
    if (a1, a2) == (r1, r2) or (comp(a1), comp(a2)) == (r1, r2):
        return 1
    if (a1, a2) == (r2, r1) or (comp(a1), comp(a2)) == (r2, r1):
        return -1
    return 0


def close(x, y):
    return abs(x - y) <= 2e-5 * max(1.0, abs(y)) + 1e-6


Rec = namedtuple("Rec", "seq multi key chr pos a1 a2 z freq n")   # key = the identifier the pipeline matches on


def reproduces(a1, a2, z, r):
    """does an input record (alleles, Z) reproduce source record r?"""
    if (a1, a2) == (r.a1, r.a2):
        return close(z, r.z)
    if (a1, a2) == (r.a2, r.a1):
        return close(z, -r.z)
    return False


def same_record(r, s):
    return reproduces(s.a1, s.a2, s.z, r)


def fnum(x):
    try:
        v = float(x)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


# ---------------------------------------------------------------- source readers
# each yields Rec for valid autosomal records with the pipeline's key(s) (rsIDs); counters in cnt
def read_gwascat(path, q, cnt):
    """GWAS Catalog FECD files: Z = ln(OR) / SE, SE from the 95% CI (the SE column is #NA)"""
    with gzip.open(path, "rt") as fh:
        col = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
        iC, iP, iE, iO, iOR, iL, iU, iR, iF, iN = (col[k] for k in (
            "chromosome", "base_pair_location", "effect_allele", "other_allele", "odds_ratio", "ci_lower",
            "ci_upper", "rsid", "effect_allele_frequency", "n"))
        for seq, line in enumerate(fh):
            p = line.rstrip("\n").split("\t")
            cnt["source_rows"] += 1
            if p[iC] not in AUTO:
                cnt["source_non_autosomal"] += 1
                continue
            o, lo, up = fnum(p[iOR]), fnum(p[iL]), fnum(p[iU])
            if not (o and lo and up and o > 0 and lo > 0 and up > lo) or not p[iP].isdigit():
                cnt["source_invalid_statistics"] += 1
                continue
            z = math.log(o) / ((math.log(up) - math.log(lo)) / (2 * q))
            rs = p[iR]
            keys = [rs] if rs.startswith("rs") else []
            yield keys, Rec(seq, False, rs, p[iC], int(p[iP]), p[iE].upper(), p[iO].upper(), z, p[iF], p[iN])


def read_finngen(path, cnt):
    with gzip.open(path, "rt") as fh:
        col = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
        iC, iP, iR, iA, iS, iPv, iB, iSE, iF = (col[k] for k in (
            "#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"))
        for seq, line in enumerate(fh):
            p = line.rstrip("\n").split("\t")
            cnt["source_rows"] += 1
            c = p[iC].replace("chr", "")
            if c not in AUTO:
                cnt["source_non_autosomal"] += 1
                continue
            b, se, pv = fnum(p[iB]), fnum(p[iSE]), fnum(p[iPv])
            if b is None or se is None or pv is None or se <= 0:
                cnt["source_invalid_statistics"] += 1
                continue
            rs = p[iS]
            keys = rs.split(",") if rs.startswith("rs") else []
            # the pipeline's key is the whole rsids field: a record with several rsIDs never matches an rsID
            yield keys, Rec(seq, len(keys) > 1, rs, c, int(p[iP]), p[iA].upper(), p[iR].upper(), b / se, p[iF], "")


def read_iop(path, cnt):
    with gzip.open(path, "rt") as fh:
        col = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
        iC, iP, iS, iE, iO, iF, iB, iSE, iN = (col[k] for k in (
            "CHR", "POS", "SNP", "EffectAllele", "NonEffectAllele", "MAF", "Effect", "StdErr", "Ntot"))
        for seq, line in enumerate(fh):
            p = line.rstrip("\n").split("\t")
            cnt["source_rows"] += 1
            if p[iC] not in AUTO:
                cnt["source_non_autosomal"] += 1
                continue
            b, se = fnum(p[iB]), fnum(p[iSE])
            if b is None or se is None or se <= 0:
                cnt["source_invalid_statistics"] += 1
                continue
            rs = p[iS]
            keys = [rs] if rs.startswith("rs") else []
            # POS is GRCh37 (step 50 lifts it to GRCh38 for the analysis files)
            yield keys, Rec(seq, False, rs, p[iC], int(p[iP]), p[iE].upper(), p[iO].upper(), b / se, p[iF], p[iN])


def read_myopia(path, cnt):
    with open(path) as fh:
        col = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
        iS, iC, iP, i1, i2, iB, iF, iN = (col[k] for k in ("SNP", "CHR", "BP", "A1", "A2", "BETA", "FREQ", "N"))
        for seq, line in enumerate(fh):
            p = line.rstrip("\n").split("\t")
            cnt["source_rows"] += 1
            z = fnum(p[iB])
            if z is None:
                cnt["source_invalid_statistics"] += 1
                continue
            yield [p[iS]], Rec(seq, False, p[iS], p[iC], int(p[iP]), p[i1], p[i2], z, p[iF], p[iN])


# ---------------------------------------------------------------- trait definitions
WB = "results/twas_wb_conditional"
TRAITS = [
    dict(trait="FECD_multiancestry", label="FECD, multi-ancestry GWAS (GCST90480062)",
         source="data/gwas/raw/GCST90480062.tsv.gz", reader=lambda c: read_gwascat(
             "data/gwas/raw/GCST90480062.tsv.gz", 1.96, c),
         sprediXcan=("data/gwas/fecd_alltissue_sprediXcan.txt.gz", "position", MASHR49,
                     "results/sprediXcan_alltissue_remap/fecd_{T}.csv"),
         mr=("data/gwas/fecd_nr1h3_eqtlgen_region.tsv", False), mr_reported=None, mr_label="FECD_multiancestry"),
    dict(trait="FECD_European", label="FECD, European GWAS (GCST90475865)",
         source="data/gwas/fecd_eur_GCST90475865.tsv.gz", reader=lambda c: read_gwascat(
             "data/gwas/fecd_eur_GCST90475865.tsv.gz", 1.959964, c),
         fusion="data/gwas/fecd_eur_fusion.tsv",
         fusion_subsets=["results/fusion_11p11_conditional/fecd_eur_chr11.sumstats"],
         conditional=[("Cells_Cultured_fibroblasts", "results/fusion_11p11_conditional/Cells_Cultured_fibroblasts_11p11.dat",
                       "results/fusion_11p11_conditional/Cells_Cultured_fibroblasts_11p11_minp1"),
                      ("Skin_Not_Sun_Exposed_Suprapubic", "results/fusion_11p11_conditional/Skin_Not_Sun_Exposed_Suprapubic_11p11.dat",
                       "results/fusion_11p11_conditional/Skin_Not_Sun_Exposed_Suprapubic_11p11_minp1")],
         fusion_nr1h3=[f"results/twas_fecd_eurGWAS/{t}/fecd_eurGWAS_{t}_all_chr.dat" for t in
                       ("Cells_Cultured_fibroblasts", "Esophagus_Mucosa", "Skin_Not_Sun_Exposed_Suprapubic")],
         mr=("data/gwas/fecd_eur_nr1h3_eqtlgen_region.tsv", False), mr_reported=None, mr_label="FECD_European"),
    dict(trait="Glaucoma", label="Glaucoma, FinnGen R12 H7_GLAUCOMA",
         source="data/gwas/finngen_R12_H7_GLAUCOMA.gz",
         reader=lambda c: read_finngen("data/gwas/finngen_R12_H7_GLAUCOMA.gz", c),
         fusion="data/gwas/glaucoma_gwas_fusion.tsv", fusion_subsets=["data/gwas/glaucoma_gwas_fusion_chr11.tsv"],
         conditional=[("Whole_Blood", f"{WB}/glaucoma_wb_chr11.dat", f"{WB}/glaucoma_wb_chr11.joint")],
         sprediXcan=("data/gwas/glaucoma_gwas_sprediXcan.txt.gz", "rsid",
                     [f"{MASHR5}/mashr_{t}.db" for t in ("Cells_Cultured_fibroblasts", "Esophagus_Mucosa",
                      "Skin_Not_Sun_Exposed_Suprapubic", "Brain_Cortex", "Whole_Blood")],
                     "results/sprediXcan/glaucoma_{t}.csv"),
         mr=("data/gwas/glaucoma_gwas_full.tsv", False), mr_reported="Glaucoma", mr_label="Glaucoma"),
    dict(trait="IOP", label="Intraocular pressure (GCST009413)",
         source="data/gwas/iop_gwas_raw.txt.gz", reader=lambda c: read_iop("data/gwas/iop_gwas_raw.txt.gz", c),
         fusion="data/gwas/iop_gwas_fusion.tsv", fusion_subsets=["data/gwas/iop_gwas_fusion_chr11.tsv"],
         conditional=[("Whole_Blood", f"{WB}/iop_wb_chr11.dat", f"{WB}/iop_wb_chr11.joint")],
         sprediXcan=("data/gwas/iop_gwas_sprediXcan.txt.gz", "rsid",
                     [f"{MASHR5}/mashr_{t}.db" for t in ("Cells_Cultured_fibroblasts", "Esophagus_Mucosa",
                      "Skin_Not_Sun_Exposed_Suprapubic")], "results/sprediXcan/iop_{t}.csv"),
         mr=("data/gwas/iop_gwas_full.tsv", False), mr_reported="IOP", mr_label="IOP"),
    dict(trait="Refractive_error", label="Refractive error (Hysi 2020), identity-merged file of step 55",
         source="data/gwas/myopia_gwas_full.tsv", reader=lambda c: read_myopia("data/gwas/myopia_gwas_full.tsv", c),
         fusion="data/gwas/myopia_gwas_fusion.tsv", fusion_subsets=["data/gwas/myopia_gwas_fusion_chr11.tsv"],
         conditional=[("Whole_Blood", f"{WB}/myopia_wb_chr11.dat", f"{WB}/myopia_wb_chr11.joint")],
         sprediXcan=("data/gwas/myopia_gwas_sprediXcan.txt.gz", "rsid",
                     [f"{MASHR5}/mashr_{t}.db" for t in ("Cells_Cultured_fibroblasts", "Esophagus_Mucosa",
                      "Skin_Not_Sun_Exposed_Suprapubic")], "results/sprediXcan/myopia_{t}.csv"),
         mr=("data/gwas/myopia_gwas_full.tsv", True), mr_reported="Myopia", mr_label="Refractive_error"),
]

summary, records, model_cov, regions, funnel, instr, sens = [], [], [], [], [], [], []


def S(trait, analysis, inp, ref, cat, n):
    summary.append((trait, analysis, inp, ref, cat, n))


# ---------------------------------------------------------------- shared references
log("loading 1000G EUR LD reference")
BIM, BIM_CHR = {}, defaultdict(list)
for c in range(1, 23):
    for line in open(f"{LDREF}.{c}.bim"):
        ch, rs, _, pos, a1, a2 = line.split()
        BIM[rs] = (ch, int(pos), a1.upper(), a2.upper())
        BIM_CHR[ch].append(rs)

E = []          # eQTLGen NR1H3 rows
with open(EQTLGEN) as fh:
    col = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        E.append(dict(SNP=p[col["SNP"]], pos=int(p[col["SNPPos"]]), A1e=p[col["AssessedAllele"]].upper(),
                      A2e=p[col["OtherAllele"]].upper(), Z=float(p[col["Zscore"]]), N=int(p[col["NrSamples"]])))
# eQTLGen allele frequencies (step 83a): the cis-MR steps define the MAF filter on the exposure (eQTLGen) MAF
# from release v14 (10 Nov 2026; up to v13 on the outcome GWAS MAF)
MAF_E = {}
with open("data/eqtl_eqtlgen/nr1h3_eqtlgen_af.tsv") as fh:
    col = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        f = float(p[col["AlleleB_all"]])
        MAF_E[p[col["SNP"]]] = min(f, 1 - f)
for e in E:
    e["maf_e"] = MAF_E.get(e["SNP"], float("nan"))
E_RS = {e["SNP"] for e in E}
log(f"eQTLGen NR1H3 rows {len(E)}, unique SNPs {len(E_RS)}")

SYMBOL = {}
for db in MASHR49:
    con = sqlite3.connect(db)
    for g, s in con.execute("SELECT gene, genename FROM extra"):
        SYMBOL.setdefault(g.split(".")[0], s)
    con.close()
for line in open("data/gene_symbol_display.tsv"):
    p = line.rstrip("\n").split("\t")
    SYMBOL.setdefault(p[0], p[1])


def model_variants(dbs):
    """varID -> (chr, pos, ref, alt, rsid); per (tissue, gene) varIDs for NR1H3/ZNF408"""
    mv, gene_v = {}, defaultdict(list)
    rs_to_v = defaultdict(set)
    for db in dbs:
        t = os.path.basename(db)[6:-3]
        con = sqlite3.connect(db)
        for g, rs, v, ref, alt in con.execute("SELECT gene, rsid, varID, ref_allele, eff_allele FROM weights"):
            c, p, r1, r2, _ = v.split("_")
            mv[v] = (c.replace("chr", ""), int(p), r1, r2, rs)
            rs_to_v[rs].add(v)
            if g.split(".")[0] in (NR1H3, ZNF408):
                gene_v[(t, g.split(".")[0])].append((v, ref.upper(), alt.upper()))
        con.close()
    return mv, gene_v, rs_to_v


# ---------------------------------------------------------------- FUSION weights (R)
R_DUMP = r'''
args <- commandArgs(TRUE)
lst <- read.table(args[1], header=TRUE, sep="\t", as.is=TRUE, quote="", comment.char="")
con <- file(args[2], "w"); writeLines("FILE\tMODEL\tSNP\tA1\tA2\tW", con)
for (i in seq_len(nrow(lst))) {
  load(lst$FILE[i]); wgt.matrix[is.na(wgt.matrix)] <- 0
  m <- lst$MODEL[i]
  w <- if (!is.na(m) && m %in% colnames(wgt.matrix)) wgt.matrix[, m] else rep(NA_real_, nrow(snps))
  writeLines(paste(lst$FILE[i], m, snps[, 2], toupper(snps[, 5]), toupper(snps[, 6]), w, sep="\t"), con)
}
close(con)
'''


def read_dat(path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().split("\t")
        hdr = [h.strip() for h in hdr]
        for line in fh:
            p = [x.strip() for x in line.rstrip("\n").split("\t")]
            rows.append(dict(zip(hdr, p)))
    return rows


def dump_weights(pairs):
    """{(FILE, MODEL)} -> {(FILE, MODEL): [(SNP, A1, A2, W)]}"""
    with tempfile.TemporaryDirectory() as td:
        lst, out, rs = f"{td}/list.tsv", f"{td}/w.tsv", f"{td}/dump.R"
        with open(lst, "w") as f:
            f.write("FILE\tMODEL\n")
            for fl, m in sorted(pairs):
                f.write(f"{fl}\t{m}\n")
        open(rs, "w").write(R_DUMP)
        subprocess.run(["Rscript", rs, lst, out], check=True)
        W = defaultdict(list)
        with open(out) as fh:
            next(fh)
            for line in fh:
                fl, m, s, a1, a2, w = line.rstrip("\n").split("\t")
                W[(fl, m)].append((s, a1, a2, float(w) if w not in ("NA", "") else float("nan")))
    return W


def fusion_loci(dat_rows, W, win=100000):
    """FUSION.post_process.R: features = models with non-NA TWAS.P < 1 (minp_input 1); every LD-reference SNP
    in any feature's weights is marked; contiguous runs (in .bim order), +/- locus_win, overlaps merged"""
    feats = [r for r in dat_rows if r["TWAS.P"] not in ("NA", "") and float(r["TWAS.P"]) < 1.0]
    rs11 = BIM_CHR["11"]
    idx = {rs: i for i, rs in enumerate(rs11)}
    keep = np.zeros(len(rs11), bool)
    for r in feats:
        for s, *_ in W[(r["FILE"], r["MODEL"])]:
            if s in idx:
                keep[idx[s]] = True
    pos = np.array([BIM[rs][1] for rs in rs11])
    loci, i, n = [], 0, len(rs11)
    while i < n:
        if keep[i]:
            j = i
            while j + 1 < n and keep[j + 1]:
                j += 1
            loci.append([pos[i] - win, pos[j] + win])
            i = j + 1
        else:
            i += 1
    cons = [loci[0]]
    for s, e in loci[1:]:
        if s < cons[-1][1]:
            cons[-1][1] = max(cons[-1][1], e)
        else:
            cons.append([s, e])
    return feats, cons


# ---------------------------------------------------------------- per-trait audit
def classify_dups(keyed):
    """ambiguity of a key with >= 2 source records"""
    if len(keyed) < 2:
        return None
    return "identical" if all(same_record(keyed[0], r) for r in keyed[1:]) else "conflicting"


def audit(T):
    tr = T["trait"]
    log(f"===== {tr}: {T['label']}")
    want_rs = set(E_RS)
    if T.get("fusion"):
        want_rs |= set(BIM)
    sp = T.get("sprediXcan")
    if sp:
        MV, GENE_V, RS2V = model_variants(sp[2])
        want_pos = {(c, p) for (c, p, *_r) in MV.values()} if sp[1] == "position" else set()
        if sp[1] == "rsid":
            want_rs |= set(RS2V)
        # the S-PrediXcan input rows written under a model varID (read first: their own rsIDs are needed
        # to check each row against its source record)
        SPI, csp = {}, Counter()
        with gzip.open(sp[0], "rt") as fh:
            h = fh.readline().rstrip("\n").split("\t")
            iV, iR, iC, iP, iE, iN, iZ = (h.index(k) for k in ("panel_variant_id", "rsid", "chromosome", "position",
                                                              "effect_allele", "non_effect_allele", "zscore"))
            for line in fh:
                p = line.rstrip("\n").split("\t")
                csp["input_rows"] += 1
                if p[iV] in MV:
                    if p[iV] in SPI:
                        csp["input_duplicate_varid_rows"] += 1
                    else:
                        SPI[p[iV]] = (p[iR], p[iC].replace("chr", ""), int(p[iP]), p[iE].upper(), p[iN].upper(),
                                      float(p[iZ]))
        if sp[1] == "rsid":
            want_rs |= {x[0] for x in SPI.values()}
    else:
        want_pos = set()

    # ---- 1. source records
    cnt = Counter()
    by_rs, by_pos = defaultdict(list), defaultdict(list)
    rsnum = array("q")
    for keys, r in T["reader"](cnt):
        cnt["source_valid_autosomal"] += 1
        if not keys:
            cnt["source_without_rsid"] += 1
        elif r.multi:
            cnt["source_with_several_rsids_in_one_field"] += 1
        else:
            k = keys[0]
            if k[2:].isdigit():
                rsnum.append(int(k[2:]))
        for k in set(keys) | {r.key}:          # every rsID of the record, and the pipeline's own key
            if k in want_rs:
                by_rs[k].append(r)
        if want_pos and (r.chr, r.pos) in want_pos:
            by_pos[(r.chr, r.pos)].append(r)
    u, n = np.unique(np.frombuffer(rsnum, dtype=np.int64), return_counts=True)
    cnt["source_rsids_with_more_than_one_record"] = int((n > 1).sum())
    cnt["source_records_sharing_an_rsid"] = int(n[n > 1].sum())
    del rsnum, u, n
    for k in ("source_rows", "source_non_autosomal", "source_invalid_statistics", "source_valid_autosomal",
              "source_without_rsid", "source_with_several_rsids_in_one_field",
              "source_rsids_with_more_than_one_record", "source_records_sharing_an_rsid"):
        S(tr, "source", T["source"], "", k, cnt[k])
    log(f"source read: {dict(cnt)}")

    def keyed(k):
        """source records the pipeline would match on key k (its rsID, or a FinnGen multi-rsID field)"""
        return sorted((r for r in by_rs.get(k, ()) if r.key == k), key=lambda r: r.seq)

    # ---- 2. FUSION input vs 1000G EUR LD reference
    FIN = {}
    if T.get("fusion"):
        fin_file = T["fusion"]
        c = Counter()
        with open(fin_file) as fh:
            hdr = fh.readline().split()
            iS, i1, i2, iZ = (hdr.index(k) for k in ("SNP", "A1", "A2", "Z"))
            for line in fh:
                p = line.split()
                c["input_rows"] += 1
                if p[iS] in BIM:
                    if p[iS] in FIN:
                        c["input_duplicate_rsid_rows"] += 1
                        continue
                    z = fnum(p[iZ])
                    if z is None:
                        c["input_non_finite_z"] += 1
                        continue
                    FIN[p[iS]] = (p[i1].upper(), p[i2].upper(), z)
        for sub in T.get("fusion_subsets", []):
            bad = nsub = 0
            with open(sub) as fh:
                hdr = fh.readline().split()
                iS, i1, i2, iZ = (hdr.index(k) for k in ("SNP", "A1", "A2", "Z"))
                for line in fh:
                    p = line.split()
                    nsub += 1
                    x = FIN.get(p[iS])
                    if x is None or (p[i1].upper(), p[i2].upper()) != x[:2] or not close(float(p[iZ]), x[2]):
                        bad += 1
            S(tr, "FUSION", sub, "rows of the main FUSION input", "subset_rows", nsub)
            S(tr, "FUSION", sub, "rows of the main FUSION input", "subset_rows_not_identical_to_main_input", bad)
        for rs, (ch, pos, b1, b2) in BIM.items():
            kk = keyed(rs)
            x = FIN.get(rs)
            if x is None:
                cat = ("missing_source_record_not_carried" if kk else
                       "missing_source_record_under_multi_rsid_key" if by_rs.get(rs) else "missing_no_source_record")
            else:
                a1, a2, z = x
                if fusion_drops(a1, a2):
                    cat = "removed_by_fusion_strand_ambiguous_or_not_snv"
                elif orient(a1, a2, b1, b2):
                    cat = "matched"
                else:
                    cat = "allele_inconsistent_used_unaligned_by_fusion"
                ver = ("reproduces_first_source_record" if kk and reproduces(a1, a2, z, kk[0]) else
                       "reproduces_later_source_record" if any(reproduces(a1, a2, z, r) for r in kk[1:]) else
                       "reproduces_no_source_record")
                c[ver] += 1
                if ver == "reproduces_no_source_record":
                    records.append((tr, "FUSION", "LD reference SNP (any model)", rs, rs, "input_not_reproducing_source",
                                    "", len(kk), f"{a1}/{a2} Z={z:.4g}", "; ".join(f"{r.a1}/{r.a2} Z={r.z:.4g}" for r in kk)))
            c[cat] += 1
            d = classify_dups(kk)
            if d:
                c[f"ambiguous_{d}_source_records"] += 1
                if d == "conflicting" and x is not None:
                    ok = [r for r in kk if orient(r.a1, r.a2, b1, b2)]
                    if not orient(x[0], x[1], b1, b2) and ok:
                        c["ambiguous_conflicting_kept_record_inconsistent_other_consistent"] += 1
        S(tr, "FUSION", fin_file, f"1000G EUR LD reference SNPs ({len(BIM)})", "reference_variants", len(BIM))
        for k in ("input_rows", "input_duplicate_rsid_rows", "input_non_finite_z", "matched",
                  "allele_inconsistent_used_unaligned_by_fusion",
                  "removed_by_fusion_strand_ambiguous_or_not_snv", "missing_no_source_record",
                  "missing_source_record_not_carried", "missing_source_record_under_multi_rsid_key",
                  "ambiguous_identical_source_records", "ambiguous_conflicting_source_records",
                  "ambiguous_conflicting_kept_record_inconsistent_other_consistent",
                  "reproduces_first_source_record", "reproduces_later_source_record", "reproduces_no_source_record"):
            S(tr, "FUSION", fin_file, "1000G EUR LD reference SNPs", k, c[k])
        log(f"FUSION: {dict(c)}")

    # ---- 3. S-PrediXcan input vs model variants
    if sp:
        sp_file, link, dbs, res_pat = sp
        c = csp
        if link == "rsid":
            c["model_rsids_with_more_than_one_varid"] = sum(1 for v in RS2V.values() if len(v) > 1)
        status = {}
        for v, (ch, pos, r1, r2, mrs) in MV.items():
            x = SPI.get(v)
            src = (sorted((r for r in by_pos.get((ch, pos), ()) if {r.a1, r.a2} == {r1, r2}), key=lambda r: r.seq)
                   if link == "position" else keyed(mrs))
            if x is None:
                cat = "missing_source_record_not_carried" if src else "missing_no_source_record"
            else:
                rs, xc, xp, ea, nea, z = x
                if (xc, xp) != (ch, pos):
                    cat = "position_inconsistent"
                elif {ea, nea} != {r1, r2}:
                    cat = "allele_inconsistent_dropped_by_sprediXcan"
                else:
                    cat = "matched"
                if link == "rsid":     # rows reach a model varID through their rsID, or by position (chr_pos_A2_A1)
                    c["linked_by_model_rsid" if rs == mrs else "linked_by_position_under_another_rsid"] += 1
                vs = keyed(rs) if link == "rsid" else src
                ver = ("reproduces_first_source_record" if vs and reproduces(ea, nea, z, vs[0]) else
                       "reproduces_later_source_record" if any(reproduces(ea, nea, z, r) for r in vs[1:]) else
                       "reproduces_no_source_record")
                c[ver] += 1
                if ver == "reproduces_no_source_record":
                    records.append((tr, "S-PrediXcan", "model variant (any gene)", v, rs, "input_not_reproducing_source",
                                    "", len(vs), f"{ea}/{nea} Z={z:.4g}", "; ".join(f"{r.a1}/{r.a2} Z={r.z:.4g}" for r in vs)))
            c[cat] += 1
            d = classify_dups(src)
            if d:
                c[f"ambiguous_{d}_source_records"] += 1
            status[v] = (cat, d, x, src)
            if cat == "position_inconsistent":      # an rsID linked to a variant at another position: list every case
                records.append((tr, "S-PrediXcan", "model variant (any gene)", v, mrs, cat, d or "", len(src),
                                f"{x[0]} chr{x[1]}:{x[2]} {x[3]}/{x[4]} Z={x[5]:.4g}", ""))
        ref = f"GTEx v8 MASHR model variants, {len(dbs)} tissues ({len(MV)})"
        S(tr, "S-PrediXcan", sp_file, ref, "reference_variants", len(MV))
        for k in ("input_rows", "input_duplicate_varid_rows", "model_rsids_with_more_than_one_varid",
                  "linked_by_model_rsid", "linked_by_position_under_another_rsid", "matched",
                  "allele_inconsistent_dropped_by_sprediXcan", "position_inconsistent", "missing_no_source_record",
                  "missing_source_record_not_carried", "ambiguous_identical_source_records",
                  "ambiguous_conflicting_source_records", "reproduces_first_source_record",
                  "reproduces_later_source_record", "reproduces_no_source_record"):
            if k in ("model_rsids_with_more_than_one_varid", "linked_by_model_rsid",
                     "linked_by_position_under_another_rsid") and link != "rsid":
                continue
            S(tr, "S-PrediXcan", sp_file, ref, k, c[k])
        log(f"S-PrediXcan: {dict(c)}")
        # NR1H3 / ZNF408 model coverage, each tissue
        for (t, g), vs in sorted(GENE_V.items()):
            cc = Counter(status[v][0] for v, _, _ in vs)
            amb = sum(1 for v, _, _ in vs if status[v][1] == "conflicting")
            used = sum(1 for v, _, _ in vs if status[v][2] is not None and status[v][0] in ("matched", "position_inconsistent"))
            rep = ""
            f = res_pat.format(T=t, t=t.lower())
            if os.path.exists(f):
                with open(f) as fh:
                    h = fh.readline().rstrip("\n").split(",")
                    for line in fh:
                        p = line.rstrip("\n").split(",")
                        if p[0].split(".")[0] == g:
                            rep = p[h.index("n_snps_used")]
            model_cov.append(dict(trait=tr, method="S-PrediXcan", tissue=t, gene=SYMBOL.get(g, g), model_id=g,
                                  model="MASHR", n_weight_variants=len(vs), matched=cc["matched"],
                                  allele_inconsistent=cc["allele_inconsistent_dropped_by_sprediXcan"],
                                  position_inconsistent=cc["position_inconsistent"],
                                  missing=cc["missing_no_source_record"] + cc["missing_source_record_not_carried"],
                                  ambiguous_conflicting=amb, used_audit=used, used_reported=rep,
                                  agrees=str((int(float(rep)) if rep != "" else 0) == used)))   # absent = no SNP used
            for v, _, _ in vs:
                cat, d, x, src = status[v]
                if cat not in ("matched", "missing_no_source_record") or d == "conflicting":
                    records.append((tr, "S-PrediXcan", f"{t} {SYMBOL.get(g, g)}", v, MV[v][4], cat, d or "",
                                    len(src), "" if x is None else f"{x[3]}/{x[4]} Z={x[5]:.4g}",
                                    "; ".join(f"{r.a1}/{r.a2} Z={r.z:.4g}" for r in src)))
        del SPI, status

    # ---- 4. cis-MR input vs eQTLGen NR1H3 cis-variants
    mr_file, z_only = T["mr"]
    rows = defaultdict(list)
    with open(mr_file) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        iS, i1, i2, iB, iSE, iF = (h.index(k) for k in ("SNP", "A1", "A2", "BETA", "SE", "FREQ"))
        iN = h.index("N") if "N" in h else None
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[iS] in E_RS:
                b, se = float(p[iB]), float(p[iSE])
                rows[p[iS]].append(dict(A1=p[i1].upper(), A2=p[i2].upper(), Z=b if z_only else b / se,
                                        FREQ=float(p[iF]) if p[iF] not in ("NA", "") else float("nan"),
                                        N=p[iN] if iN is not None else ""))
    c = Counter()
    ed = {}
    for e in E:
        ed.setdefault(e["SNP"], e)
    for rs, e in ed.items():
        rr, kk = rows.get(rs, []), keyed(rs)
        if not rr:
            cat = "missing_source_record_not_carried" if kk else "missing_no_source_record"
        else:
            o = [orient(r["A1"], r["A2"], e["A1e"], e["A2e"]) for r in rr]
            cat = "matched" if all(o) else "allele_inconsistent_dropped_by_harmonisation"
            if len(rr) > 1:
                c["ambiguous_several_input_rows"] += 1
            for r in rr:
                ok = any(reproduces(r["A1"], r["A2"], r["Z"], s) for s in kk)
                c["reproduces_source_record" if ok else "reproduces_no_source_record"] += 1
                if not ok or cat != "matched":
                    records.append((tr, "cis-MR", "eQTLGen NR1H3", rs, rs, cat if ok else "input_not_reproducing_source",
                                    classify_dups(kk) or "", len(kk), f"{r['A1']}/{r['A2']} Z={r['Z']:.4g}",
                                    "; ".join(f"{s.a1}/{s.a2} Z={s.z:.4g}" for s in kk)))
        d = classify_dups(kk)
        if d:
            c[f"ambiguous_{d}_source_records"] += 1
            if d == "conflicting":
                records.append((tr, "cis-MR", "eQTLGen NR1H3", rs, rs, "ambiguous_conflicting_source_records", d,
                                len(kk), "; ".join(f"{r['A1']}/{r['A2']} Z={r['Z']:.4g}" for r in rr),
                                "; ".join(f"{s.a1}/{s.a2} Z={s.z:.4g}" for s in kk)))
        c[cat] += 1
    ref = f"eQTLGen NR1H3 cis-eQTL variants ({len(ed)})"
    S(tr, "cis-MR", mr_file, ref, "reference_variants", len(ed))
    for k in ("matched", "allele_inconsistent_dropped_by_harmonisation", "missing_no_source_record",
              "missing_source_record_not_carried", "ambiguous_several_input_rows", "ambiguous_identical_source_records",
              "ambiguous_conflicting_source_records", "reproduces_source_record", "reproduces_no_source_record"):
        S(tr, "cis-MR", mr_file, ref, k, c[k])
    log(f"cis-MR: {dict(c)}")
    mr_funnel(T, E, rows, keyed)

    # ---- 5. FUSION models: 11p11.2 conditional region(s) and NR1H3
    if T.get("conditional"):
        fusion_models(T, FIN, keyed)
    del by_rs, by_pos, FIN


def mr_funnel(T, E, rows, keyed):
    tr = T["trait"]
    # the cis-MR steps merge eQTLGen with the GWAS by rsID (all rows), keep eQTLGen MAF > 0.01 and |Z| > 5.45,
    # clump (r2 < 0.01, 5,000 kb, 1000G EUR chr11), then harmonise alleles (TwoSampleMR action 2)
    m = []
    for e in E:
        for r in rows.get(e["SNP"], []):
            m.append((e, r))
    in_gwas = {e["SNP"] for e, _ in m}
    m2 = [(e, r) for e, r in m if math.isfinite(e["maf_e"]) and e["maf_e"] > 0.01
          and math.isfinite(e["Z"]) and math.isfinite(r["Z"])]
    cands = [(e, r) for e, r in m2 if abs(e["Z"]) > 5.45]
    snps = sorted({e["SNP"] for e, _ in cands})
    pz = {e["SNP"]: e["Z"] for e, _ in cands}
    with tempfile.TemporaryDirectory() as td:
        with open(f"{td}/c.tsv", "w") as f:
            f.write("SNP\tP\n")
            for s in snps:
                f.write(f"{s}\t{2 * norm.sf(abs(pz[s])):.15g}\n")
        subprocess.run([PLINK, "--bfile", f"{LDREF}.11", "--clump", f"{td}/c.tsv", "--clump-p1", "5e-8",
                        "--clump-p2", "1", "--clump-r2", "0.01", "--clump-kb", "5000", "--out", f"{td}/o",
                        "--memory", "4000", "--silent"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        clumped = set()
        with open(f"{td}/o.clumped") as fh:
            h = fh.readline().split()
            for line in fh:
                p = line.split()
                if p:
                    clumped.add(p[h.index("SNP")])
    ins = [(e, r) for e, r in cands if e["SNP"] in clumped]
    harm = [(e, r) for e, r in ins if orient(r["A1"], r["A2"], e["A1e"], e["A2e"])]
    got = sorted({e["SNP"] for e, _ in harm})
    rep = ""
    if T["mr_reported"]:
        rep_set = set()
        with open("analysis/nr1h3_cismr_instruments_exposureN.tsv") as fh:
            h = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if p[h.index("outcome")] == T["mr_reported"]:
                    rep_set.add(p[h.index("SNP")])
        rep = ",".join(sorted(rep_set))
        same = str(set(got) == rep_set)
    else:   # FECD steps 83/91 report k only
        f = ("analysis/nr1h3_eqtlgen_fecd_mr_EUR_diagnostics.tsv" if tr == "FECD_European"
             else "analysis/nr1h3_fecd_cismr_ancestry_sensitivity.tsv")
        with open(f) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            p = fh.readline().rstrip("\n").split("\t")
        k = p[h.index("k")] if "k" in h else p[h.index("nsnp")]
        rep = f"k = {k}"
        same = str(int(k) == len(got))
    funnel.append(dict(outcome=T["mr_label"], mr_input=T["mr"][0], eqtlgen_rows=len(E),
                       eqtlgen_variants=len({e["SNP"] for e in E}), in_gwas=len(in_gwas),
                       input_rows_merged=len(m), maf_gt_0_01=len({e["SNP"] for e, _ in m2}), maf_gt_0_01_rows=len(m2),
                       genome_wide_significant=len(snps), genome_wide_significant_rows=len(cands), clumped=len(clumped),
                       after_harmonisation=len(got), instruments=",".join(got), reported=rep, same_as_reported=same))
    log(f"cis-MR funnel: {funnel[-1]}")
    for e, r in ins:
        kk = keyed(e["SNP"])
        o = orient(r["A1"], r["A2"], e["A1e"], e["A2e"])
        ok = len(kk) == 1 and reproduces(r["A1"], r["A2"], r["Z"], kk[0]) and o != 0 and \
            (not math.isfinite(r["FREQ"]) or fnum(kk[0].freq) is None or
             abs((fnum(kk[0].freq) if (r["A1"], r["A2"]) == (kk[0].a1, kk[0].a2) else 1 - fnum(kk[0].freq)) - r["FREQ"]) < 1.5e-4)
        s = kk[0] if kk else None
        instr.append(dict(outcome=T["mr_label"], SNP=e["SNP"], eqtl_assessed=e["A1e"], eqtl_other=e["A2e"],
                          eqtl_Z=e["Z"], n_source_records=len(kk),
                          source_A1="" if s is None else s.a1, source_A2="" if s is None else s.a2,
                          source_Z="" if s is None else f"{s.z:.6g}", source_FREQ="" if s is None else s.freq,
                          source_N="" if s is None else s.n,
                          input_A1=r["A1"], input_A2=r["A2"], input_Z=f"{r['Z']:.6g}", input_FREQ=f"{r['FREQ']:.6g}",
                          input_N=r["N"], orientation_vs_eqtl={1: "same", -1: "swapped", 0: "inconsistent"}[o],
                          verified="OK" if ok else "CHECK"))


def fusion_models(T, FIN, keyed):
    tr = T["trait"]
    dats = [(t, read_dat(d), d, pre) for t, d, pre in T["conditional"]]
    extra = []
    for d in T.get("fusion_nr1h3", []):
        t = d.split("/")[2]
        extra += [(t, r, d) for r in read_dat(d) if r["ID"].split(".")[0] == NR1H3]
    pairs = {(r["FILE"], r["MODEL"]) for _, rows, _, _ in dats for r in rows} | {(r["FILE"], r["MODEL"]) for _, r, _ in extra}
    W = dump_weights(pairs)
    log(f"FUSION weights read: {len(W)} models")

    def cover(t, r, region, joint):
        ws = W[(r["FILE"], r["MODEL"])]
        cc = Counter()
        used = []
        for s, a1, a2, w in ws:
            b = BIM.get(s)
            if b is None or b[0] != "11" and r["CHR"] == "11":
                cc["not_in_ld_reference"] += 1
                continue
            x = FIN.get(s)
            al = (x[0], x[1]) if x else (b[2], b[3])           # FUSION fills missing SNPs with the reference alleles
            if fusion_drops(*al):
                cc["removed_strand_ambiguous"] += 1
                continue
            used.append((s, w, x, b))
        if r["MODEL"] == "top1" and used:                     # FUSION keeps only the largest weight of a top1 model
            k = max(range(len(used)), key=lambda i: used[i][1] ** 2)
            used = [(s, w if i == k else 0.0, x, b) for i, (s, w, x, b) in enumerate(used)]
        for s, w, x, b in used:
            nz = w != 0 and math.isfinite(w)
            if x is None:
                st = "missing_imputed"
            elif orient(x[0], x[1], b[2], b[3]):
                st = "matched"
            else:
                st = "allele_inconsistent"
            d = classify_dups(keyed(s))
            cc[st] += 1
            if nz:
                cc[f"nonzero_{st}"] += 1
                cc["nonzero"] += 1
            if d == "conflicting":
                cc["ambiguous_conflicting"] += 1
                if nz:
                    cc["nonzero_ambiguous_conflicting"] += 1
            if st == "allele_inconsistent" or d == "conflicting":
                kk = keyed(s)
                records.append((tr, "FUSION", f"{t} {SYMBOL.get(r['ID'].split('.')[0], r['ID'])} ({r['MODEL']})", s, s,
                                st, d or "", len(kk), "" if x is None else f"{x[0]}/{x[1]} Z={x[2]:.4g}",
                                "; ".join(f"{q.a1}/{q.a2} Z={q.z:.4g}" for q in kk) + f" | LD ref {b[2]}/{b[3]}"))
        nsnp = len(used)
        g = r["ID"].split(".")[0]
        model_cov.append(dict(trait=tr, method="FUSION", tissue=t, gene=SYMBOL.get(g, g), model_id=r["ID"],
                              model=r["MODEL"], n_weight_variants=len(ws), matched=cc["matched"],
                              allele_inconsistent=cc["allele_inconsistent"], position_inconsistent="",
                              missing=cc["missing_imputed"], ambiguous_conflicting=cc["ambiguous_conflicting"],
                              used_audit=nsnp, used_reported=r["NSNP"],
                              agrees=str(r["NSNP"] not in ("NA", "") and int(float(r["NSNP"])) == nsnp),
                              impute_fraction=f"{cc['missing_imputed'] / nsnp:.3f}" if nsnp else "",
                              nonzero_weights=cc["nonzero"], nonzero_missing=cc["nonzero_missing_imputed"],
                              nonzero_allele_inconsistent=cc["nonzero_allele_inconsistent"],
                              nonzero_ambiguous_conflicting=cc["nonzero_ambiguous_conflicting"],
                              TWAS_Z=r["TWAS.Z"], TWAS_P=r["TWAS.P"], region=region, joint=joint))

    for t, rows, d, pre in dats:
        feats, loci = fusion_loci(rows, W)
        nr = [r for r in feats if r["ID"].split(".")[0] == NR1H3]
        if not nr:
            continue
        p0, p1 = int(nr[0]["P0"]), int(nr[0]["P1"])
        lo, hi = next((s, e) for s, e in loci if int(p0) < e and int(p1) > s)
        inloc = [r for r in rows if r["CHR"] == "11" and int(r["P0"]) < hi and int(r["P1"]) > lo]
        tested = [r for r in inloc if r in feats]
        joint = {r["ID"] for r in read_dat(pre + ".joint_included.dat")}   # FUSION.post_process --out prefix
        rep = ""
        if os.path.exists(pre + ".report"):
            rr = read_dat(pre + ".report")
            if rr:
                rep = f"report P0-P1 {rr[0]['P0']}-{rr[0]['P1']}, hit genes {rr[0]['HIT.GENES']}"
        region = f"chr11:{lo}-{hi}"
        uw = set()
        for r in inloc:
            for s, *_ in W[(r["FILE"], r["MODEL"])]:
                if s in BIM and BIM[s][0] == "11":
                    uw.add(s)
        cs = Counter()
        for s in uw:
            x = FIN.get(s)
            b = BIM[s]
            al = (x[0], x[1]) if x else (b[2], b[3])
            cs["strand_ambiguous" if fusion_drops(*al) else "missing" if x is None else
               "matched" if orient(x[0], x[1], b[2], b[3]) else "allele_inconsistent"] += 1
        regions.append(dict(trait=tr, tissue=t, conditional_input=d, region_GRCh37=region,
                            models_in_region=len(inloc), models_tested=len(tested),
                            models_in_joint_model=",".join(sorted(SYMBOL.get(r["ID"].split(".")[0], r["ID"])
                                                                  for r in inloc if r["ID"] in joint)),
                            weight_snps_in_ld_reference=len(uw), matched=cs["matched"], missing_imputed=cs["missing"],
                            allele_inconsistent=cs["allele_inconsistent"], strand_ambiguous_removed=cs["strand_ambiguous"],
                            fusion_report=rep))
        log(f"conditional region {t}: {regions[-1]}")
        for r in inloc:
            cover(t, r, region, str(r["ID"] in joint))
    for t, r, d in extra:
        cover(t, r, "genome-wide TWAS (S3, Figure 1)", "")
    for reg in [x for x in regions if x["trait"] == tr and x["tissue"] == "Whole_Blood"]:
        if reg["allele_inconsistent"] or reg["strand_ambiguous_removed"]:
            fusion_sensitivity(T, FIN, keyed)


def fusion_sensitivity(T, FIN, keyed):
    """step 76d rerun on chr11 with allele-consistent source records where the kept record disagrees with the
    LD reference (multi-allelic sites listed once per allele under one rsID; the pipeline keeps the first)"""
    tr = T["trait"]
    sub = T["fusion_subsets"][0]
    od = f"results/input_audit_sensitivity/{tr}"
    os.makedirs(od, exist_ok=True)
    n = 0
    with open(sub) as fi, open(f"{od}/fusion_chr11_allele_consistent.tsv", "w") as fo:
        fo.write(fi.readline())
        for line in fi:
            s, a1, a2, z, nn = line.split()
            b = BIM.get(s)
            if b and (fusion_drops(a1, a2) or not orient(a1, a2, b[2], b[3])):
                ok = [r for r in keyed(s) if not fusion_drops(r.a1, r.a2) and orient(r.a1, r.a2, b[2], b[3])]
                if ok:
                    a1, a2, z = ok[0].a1, ok[0].a2, repr(ok[0].z)
                    n += 1
            fo.write(f"{s}\t{a1}\t{a2}\t{z}\t{nn}\n")
    with open(f"{BASE}/data/weights_eur/GTExv8.EUR.Whole_Blood.pos") as fi, open(f"{od}/wb_chr11.pos", "w") as fo:
        fo.write(fi.readline())
        for line in fi:
            if line.split("\t")[3] == "11":
                fo.write(line)
    ld = f"{LDREF}."
    with open(f"{od}/fusion.log", "w") as lg:
        subprocess.run(["Rscript", f"{FUSION_DIR}/FUSION.assoc_test.R", "--sumstats", f"{od}/fusion_chr11_allele_consistent.tsv",
                        "--weights", f"{od}/wb_chr11.pos", "--weights_dir", "data/weights_eur", "--ref_ld_chr", ld,
                        "--chr", "11", "--out", f"{od}/wb_chr11.dat"], check=True, stdout=lg, stderr=lg)
        subprocess.run(["Rscript", f"{FUSION_DIR}/FUSION.post_process.R", "--sumstats", f"{od}/fusion_chr11_allele_consistent.tsv",
                        "--input", f"{od}/wb_chr11.dat", "--out", f"{od}/wb_chr11.joint", "--ref_ld_chr", ld, "--chr", "11",
                        "--locus_win", "100000"], check=True, stdout=lg, stderr=lg)
    t0 = T["conditional"][0]
    a = {r["ID"]: r for r in read_dat(t0[1])}
    b = {r["ID"]: r for r in read_dat(f"{od}/wb_chr11.dat")}
    dz = []
    for k, r in a.items():
        za, zb = fnum(r["TWAS.Z"]), fnum(b[k]["TWAS.Z"])
        if za is not None and zb is not None and za != zb:
            dz.append((abs(za - zb), k))
    dz.sort(reverse=True)
    ja = {r["ID"]: r for r in read_dat(t0[2] + ".joint_included.dat")}
    jb = {r["ID"]: r for r in read_dat(f"{od}/wb_chr11.joint.joint_included.dat")}
    nr = next(k for k in a if k.split(".")[0] == NR1H3)
    big = dz[0] if dz else (0.0, "")
    sens.append(dict(trait=tr, records_replaced_chr11=n, models_chr11=len(a), models_with_changed_TWAS_Z=len(dz),
                     largest_abs_change_TWAS_Z=f"{big[0]:.4f}",
                     largest_change_gene=SYMBOL.get(big[1].split(".")[0], big[1]),
                     largest_change_position=a[big[1]]["P0"] if big[1] else "",
                     NR1H3_TWAS_Z_reported=a[nr]["TWAS.Z"], NR1H3_TWAS_Z_corrected=b[nr]["TWAS.Z"],
                     joint_genes_reported=",".join(sorted(SYMBOL.get(k.split(".")[0], k) for k in ja)),
                     joint_genes_corrected=",".join(sorted(SYMBOL.get(k.split(".")[0], k) for k in jb)),
                     NR1H3_joint_reported="" if nr not in ja else f"Z {ja[nr]['JOINT.Z']} P {ja[nr]['JOINT.P']}",
                     NR1H3_joint_corrected="" if nr not in jb else f"Z {jb[nr]['JOINT.Z']} P {jb[nr]['JOINT.P']}"))
    log(f"allele-consistent FUSION sensitivity: {sens[-1]}")


ONLY = os.environ.get("AUDIT_TRAITS")      # debugging: a comma-separated subset of traits
for T in TRAITS:
    if not ONLY or T["trait"] in ONLY.split(","):
        audit(T)


# ---------------------------------------------------------------- write
def wtsv(path, hdr, rows):
    with open(path, "w") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows:
            vals = [r.get(h, "") for h in hdr] if isinstance(r, dict) else list(r)
            f.write("\t".join(map(str, vals)) + "\n")


wtsv(f"{OUT}/input_audit_summary.tsv", ["trait", "analysis", "input", "reference", "category", "n"], summary)
wtsv(f"{OUT}/input_audit_records.tsv", ["trait", "analysis", "used_by", "reference_variant", "rsid", "status",
                                         "source_duplicates", "n_source_records", "input_record", "source_records"],
     records)
mc_hdr = ["trait", "method", "tissue", "gene", "model_id", "model", "n_weight_variants", "used_audit", "used_reported",
          "agrees", "matched", "missing", "allele_inconsistent", "position_inconsistent", "ambiguous_conflicting",
          "impute_fraction", "nonzero_weights", "nonzero_missing", "nonzero_allele_inconsistent",
          "nonzero_ambiguous_conflicting", "TWAS_Z", "TWAS_P", "region", "joint"]
wtsv(f"{OUT}/twas_model_coverage.tsv", mc_hdr, model_cov)
wtsv(f"{OUT}/fusion_conditional_regions.tsv", list(regions[0].keys()) if regions else ["trait"], regions)
wtsv(f"{OUT}/cis_instrument_coverage.tsv", list(funnel[0].keys()) if funnel else ["trait"], funnel)
wtsv(f"{OUT}/mr_instruments_verified.tsv", list(instr[0].keys()) if instr else ["trait"], instr)
if sens:
    wtsv(f"{OUT}/fusion_allele_sensitivity.tsv", list(sens[0].keys()), sens)
elif os.path.exists(f"{OUT}/fusion_allele_sensitivity.tsv"):   # left by an earlier run whose inputs needed it
    os.remove(f"{OUT}/fusion_allele_sensitivity.tsv")
log("written:", ", ".join(sorted(os.listdir(OUT))))
