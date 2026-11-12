#!/usr/bin/env python3
# ============================================================
# Step 140: check the reported tables against the pipeline outputs.
#   (1) Writes the supplementary tables that are pure exports of pipeline outputs, exactly as they appear in the
#       workbook (row order, rounding, derived columns): S1, S2, S3, S10, S11 and the S4 power column ->
#       analysis/supplementary/<sheet>.tsv
#   (2) Compares every data cell of the supplementary workbook (S1-S14) and of manuscript Tables 1-3 with the
#       outputs of this run. S1, S2, S3, S10 and S11 must match the exports exactly. Elsewhere a reported number must
#       equal the output value or a rounding of it to >= 2 significant figures, and text must match exactly;
#       captions, notes and the S9 Panel D narrative summary are not checked.
#   Reference copies of the reported tables are in reference/: the supplementary workbook and
#   manuscript_tables.tsv (Tables 1-3 as printed in the manuscript).
#   Output: analysis/reported_tables_check.tsv (one row per table; mismatches listed in
#   analysis/reported_tables_mismatches.tsv). Exits 1 if any reported value disagrees with the outputs.
# Needs pandas and openpyxl.
# ============================================================
import csv, glob, gzip, math, os, re, sys
import numpy as np, pandas as pd, openpyxl

B = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
XLSX = os.environ.get("SUPP_XLSX", f"{B}/reference/Supplementary_Tables_NR1H3_FECD.xlsx")
MTAB = os.environ.get("MANUSCRIPT_TABLES", f"{B}/reference/manuscript_tables.tsv")
OUTD = f"{B}/analysis/supplementary"; os.makedirs(OUTD, exist_ok=True)
T3 = ["Cells_Cultured_fibroblasts", "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Mucosa"]

# ---------------- helpers ----------------
NONE_LIKE = {"", "-", "—", "n/a", "NA", "nan", "None", "NaN"}
def num(x):
    if x is None: return None
    if isinstance(x, str) and x.strip() in NONE_LIKE: return None
    if isinstance(x, (int, float, np.integer, np.floating)):
        return None if (isinstance(x, float) and math.isnan(x)) else float(x)
    s = str(x).strip().replace("−", "-").replace("–", "-").replace(" ", "").replace(",", "")
    m = re.fullmatch(r"([+-]?[\d.]+)\s*[×x]\s*10([⁻⁺]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+)", s)
    if m:
        sup = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789-+")
        return float(m.group(1)) * 10 ** int(m.group(2).translate(sup))
    try: return float(s.lstrip("+"))
    except ValueError: return None

def sig(x, k):
    return float(f"{x:.{k}g}")

def rounding_of(v, s):
    """sheet value v is the output value s or s rounded to >= 2 significant figures (or integer-equal)."""
    if s is None: return v is None
    if v is None: return False
    if v == s or math.isclose(v, s, rel_tol=1e-12, abs_tol=0): return True
    if s == 0: return v == 0
    from decimal import Decimal, ROUND_HALF_UP
    for k in range(2, 16):
        if math.isclose(sig(s, k), v, rel_tol=1e-12): return True
        d = Decimal(repr(s))                    # half-up rounding of the printed value (e.g. 1.85e-04 -> 1.9e-04)
        e = d.adjusted() - (k - 1)
        if math.isclose(float(d.quantize(Decimal(1).scaleb(e), rounding=ROUND_HALF_UP)), v, rel_tol=1e-12): return True
    for d in range(0, 7):                       # fixed decimals (e.g. 0.07, 1.0)
        if math.isclose(round(s, d), v, rel_tol=1e-12, abs_tol=0) and abs(s) >= 0.5 * 10 ** (-d) * 1.0001: return True
    return False

def sig4(x):
    if x is None or (isinstance(x, float) and math.isnan(x)) or x == 0: return x
    return sig(float(x), 4)

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p, kind="stable")
    q = p[o] * n / np.arange(1, n + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out

def ensg(g): return re.sub(r"\..*", "", str(g))

WB = openpyxl.load_workbook(XLSX, read_only=True) if os.path.exists(XLSX) else None
def sheet(name, header_row=2, first=3, last=None, ncol=None):
    ws = WB[name]
    rows = list(ws.iter_rows(min_row=header_row, max_row=last, values_only=True))
    hdr = [h for h in rows[0] if h is not None] if ncol is None else list(rows[0][:ncol])
    data = [list(r[:len(hdr)]) for r in rows[first - header_row:] if any(v is not None for v in r[:len(hdr)])]
    return pd.DataFrame(data, columns=hdr)

RESULTS, MISMATCH = [], []
def record(table, n_rows, n_cells, bad, note=""):
    RESULTS.append((table, n_rows, n_cells, len(bad), note)); MISMATCH.extend((table,) + b for b in bad[:200])
    print(f"{table:34s} rows {n_rows:6d}  cells {n_cells:8d}  mismatches {len(bad):5d}  {note}")

def compare(table, rep, out, key, cols, exact=False, note="", subset=False):
    """rep = reported table, out = expected (from outputs); key columns join the two. subset=True: the table shows
    selected rows only (e.g. top 25), so output rows absent from the table are not flagged."""
    rep = rep.copy(); out = out.copy()
    for k in key:
        for df in (rep, out):
            df[k] = df[k].map(lambda x: "" if x is None or (isinstance(x, float) and math.isnan(x)) or str(x).strip() in NONE_LIKE else str(x).strip())
    m = rep.merge(out, on=key, how="outer", suffixes=("_r", "_o"), indicator=True)
    bad = [(" / ".join(str(r[k]) for k in key), "row", "missing in " + ("outputs" if r["_merge"] == "left_only" else "table"), "")
           for _, r in m[(m["_merge"] == "left_only") | ((m["_merge"] == "right_only") & (not subset))].iterrows()]
    m = m[m["_merge"] == "both"]; cells = 0
    for c in cols:
        for _, r in m.iterrows():
            v, s = r[c + "_r"], r[c + "_o"]; cells += 1
            nv, ns = num(v), num(s)
            if isinstance(s, str) and ns is None or isinstance(v, str) and nv is None and not (v in ("", "NA", "n/a")):
                norm = lambda x: "" if x is None or (isinstance(x, float) and math.isnan(x)) or str(x).strip() in NONE_LIKE else str(x).strip()
                ok = norm(v) == norm(s)
            elif exact:
                ok = (nv is None and ns is None) or (nv is not None and ns is not None and math.isclose(nv, ns, rel_tol=1e-12, abs_tol=0))
            else:
                ok = (nv is None and (ns is None or v in ("n/a", "NA", ""))) or rounding_of(nv, ns)
            if not ok: bad.append((" / ".join(str(r[k]) for k in key), c, str(v), str(s)))
    record(table, len(m), cells, bad, note)

# ================= exports: S1, S2, S10, S11 (exact) =================
SPX = f"{B}/results/sprediXcan_alltissue_remap"
def s1():
    parts = []
    for t in T3:
        d = pd.read_csv(f"{SPX}/fecd_{t}.csv"); d = d[d.pvalue.notna()].copy(); n = len(d)
        d["TISSUE"] = t; d["ENSG"] = d.gene.map(ensg); d["SYMBOL"] = d.gene_name; d["bonf"] = 0.05 / n
        d["FDR"] = bh(d.pvalue.values); d["sig"] = np.where(d.pvalue < d.bonf, "TRUE", "FALSE"); parts.append(d)
    a = pd.concat(parts, ignore_index=True)
    a["_p"] = a.pvalue.map(sig4); a = a.sort_values("_p", kind="stable").drop(columns="_p")
    for c in ["zscore", "effect_size", "pvalue", "var_g", "pred_perf_r2", "pred_perf_pval", "pred_perf_qval", "bonf", "FDR"]:
        a[c] = a[c].map(sig4)
    return a[["gene", "gene_name", "zscore", "effect_size", "pvalue", "var_g", "pred_perf_r2", "pred_perf_pval",
              "pred_perf_qval", "n_snps_used", "n_snps_in_cov", "n_snps_in_model", "TISSUE", "ENSG", "SYMBOL", "bonf",
              "FDR", "sig"]]

def s2():
    s = pd.read_csv(f"{B}/results/multixcan/fecd_smultixcan_remapped.txt", sep="\t"); s = s[s.pvalue.notna()].copy()
    s["ENSG"] = s.gene.map(ensg); s["SYMBOL"] = s.gene_name; s["FDR"] = bh(s.pvalue.values)
    s["_p"] = s.pvalue.map(sig4); s = s.sort_values("_p", kind="stable").drop(columns="_p")
    for c in ["pvalue", "n", "n_indep", "p_i_best", "p_i_worst", "eigen_max", "eigen_min", "eigen_min_kept", "z_min",
              "z_max", "z_mean", "z_sd", "tmi", "status", "FDR"]:
        s[c] = s[c].map(sig4)
    return s[["gene", "gene_name", "pvalue", "n", "n_indep", "p_i_best", "t_i_best", "p_i_worst", "t_i_worst", "eigen_max",
              "eigen_min", "eigen_min_kept", "z_min", "z_max", "z_mean", "z_sd", "tmi", "status", "ENSG", "SYMBOL", "FDR"]]

def s10():
    rows = []
    for f in sorted(glob.glob(f"{SPX}/fecd_*.csv")):
        t = os.path.basename(f)[5:-4]; d = pd.read_csv(f); d = d[d.pvalue.notna()].copy()
        d["rank_in_tissue"] = d.pvalue.rank(method="first").astype(int); d["ENSG"] = d.gene.map(ensg)
        for g, e in (("NR1H3", "ENSG00000025434"), ("ZNF408", "ENSG00000175213")):
            x = d[d.ENSG == e]
            if len(x):
                r = x.iloc[0]
                rows.append(dict(SYMBOL=g, ENSG=e, tissue=t, zscore=sig4(r.zscore), effect_size=sig4(r.effect_size),
                                 pvalue=sig(r.pvalue, 3), rank_in_tissue=int(r.rank_in_tissue), n_genes_tested=len(d),
                                 n_snps_used=int(r.n_snps_used), n_snps_in_model=int(r.n_snps_in_model), _p=r.pvalue))
    b = pd.DataFrame(rows); b["_g"] = b.SYMBOL.map({"NR1H3": 0, "ZNF408": 1})
    return b.sort_values(["_g", "_p"], kind="stable").drop(columns=["_g", "_p"])

def s11():
    s = pd.read_csv(f"{B}/results/multixcan/fecd_smultixcan_alltissue49.txt", sep="\t"); s = s[s.pvalue.notna()].copy()
    s["ENSG"] = s.gene.map(ensg); s["SYMBOL"] = s.gene_name; s["n_tissues"] = s.n.astype(int)
    s["n_indep"] = s.n_indep.astype(int); s["best_tissue"] = s.t_i_best.str.replace("^mashr_", "", regex=True)
    s["bonferroni_sig"] = np.where(s.pvalue < 0.05 / len(s), "yes", "no")
    s["_p"] = s.pvalue.map(sig4); s = s.sort_values("_p", kind="stable").drop(columns="_p")
    s["rank"] = np.arange(1, len(s) + 1)
    for c in ["pvalue", "z_min", "z_max", "z_mean"]: s[c] = s[c].map(sig4)
    return s[["rank", "ENSG", "SYMBOL", "pvalue", "n_tissues", "n_indep", "best_tissue", "z_min", "z_max", "z_mean",
              "bonferroni_sig"]]

def fusion_eur():
    """FUSION, European GWAS: the three tissues' outputs in file order (fibroblast, skin, oesophagus), models without a
    result (TWAS.P = NA) omitted."""
    fu = []
    for t in T3:
        d = pd.read_csv(f"{B}/results/twas_fecd_eurGWAS/{t}/fecd_eurGWAS_{t}_all_chr.dat", sep="\t",
                        skipinitialspace=True)          # FUSION right-aligns numbers and NA with spaces
        d.insert(0, "TISSUE", t); fu.append(d)
    fu = pd.concat(fu, ignore_index=True)
    return fu[fu["TWAS.P"].notna()].reset_index(drop=True)

def s3():
    """S3 = fusion_eur() with positions and counts exact and the other numbers to 4 significant figures."""
    fu = fusion_eur()
    for c in ["HSQ", "BEST.GWAS.Z", "EQTL.R2", "EQTL.Z", "EQTL.GWAS.Z", "MODELCV.R2", "MODELCV.PV", "TWAS.Z", "TWAS.P"]:
        fu[c] = fu[c].map(sig4)
    for c in ["CHR", "P0", "P1", "NSNP", "NWGT"]:
        fu[c] = fu[c].astype(int)
    return fu

EXPORTS = {"S1_FECD_S-PrediXcan_all": s1(), "S2_FECD_S-MultiXcan_all": s2(), "S3_FUSION_all": s3(),
           "S10_NR1H3_ZNF408_49tissues": s10(), "S11_S-MultiXcan_49tissues": s11()}
for name, df in EXPORTS.items():
    df.to_csv(f"{OUTD}/{name}.tsv", sep="\t", index=False, na_rep="")
print(f"exports written to {OUTD}: " + ", ".join(f"{k} ({len(v)} rows)" for k, v in EXPORTS.items()))

if WB is None:
    print(f"reference workbook {XLSX} not found: exports only"); sys.exit(0)

def exact_sheet(name, key):
    rep = sheet(name); out = EXPORTS[name]
    order_ok = len(rep) == len(out) and all(str(a) == str(b) for a, b in zip(rep[key[0]], out[key[0]]))
    compare(name, rep, out, key, [c for c in out.columns if c not in key], exact=True,
            note="exact; row order " + ("identical" if order_ok else "DIFFERS"))
    if not order_ok: MISMATCH.append((name, "row order", "differs", ""))
exact_sheet("S1_FECD_S-PrediXcan_all", ["gene", "TISSUE"])
exact_sheet("S2_FECD_S-MultiXcan_all", ["gene"])
exact_sheet("S3_FUSION_all", ["TISSUE", "ID"])
exact_sheet("S10_NR1H3_ZNF408_49tissues", ["SYMBOL", "tissue"])
exact_sheet("S11_S-MultiXcan_49tissues", ["ENSG"])

fu = fusion_eur()     # unrounded FUSION outputs (European GWAS) for the S9 Panel A check below

# ================= S4 FinnGen look-up =================
fg = pd.read_csv(f"{B}/analysis/nr1h3_finngen_replication_H184_vs_H185.tsv", sep="\t")
pw = pd.read_csv(f"{B}/analysis/nr1h3_finngen_lookup_power.tsv", sep="\t")
disc = {}
for t in T3:
    d = pd.read_csv(f"{SPX}/fecd_{t}.csv"); d["ENSG"] = d.gene.map(ensg)
    for g, e in (("NR1H3", "ENSG00000025434"), ("ZNF408", "ENSG00000175213")):
        x = d[(d.ENSG == e) & d.pvalue.notna()]
        if len(x): disc[(g, t)] = (x.zscore.iloc[0], x.pvalue.iloc[0])
o = fg.rename(columns={"zscore": "finngen_z", "pvalue": "finngen_p"}).copy()
o["discovery_z"] = [disc.get((g, t), (None, None))[0] for g, t in zip(o.gene, o.tissue)]
o["discovery_p"] = [disc.get((g, t), (None, None))[1] for g, t in zip(o.gene, o.tissue)]
o["concordant_direction"] = [("n/a" if dz is None or pd.isna(dz) else ("yes" if np.sign(dz) == np.sign(fz) else "no"))
                             for dz, fz in zip(o.discovery_z, o.finngen_z)]
o["P_below_0.025"] = np.where(o.finngen_p < 0.025, "yes", "no")
o = o.merge(pw[["gene", "endpoint", "tissue", "power"]], on=["gene", "endpoint", "tissue"], how="left")
o["power"] = o.power.map(lambda x: None if pd.isna(x) else round(x, 2))
o.to_csv(f"{OUTD}/S4_FinnGen_lookup.tsv", sep="\t", index=False, na_rep="")
rep = sheet("S4_FinnGen_lookup")
compare("S4_FinnGen_lookup", rep, o, ["gene", "endpoint", "tissue"],
        ["ICD10", "cases", "controls", "discovery_z", "discovery_p", "finngen_z", "finngen_p", "concordant_direction",
         "P_below_0.025", "power"])

# ================= S5 conditional 11p11.2 =================
s5 = []
for t in ("Cells_Cultured_fibroblasts", "Skin_Not_Sun_Exposed_Suprapubic"):
    base = f"{B}/results/fusion_11p11_conditional/{t}_11p11_minp1"
    inc = pd.read_csv(base + ".joint_included.dat", sep=r"\s+"); dro = pd.read_csv(base + ".joint_dropped.dat", sep=r"\s+")
    for d, st, zc, pc in ((inc, "included (joint)", "JOINT.Z", "JOINT.P"), (dro, "dropped (conditioned)", "COND.Z", "COND.P")):
        for _, r in d.iterrows():
            s5.append(dict(tissue=t, ID=r.ID, **{"TWAS.Z": r["TWAS.Z"], "TWAS.P": r["TWAS.P"], "status": st,
                                                 "joint_or_cond_Z": r[zc], "joint_or_cond_P": r[pc]}))
rep = sheet("S5_Conditional_11p11.2")
compare("S5_Conditional_11p11.2", rep, pd.DataFrame(s5), ["tissue", "ID"],
        ["TWAS.Z", "TWAS.P", "status", "joint_or_cond_Z", "joint_or_cond_P"])

# ================= S6 colocalisation =================
s6 = pd.read_csv(f"{B}/analysis/fecd_coloc_S6_fullcis.tsv", sep="\t")
rep = sheet("S6_Colocalisation_FECD")
compare("S6_Colocalisation_FECD", rep, s6, ["ENSG", "tissue", "GWAS"],
        [c for c in rep.columns if c not in ("ENSG", "tissue", "GWAS") and c in s6.columns])

# ================= S7 FECD cis-MR =================
s7 = []
for gene, (lab, acc, nca, nco, f) in [(g, x) for g in ("nr1h3", "tspan10")       # TSPAN10 rows from v2.21 (13 Nov 2026)
                                      for x in (("European subset (primary)", "GCST90475865", 3473, 445239, "EUR"),
                                                ("Multi-ancestry (sensitivity)", "GCST90480062", 3928, 566393, "multiancestry"))]:
    m = pd.read_csv(f"{B}/analysis/{gene}_eqtlgen_fecd_mr_{f}.tsv", sep="\t")
    dg = pd.read_csv(f"{B}/analysis/{gene}_eqtlgen_fecd_mr_{f}_diagnostics.tsv", sep="\t").iloc[0]
    for _, r in m.iterrows():
        lo, hi = r.OR_lo, r.OR_hi
        if r.method == "MR Egger":                # CI from t with k - 2 df, consistent with the MR-Egger P-value
            from scipy.stats import t as tdist
            c = tdist.ppf(0.975, int(r.nsnp) - 2); lo, hi = math.exp(r.b - c * r.se), math.exp(r.b + c * r.se)
        s7.append(dict(outcome_gwas=lab, gwas_accession=acc, n_cases=nca, n_controls=nco, exposure=r.exposure,
                       method=r.method, nsnp=r.nsnp,
                       b=r.b, se=r.se, pval=r.pval, OR=r.OR, OR_lo=lo, OR_hi=hi, F_set=dg.F_set,
                       clump_kb=dg.clump_kb, Q_pval=dg.Q_pval, egger_intercept=dg.egger_int,
                       egger_intercept_p=dg.egger_int_p))
rep = sheet("S7_cisMR_FECD")
compare("S7_cisMR_FECD", rep, pd.DataFrame(s7), ["outcome_gwas", "exposure", "method"],
        ["gwas_accession", "n_cases", "n_controls", "nsnp", "b", "se", "pval", "OR", "OR_lo", "OR_hi", "F_set", "clump_kb",
         "Q_pval", "egger_intercept", "egger_intercept_p"])

# ================= S8 cross-trait FUSION whole blood =================
s8 = pd.read_csv(f"{B}/analysis/suppl_crosstrait_fusion_wb_full.tsv", sep="\t")
s8 = s8.rename(columns={"cond_status": "conditional_status"})
rep = sheet("S8_crosstrait_FUSION_WB")
compare("S8_crosstrait_FUSION_WB", rep, s8, ["trait", "ENSG"],
        ["TWAS.Z", "TWAS.P", "conditional_status", "joint_or_cond_Z", "joint_or_cond_P"])

# ================= S9 Panels A-C =================
ws = WB["S9_crosstrait_analyses"]; rows = list(ws.iter_rows(values_only=True))
def panel(title_start, ncol):
    i = next(k for k, r in enumerate(rows) if isinstance(r[0], str) and r[0].startswith(title_start))
    hdr = list(rows[i + 1][:ncol]); out = []
    for r in rows[i + 2:]:
        if r[0] is None or (isinstance(r[0], str) and (r[0].startswith("Note") or r[0].startswith("Panel"))): break
        out.append(list(r[:ncol]))
    return pd.DataFrame(out, columns=hdr)
# Panel A: whole-blood FUSION (S8 outputs) for the cross-trait rows; proxy-tissue FUSION (S3 outputs) for FECD
pa = panel("Panel A", 7)
oa = []
for lab, tr in (("Glaucoma", "glaucoma"), ("Refractive error", "myopia"), ("IOP", "iop")):
    x = s8[(s8.trait == tr) & (s8.ENSG.map(ensg) == "ENSG00000025434")].iloc[0]
    oa.append({"row": lab, "TWAS.Z (marginal)": x["TWAS.Z"], "TWAS.P (marginal)": x["TWAS.P"],
               "Conditional Z": x["joint_or_cond_Z"], "Conditional P": x["joint_or_cond_P"]})
for lab, t in (("FECD (fibroblast)", "Cells_Cultured_fibroblasts"), ("FECD (skin)", "Skin_Not_Sun_Exposed_Suprapubic")):
    x = fu[(fu.TISSUE == t) & (fu.ID.map(ensg) == "ENSG00000025434")].iloc[0]
    oa.append({"row": lab, "TWAS.Z (marginal)": x["TWAS.Z"], "TWAS.P (marginal)": x["TWAS.P"]})
pa["row"] = [next((k for k in ("Glaucoma", "Refractive error", "IOP", "FECD (fibroblast)", "FECD (skin)")
                   if str(ph).startswith(k)), None) for ph in pa.Phenotype]
compare("S9 Panel A", pa, pd.DataFrame(oa), ["row"],
        ["TWAS.Z (marginal)", "TWAS.P (marginal)", "Conditional Z", "Conditional P"])
# Panel B: glaucoma per-tissue S-PrediXcan and S-MultiXcan P for NR1H3
pb = panel("Panel B", 6); ob = []
for lab, d in (("Cultured fibroblasts", "sprediXcan_glaucoma_brain_staged/mashr_Cells_Cultured_fibroblasts.csv"),
               ("Skin (not sun exposed)", "sprediXcan_glaucoma_brain_staged/mashr_Skin_Not_Sun_Exposed_Suprapubic.csv"),
               ("Brain cortex", "sprediXcan_glaucoma_brain_staged/mashr_Brain_Cortex.csv"),
               ("Whole blood", "sprediXcan_glaucoma_wb_staged/mashr_Whole_Blood.csv")):
    d = pd.read_csv(f"{B}/results/{d}"); x = d[d.gene.map(ensg) == "ENSG00000025434"].iloc[0]
    ob.append({"Tissue": lab, "S-PrediXcan z": x.zscore, "S-PrediXcan P": x.pvalue})
for lab, f in (("S-MultiXcan (fibroblast + brain + skin)", "glaucoma_smultixcan_brain.txt"),
               ("S-MultiXcan (fibroblast + skin, no brain cortex)", "glaucoma_smultixcan_nobrain.txt")):
    d = pd.read_csv(f"{B}/results/multixcan/{f}", sep="\t"); x = d[d.gene.map(ensg) == "ENSG00000025434"].iloc[0]
    ob.append({"Tissue": lab, "S-PrediXcan z": None, "S-PrediXcan P": x.pvalue})
compare("S9 Panel B", pb[pb.Tissue.isin([o["Tissue"] for o in ob])], pd.DataFrame(ob), ["Tissue"],
        ["S-PrediXcan z", "S-PrediXcan P"], subset=True)
# Panel C: cross-trait cis-MR summary (Table 3 output)
pc = panel("Panel C", 16)
t3 = pd.read_csv(f"{B}/analysis/table3_nr1h3_crosstrait_cismr.tsv", sep="\t")
t3["trait"] = t3.outcome.replace({"Myopia": "Refractive error"})
oc = []
for tr, d in t3.groupby("trait", sort=False):
    g = lambda m, c: d.loc[d.method == m, c].iloc[0]
    oc.append({"trait": tr, "k": g("Inverse variance weighted", "nsnp"), "I2GX": g("Inverse variance weighted", "Isq_GX"),
               "ivw_b": g("Inverse variance weighted", "b"), "ivw_p": g("Inverse variance weighted", "pval"),
               "wm_b": g("Weighted median", "b"), "wm_p": g("Weighted median", "pval"),
               "wmo_b": g("Weighted mode", "b"), "wmo_p": g("Weighted mode", "pval"),
               "egg_b": g("MR Egger", "b"), "egg_p": g("MR Egger", "pval"),
               "egg_int_p": g("Inverse variance weighted", "egger_intercept_pval"),
               "Q_p": g("Inverse variance weighted", "cochran_Q_pval")})
compare("S9 Panel C", pc, pd.DataFrame(oc), ["trait"],
        [c for c in ("k", "I2GX", "ivw_b", "ivw_p", "wm_b", "wm_p", "wmo_b", "wmo_p", "egg_b", "egg_p", "egg_int_p", "Q_p")
         if c in pc.columns])

# ================= S12 reference loci + gene classification =================
ws = WB["S12_FECD_loci_classification"]; rows = list(ws.iter_rows(values_only=True))
ia = next(k for k, r in enumerate(rows) if r[0] == "Panel A - Reported FECD loci")
ib = next(k for k, r in enumerate(rows) if r[0] == "Panel B - TWAS gene classification")
ha = list(rows[ia + 1][:9]); A = pd.DataFrame([list(r[:9]) for r in rows[ia + 2:ib] if r[0]], columns=ha)
hb = [h for h in rows[ib + 1] if h is not None]; Bp = pd.DataFrame([list(r[:len(hb)]) for r in rows[ib + 2:] if r[0]], columns=hb)
# S12 includes loci re-derived from the Insawang et al. summary statistics, which the authors shared and which are not
# redistributed (10 Nov 2026: nor is the lead-variant table derived from them). Without either file step 109 leaves
# those loci out, so S12 cannot be reproduced and is reported as not checked rather than as mismatches.
H185 = os.environ.get("H185_FILE", f"{B}/data/private/H185_MVP_FinnGen_UKB+AoU_SNPsinatleast2studies.gz")
H185_LEADS = os.environ.get("H185_LEADS_FILE", f"{B}/data/derived/insawang2026_clumped_lead_variants.tsv")
if os.path.exists(H185) or os.path.exists(H185_LEADS):
    la = pd.read_csv(f"{B}/analysis/fecd_reported_loci.tsv", sep="\t", dtype={"chr": str})
    compare("S12 Panel A", A, la, ["locus_id"], [c for c in ha if c != "locus_id"])
    lb = pd.read_csv(f"{B}/analysis/fecd_twas_gene_locus_class.tsv", sep="\t", dtype={"chr": str})
    disp = pd.read_csv(f"{B}/data/gene_symbol_display.tsv", sep="\t").set_index("ENSG").display_symbol.to_dict()
    lb["symbol"] = [disp.get(e, s_) for e, s_ in zip(lb.ENSG, lb.symbol)]      # current symbols, as in Table 1
    compare("S12 Panel B", Bp, lb, ["ENSG"], [c for c in hb if c != "ENSG" and c in lb.columns])
else:
    for t in ("S12 Panel A", "S12 Panel B"):
        record(t, 0, 0, [], "NOT CHECKED: needs the Insawang et al. summary statistics (restricted; set H185_FILE)")

# ================= S13 fine-mapping (tables; Panel A text items not checked) =================
ws = WB["S13_FineMapping_11p11.2"]; rows = list(ws.iter_rows(values_only=True))
def panel13(start):
    i = next(k for k, r in enumerate(rows) if isinstance(r[0], str) and r[0].startswith(start))
    hdr = [h for h in rows[i + 1] if h is not None]; out = []
    for r in rows[i + 2:]:
        if r[0] is None or (isinstance(r[0], str) and r[0].startswith("Panel")): break
        out.append(list(r[:len(hdr)]))
    return pd.DataFrame(out, columns=hdr)
sp = pd.read_csv(f"{B}/results/finemap_11p11/susie_pip.tsv", sep="\t")
sp = sp.rename(columns={"pos": "pos_GRCh38", "ea": "effect_allele", "oa": "other_allele", "z": "z_EUR", "p": "P_EUR",
                        "PIP": "SuSiE_PIP"}).assign(chr="11")
pbb = panel13("Panel B")
compare("S13 Panel B", pbb, sp, ["rsid"], [c for c in pbb.columns if c != "rsid"], subset=True)
ab = pd.read_csv(f"{B}/results/finemap_11p11/abf_credible_sets.tsv", sep="\t")
ab = ab.rename(columns={"id": "variant_id", "pos": "pos_GRCh38", "p": "P", "PIP": "ABF_PIP"})
for start, st in (("Panel C", "European"), ("Panel D", "Multi-ancestry")):
    pp = panel13(start)
    compare(f"S13 {start}", pp, ab[ab.set == st], ["variant_id"], [c for c in pp.columns if c != "variant_id"], subset=True)
ec = pd.read_csv(f"{B}/analysis/eqtlcat_susie_11p11.tsv", sep="\t")
pe = panel13("Panel E")
compare("S13 Panel E", pe, ec, ["tissue", "ENSG", "credible_set"], [c for c in pe.columns if c not in ("tissue", "ENSG", "credible_set")])
fo = pd.read_csv(f"{B}/results/focus_11p11/focus_11p11_summary.tsv", sep="\t", keep_default_na=False, na_values=["NA", ""])   # gene "NULL" = null model
fo = fo.rename(columns={"twas_z": "TWAS_z", "PIP": "model_PIP", "gene_summed_PIP": "gene_sum_of_model_PIPs"})
for start in ("Panel F", "Panel G"):
    pp = panel13(start)
    compare(f"S13 {start}", pp, fo, ["region", "gene", "tissue"], [c for c in pp.columns if c not in ("region", "gene", "tissue")], subset=True)

# ================= S14 RNA-seq samples =================
st = pd.read_csv(f"{B}/analysis/nr1h3_endothelium_samples.tsv", sep="\t") if os.path.exists(f"{B}/analysis/nr1h3_endothelium_samples.tsv") else None
if st is not None:
    rep = sheet("S14_Endothelium_RNAseq_samples", header_row=3, first=4)
    compare("S14_Endothelium_RNAseq_samples", rep[rep.included == "yes"], st, ["sample_in_analysis"], ["group", "batch_in_model"])

# ================= manuscript Tables 1-3 =================
mt = pd.read_csv(MTAB, sep="\t", dtype=str, keep_default_na=False)
t1 = mt[mt.table == "Table 1"].pivot(index="row", columns="column", values="value").reset_index()
r1 = pd.read_csv(f"{B}/analysis/table2_smultixcan_ranking.tsv", sep="\t", dtype={"rank": str})
o1 = r1.rename(columns={"rank": "Rank", "symbol": "Symbol", "smultixcan_P": "S-MultiXcan P", "cytoband": "Cytoband",
                        "FECD_locus": "Gene region", "locus_min_GWAS_P_1Mb": "P value of locus"})
compare("Table 1", t1, o1, ["Rank"], ["Symbol", "S-MultiXcan P", "Cytoband", "Gene region", "P value of locus"])
t2 = mt[mt.table == "Table 2"].pivot(index="row", columns="column", values="value").reset_index()
r2 = pd.read_csv(f"{B}/analysis/table1_nr1h3_znf408_sprediXcan.tsv", sep="\t")
r2 = r2[r2.gene == "NR1H3"].copy()
r2["Tissue"] = r2.tissue.map({"Cells_Cultured_fibroblasts": "Cultured fibroblasts",
                              "Skin_Not_Sun_Exposed_Suprapubic": "Skin not sun exposed suprapubic",
                              "Esophagus_Mucosa": "Oesophagus mucosa"})
o2 = r2.rename(columns={"n_snps": "n SNPs", "beta": "β"})
compare("Table 2", t2, o2, ["Tissue"], ["n SNPs", "Z", "β", "SE", "P"])
t3m = mt[mt.table == "Table 3"].pivot(index="row", columns="column", values="value").reset_index()
o3 = t3.rename(columns={"trait": "Outcome phenotype", "method": "MR method", "b": "Effect estimate", "se": "Standard error",
                        "pval": "P", "Isq_GX": "I2", "egger_intercept_pval": "Egger intercept P", "cochran_Q_pval": "Cochran's Q P"})
compare("Table 3", t3m, o3, ["Outcome phenotype", "MR method"],
        ["Effect estimate", "Standard error", "P", "I2", "Egger intercept P", "Cochran's Q P"])

# ================= report =================
pd.DataFrame(RESULTS, columns=["table", "rows", "cells", "mismatches", "note"]).to_csv(
    f"{B}/analysis/reported_tables_check.tsv", sep="\t", index=False)
pd.DataFrame(MISMATCH, columns=["table", "row", "column", "reported", "output"]).to_csv(
    f"{B}/analysis/reported_tables_mismatches.tsv", sep="\t", index=False)
n_bad = sum(r[3] for r in RESULTS)
print(f"\n{'ALL REPORTED VALUES MATCH THE OUTPUTS' if n_bad == 0 else f'{n_bad} MISMATCHES'} "
      f"(analysis/reported_tables_check.tsv)")
sys.exit(1 if n_bad else 0)
