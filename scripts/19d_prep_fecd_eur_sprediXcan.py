#!/usr/bin/env python3
# ============================================================
# Step 19d (30 Oct 2026): European-ancestry FECD GWAS (GCST90475865; 3,473 cases / 445,239 controls; GRCh38)
#   in S-PrediXcan format, prepared like the multi-ancestry GWAS (steps 11 and 17), then matched to the GTEx v8
#   MASHR model variants of all 49 tissues by chromosome, GRCh38 position and allele pair with step 95a (as the
#   multi-ancestry input is), for steps 19e/19f.
#   Output: data/gwas/fecd_eur_gwas_sprediXcan_unmapped.txt.gz (all variants)
#           data/gwas/fecd_eur_alltissue_sprediXcan.txt.gz (one row per matched model variant; step 95a)
#   - beta = ln(OR); SE = (ln CI upper - ln CI lower) / (2 x qnorm(0.975)) (the standard_error column is empty;
#     the European convention of steps 83a, 102c and 113); Z = beta / SE
#   - kept: autosomes, valid OR and CI, SE > 0, 0 < P <= 1, both alleles given (step 11's filters)
#   - panel_variant_id = chr{CHR}_{BP}_{other allele}_{effect allele}_b38; duplicates of it keep the first
#     record (step 17)
#   Replaces the 4 Jul 2026 input, which matched the model variants by rsID (model variants without an
#   rsID were not used); that input and its outputs are kept for reference (results/multixcan/fecd_eur_smultixcan.txt).
# Usage: python3 19d_prep_fecd_eur_sprediXcan.py
# ============================================================
import gzip, math, os, subprocess, sys

BASE = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
SRC = f"{BASE}/data/gwas/fecd_eur_GCST90475865.tsv.gz"
OUT = f"{BASE}/data/gwas/fecd_eur_gwas_sprediXcan_unmapped.txt.gz"
MAPPED = f"{BASE}/data/gwas/fecd_eur_alltissue_sprediXcan.txt.gz"
Z975 = 1.959963984540054
AUTO = {str(c) for c in range(1, 23)}


def num(x):
    try:
        v = float(x)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


n_in = n_out = 0
seen = set()
with gzip.open(SRC, "rt") as fi, gzip.open(OUT, "wt") as fo:
    col = {c: i for i, c in enumerate(fi.readline().rstrip("\n").split("\t"))}
    iC, iP, iE, iO, iOR, iF, iPv, iS, iU, iL, iN = (col[k] for k in (
        "chromosome", "base_pair_location", "effect_allele", "other_allele", "odds_ratio",
        "effect_allele_frequency", "p_value", "rsid", "ci_upper", "ci_lower", "n"))
    fo.write("panel_variant_id\trsid\tchromosome\tposition\teffect_allele\tnon_effect_allele\tzscore\tpvalue\t"
             "sample_size\teffect_size\tstandard_error\tfrequency\n")
    for line in fi:
        n_in += 1
        f = line.rstrip("\n").split("\t")
        c, bp, ea, oa = f[iC], f[iP], f[iE], f[iO]
        orr, up, lo, p = num(f[iOR]), num(f[iU]), num(f[iL]), num(f[iPv])
        if c not in AUTO or not bp.isdigit() or not ea or not oa or ea == "#NA" or oa == "#NA":
            continue
        if None in (orr, up, lo, p) or orr <= 0 or up <= 0 or lo <= 0 or not 0 < p <= 1:
            continue
        b = math.log(orr)
        se = (math.log(up) - math.log(lo)) / (2 * Z975)
        if not se > 0:
            continue
        pvid = f"chr{c}_{bp}_{oa}_{ea}_b38"
        if pvid in seen:
            continue
        seen.add(pvid)
        fo.write(f"{pvid}\t{f[iS]}\t{c}\t{bp}\t{ea}\t{oa}\t{b / se!r}\t{f[iPv]}\t{f[iN]}\t{b!r}\t{se!r}\t{f[iF]}\n")
        n_out += 1
print(f"GWAS rows read: {n_in:,} | written: {n_out:,} -> {OUT}", flush=True)
del seen
# match to the model variants by chromosome, position and allele pair (step 95a)
subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), "95a_remap_fecd_gwas_alltissue.py"),
                OUT, MAPPED], check=True)
