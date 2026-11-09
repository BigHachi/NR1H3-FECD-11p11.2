#!/usr/bin/env python3
# ============================================================
# Step 55b: Independent verification of the refractive-error merge (step 55).
#   1. Genome-wide: every 1000G EUR LDREF variant (GRCh37 position + rsID, the
#      FUSION/clumping reference) is looked up in the ORIGINAL METAL file by its
#      GRCh37 MarkerName and alleles; its Zscore, aligned to the bim A1, is compared
#      with data/gwas/myopia_gwas_full.tsv. This uses no liftOver and no GWAS
#      Catalog fields, so it checks the merge independently.
#   2. cis-MR instruments: the original record (Allele1/Allele2/Freq1/Weight/Zscore/
#      P-value), the harmonised record and the merged record side by side.
#   Outputs: analysis/myopia_gwas_merge_verification.tsv
#            analysis/myopia_mr_instruments_verified.tsv
# ============================================================
import glob, gzip, os

BASE = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
os.chdir(BASE)
LDREF = "data/ld_ref/LDREF/1000G.EUR"
ORIG = "data/gwas/myopia_gwas_original.txt.gz"
HARM = "data/gwas/myopia_gwas_harmonised.tsv.gz"
MERGED = "data/gwas/myopia_gwas_full.tsv"
INSTRUMENTS = ["rs7395581", "rs11038863", "rs3740694", "rs4647760"]   # NR1H3 eQTLGen cis-MR (k = 4)

bim, bim_rs = {}, {}
for f in sorted(glob.glob(f"{LDREF}.*.bim")):
    for line in open(f):
        c, rs, _, pos, a1, a2 = line.split()
        bim[f"chr:{c}:{pos}"] = (rs, a1.upper(), a2.upper())
        bim_rs[rs] = f"chr:{c}:{pos}"

truth, n_dup, orig_inst = {}, set(), {}
want = {bim_rs[r]: r for r in INSTRUMENTS if r in bim_rs}
with gzip.open(ORIG, "rt") as fh:
    col = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        m = p[col["MarkerName"]]
        if m in want:
            orig_inst[want[m]] = p
        b = bim.get(m)
        if b is None:
            continue
        rs, a1, a2 = b
        o1, o2, z = p[col["Allele1"]].upper(), p[col["Allele2"]].upper(), float(p[col["Zscore"]])
        if (o1, o2) == (a1, a2):
            s = 1
        elif (o2, o1) == (a1, a2):
            s = -1
        else:
            continue
        if rs in truth:
            n_dup.add(rs)
        truth[rs] = (a1, a2, s * z)
for rs in n_dup:
    truth.pop(rs)

stats = dict(checkable=len(truth), compared=0, agree=0, sign_flipped=0, other_mismatch=0)
merged_inst = {}
with open(MERGED) as fh:
    col = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if p[0] in INSTRUMENTS:
            merged_inst[p[0]] = p
        t = truth.get(p[0])
        if t is None:
            continue
        a1, a2, z = p[col["A1"]], p[col["A2"]], float(p[col["BETA"]])
        stats["compared"] += 1
        if (a1, a2) == (t[0], t[1]):
            pass
        elif (a2, a1) == (t[0], t[1]):
            z = -z
        else:
            stats["other_mismatch"] += 1
            continue
        if abs(z - t[2]) < 5e-4:
            stats["agree"] += 1
        elif abs(z + t[2]) < 5e-4:
            stats["sign_flipped"] += 1
        else:
            stats["other_mismatch"] += 1
stats["coverage_of_checkable"] = round(stats["compared"] / stats["checkable"], 4)
with open("analysis/myopia_gwas_merge_verification.tsv", "w") as fo:
    fo.write("metric\tvalue\n")
    for k, v in stats.items():
        fo.write(f"{k}\t{v}\n")
print("1000G EUR LDREF independent check:", stats)

harm_inst = {}
with gzip.open(HARM, "rt") as fh:
    col = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if p[col["hm_rsid"]] in INSTRUMENTS:
            harm_inst[p[col["hm_rsid"]]] = {k: p[col[k]] for k in
                ("hm_pos", "hm_effect_allele", "hm_other_allele", "hm_code", "effect_allele", "eaf", "p_value")}

hdr = ["SNP", "MARKER_GRCh37", "orig_Allele1", "orig_Allele2", "orig_Freq1", "orig_Weight", "orig_Zscore", "orig_P",
       "harm_pos_GRCh38", "harm_hm_effect", "harm_hm_other", "harm_hm_code",
       "merged_BP", "merged_A1", "merged_A2", "merged_Z", "merged_FREQ", "merged_N", "merged_P", "verified"]
with open("analysis/myopia_mr_instruments_verified.tsv", "w") as fo:
    fo.write("\t".join(hdr) + "\n")
    for rs in INSTRUMENTS:
        o, h, m = orig_inst.get(rs), harm_inst.get(rs), merged_inst.get(rs)
        if not (o and h and m):
            fo.write(f"{rs}\t" + "\tNA" * (len(hdr) - 2) + "\tMISSING\n")
            continue
        a1o, a2o, fr, wt, z, pv = o[1].upper(), o[2].upper(), float(o[3]), o[5], float(o[6]), o[7]
        s = 1 if (m[3], m[4]) == (a1o, a2o) else -1 if (m[4], m[3]) == (a1o, a2o) else 0
        ok = (s != 0 and abs(float(m[5]) - s * z) < 5e-4 and abs(float(m[8]) - (fr if s == 1 else 1 - fr)) < 1.5e-4
              and float(m[9]) == float(wt) and float(m[7]) == float(pv) and m[2] == h["hm_pos"]
              and (m[3], m[4]) == (h["hm_effect_allele"], h["hm_other_allele"]))
        fo.write("\t".join(map(str, [rs, o[0], a1o, a2o, o[3], wt, o[6], pv, h["hm_pos"], h["hm_effect_allele"],
                 h["hm_other_allele"], h["hm_code"], m[2], m[3], m[4], m[5], m[8], m[9], m[7],
                 "OK" if ok else "MISMATCH"])) + "\n")
print(open("analysis/myopia_mr_instruments_verified.tsv").read())
