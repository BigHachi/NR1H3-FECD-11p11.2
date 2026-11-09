#!/usr/bin/env python3
# ============================================================
# Step 60b (10 Nov 2026, review): allele-matched glaucoma GWAS for FUSION.
#   FinnGen lists a multi-allelic site once per alternative allele under one rsID, and step 60 keeps
#   the first record per rsID (data/gwas/glaucoma_gwas_fusion_firstrecord.tsv.gz). Where that record is
#   strand-ambiguous, not a single A/C/G/T, or allele-inconsistent with the 1000G EUR LD reference,
#   FUSION either drops it or uses it unaligned (allele.qc only removes ambiguous/non-ACGT pairs).
#   For every rsID on the LD reference this step keeps the first FinnGen record whose alleles match
#   the reference (step 55c's rule, which steps 76d-76h applied as a sensitivity analysis up to
#   release v13). A record that is still allele-inconsistent after the search is dropped, so FUSION
#   imputes it from LD instead of using it unaligned; a strand-ambiguous record without a
#   replacement is left for FUSION to drop. Rows not on the LD reference are written unchanged.
#   Output: data/gwas/glaucoma_gwas_fusion.tsv (the FUSION input of steps 76d and 76h)
#           analysis/glaucoma_fusion_allele_matching.tsv (one row per record replaced or dropped)
# ============================================================
import gzip, math, os
from collections import Counter, defaultdict

BASE = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
os.chdir(BASE)
LDREF = "data/ld_ref/LDREF/1000G.EUR"
IN = "data/gwas/glaucoma_gwas_fusion_firstrecord.tsv.gz"
OUT = "data/gwas/glaucoma_gwas_fusion.tsv"
LOG = "analysis/glaucoma_fusion_allele_matching.tsv"
SOURCE = "data/gwas/finngen_R12_H7_GLAUCOMA.gz"
ACGT = set("ACGT")
_COMP = str.maketrans("ACGT", "TGCA")


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


BIM = {}
for c in range(1, 23):
    for line in open(f"{LDREF}.{c}.bim"):
        ch, rs, _, pos, a1, a2 = line.split()
        BIM[rs] = (ch, int(pos), a1.upper(), a2.upper())

# ---- pass 1: rows on the LD reference whose record has to be checked
rows, need = [], set()
with gzip.open(IN, "rt") as fh:
    hdr = fh.readline()
    for line in fh:
        s, a1, a2, z, n = line.split()
        rows.append((s, a1, a2, z, n))
        b = BIM.get(s)
        if b is not None and (fusion_drops(a1, a2) or not orient(a1, a2, b[2], b[3])):
            need.add(s)

# ---- pass 2: FinnGen records whose rsids field is one of those rsIDs (the pipeline's key)
cand = defaultdict(list)
with gzip.open(SOURCE, "rt") as fh:
    col = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    iC, iR, iA, iS, iB, iSE, iPv = (col[k] for k in ("#chrom", "ref", "alt", "rsids", "beta", "sebeta", "pval"))
    autos = {str(c) for c in range(1, 23)}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if p[iS] not in need or p[iC].replace("chr", "") not in autos:
            continue
        try:
            b, se, pv = float(p[iB]), float(p[iSE]), float(p[iPv])
        except ValueError:
            continue
        if not (math.isfinite(b) and math.isfinite(se) and math.isfinite(pv)) or se <= 0:
            continue
        cand[p[iS]].append((p[iA].upper(), p[iR].upper(), b / se))

# ---- write
cnt = Counter()
with open(OUT, "w") as out, open(LOG, "w") as log:
    out.write(hdr)
    log.write("SNP\tCHR\tPOS\tLDREF_A1\tLDREF_A2\tA1_first\tA2_first\tZ_first\tfirst_record_status\t"
              "action\tA1_used\tA2_used\tZ_used\n")
    for s, a1, a2, z, n in rows:
        if s in need:
            b = BIM[s]
            st = "strand_ambiguous_or_not_snv" if fusion_drops(a1, a2) else "allele_inconsistent"
            ok = [r for r in cand.get(s, ()) if not fusion_drops(r[0], r[1]) and orient(r[0], r[1], b[2], b[3])]
            if ok:
                log.write(f"{s}\t{b[0]}\t{b[1]}\t{b[2]}\t{b[3]}\t{a1}\t{a2}\t{z}\t{st}\treplaced\t"
                          f"{ok[0][0]}\t{ok[0][1]}\t{ok[0][2]!r}\n")
                a1, a2, z = ok[0][0], ok[0][1], repr(ok[0][2])
                cnt[f"replaced_{st}"] += 1
            elif st == "allele_inconsistent":
                log.write(f"{s}\t{b[0]}\t{b[1]}\t{b[2]}\t{b[3]}\t{a1}\t{a2}\t{z}\t{st}\tdropped\t\t\t\n")
                cnt["dropped_allele_inconsistent"] += 1
                continue
            else:
                cnt["kept_for_fusion_to_drop_strand_ambiguous_or_not_snv"] += 1
        out.write(f"{s}\t{a1}\t{a2}\t{z}\t{n}\n")
print(f"input rows: {len(rows)}; records checked (on the LD reference): {len(need)}")
print(dict(cnt))
print(f"wrote {OUT} and {LOG}")
