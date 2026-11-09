#!/usr/bin/env python3
# ============================================================
# Step 95a: Remap FECD GWAS to the varIDs of ALL 49 GTEx v8 MASHR models
#   Generalises the June 2026 fibroblast fix (fecd_fib_sprediXcan.txt.gz).
#   The original prep built panel_variant_id as chr_pos_NEA_EA_b38, which only
#   matches a model varID (chr_pos_REF_ALT_b38) when the GWAS effect allele
#   happens to be ALT -> 29-46% of model SNPs used in 47/49 tissues.
#   Here every model varID is matched to the GWAS by chr+pos+allele SET, and
#   the row is written under the model's own varID. effect/non_effect alleles
#   are kept as in the GWAS, so S-PrediXcan does the allele flip itself.
#   Usage: python 95a_remap_fecd_gwas_alltissue.py [in.txt.gz] [out.txt.gz]
# ============================================================
import gzip, glob, os, sqlite3, sys

BASE = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
IN   = sys.argv[1] if len(sys.argv) > 1 else f"{BASE}/data/gwas/fecd_gwas_sprediXcan.txt.gz"
OUT  = sys.argv[2] if len(sys.argv) > 2 else f"{BASE}/data/gwas/fecd_alltissue_sprediXcan.txt.gz"
DBS  = sorted(glob.glob(f"{BASE}/data/predictdb_all/eqtl/mashr/mashr_*.db"))
assert len(DBS) == 49, f"expected 49 MASHR dbs, found {len(DBS)}"

# (chr, pos, frozenset(alleles)) -> set of model varIDs
want = {}
for db in DBS:
    con = sqlite3.connect(db)
    for (v,) in con.execute("SELECT DISTINCT varID FROM weights"):
        c, p, a1, a2, _ = v.split("_")
        want.setdefault((c, p, frozenset((a1, a2))), set()).add(v)
    con.close()
n_model = sum(len(s) for s in want.values())
print(f"model varIDs (union of 49 tissues): {n_model:,}")

done, n_in = set(), 0
with gzip.open(IN, "rt") as fi, gzip.open(OUT, "wt") as fo:
    hdr = fi.readline(); fo.write(hdr)
    h = hdr.rstrip("\n").split("\t")
    iV, iC, iP, iE, iN = (h.index(k) for k in
        ("panel_variant_id", "chromosome", "position", "effect_allele", "non_effect_allele"))
    for line in fi:
        n_in += 1
        f = line.rstrip("\n").split("\t")
        key = (f"chr{f[iC]}" if not f[iC].startswith("chr") else f[iC], f[iP],
               frozenset((f[iE], f[iN])))
        for v in want.get(key, ()):
            if v in done:        # duplicate GWAS rows at one site: keep the first
                continue
            done.add(v); f[iV] = v
            fo.write("\t".join(f) + "\n")
print(f"GWAS rows read: {n_in:,} | model varIDs matched: {len(done):,} "
      f"({100*len(done)/n_model:.1f}%)")
print(f"written: {OUT}")
