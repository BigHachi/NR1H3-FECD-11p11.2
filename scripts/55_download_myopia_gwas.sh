#!/bin/bash
# ============================================================
# Myopia Analysis - Step 55: Merge Hysi 2020 refractive-error GWAS
# GCST010003 | PMID 32231278 | public subset (excl. 23andMe), European
# Original METAL file (GRCh37 MarkerName, Allele1/Allele2, Freq1, Weight, Zscore)
# linked to the GWAS Catalog harmonised file (GRCh38, rsID, reference-aligned
# alleles) by VARIANT IDENTITY: chromosome + GRCh38 position (UCSC liftOver of
# the GRCh37 MarkerName) + submitted alleles. p-value and allele frequency are
# checked for agreement only. Z is aligned explicitly to hm_effect_allele;
# missing, ambiguous and allele-inconsistent records are rejected and counted.
#   Output : data/gwas/myopia_gwas_full.tsv
#            SNP CHR BP(GRCh38) A1 A2 BETA(=Z for A1) SE(=1) P FREQ(of A1)
#            N(METAL Weight) MARKER_GRCh37 PALINDROMIC LIFT_STRAND
#   Audit  : analysis/myopia_gwas_merge_audit.tsv
# Replaces (Oct 19 2026) the earlier merge on (round(p,8), round(eaf,4)), which
# mis-assigned Z-scores between variants sharing that key and dropped ~8% of
# variants (old script kept as 55_download_myopia_gwas.sh.bak_20261019).
# Needs pyliftover + numpy in $METAXCAN_PY (as step 50).
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -euo pipefail
cd "$FECD_BASE"

GWAS_ORIG=data/gwas/myopia_gwas_original.txt.gz
GWAS_HARM=data/gwas/myopia_gwas_harmonised.tsv.gz
GWAS_FULL=data/gwas/myopia_gwas_full.tsv
CHAIN=data/liftover/hg19ToHg38.over.chain.gz
AUDIT=analysis/myopia_gwas_merge_audit.tsv
mkdir -p data/liftover analysis

if [ ! -f "$GWAS_ORIG" ]; then
  echo "Downloading original file (Z-scores)..."
  curl -L --fail --retry 3 -o "$GWAS_ORIG" \
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010003/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020.txt.gz"
fi

if [ ! -f "$GWAS_HARM" ]; then
  echo "Downloading harmonised file (GRCh38 + rsids)..."
  curl -L --fail --retry 3 -o "$GWAS_HARM" \
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010003/harmonised/32231278-GCST010003-HP_0000539.h.tsv.gz"
fi

if [ ! -f "$CHAIN" ]; then
  echo "Downloading hg19->hg38 chain file..."
  curl -L --fail --retry 3 -o "$CHAIN" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
fi

echo "=== Merging refractive-error GWAS by variant identity ==="
WORK=$(mktemp -d "${TMPDIR:-/tmp}/myopia_idmap.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

cat > "$WORK/myopia_idmap.py" << 'PYEOF'
"""Identity-based merge of the Hysi 2020 refractive-error GWAS (GCST010003).

The original METAL file carries the Z-scores (MarkerName = chr:CHR:POS on GRCh37,
Allele1/Allele2/Freq1/Weight/Zscore/P-value; Zscore and Freq1 refer to Allele1).
The GWAS Catalog harmonised file carries rsIDs and GRCh38 reference-aligned alleles,
plus the as-submitted alleles (effect_allele = Allele1), eaf and p_value, but no Z.

Records are linked by variant identity:
    (chromosome, GRCh38 position from UCSC liftOver of the GRCh37 MarkerName,
     submitted Allele1, submitted Allele2)
  = (chromosome, base_pair_location, effect_allele, other_allele) in the harmonised file.
p-value and allele frequency are checked for agreement but are not part of the key.
The Z-score is then aligned explicitly to hm_effect_allele (complementing the
submitted alleles when the GRCh37->GRCh38 liftOver maps to the minus strand).
Duplicated keys, missing records, allele-inconsistent records and duplicated
rsIDs are rejected and counted in the audit table.

usage: python -I myopia_idmap.py ORIG HARM CHAIN OUT AUDIT TMPDIR
"""
import gzip, math, os, sys
from array import array
from collections import Counter, defaultdict
import numpy as np
from pyliftover import LiftOver

ORIG, HARM, CHAIN, OUT, AUDIT, TMP = sys.argv[1:7]
os.makedirs(TMP, exist_ok=True)
AUTOSOMES = [str(c) for c in range(1, 23)]
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(a):
    return a.translate(_COMP)[::-1]


def palindromic(a, b):
    return len(a) == 1 and len(b) == 1 and revcomp(a) == b


audit = Counter()

# ---- 1. bucket the harmonised file by chromosome (as-submitted key + hm fields) ----
hfh = {c: open(os.path.join(TMP, f"harm_{c}.tsv"), "w") for c in AUTOSOMES}
with gzip.open(HARM, "rt") as fh:
    col = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        audit["harm_rows"] += 1
        c = p[col["chromosome"]]
        if c not in hfh:
            audit["harm_non_autosomal"] += 1
            continue
        hfh[c].write("\t".join([
            p[col["base_pair_location"]], p[col["effect_allele"]].upper(), p[col["other_allele"]].upper(),
            p[col["hm_chrom"]], p[col["hm_pos"]], p[col["hm_effect_allele"]].upper(),
            p[col["hm_other_allele"]].upper(), p[col["hm_code"]], p[col["hm_rsid"]],
            p[col["eaf"]], p[col["p_value"]]]) + "\n")
for f in hfh.values():
    f.close()

# ---- 2. lift the original GRCh37 MarkerName to GRCh38 and bucket by chromosome ----
lo = LiftOver(CHAIN)
ofh = {c: open(os.path.join(TMP, f"orig_{c}.tsv"), "w") for c in AUTOSOMES}
with gzip.open(ORIG, "rt") as fh:
    col = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        audit["orig_rows"] += 1
        marker = p[col["MarkerName"]]
        try:
            _, c37, pos37 = marker.split(":")
            pos37 = int(pos37)
            z = float(p[col["Zscore"]]); pv = float(p[col["P-value"]]); fr = float(p[col["Freq1"]])
            float(p[col["Weight"]])
        except ValueError:
            audit["orig_unparseable"] += 1
            continue
        if c37 not in ofh:
            audit["orig_non_autosomal"] += 1
            continue
        if not (math.isfinite(z) and 0 < pv <= 1 and 0 <= fr <= 1):
            audit["orig_invalid_values"] += 1
            continue
        hits = lo.convert_coordinate(f"chr{c37}", pos37 - 1)  # pyliftover is 0-based
        if not hits:
            audit["liftover_failed"] += 1
            continue
        if len(hits) > 1:
            audit["liftover_multiple"] += 1
            continue
        c38, pos38, strand = hits[0][0].replace("chr", ""), hits[0][1] + 1, hits[0][2]
        if c38 != c37:
            audit["liftover_other_chrom"] += 1
            continue
        ofh[c37].write("\t".join([str(pos38), p[col["Allele1"]].upper(), p[col["Allele2"]].upper(),
                                  marker, strand, p[col["Zscore"]], p[col["P-value"]],
                                  p[col["Freq1"]], p[col["Weight"]]]) + "\n")
for f in ofh.values():
    f.close()

# ---- 3. per-chromosome identity join + explicit allele alignment ----
codecross = Counter()     # (hm_code, liftover strand, palindromic, sign) for the audit
tmp_out = os.path.join(TMP, "matched.tsv")
rs_nums = array("q")
with open(tmp_out, "w") as fout:
    for c in AUTOSOMES:
        harm = defaultdict(list)
        harm_pos = set()
        with open(os.path.join(TMP, f"harm_{c}.tsv")) as fh:
            for line in fh:
                h = line.rstrip("\n").split("\t")
                harm[(h[0], h[1], h[2])].append(h)
                harm_pos.add(h[0])
        orig = defaultdict(list)
        with open(os.path.join(TMP, f"orig_{c}.tsv")) as fh:
            for line in fh:
                o = line.rstrip("\n").split("\t")
                orig[(o[0], o[1], o[2])].append(o)
        for key, olist in orig.items():
            if len(olist) > 1:
                audit["ambiguous_source_duplicate_key"] += len(olist)
                continue
            o = olist[0]
            hlist = harm.get(key)
            if hlist is None:
                if (key[0], key[2], key[1]) in harm:
                    audit["allele_inconsistent_submitted_order"] += 1
                elif key[0] in harm_pos:
                    audit["allele_inconsistent_other_alleles_at_pos"] += 1
                else:
                    audit["missing_no_harmonised_record"] += 1
                continue
            if len(hlist) > 1:
                audit["ambiguous_harmonised_duplicate_key"] += 1
                continue
            h = hlist[0]
            pos38, a1, a2, marker, strand, z, pv, fr, wt = o
            _, _, _, hm_chr, hm_pos, hm_ea, hm_oa, hm_code, rsid, eaf_h, p_h = h
            # identity carried through harmonisation unchanged?
            if hm_chr != c or hm_pos != pos38:
                audit["inconsistent_hm_position"] += 1
                continue
            # p-value / frequency agreement (consistency check, not part of the key)
            pv_f, fr_f = float(pv), float(fr)
            try:
                ph_f, eh_f = float(p_h), float(eaf_h)
            except ValueError:
                audit["inconsistent_p_or_eaf"] += 1
                continue
            if abs(math.log10(pv_f) - math.log10(ph_f)) > 1e-3 or abs(fr_f - eh_f) > 1e-4:
                audit["inconsistent_p_or_eaf"] += 1
                continue
            if not (rsid.startswith("rs") and rsid[2:].isdigit()):
                audit["no_rsid"] += 1
                continue
            # explicit alignment to hm_effect_allele
            b1, b2 = (revcomp(a1), revcomp(a2)) if strand == "-" else (a1, a2)
            if (b1, b2) == (hm_ea, hm_oa):
                sign = 1
            elif (b1, b2) == (hm_oa, hm_ea):
                sign = -1
            else:
                audit["allele_inconsistent_with_harmonised"] += 1
                continue
            pal = palindromic(a1, a2)
            codecross[(hm_code, strand, pal, sign)] += 1
            zz = sign * float(z)
            f1 = fr_f if sign == 1 else 1 - fr_f
            fout.write(f"{rsid}\t{c}\t{pos38}\t{hm_ea}\t{hm_oa}\t{zz:.6f}\t1\t{pv}\t{f1:.4f}\t{wt}\t{marker}\t{int(pal)}\t{strand}\n")
            rs_nums.append(int(rsid[2:]))
        del harm, orig, harm_pos

# ---- 4. reject rsIDs that map to more than one accepted record ----
u, n = np.unique(np.frombuffer(rs_nums, dtype=np.int64), return_counts=True)
dup_rs = {f"rs{r}" for r in u[n > 1]}
audit["ambiguous_duplicate_rsid"] = int(n[n > 1].sum())
del rs_nums, u, n
with open(tmp_out) as fin, open(OUT, "w") as fout:
    fout.write("SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\tFREQ\tN\tMARKER_GRCh37\tPALINDROMIC\tLIFT_STRAND\n")
    for line in fin:
        if line.split("\t", 1)[0] in dup_rs:
            continue
        fout.write(line)
        audit["written"] += 1
        if line.split("\t")[11] == "1":
            audit["written_palindromic"] += 1

# the GWAS Catalog harmoniser's own orientation call (hm_code) vs the explicit alignment
FLIPPED = {"2", "4", "6", "8", "11", "13"}
agree = sum(n for (code, s, pal, sign), n in codecross.items() if (code in FLIPPED) == (sign == -1))
audit["hm_code_orientation_agrees"] = agree
audit["hm_code_orientation_disagrees"] = sum(codecross.values()) - agree

with open(AUDIT, "w") as fa:
    fa.write("category\tn\n")
    for k in ["orig_rows", "orig_unparseable", "orig_non_autosomal", "orig_invalid_values",
              "liftover_failed", "liftover_multiple", "liftover_other_chrom",
              "harm_rows", "harm_non_autosomal",
              "missing_no_harmonised_record", "allele_inconsistent_submitted_order",
              "allele_inconsistent_other_alleles_at_pos", "allele_inconsistent_with_harmonised",
              "ambiguous_source_duplicate_key", "ambiguous_harmonised_duplicate_key",
              "inconsistent_hm_position", "inconsistent_p_or_eaf", "no_rsid",
              "ambiguous_duplicate_rsid", "written", "written_palindromic",
              "hm_code_orientation_agrees", "hm_code_orientation_disagrees"]:
        fa.write(f"{k}\t{audit[k]}\n")
    fa.write("# hm_code x liftover strand x palindromic x applied sign (accepted records before rsID filter)\n")
    for (code, s, pal, sign), n in sorted(codecross.items(), key=lambda kv: -kv[1]):
        fa.write(f"#  hm_code={code}\tstrand={s}\tpalindromic={int(pal)}\tsign={sign:+d}\t{n}\n")

for c in AUTOSOMES:
    os.remove(os.path.join(TMP, f"harm_{c}.tsv")); os.remove(os.path.join(TMP, f"orig_{c}.tsv"))
os.remove(tmp_out)
print(open(AUDIT).read())
PYEOF

"$METAXCAN_PY" -I "$WORK/myopia_idmap.py" "$GWAS_ORIG" "$GWAS_HARM" "$CHAIN" \
  "$GWAS_FULL.tmp" "$AUDIT" "$WORK/buckets"
mv "$GWAS_FULL.tmp" "$GWAS_FULL"

echo ""
echo "=== Myopia merge complete ==="
echo "Output: $GWAS_FULL  ($(wc -l < $GWAS_FULL) lines); audit: $AUDIT"
