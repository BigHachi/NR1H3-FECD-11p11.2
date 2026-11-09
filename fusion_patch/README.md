# FUSION local modification (required to reproduce any --coloc_P result)

Base: gusevlab/fusion_twas master ZIP at commit e1ba5f7 (2022-05-24).

1. `FUSION.assoc_test.coloc.patch` — apply to that release:
       patch -p0 FUSION.assoc_test.R < FUSION.assoc_test.coloc.patch
   It names the beta vectors by SNP, drops variants with a missing GWAS Z, and passes
   `snp=` to coloc.abf. Required for coloc >= 5.x; varbeta is rebuilt after the filter
   so beta/varbeta lengths stay aligned.

2. `FUSION.post_process.R` in the analysis tree is byte-identical to CURRENT upstream
   master (2026), not the 2022 release. It needs the `here` package and
   `utils/plink_utils.R` from upstream master. The install therefore mixes two upstream
   versions; both files must be taken as described to reproduce the published output.
