# Code availability: NR1H3 (LXRα) at 11p11.2 in Fuchs endothelial corneal dystrophy

Analysis code for the manuscript *"Integrative post-GWAS analyses nominate 11p11.2 as a candidate locus for Fuchs
endothelial corneal dystrophy"* (release **v15.2**, `NR1H3_FECD_code_bundle_20261113_v15.2`, 13 November 2026; the
version used for the article). It is provided so that editors, reviewers and readers can inspect and reproduce every
reported result. The code is maintained in a version-controlled GitHub repository, where later corrections or
compatibility updates are recorded as new releases; release v15.2 is the one archived on Zenodo for the article (links
in the manuscript's Code Availability statement).

| Path | Contents |
|---|---|
| `scripts/` | every script that produces a reported result (R, shell, Python), in run order below |
| `scripts/run_full_pipeline.sh` | master script; runs the whole pipeline in order |
| `scripts/config.sh`, `scripts/config.local.sh.example` | all paths and tool locations in one place |
| `scripts/fusion_genomewide.sh` | the genome-wide FUSION TWAS + conditional/joint procedure used for every trait (FECD step 82; glaucoma 76h; IOP and refractive error 76i) |
| `scripts/fig_common.R` | shared figure style and export (italic gene symbols; vector PDF, SVG and 600-dpi TIFF) used by every plotting step |
| `scripts/00_download_inputs.sh` | downloads the public inputs and the pinned tools; checksums are verified where recorded and supplied to the download workflow |
| `scripts/00_install_R_packages.R`, `scripts/environment_metaxcan.yml` | R packages and the MetaXcan Python environment, at the versions used |
| `data/gene_symbol_display.tsv` | the gene-name convention used in Figures 1–3, Table 1 and the supplementary tables |
| `data/derived/` | the S8 gene symbols, the gene positions used by step 109 (`ensembl_gene_tss.tsv`) and a derived table standing in for one non-public input (see below) |
| `fusion_patch/` | the local patch to FUSION required to reproduce the FUSION results (see below) |
| `analysis_outputs/` | the small result tables behind the main-text numbers, and the input-audit tables (steps 55, 55b, 55c), for checking without re-running |
| `additional_outputs/` | the outputs requested in peer review on 29 October 2026, with the tables derived for the request (step 139; see *Additional outputs requested in review*) |
| `reference/` | the reported tables (Supplementary Tables S1–S14 workbook v2.22; manuscript Tables 1–3 as printed), which step 140 compares with a run's outputs |
| `pipeline_reference/full_pipeline_all_scripts.txt` | all bundled scripts concatenated, for reading in one place |

## Licence

The code in this repository is released under the MIT License (`LICENSE`). The licence covers the scripts and
documentation written for this study. It does not cover the third-party software the pipeline downloads (MetaXcan,
FUSION, PLINK, FOCUS and the R and Python packages, each under its own licence) or the input datasets, which are
downloaded from their providers and are not redistributed here.

## Quick start

1. Unzip the bundle; the unzipped folder is the project folder (`FECD_BASE`).
2. Copy `scripts/config.local.sh.example` to `scripts/config.local.sh` and set `FECD_BASE`
   (default `~/FECD_GWAS`) and, if needed, `METAXCAN_PY`, `FOCUS_PY` and `TABIX`.
3. Install the software: `Rscript scripts/00_install_R_packages.R` and
   `conda env create -f scripts/environment_metaxcan.yml`. FOCUS (step 113) gets its own Python 3.11 virtual
   environment, which the master script creates if it is missing. `FOCUS_PY` (default `python3.11`) must point to a
   Python 3.11 interpreter: FOCUS 0.802 needs pandas 1.5.3, which has no wheels for Python 3.12 or later, and the
   installer stops with a message if another version is given.
4. Run `bash scripts/run_full_pipeline.sh`. Step 00 downloads every input; each step
   writes `logs/<step>.log`, and the run stops at the first failed step.

Requirements: about 12 GB of downloads and up to 15 GB of
intermediate files; 16 GB RAM recommended (the full run completes on 8 GB, with heavy swapping in
step 11). Run time: about 9 hours on a laptop, mostly S-PrediXcan over 49 tissues
(step 95) and FUSION (step 82).

---

## Data sources (inputs, not redistributed)

| Dataset | Accession / source | Used for |
|---|---|---|
| FECD GWAS, multi-ancestry | **GCST90480062**, Verma et al., *Science* 2024 (PMID 39024449); MVP PheCode 364.51; 3,928 cases / 566,393 controls | S-PrediXcan, S-MultiXcan, colocalisation |
| FECD GWAS, European | **GCST90475865**, Verma et al. 2024; 3,473 cases / 445,239 controls (SE derived from the 95% CI) | FUSION TWAS, cis-MR (primary) |
| FinnGen R12 | `H7_CORNEALDYSTROPHY` (ICD-10 H18.4), `H7_KERATOCONUS` (H18.5), `H7_GLAUCOMA` (26,591 cases / 473,757 controls; R12 manifest); [finngen_R12_H7_GLAUCOMA.gz](https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_H7_GLAUCOMA.gz) | FECD look-up; glaucoma cross-trait |
| Refractive error | Hysi et al. 2020, **GCST010003**, public subset without 23andMe (UKB + GERA + CREAM; Z-score meta-analysis); n = 351,091 (per-variant METAL weight up to 260,974); [GWAS Catalog files](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010003/) | cross-trait |
| IOP | IGGC meta-analysis released with Bonnemaijer et al. 2019, **GCST009413** (n = 31,269 Europeans); [GWAS Catalog files](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009413/) | cross-trait |
| GTEx v8 | PredictDB MASHR models (Zenodo 3518299); FUSION EUR weights; GTEx Portal API v2 (significant eQTLs) | TWAS, FOCUS, tissue MR, eGene status |
| eQTLGen phase I | cis-eQTLs FDR < 0.05 (N ≤ 31,684); per-SNP allele frequencies (`2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz`) | cis-MR exposure (effects standardised with the eQTLGen allele frequency, step 83a) |
| eQTL Catalogue release 6 | GTEx v8 full cis summary statistics (all variant–gene pairs) for every S6 gene in fibroblasts, skin and oesophagus mucosa, and for NR1H3 including whole blood; SuSiE eQTL credible sets | all colocalisation (steps 102c, 79b, 108); eQTL fine-mapping (112b) |
| 1000 Genomes phase 3, GRCh38 | biallelic SNV/indel release 20190312, chr11:46.2–48.3 Mb, 503 Europeans, MAF ≥ 1% (step 112a, remote tabix) | LD for SuSiE-RSS and FOCUS (steps 112, 113) |
| Corneal endothelium RNA-seq | GEO GSE112201, GSE101872 | NR1H3 expression, FECD vs control |
| GWAS Catalog | REST API v2, nine FECD analyses (step 109a) | reported FECD loci (S12) |
| LD reference | 1000 Genomes Phase 3 EUR, HapMap3 variants (FUSION LDREF) | FUSION, clumping |

Two inputs are not publicly downloadable:

- **Insawang et al. 2026 meta-analysis summary statistics** (`H185_MVP_FinnGen_UKB+AoU_SNPsinatleast2studies.gz`), shared by the authors and not redistributed, nor are tables derived from them. Point `H185_FILE` to the file if you have it. Without it, step 109 leaves out the loci re-derived from these statistics (the Insawang et al. loci in the GWAS Catalog are still used), step 109b (the 11p11.2 look-up in the Discussion) is skipped, and step 140 reports Supplementary Table S12 as not checked; no other table or figure depends on the file.
- **Genome-wide FUSION whole-blood conditional TWAS for FinnGen glaucoma, run 12 June 2025** (`Eye Conditional results 12_6_25 FinnGEN GTEV8 whole blood.xlsx`). Step 76f, for the genome-wide rank of NR1H3; fallback `glaucoma_wb_conditional_joint_included.tsv` (the jointly included genes). Step 76d reruns chr11 and reproduces the NR1H3 values (Z −7.7, joint P 1.2 × 10⁻¹⁴).

## Software

Versions are pinned in `scripts/00_install_R_packages.R` and `scripts/environment_metaxcan.yml`;
step 00 fetches MetaXcan, FUSION and PLINK at the versions below.

- **R 4.5.0**, Bioconductor 3.21: TwoSampleMR 0.7.0 (MRCIEU @ b34cc5e), ieugwasr 1.1.0.9000 (@ 8d34780), coloc 5.2.3, susieR 0.14.2, data.table 1.18.2.1, ggplot2 4.0.2, DESeq2 1.48.2, ggrepel 0.9.8, patchwork 1.3.2, ggtext 0.2.0, svglite 2.2.2, ragg 1.5.2, ggrastr 1.0.2
- **Google Chrome or Chromium** (optional; step 139b prints each figure's SVG to a vector PDF)
- **MetaXcan** (S-PrediXcan / S-MultiXcan), `hakyimlab/MetaXcan` @ `v0.8.1-1-g9b9f76b`, unmodified
- **FUSION** TWAS, `gusevlab/fusion_twas` @ `e1ba5f7` (2022-05-24), **locally patched** (see `fusion_patch/README.md`):
  `FUSION.assoc_test.R` needs the coloc patch for coloc ≥ 5, and `FUSION.post_process.R` and `utils/plink_utils.R` are from commit `9346c12`.
- **PLINK v1.90b7.2**
- Python 3.11 for MetaXcan (numpy 1.26.4, scipy 1.15.3, pandas 2.2.3, statsmodels 0.14.6, pyliftover 0.4.1; htslib 1.24 for tabix); step 55c also runs in this interpreter (`METAXCAN_PY`, for numpy and scipy); other Python steps use only the standard library
- **FOCUS** (pyfocus 0.802), in its own virtual environment `software/focus_venv` made by `00_install_focus.sh` (Python 3.11).
  FOCUS predates pandas 2, numpy 1.24 and SQLAlchemy 2, so the script pins pandas 1.5.3, numpy 1.26.4, scipy 1.11.4,
  xarray 2023.1.0 and SQLAlchemy 1.4.54, and patches pyfocus twice: `np.warnings` → `warnings` (`data/exprref.py`,
  `data/ldref.py`) and keyword arguments to `DataFrame.pivot` (`finemap.py`). No change affects the statistics.

---

## Provenance of reported outputs

Each reported output has one of three routes:

1. **Reanalysis of the original inputs by this release.** Every table and figure is regenerated from the public
   inputs by the scripts in this release, and step 140 checks every value in manuscript Tables 1–3 and
   Supplementary Tables S1–S14 against the outputs of the same run (`analysis/reported_tables_check.tsv`). Until
   21 October 2026, S1, S2, S11 and Figure 2 were retained outputs of an earlier version of step 17 that differed
   for two genes (P2RX5, TXN2), and S3 came from a 28 June build of the European-GWAS FUSION input that stored Z to
   6 significant figures; all have since been regenerated with this release (see *Changes for Final Draft v8*).
2. **Restricted or retained inputs.** Step 76f reads `data/derived/glaucoma_wb_conditional_joint_included.tsv` when the
   12 June 2025 genome-wide glaucoma conditional TWAS file is absent (step 76h reruns that analysis from public data).
   Supplementary Table S12 and the 11p11.2 look-up in the Discussion need the Insawang et al. summary statistics,
   which are not redistributed (steps 109, 109b; see *Data sources*).
3. **Historical genome-wide glaucoma ranking.** NR1H3's rank as the third-strongest conditionally independent
   whole-blood FUSION signal for glaucoma (after LOXL1-AS1 and AFAP1) comes from the 12 June 2025 genome-wide
   conditional analysis (`Eye Conditional results 12_6_25 FinnGEN GTEV8 whole blood.xlsx`, read by step 76f). Step
   76h (since release v11) reruns that analysis on all 22 chromosomes, with the allele-matched input reported from
   release v14 and the first-record input reported up to v13; both place NR1H3 third (see *Changes for release v14*).
4. **One genome-wide FUSION procedure for every trait (release v14.6).** FECD (step 82), glaucoma (76h), IOP and
   refractive error (76i) all run `scripts/fusion_genomewide.sh`: per-chromosome inputs restricted to the 1000G EUR
   LD-reference SNPs, FUSION.assoc_test, models in the MHC kept aside (FUSION writes them to a separate file), a
   genome-wide Bonferroni threshold (0.05 / models with a TWAS P), and FUSION.post_process (100-kb window) on the
   models below it. The chromosome-11 results in S8 and Figure 5 are those of the genome-wide runs (step 76d). Before
   v14.6, S8 came from a chromosome-11-only run in which every chr11 model entered the joint model (joint threshold
   0.05/463); with the genome-wide threshold NR1H3 is conditionally independent for refractive error (JOINT Z 5.2,
   P 2.4 × 10⁻⁷; it was conditioned away, P 4.0 × 10⁻⁴), although the joint model also keeps the correlated
   neighbours SLC39A13 and C1QTNF4. All three have inflated joint effects (JOINT.BETA: NR1H3 16.3, SE 3.2, marginal
   Z 6.1; C1QTNF4 −16.1, SE 3.0, marginal Z −6.4), and that of SLC39A13 also reverses in sign (JOINT.BETA −20.9,
   SE 5.0; marginal Z 6.3), so the independence of NR1H3 from them is not robust; for IOP, NR1H3 (TWAS P 1.3 × 10⁻⁴) is
   below the genome-wide threshold (6.3 × 10⁻⁶) and does not enter the conditional analysis. Glaucoma is unchanged
   (conditionally independent, third genome-wide).

## Run order and manuscript mapping (release v14)

### 0 · Inputs
| Script | Purpose |
|---|---|
| `11_prep_fecd_gwas.R` | harmonise the multi-ancestry FECD GWAS |
| `00_download_inputs.sh` | public inputs and pinned tools; checksums verified where recorded and supplied |
| `00_install_focus.sh` | FOCUS virtual environment (run only if `software/focus_venv` is missing) |
| `15_prep_fecd_fusion.sh` | multi-ancestry FECD GWAS → FUSION format |
| `17_prep_fecd_sprediXcan.sh` | FECD GWAS → S-PrediXcan format |
| `95a_remap_fecd_gwas_alltissue.py` | map each GWAS variant to the varIDs of all 49 MASHR models (fixes allele-orientation loss in fibroblasts) |
| `39_download_eur_weights.sh` | GTEx v8 EUR FUSION weights |
| `45_…`, `50_…`, `55_download_myopia_gwas.sh` | glaucoma, IOP and refractive-error GWAS |
| `55b_verify_myopia_gwas.py` | refractive-error merge: independent 1000G check, MR instrument records |
| (`55c_audit_gwas_inputs.py`) | input audit of every GWAS; runs at the end, in section 7b, because it reads every analysis input |
| `60_prep_fusion_gwas.sh` | cross-trait GWAS → FUSION format |
| `83a_extract_eqtlgen_nr1h3.sh` | eQTLGen NR1H3 cis-eQTLs |
| `100_fetch_eqtl_catalogue_region.sh` | full cis GTEx summary statistics for NR1H3, incl. whole blood (steps 79b, 108) |
| `109a_fetch_gwascat_fecd.py` | GWAS Catalog FECD associations (step 109) |

### 1 · FECD TWAS → Figures 1–2, Tables 1–2, S1–S3, S5, S10–S12
| Script | Output |
|---|---|
| `95_alltissue_sprediXcan_fecd.sh` (run with `GWAS_FILE=data/gwas/fecd_alltissue_sprediXcan.txt.gz OUT_DIR=results/sprediXcan_alltissue_remap TAG=fecd`) | S-PrediXcan, 49 tissues → **S1**, **S10** |
| `19b_rerun_fecd_smultixcan_remapped.sh` | 3-tissue S-MultiXcan (17,285 genes) → **S2** |
| `19d_prep_fecd_eur_sprediXcan.py` | European FECD GWAS (GCST90475865) → S-PrediXcan format as for the multi-ancestry GWAS (steps 11, 17): beta = ln(OR), SE from the 95% CI; then step 95a matches it to the MASHR model variants by chromosome, GRCh38 position and allele pair (`data/gwas/fecd_eur_alltissue_sprediXcan.txt.gz`) |
| `19e_fecd_eur_smultixcan.sh` | European 3-tissue S-PrediXcan + S-MultiXcan (`results/multixcan/fecd_eur_smultixcan_remapped.txt`; 17,303 genes with a P; 93% of model SNPs used per tissue) |
| `19f_fecd_eur_smultixcan_rank.py` | European ranking (`analysis/fecd_eur_smultixcan_ranking.tsv`); checks that every gene ranked down to NR1H3 (23rd) is in the S6 gene list and writes each S6 gene's European rank into `data/tophits21/genes_coloc.tsv` (→ 102c → 81c → **S6** `S_MultiXcan_rank`; `analysis/fecd_eur_smultixcan_S6_genes.tsv`). Table 1 and S2 use the multi-ancestry run (19b) |
| `95b_alltissue_smultixcan_fecd.sh` | 49-tissue S-MultiXcan → **S11** |
| `82_fecd_eur_matched_twas.sh` | FUSION, European GWAS, 3 tissues, with `fusion_genomewide.sh` → **S3**, Figure 1, FECD values in Figure 5 (plus a genome-wide conditional/joint analysis of the Bonferroni-significant models, not reported) |
| `19c_tables1_2_twas_summary.py` | **Table 1** (S-MultiXcan ranking, `analysis/table2_smultixcan_ranking.tsv`), **Table 2** (NR1H3 per tissue, `analysis/table1_nr1h3_znf408_sprediXcan.tsv`), per-tissue FDR. The file names predate the swap of Tables 1 and 2 in manuscript v3.9. |
| `69b_fig1_fig2_v2.R` | **Figures 1–2** |
| `24b_fusion_11p11_conditional.sh` | conditional analysis within 11p11.2, European GWAS (NR1H3 independent in fibroblasts, not skin); runs minp 1 (→ **S5**) and minp 0.05 |
| `109_fecd_reported_loci.py` | 51 reported FECD loci and classification of TWAS genes → **S12** |
| `109b_insawang_11p11_lookup.py` | 11p11.2 in the Insawang et al. meta-analysis (Discussion) |
| `97b_egene_status_nr1h3_znf408.py` | GTEx eGene status of NR1H3 and ZNF408; skin genome-wide eQTLs absent from the LD reference |
| `110_nr1h3_endothelium_fecd_vs_control.R` | NR1H3 / ZNF408 (and TSPAN10, from release v14.9) in FECD vs control corneal endothelium (DESeq2 `~ batch + status` only; the sample sheet for **S14**) |

### 2 · FinnGen look-up → S4
| Script | Output |
|---|---|
| `84b_finngen_replication_H184_H185.sh` | NR1H3 / ZNF408 S-PrediXcan in H18.4 and H18.5 → **S4** |

### 3 · Colocalisation → Figure 3, S6
| Script | Output |
|---|---|
| `101c_fetch_coloc_allpairs.sh` (run with `GENES=data/tophits21/genes_coloc.tsv`) | all variant–gene pairs in each gene's cis window, eQTL Catalogue r6 (GTEx v8), 27 genes × 3 tissues |
| `102c_coloc_abf_allgenes.R` | coloc.abf, European GWAS (primary) and multi-ancestry GWAS (sensitivity); priors p1 = p2 = 10⁻⁴, p12 = 10⁻⁵, PP4 also at p12 = 5 × 10⁻⁶ |
| `81c_fecd_coloc_heatmap_fullcis.R` | **Figure 3** and **S6** (gene list: European S-MultiXcan ranks 1–23, down to NR1H3, plus four genes from earlier versions of S6; `S_MultiXcan_rank` from step 19f) (NR1H3 PP4 0.08 fibroblasts / 0.23 skin / 0.14 oesophagus; ZNF408 ≤ 0.60; PP4 ≥ 0.8 for PIDD1, EPB41L2, PANO1, TNS3, TSPAN10) |

### 3b · Fine-mapping at 11p11.2 → S13
| Script | Output |
|---|---|
| `112a_ld_1000g_11p11.sh` | 1000 Genomes GRCh38 European LD panel for chr11:46.2–48.3 Mb, and the European GWAS for the region |
| `112_finemap_11p11.R` | SuSiE-RSS (L = 10; LD-mismatch and allele-flip checks) and single-variant ABF credible sets: no 95% SuSiE credible set (max PIP 0.016, rs3136447); ABF 95% sets of 6,902 (European) and 5,118 (multi-ancestry) variants |
| `112b_eqtlcat_susie_11p11.py` | eQTL Catalogue SuSiE credible sets in the region (none for NR1H3 or ZNF408) |
| `113_focus_11p11.sh` | FOCUS gene-level fine-mapping: summed PIP ZNF408 0.39, ARHGAP1 0.28, NR1H3 0.21 (null 0.026); default-LD-block sensitivity ZNF408 0.47, NR1H3 0.43, ARHGAP1 0.33 |

### 4 · Causal inference for FECD → Figure 4, S7
| Script | Output |
|---|---|
| `91_mr_nr1h3_eur.R` | primary cis-MR, European GWAS (IVW OR 1.06, P 0.37; F ≈ 597); exposure effects standardised from the eQTLGen Z with each instrument's eQTLGen N and allele frequency (`data/eqtl_eqtlgen/nr1h3_eqtlgen_af.tsv`, step 83a) |
| `83_mr_nr1h3_eqtlgen_fecd.R` | multi-ancestry sensitivity → **S7** |
| `83a_extract_eqtlgen_nr1h3.sh TSPAN10 ENSG00000182612`, `91_mr_nr1h3_eur.R TSPAN10 17`, `83_mr_nr1h3_eqtlgen_fecd.R TSPAN10 17` | the same steps run for TSPAN10 (17q25.3, the colocalising gene outside established loci; LD reference chromosome 17): 2 instruments (F ≈ 171), IVW only; European OR 0.67 (0.47–0.95), P 0.025; multi-ancestry OR 0.81, P 0.21 → **S7** TSPAN10 rows (release v14.8) |
| `94_mr_ancestry_sidebyside.R`, `88_verify_nr1h3_fecd_mr.R` | ancestry comparison; independent check |
| `87_nr1h3_fecd_mr_forest.R EUR` | **Figure 4** |

### 5 · Cross-trait TWAS → Figure 5, S8, S9 (Panels A, B)
| Script | Output |
|---|---|
| `46_…` → `47_…` → `48_…` → `49_analyse_glaucoma.R` | glaucoma S-PrediXcan / S-MultiXcan |
| `75_glaucoma_brain_panel.sh` | fibroblast + brain cortex + skin panel (S-MultiXcan P 3.0 × 10⁻⁸) |
| `75b_glaucoma_smultixcan_nobrain.sh` | without brain cortex (P 8.4 × 10⁻⁶, not Bonferroni) |
| `77_glaucoma_wholeblood_panel.sh` | whole-blood S-PrediXcan (Z −7.63) |
| `51_…`–`54_analyse_iop.R`, `56_…`–`59_analyse_myopia.R` | IOP and refractive-error S-PrediXcan / S-MultiXcan (S9 Panel D) |
| `60b_glaucoma_fusion_allele_match.py` | glaucoma FUSION input matched to the LD-reference alleles at multi-allelic FinnGen sites → `data/gwas/glaucoma_gwas_fusion.tsv`, `analysis/glaucoma_fusion_allele_matching.tsv` (1,056 records replaced, 8 dropped) |
| `76h_glaucoma_wb_genomewide.sh` (+ `76h_glaucoma_wb_genomewide_rank.py`) | genome-wide whole-blood FUSION TWAS + conditional/joint with `fusion_genomewide.sh`, glaucoma, with the allele-matched input (`allele_corrected/`, reported) and the first-record input of releases up to v13 (`as_reported/`) → `analysis/glaucoma_wb_genomewide_joint_ranking.tsv`, `analysis/glaucoma_wb_genomewide_rerun_summary.tsv`; NR1H3 third in both (7,966 models with a P outside the MHC; Bonferroni 6.3 × 10⁻⁶; 20 jointly included genes) |
| `76i_crosstrait_wb_genomewide.sh` (+ `76i_crosstrait_wb_genomewide_summary.py`) | the same for IOP and refractive error → `analysis/crosstrait_wb_genomewide_summary.tsv` (NR1H3 status, comparison with the June 2025 run when its workbook is present), `analysis/crosstrait_wb_genomewide_joint_ranking.tsv` |
| `76d_fusion_wb_conditional.sh` | chromosome-11 whole-blood FUSION results of the genome-wide runs (steps 76h, 76i; conditional status with the genome-wide threshold) → **S8**, S9 Panel A, `analysis/nr1h3_crosstrait_fusion_conditional.tsv`; glaucoma also with the first-record input of releases up to v13, for comparison |
| `76e_fig5_glaucoma_label.R` | **Figure 5** |
| `76g_s8_crosstrait_fusion_wb_table.R` | **S8** table from the step-76d outputs |
| `76f_glaucoma_wb_conditional_rank.py` | genome-wide rank of NR1H3 for glaucoma (third, after LOXL1-AS1 and AFAP1), from the June 2025 workstation run |

### 6 · Cross-trait causal inference → Table 3, Figure 6, S9 Panels C–D
| Script | Output |
|---|---|
| `89c_table3_crosstrait_cismr.R` | **Table 3** (all estimators; bootstrap seed fixed) and S9 Panel C; exposure standardised with eQTLGen N and allele frequency |
| `89b_nr1h3_cismr_standardised_perSNP_N.R` | **Figure 6**, drawn from the step 89c output: glaucoma odds ratio, IOP mmHg, refractive error SD units (`analysis/nr1h3_cismr_fig6_values.tsv`) |
| `79b_coloc_nr1h3_glaucoma_fullcis.R` | NR1H3 eQTL × glaucoma colocalisation on full cis summary statistics (fibroblast PP4 0.41; whole blood 0.74, 0.59 at p12 = 5 × 10⁻⁶; skin 0.04; oesophagus 0.02) → S9 Panel D |
| `108b_coloc_fecd_refractive_17q25.R` | FECD x refractive error colocalisation at 17q25.3 (TSPAN10) and the refraction direction of the FECD risk allele (release v15.0) |
| `108_coloc_nr1h3_glaucoma_fullcis.R` | FECD × glaucoma GWAS colocalisation at 11p11.2 (PP4 0.01; PP3 0.75) |

### 7 · Figure 7
`70b_flowchart_NR1H3_v2.R` → **Figure 7**.

### 7b · Input audit (all GWAS)
`55c_audit_gwas_inputs.py` (run with `METAXCAN_PY`) applies one set of checks to every GWAS input of the TWAS,
conditional TWAS and cis-MR analyses: FECD multi-ancestry (S-PrediXcan/S-MultiXcan, cis-MR sensitivity), FECD
European (FUSION, 11p11.2 conditional FUSION, primary cis-MR), glaucoma, IOP and refractive error (S-PrediXcan,
whole-blood FUSION and its chr11 conditional analysis, cis-MR). Each input the pipeline wrote is linked back to the
source GWAS records and to the variants it is matched against (FUSION: the 1000G EUR LD reference; S-PrediXcan: the
GTEx v8 MASHR model variants; cis-MR: the eQTLGen NR1H3 cis-variants). Every reference variant is counted as matched,
missing, ambiguous (more than one source record for the key; the pipeline keeps the first) or allele-inconsistent, and
every input record is checked against its source record. The step then re-derives, from the audited inputs, the model
variants of NR1H3 and ZNF408 in every S-PrediXcan tissue, the weights of every FUSION model in the 11p11.2 conditional
regions (rebuilt as `FUSION.post_process.R` builds them) and the cis-MR instrument selection. Outputs in
`analysis/input_audit/`: `input_audit_summary.tsv` (counts), `input_audit_records.tsv` (each ambiguous-conflicting,
allele- or position-inconsistent or unverified record that a reported analysis uses), `twas_model_coverage.tsv`,
`fusion_conditional_regions.tsv`, `cis_instrument_coverage.tsv`, `mr_instruments_verified.tsv` and, when a conditional
region holds an allele-inconsistent weight SNP, `fusion_allele_sensitivity.tsv` (chr11 FUSION rerun with allele-consistent
records, joint model of the models below the genome-wide threshold; FUSION files in `results/input_audit_sensitivity/`). Since release v14 no conditional region does, so that
file is not written. The cis-MR MAF filter is applied to the eQTLGen allele frequency, as in the cis-MR steps.
Results: *Notes for reviewers*.

### 7d · Figure files
`139b_export_figures.sh` re-plots Figures 1–7 into `figures/` as `FigureN.pdf` (vector; the SVG is printed to PDF by
headless Google Chrome or Chromium, with the fonts embedded), `FigureN.svg` and `FigureN.tif` (600 dpi, LZW). R's
cairo devices need XQuartz on macOS and the base `pdf()` device cannot set Unicode text in Arial, hence the browser
step; without Chrome the PDF is skipped with a message. Set `CHROME` to the browser's path if it is not found.

### 7c · Additional outputs requested in review
`139_additional_outputs.py` collects the outputs requested in peer review (29 October 2026) in
`analysis/additional_outputs/` (bundled as `additional_outputs/`) and writes the tables derived for the request. It
re-estimates nothing. See *Additional outputs requested in review*.

### 8 · Check of the reported tables
`140_check_reported_tables.py` writes S1, S2, S3, S10 and S11 exactly as reported (`analysis/supplementary/`) and checks every value in manuscript Tables 1–3 and Supplementary Tables S1–S14 (copies in `reference/`) against the outputs of the run (`analysis/reported_tables_check.tsv`, mismatches in `analysis/reported_tables_mismatches.tsv`); the run fails on any mismatch.

---

## Notes for reviewers

- **Colocalisation input.** Every colocalisation now uses full cis summary statistics (all variant–gene pairs,
  eQTL Catalogue release 6 = GTEx v8), as coloc.abf assumes. Bundles up to v3.9 used only the significant pairs
  returned by the GTEx Portal (steps 83b, 89, 81b and 79), which these steps replace. The main changes: NR1H3 × FECD
  PP4 0.81 → 0.08 in fibroblasts (no colocalisation in any tissue); NR1H3 × glaucoma 0.68 → 0.41 (fibroblasts)
  and 0.61 → 0.74 (whole blood).
- **Bootstrap estimators.** Weighted-median and weighted-mode standard errors are bootstrapped; steps 89c
  (`set.seed(1)`) and 83, 88 and 91 (`set.seed(2024)`) fix the seed so Table 3 and S7 reproduce exactly. Steps 83 and
  91 also fix the TwoSampleMR exposure and outcome IDs (otherwise random strings), so their output files are
  byte-identical between runs.
- **Figure labels.** Every ggrepel label layer (steps 69b, 49, 54, 59) has a fixed seed and no time limit
  (`max.time = Inf`; by default ggrepel stops arranging labels after 0.5 s, so a crowded plot depended on the speed
  of the computer), so label positions in the PNG and PDF figures are the same in every run.
- **Input audit (step 55c; `analysis/input_audit/`).** The same checks are applied to every GWAS input
  (counts are of reference variants; "matched" = one source record, alleles consistent with the reference, and the
  input record reproduces the source record's alleles and Z):
  - *Every input record reproduces its source record*, for all five GWAS inputs and all three analysis types. No
    S-PrediXcan model variant is matched to a GWAS record at another position except two that S-PrediXcan does not
    use or uses correctly (glaucoma rs9282516, an indel written at a shifted position with a different allele set;
    IOP rs4993392, lifted by step 50 onto an alternate chr22 contig, with matching alleles).
  - *cis-MR.* For each outcome (FECD multi-ancestry and European, glaucoma, IOP, refractive error) the instrument
    selection rerun on the audited input returns the four reported instruments (rs7395581, rs11038863, rs3740694,
    rs4647760); each has exactly one source record whose alleles, Z, frequency and N match the cis-MR input
    (`mr_instruments_verified.tsv`). Of 2,444 eQTLGen NR1H3 cis-variants, 2,182–2,342 are matched per outcome; the
    allele-inconsistent ones (0–4 per outcome) are removed by TwoSampleMR harmonisation and none is an instrument. The
    selection counts equal those printed by steps 83 and 91 (rows after the rsID merge, `cis_instrument_coverage.tsv`).
  - *TWAS model variants.* Every NR1H3 and ZNF408 S-PrediXcan model, in every tissue and trait analysed, uses all of
    its model variants that the GWAS contains, and the counts equal `n_snps_used` in the S-PrediXcan outputs
    (117 gene–tissue models). For every FUSION model in the 11p11.2 conditional regions (FECD fibroblasts 16 models
    and skin 22, chr11:45.76–48.39 and 45.80–48.60 Mb, matching FUSION's own report; whole blood 17 models,
    chr11:45.09–48.60 Mb, for glaucoma, IOP and refractive error) the number of SNPs FUSION used equals the audit's,
    and at most 4 weight SNPs per region are missing and imputed (impute fraction ≤ 1.1%). NR1H3's weights are complete
    except one SNP imputed for FECD and two for glaucoma.
  - *FUSION inputs genome-wide* (1,190,321 LD-reference SNPs): FECD European 1,183,054 matched, 0 allele-inconsistent;
    IOP 1,189,381 and 13; refractive error 1,179,368 and 3. **Glaucoma** (up to release v13; corrected in v14 by step
    60b, after which the audit finds no allele-inconsistent glaucoma record in the FUSION input): FinnGen lists each
    allele of a multi-allelic site as a separate record under one rsID, and step 45 keeps the first. For 1,057
    LD-reference SNPs the kept record
    is a different allele from the reference SNP although FinnGen also has the matching record; as a result 575 SNPs
    reach FUSION with alleles that do not match the reference (FUSION uses such Z-scores without alignment) and 496 are
    dropped as strand-ambiguous (other GWAS: 0–13 and 0–4). One such SNP
    lay in the whole-blood conditional region of releases up to v14.5 (rs12798703, a weight SNP of PTPRJ and FNBP4, not
    NR1H3; since v14.6 neither gene enters the glaucoma joint model). Rerunning the chromosome-11 analysis of those releases with the allele-consistent FinnGen record for all 53 affected chr11 SNPs (`fusion_allele_sensitivity.tsv`)
    leaves NR1H3 unchanged (TWAS Z −7.72; joint Z −7.7, P 1.2 × 10⁻¹⁴) and the jointly included genes unchanged
    (BANF1, ME3, NR1H3); the largest change on chr11 is |ΔZ| 0.34 (RNASEH2C, 65.7 Mb). The glaucoma S-PrediXcan input
    is unaffected (S-PrediXcan drops records whose alleles differ from the model's).
  - *Other notes.* The IOP file's frequency column is the minor-allele frequency, not the effect-allele frequency;
    it is used only for the MAF filter and for TwoSampleMR's strand check of palindromic SNPs, and none of the
    instruments is palindromic. 476 LD-reference SNPs are missing from the IOP FUSION input because step 50 drops
    records whose liftOver fails, although FUSION does not use positions; none is a weight SNP in the conditional
    region.
- **GTEx-instrumented MR and SMR** (steps 31, 32, 78 and 90–93) are not part of this release (see *Changes for
  Final Draft v7*).
- **Gene names.** Current HGNC symbols are used where one exists (GATD1 = GENCODE v26 PDDC1; GATD1-DT = AP006621.5);
  otherwise the GENCODE v26 name (AP006621.6) or the Ensembl ID. The full supplementary tables S1, S2 and S11 keep
  the GENCODE v26 names distributed with the GTEx v8 models.

## Additional outputs requested in review (`additional_outputs/`)

Written by step 139 from the outputs of earlier steps and the downloaded inputs. Related files already in
`analysis_outputs/` are named below rather than duplicated.

**`glaucoma_allele_matching/`.** Since release v14 the glaucoma FUSION input is matched to the LD-reference alleles
(step 60b): where FinnGen lists a multi-allelic site as several records under one rsID and the first record (kept by
step 45, and the input of releases up to v13) does not match the 1000G EUR reference, the first FinnGen record whose
alleles match is used (1,056 records); 8 records with no matching record are dropped so that FUSION imputes them
instead of using them unaligned. This folder compares the reported (allele-matched) results with the first-record
results of releases up to v13. `allele_matching_records.tsv` lists every record replaced or dropped.
- `chr11/`: the chromosome-11 rows of the genome-wide whole-blood runs (step 76h; GTEx v8 EUR weights; conditional/joint
  analysis of the models below the genome-wide threshold, 1000G EUR LD, 100-kb window), all 463 chr11 models,
  `reported_allele_matched/` (behind S8, step 76d) and `first_record_v13/`; `chr11_first_record_vs_reported.tsv` puts
  all 463 models side by side. 41 models differ (largest |ΔZ| 0.34, RNASEH2C: Z 1.17 → 1.50); NR1H3 is unchanged
  (TWAS Z −7.72; JOINT Z −7.7, P 1.2 × 10⁻¹⁴), the same two genes are jointly included (BANF1, NR1H3) and no gene
  changes joint/conditional status. (Up to release v14.5 this folder held a chromosome-11-only run in which every chr11
  model entered the joint model; BANF1, ME3 and NR1H3 were jointly included.)
- `genomewide/`: step 76h on all 22 chromosomes with both inputs: all marginal results outside the MHC (7,980 models,
  7,966 with a P-value; the 145 MHC models are in `glaucoma_wb_MHC.dat` upstream and not analysed further), the jointly
  included and dropped genes of every chromosome, the Bonferroni threshold (0.05/7,966 = 6.3 × 10⁻⁶) and the
  post_process logs. 1,058 models differ (largest |ΔZ| 0.52, GTF2H5); NR1H3 is third with both inputs (JOINT Z −7.7,
  P 1.2 × 10⁻¹⁴), behind LOXL1-AS1 and AFAP1, and the 20 jointly included genes are the same. Up to release v14.5 the
  MHC models were analysed too and four MHC genes (DDX39B, RPL32P1, HLA-DQB1, AGPAT1, all ranked below NR1H3) were
  jointly included (24 genes).
  Summaries in `analysis_outputs/`: `glaucoma_wb_genomewide_joint_ranking.tsv`, `glaucoma_wb_genomewide_rerun_summary.tsv`.

**cis-MR instruments** (`analysis_outputs/nr1h3_cismr_instruments_exposureN.tsv`, step 89c). For each outcome
(glaucoma, IOP, refractive error = "Myopia") and each of the four eQTLGen whole-blood NR1H3 instruments (rs11038863,
rs3740694, rs4647760, rs7395581): the eQTLGen Z-score (`exposure_Z`), the eQTLGen sample size for that variant
(`NrSamples`), the minor-allele frequency in eQTLGen (`maf_eqtlgen`, used since release v14) and in the outcome GWAS
(`maf_outcome_gwas`, used up to v13), and the standardised exposure beta and SE under a fixed N of 31,684
(`*_fixedN`), under the variant's own N (`*_variantN`, used for Table 3 and Figure 6; both with the eQTLGen frequency)
and as computed up to v13 (`*_v13_outcome_maf`). The eQTLGen and FinnGen frequencies differ most (rs7395581 0.32 vs
0.40; rs4647760 0.026 vs 0.015), so the glaucoma estimates changed most (IVW OR 0.83 → 0.84, P 2.0 × 10⁻⁵ →
6.9 × 10⁻⁶).

**`fecd_european_smultixcan/`.** The European FECD S-MultiXcan behind the `S_MultiXcan_rank` column and the gene list
of S6 (steps 19d, 95a, 19e, 19f).
- `fecd_eur_smultixcan_remapped.txt`: MetaXcan v0.8.1-1-g9b9f76b (its internal version string reads 0.7.5); GTEx v8
  MASHR models of cultured fibroblasts, skin (not sun-exposed) and oesophagus mucosa; covariance
  `gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz`; condition-number cutoff 30. GWAS GCST90475865
  (3,473 cases / 445,239 controls, GRCh38), prepared as the multi-ancestry GWAS is (beta = ln OR, SE from the 95% CI)
  and matched to the model variants by chromosome, GRCh38 position and allele pair (93% of model SNPs used per
  tissue). 17,345 genes, 17,303 with a P-value; NR1H3 23rd (P 7.07 × 10⁻⁵), ZNF408 21st (P 5.76 × 10⁻⁵).
- `per_tissue_S-PrediXcan/`: the three per-tissue S-PrediXcan results it combines.
- `S6_rank_lookup.tsv`: each S6 gene's rank in Supplementary Tables v2.15/v2.16, in the 4 July 2026 run, in v2.17
  (this run) and in the multi-ancestry S-MultiXcan of Table 1.
- `previous_run_4Jul2026_rsID_matched/`: the run behind S6 up to v2.16, kept for reference (input matched by rsID,
  so model variants without an rsID were not used; 88% of model SNPs per tissue; 16,922 genes with a P-value). It is
  not rerun by this release; step 139 copies it if present.
- Already in `analysis_outputs/`: `fecd_eur_smultixcan_ranking.tsv` (all 17,303 ranked genes) and
  `fecd_eur_smultixcan_S6_genes.tsv` (the 27 S6 genes with their European rank and why each is included).
- What changed in v12, and how the European ranks relate to Table 1: *Changes for release v12*. The European ranks
  1–23 are Table 1's 21 genes plus TNS3 and CUL3 (multi-ancestry ranks 25 and 40).

**`rnaseq_batch_metadata/`.** Basis for assigning the two GSE101872 controls to batch 1 in step 110. GEO gives no
batch label for GSE101872, so the assignment is inferred from the GEO and SRA records:
1. GSE101872 is the sequencing run that GSE112201 calls batch 1. Its five FECD samples are the same data as five of
   the six GSE112201 batch-1 samples: the gene counts are identical at all 23,398 genes
   (`shared_donor_count_identity.tsv`), and so are the SRA read counts (e.g. 2011-024 = RNA16, 74,393,209 read
   pairs; `sra_runinfo_extract.tsv`).
2. All seven GSE101872 samples (both controls included) and all six GSE112201 batch-1 samples were sequenced on an
   Illumina HiSeq 2000 (GPL11154) with 2 × 51 bp paired reads, TruSeq RNA v2 libraries and the same processing
   (MAP-RSeq v1.2, TopHat 2.0.6, HTSeq-count 0.5.3p9, hg19). The GSE101872 samples were deposited together in July
   2017 (`geo_sample_metadata.tsv`, `sra_runinfo_extract.tsv`).
3. GSE112201 batches 2 and 3, and its four controls (which GEO labels batch 2), were sequenced on an Illumina
   HiSeq 4000 (GPL20301) with 2 × 101 bp reads.

Batch 1 is therefore the only batch that shares instrument, read length and deposition with the GSE101872 controls.
Sources: GEO series matrices GSE101872, GSE112201-GPL11154 and GSE112201-GPL20301 (downloaded by step 00); NCBI SRA
run information (fetched by step 139 from NCBI E-utilities and cached as `data/endothelium_rnaseq/sra_runinfo.xml`).
The sample-level assignments are in S14.

## Changes for release v15.2 (13 November 2026)

- **S4 without the power column.** The manuscript no longer reports power for the FinnGen look-up, so Supplementary
  Table S4 (workbook v2.22) drops its power column, step `84c_finngen_replication_power.py` is removed from the
  pipeline, and step 140 no longer checks that column. Every other table cell, script and output is unchanged.
- **Checks.** Step 140: 0 mismatches against Tables 1–3 and Supplementary Tables S1–S14 v2.22.

## Changes for release v15.1 (13 November 2026)

- **One script per analysis.** The TSPAN10 analyses no longer use copied scripts. Steps `83a_extract_eqtlgen_nr1h3.sh`,
  `91_mr_nr1h3_eur.R` and `83_mr_nr1h3_eqtlgen_fecd.R` take the gene (and LD-reference chromosome) as arguments, and
  step `108_coloc_nr1h3_glaucoma_fullcis.R` takes the second trait and the locus; the master script runs each once with
  its defaults (NR1H3; glaucoma at 11p11.2) and once for TSPAN10 (refractive error at 17q25.3). The copies
  `83a_extract_eqtlgen_tspan10.sh`, `91_mr_tspan10_eur.R`, `83_mr_tspan10_eqtlgen_fecd.R` and
  `108b_coloc_fecd_refractive_17q25.R` are removed. All NR1H3 outputs, the glaucoma colocalisation and all TSPAN10 MR
  outputs are byte-identical to release v15.0.
- **Correction, FECD x refractive error.** Step 108b of v15.0 filtered the refractive-error GWAS by position but not by
  chromosome (`CHR == CHR` compared the column with itself), so records from other chromosomes at the same position
  and allele pair could be matched to chromosome-17 FECD variants (7,868 matched variants; 7,720 with the chromosome
  filter). Corrected values: PP4 = 0.865 (European; was 0.864) and 0.835 (multi-ancestry; was 0.834); the top variant
  (rs9895741) and the direction at the FECD lead rs6565597 (risk allele T, more myopic, Z -11.1) are unchanged.
- **Input checks for the TSPAN10 instruments.** rs7217415 (A/G) and rs34232035 (T/C) are on chromosome 17 in eQTLGen,
  the FECD GWAS and the LD reference, are not strand-ambiguous, and their allele frequencies agree across sources. One
  multi-allelic site in the region (rs72858871) has two FECD records; its eQTL Z (4.5) is below the instrument
  threshold, so it does not enter the analysis.
  Of 518 eQTLGen TSPAN10 cis-eQTL variants, 513 are in the FECD GWAS; six genome-wide-significant ones are absent from
  both FECD files and from the LD reference, so they could not have been selected by clumping. The full input audit
  (step 55c: allele, position and multi-allelic checks of every TWAS, conditional and cis-MR input) was rerun and
  its six output files are byte-identical to release v15.0.
- **Checks.** Step 140: 0 mismatches. The changed steps were rerun from this bundle's scripts in the clean folder of the
  v14 reproducibility check: all 19 outputs are byte-identical to the project's.

## Changes for release v15.0 (13 November 2026)

- **Table 1 "P value of locus" uses MAF >= 1%.** Step 19c now takes the smallest multi-ancestry FECD GWAS P within 1 Mb
  of each gene's TSS among variants with minor allele frequency >= 1%, the filter of the colocalisation, fine-mapping and
  cis-MR analyses. Without it, ultra-rare variants set the value at three loci outside the established regions:
  ZNF408 3.5e-5 -> 3.6e-5, SHISA9 1.1e-7 -> 5.1e-5, TSPAN10 2.5e-7 -> 5.2e-5 (rs117831598, MAF 0.0002, had set the
  TSPAN10 value). Every other row, and Table 2, is unchanged; `reference/manuscript_tables.tsv` is updated.

- **FECD x refractive error at 17q25.3 (new step 108b).** `108b_coloc_fecd_refractive_17q25.R` colocalises the FECD
  GWAS (European primary, multi-ancestry sensitivity) with the Hysi et al. refractive-error GWAS at the TSPAN10 locus
  (TSS +/- 1.05 Mb), with the coloc.abf set-up of step 108 (refractive-error effects standardised from Z, per-variant N and
  allele frequency). PP4 = 0.86 (European) and 0.83 (multi-ancestry); the FECD risk allele at the lead variant rs6565597
  (T; OR 1.11, P 7.3e-5) is associated with more myopic refraction (Z -11.1, P 1.8e-28). Outputs:
  `analysis_outputs/coloc_17q25_fecd_refractive.tsv`, `..._direction.tsv`. Note: the smallest FECD P within 1 Mb of the
  TSPAN10 TSS (2.5e-7, Table 1 "P value of locus") is an ultra-rare variant (rs117831598, MAF 0.0002) absent from the
  refractive-error GWAS; the common-variant lead is rs6565597.
- **Checks.** Step 140: 0 mismatches against Tables 1–3 (Table 1 as above) and Supplementary Tables S1–S14 v2.21. Steps
  19c and 108b were rerun from this bundle's scripts in the clean folder of the v14 reproducibility check; their outputs
  are byte-identical to the project's. No other script or output changed.

## Changes for release v14.9 (13 November 2026)

- **TSPAN10 in step 110.** `110_nr1h3_endothelium_fecd_vs_control.R` (corneal endothelium, FECD vs control; DESeq2
  ~ batch + status, unchanged) now also reports TSPAN10, the colocalising gene outside established loci:
  median normalised count 40 in FECD and 9 in controls, log2 fold change 1.82, nominal P 0.006
  (`analysis_outputs/nr1h3_endothelium_fecd_vs_control.tsv`, as reported in the Results). The NR1H3, ZNF408 and TCF4
  rows are byte-identical to release v14.8.
- **Checks.** Step 140: 0 mismatches against Tables 1–3 and Supplementary Tables S1–S14 v2.21 (no table changed). Step
  110 was rerun from this bundle's scripts in the clean folder of the v14 reproducibility check; its outputs are
  byte-identical to the project's. No other script or output changed.

## Changes for release v14.8 (13 November 2026)

- **TSPAN10 cis-MR (new S7 rows).** Steps `83a_extract_eqtlgen_tspan10.sh`, `91_mr_tspan10_eur.R` and
  `83_mr_tspan10_eqtlgen_fecd.R` are copies of the NR1H3 steps 83a, 91 and 83 with only the gene (ENSG00000182612) and
  the LD-reference chromosome (17) changed; the master script runs them after the NR1H3 steps. Two independent
  eQTLGen whole-blood instruments remain (instrument-set F ≈ 171), so only the IVW estimate is computed: European FECD
  GWAS OR 0.67 per SD (95% CI 0.47–0.95), P 0.025; multi-ancestry OR 0.81 (0.58–1.13), P 0.21. The estimate rests
  largely on rs7217415 (Wald ratio OR 0.63, P 0.023).
- Supplementary Tables v2.21: S7 has one TSPAN10 IVW row per FECD GWAS (egger_intercept n/a), legend and Contents
  entry updated; every other cell is unchanged. Step 140 joins S7 on outcome GWAS, exposure and method.
- **Checks.** Step 140: 0 mismatches against Tables 1–3 and Supplementary Tables S1–S14 v2.21. Clean-folder check: the
  three new steps and step 140 were run from this bundle's scripts in the folder of the v14 reproducibility check; all
  eight TSPAN10 outputs (exposure, allele-frequency and outcome files; MR results and diagnostics) are byte-identical
  to the project's, and step 140 gives 0 mismatches (S12 not checked without the restricted Insawang et al. file). No
  other script or output changed.

## Changes for release v14.7 (12 November 2026)

- **Wording only: no analysis, output, table value or figure changed.** *Notes for reviewers*, item 4, said that the
  refractive-error joint model keeps SLC39A13 and C1QTNF4 "with inflated effects of opposite sign to their marginal
  ones". Only SLC39A13 reverses in sign (marginal Z 6.3, JOINT.Z −4.2); C1QTNF4 keeps its sign (marginal Z −6.4,
  JOINT.BETA −16.1). Item 4 now says that all three genes have inflated joint effects and that SLC39A13's also
  reverses in sign; the S9 Panel A note says the same.
- `reference/`: Supplementary Tables workbook v2.20 (v2.19 with only that S9 note changed; every value identical).
- **Checks.** Step 140: 0 mismatches against Tables 1–3 and Supplementary Tables S1–S14 v2.20. The scripts are those
  of v14.6 (only the release label in `run_full_pipeline.sh` changed), so v14.6's clean-folder check applies.

## Changes for release v14.6 (12 November 2026)

- **One genome-wide FUSION procedure for every trait.** New `scripts/fusion_genomewide.sh` (per-chromosome inputs
  restricted to the 1000G EUR LD-reference SNPs; FUSION.assoc_test; models in the MHC kept aside; genome-wide
  Bonferroni threshold; FUSION.post_process, 100-kb window, on the models below it) is now run by step 82 (FECD,
  3 tissues), step 76h (glaucoma, both inputs; `76h_glaucoma_wb_genomewide_inputs.py` is no longer needed) and new step
  76i (IOP and refractive error, genome-wide, as the June 2025 workstation run had them). The GTEx v8 EUR whole-blood
  weights are now downloaded by step 39. The master script runs 76h and 76i before 76d.
- **FECD unchanged.** The FUSION results of the three tissues (S3, Figure 1) and the 11p11.2 conditional analysis (S5)
  are byte-identical to release v14.5 (models in the MHC were never in them). Step 82 also writes a genome-wide
  conditional/joint analysis of its Bonferroni-significant models (not reported).
- **Glaucoma.** Marginal results identical. The MHC models are no longer in the conditional analysis: 20 jointly included
  genes instead of 24 (DDX39B, RPL32P1, HLA-DQB1 and AGPAT1 drop out, all ranked below NR1H3); NR1H3 is still third,
  behind LOXL1-AS1 and AFAP1 (JOINT Z −7.7, P 1.2 × 10⁻¹⁴).
- **S8 and the cross-trait conditional status.** Step 76d takes the chromosome-11 results from the genome-wide runs, so
  the conditional status uses each trait's genome-wide threshold (6.3 × 10⁻⁶) instead of every chr11 model entering
  the joint model (0.05/463). Refractive error: NR1H3 conditionally independent (JOINT Z 5.2, P 2.4 × 10⁻⁷; rank 78 of
  131 genome-wide), alongside the correlated neighbours SLC39A13 and C1QTNF4 (see *Notes for reviewers*, item 4); it
  was conditioned away (P 4.0 × 10⁻⁴). IOP: NR1H3 below the genome-wide threshold, not in the conditional analysis (it
  was conditioned away, P 0.26). Glaucoma: unchanged. Supplementary Tables v2.19 (S8: conditional status of every chr11
  row; S9: Panel A, its note, Panel D and the legend); every other table is unchanged.
- **Figure 5.** The glaucoma label reads Z = −7.72 (was −7.70): step 76d now reads NR1H3's marginal Z at full precision
  from the results file (P = 2Φ(−|Z|), FUSION's formula) instead of the rounded joint-model file. Figures 1–4, 6 and 7
  are unchanged.
- **Input audit (55c) and additional outputs (139)** follow the new runs: the whole-blood conditional regions are rebuilt
  from the models that entered the joint model.
- `fig_common.R`: headless Chrome PDF printing is retried up to three times (it occasionally stalled under heavy CPU
  load). Run-order table: the no-brain-cortex S-MultiXcan P is 8.4 × 10⁻⁶ (the output; was given as 8.5 × 10⁻⁶).
- **Checks.** Step 140: 0 mismatches against Tables 1–3 and Supplementary Tables S1–S14 v2.19. All seven figures
  regenerate; only Figure 5 differs from v14.5 (the label above). Clean-folder check: the steps this release changes
  (39, 82, 24b, 76h, 76i, 76d, 76g, 76e, 76f, 70b, 55c, 139, 139b, 140) were rerun from this bundle's scripts in the
  folder of the v14 reproducibility check (public inputs only): step 140 gives 0 mismatches (S12 not checked without
  the restricted Insawang et al. file), and 976 of 990 outputs are byte-identical to the project's. The others differ
  only in the run-folder path (4 FUSION .report files), the PDF creation date (7 figures; identical when rendered) or
  the non-public inputs (the June 2025 comparison in `crosstrait_wb_genomewide_*.tsv`; the 4 July 2026 S-MultiXcan run
  kept in `additional_outputs/fecd_european_smultixcan/`), plus one log left from an earlier project run.

## Changes for release v14.5 (12 November 2026)

Figure 7 only (step 70b): row (v) of the glaucoma column now has two boxes, as in the FECD column. "Colocalisation"
gives the NR1H3 eQTL x glaucoma PP4 in whole blood (0.74) and cultured fibroblasts (0.41) and the FECD x glaucoma
GWAS PP4 at 11p11.2 (0.01), read from the step 108 output `analysis/coloc_11p11_glaucoma_fullcis.tsv`; "cis-MR"
gives the glaucoma IVW estimate as before. Row (v) is 1.75 in high in both columns so the colocalisation box fits.
No analysis, table or other figure changed; v14's reproducibility check applies unchanged.

## Changes for release v14.4 (12 November 2026)

Figure 7 only (step 70b): the glaucoma column now follows the rows of the FECD column. Row (ii) shows the glaucoma
GWAS preprocessing (allele matching to model and LD-reference variants; steps 46, 60, 60b), row (iii) the glaucoma
TWAS methods (S-PrediXcan + S-MultiXcan in brain cortex, fibroblasts and skin; FUSION in whole blood) and row (iv)
the NR1H3 results (FUSION whole-blood Z and conditional status, read from the step 76d output
`analysis/nr1h3_crosstrait_fusion_conditional.tsv`; S-MultiXcan P with and without brain cortex). An arrow from
"Conditional analysis and fine-mapping" (locus prioritisation) into the glaucoma column shows that the cross-trait
analyses follow up the prioritised FECD locus. The GTEx v8 box of the glaucoma column reads "MASHR models and FUSION
weights". No analysis, table or other figure changed; v14's reproducibility check applies unchanged.

## Changes for release v14.3 (11 November 2026)

Documentation only: the download steps of the cross-trait GWAS state the sample sizes and the download links, and
*Data sources* lists them. Step 45: FinnGen R12 `H7_GLAUCOMA`, 26,591 cases / 473,757 controls (FinnGen R12 manifest).
Step 50: IOP, GCST009413, 31,269 Europeans. Step 55: refractive error, GCST010003 public subset without 23andMe,
351,091 Europeans (GWAS Catalog; the per-variant METAL weight, at most 260,974, is the per-variant N used). Every
link was checked on 11 November 2026. No code that runs, analysis, table or figure changed.

## Changes for release v14.2 (11 November 2026)

Figure 7 only (step 70b): two labels in the study-design row now give the samples the analyses used. The 1000 Genomes
box reads "n = 489" (was 503): FUSION and the cis-MR clumping use the FUSION LD reference (`data/ld_ref/LDREF`, 489
Europeans, HapMap3 variants, GRCh37); 503 is the GRCh38 1000 Genomes panel of the fine-mapping steps (112a, 112, 113).
The GTEx v8 box reads "403–430 EUR donors" (was "~500 donors each"), the European donors behind the FUSION weights
(403 cultured fibroblasts, 411 oesophagus mucosa, 430 skin not sun exposed). No analysis, table or other figure
changed; v14's reproducibility check applies unchanged (step 140 does not involve Figure 7).

## Changes for release v14.1 (10 November 2026)

Figure 7 only (step 70b): the glaucoma column now has a single arrow running down the centre line from between the
FinnGen R12 and GTEx v8 boxes to "Glaucoma TWAS" (v14 joined the two boxes with a bracket). No analysis, table or
other figure changed; v14's reproducibility check applies unchanged (step 140 does not involve Figure 7).

## Changes for release v14 (10 November 2026)

Final review (10 November 2026). Every value in Tables 1–3 and Supplementary Tables S1–S14 (v2.18) matches the
outputs of this run (step 140, 0 mismatches).
- **Glaucoma FUSION input allele-matched** (new step 60b; steps 60, 76d, 76h, 139 and 55c updated). The
  allele-corrected glaucoma outputs released as a sensitivity analysis in v11–v13 are now the reported results. S8:
  41 glaucoma rows change (TWAS Z/P, conditional Z/P); no conditional status changes. NR1H3 (Z −7.72; joint P
  1.2 × 10⁻¹⁴), the jointly included genes and the genome-wide rank (third) are unchanged; Figure 5 and S9 Panel A are
  unchanged. Comparison with the first-record input: `additional_outputs/glaucoma_allele_matching/`.
- **cis-MR exposure standardised with eQTLGen allele frequencies** (steps 83a, 83, 88, 91, 89b, 89c, 55c). Up to v13
  the Z → beta conversion and the MAF filter used the outcome GWAS frequency. The four instruments, F (597) and the
  FECD conclusions are unchanged (European IVW OR 1.06, 95% CI 0.93–1.22, P 0.37; multi-ancestry P 0.31 → 0.32).
  Cross-trait (Table 3): glaucoma IVW log-OR −0.181 → −0.176 (OR 0.84), P 2.0 × 10⁻⁵ → 6.9 × 10⁻⁶; IOP −0.205 mmHg,
  P 0.066; refractive error 0.043 SD, P 1.3 × 10⁻⁷ → 7.3 × 10⁻⁸; I²GX 0.994 for all outcomes. S7, S9 Panels C–D and
  Table 3 updated.
- **Figure 6 on each outcome's own scale** (step 89b rewritten). The common "standardised beta" axis treated
  glaucoma as an effect on the observed 0/1 scale, which depends on FinnGen's case fraction; the figure now draws the
  four estimators of Table 3 per outcome: odds ratio for glaucoma, mmHg for IOP, SD units for refractive error.
- **Figures** (new `scripts/fig_common.R`, step 139b): one style for Figures 1–7, gene symbols in italics, P-values as
  a × 10^b, broken y-axes in Figures 1 and 2 (TCF4/LINC01415), larger text in Figure 7 (12.5–14 pt on a 10-in-wide
  canvas), vector PDFs with embedded fonts plus SVG and 600-dpi TIFF in `figures/`. Figures 1–6 are drawn at the
  journal's full page width (190 mm) so text prints at its stated size; Figure 3 fits one page. The background points
  of Figures 1 and 2 (9,505 and 17,285 genes) are rasterised at 600 dpi inside the vector PDF (ggrastr), because
  headless Chrome did not reliably print SVGs with that many elements; all text, labels and highlighted genes stay
  vector.
- **Step 109 robustness**: the Ensembl REST gene lookup is retried with back-off (it returned HTTP 503 during the
  clean-folder test) and, if Ensembl stays unavailable, the gene positions this release used are read from
  `data/derived/ensembl_gene_tss.tsv` (GRCh38 TSS of the 20 genes; identical outputs either way).
- **Supplementary Tables v2.18**: S7, S8 and S9 values and legends as above; gene symbols in italics throughout.
- **Licence** (MIT, `LICENSE`) and **restricted data**: the two tables derived from the Insawang et al. summary
  statistics (`insawang2026_*.tsv`) are no longer distributed; steps 109, 109b and 140 run without them (S12 then
  reported as not checked).

## Reproducibility check (10 November 2026, v14 release)

Release v14 was run from an empty folder on an Apple M1 laptop (8 GB RAM) with only the files in this release and
the public inputs; as in the v8 check, the large raw downloads were linked from an earlier download and every other
file was regenerated. The restricted inputs were absent (Insawang et al. summary statistics; the June 2025 glaucoma
run), as for any user of this release.
- **Step 140:** every value in manuscript Tables 1–3 and Supplementary Tables S1–S14 (v2.18) matches the outputs
  (0 mismatches; S1, S2, S3, S10 and S11 cell by cell). S12 is reported as not checked, as designed, because it
  needs the Insawang et al. summary statistics.
- **Comparison with the project's outputs** (911 files in `analysis/`, `results/` and `figures/`): 880 are
  byte-identical, including every table, MR, colocalisation, fine-mapping and FUSION result and the allele-matched
  glaucoma input of step 60b. The other 31: 17 PDF figures (creation date only), 4 FUSION `.report` files
  (absolute paths), 4 logs, 1 gzip file (header only; identical content), 3 files that depend on the restricted
  Insawang file (`fecd_reported_loci.tsv`, `fecd_twas_gene_locus_class.tsv`, `reported_tables_check.tsv`; Tables 1
  and 2 are unaffected), `S6_rank_lookup.tsv` (the project also holds the unpublished 4 July 2026 run it is
  compared with) and the diagnostic `_snp_usage.tsv` of step 95, cut short when an accidental repeat of step 95 was
  stopped during the test (the 49-tissue results themselves are identical; S1, S10 and S11 match).
- **Interruptions:** the run stopped twice and was resumed at the failed step: Ensembl REST was unavailable (HTTP
  503) at step 109, which now retries and falls back to `data/derived/ensembl_gene_tss.tsv`; and a disk-space
  guard stopped it at step 53 on a nearly full disk. Step 11 needs about 10 GB of memory and swapped heavily on
  8 GB (2.5 h).

## Changes for release v13 (30 October 2026)

- **Additional outputs requested in review are now part of the release** (`additional_outputs/`, new step 139, and
  `analysis_outputs/nr1h3_cismr_instruments_exposureN.tsv`). Until v12 they were sent as a separate package. Step
  139 rebuilt 32 of that package's files from this release's outputs (the SRA run information re-fetched from NCBI),
  and all 32 are byte-identical to the package; its other data files were already in `analysis_outputs/`
  (identical). Step 139 runs before
  step 140, which remains the final step. No analysis and no reported value changed.
- **Provenance** (route 3) now records that step 76h reruns the genome-wide glaucoma ranking.
- **Manuscript** (Final Draft v5): Code Availability names this release. It is otherwise identical to Final Draft
  v4.

## Changes for release v12 (30 October 2026)

- **European S-MultiXcan input now matched by chromosome, GRCh38 position and allele pair (steps 19d–19f).** In
  release v11 the European input reproduced the 4 July 2026 file, which matched the MASHR model variants by rsID, so
  model variants without an rsID in the model were not used. 19d now prepares the European GWAS as steps 11 and 17
  prepare the multi-ancestry GWAS and runs step 95a on it, so both analyses use the matching the Methods describe.
  93% of model SNPs are used in each tissue (88% before), and 17,303 genes have a P-value (16,922 before). Most
  P-values are unchanged (NR1H3 7.07 × 10⁻⁵, ZNF408 5.76 × 10⁻⁵), but three genes gain model SNPs and rise: CCDC68
  (18th → 3rd), GATD1 (714th → 11th) and TALDO1 (566th → 12th). NR1H3 is now 23rd (was 21st) and ZNF408 21st (was
  19th). Every gene ranked down to NR1H3 was already in S6, so no colocalisation changed: rerunning 102c and 81c
  changed only the rank column (all posterior probabilities and the order of S6 identical; Figure 3 identical
  in content). The 4 July outputs (`results/multixcan/fecd_eur_smultixcan.txt`) are kept for reference.
- **S6 gene list** (`data/tophits21/genes_coloc.tsv`): its rank column is now written by 19f from the European
  ranking (the gene set is unchanged: European ranks 1–23 plus CHID1, RPLP2, LAMC1-AS1 and ENSG00000279672).
  CHID1 (39) and RPLP2 (3,401) now show their European rank instead of 0; LAMC1-AS1 and ENSG00000279672 have no
  European P-value and are blank. `data/tophits21/genes.tsv` (the old ranks 1–21 list) is no longer used or bundled.
- **Supplementary Tables v2.17** (`reference/`): S6 rank column and legend only (150 rank cells change). Step 140
  passes with v2.17.
- **RNA-seq (step 110): the CTG18.1 subgroup model was removed.** 110 no longer fits `~ batch + group` and no
  longer writes the expansion-positive and expansion-negative estimates (`log2FC_REpos`, `P_REpos`, `log2FC_REneg`,
  `P_REneg`); the FECD-versus-control results are unchanged. The expansion status stays in the S14 sample sheet as a
  sample description. No table reports a subgroup result.
- **Manuscript** (Final Draft v4): the subgroup P-values were removed from the Results, the matching Methods clause
  was removed, and Code Availability names this release.

## Changes for release v11 (29 October 2026)

- **European S-MultiXcan added (steps 19d, 19e, 19f).** The `S_MultiXcan_rank` column of S6 and the S6 gene list
  (`data/tophits21/genes.tsv`) come from the European-ancestry FECD S-MultiXcan run on 4 July 2026, which earlier
  releases did not contain. 19d rebuilds its input from GCST90475865: the rebuilt file is byte-identical to the 4 July
  input (66,928 rows; SE = (ln upper − ln lower) / (2 × 1.959963985)). Variants were matched to the five MASHR models
  in `data/predictdb` by the model's rsID (model SNPs without an rsID were not used), unlike the multi-ancestry input
  (95a: position and allele pair). 19e runs S-PrediXcan in the three proxy tissues and S-MultiXcan with the
  4 July commands; the rerun reproduced every S-PrediXcan Z and all 17,021 S-MultiXcan P-values. 19f ranks the
  16,922 genes with a P-value and confirms that S6 ranks 1–21 are the European ranks 1–21. NR1H3 is 21st here and in
  the multi-ancestry ranking of Table 1; other genes differ (e.g. ZNF408 19 European, 17 multi-ancestry).
- **Supplementary Tables v2.16** (`reference/`): only the S6 legend changes. It now defines `S_MultiXcan_rank` as the
  European rank among 16,922 genes, explains rank 0 (six genes kept from the previous S6, outside European ranks
  1–21: CHID1 41, TALDO1 566, GATD1 714; LAMC1-AS1, RPLP2 and ENSG00000279672 have no European P-value), describes
  the gene list as the European top 21 (it said "the Table 2 genes") and names the ranking files. Every other cell is
  identical to v2.15.
- **Genome-wide glaucoma rerun (step 76h).** The genome-wide rank of NR1H3 among conditionally independent glaucoma
  whole-blood signals (third, behind LOXL1-AS1 and AFAP1) came from a June 2025 workstation run (step 76f) and had
  not been repeated after the allele-matching correction of v9 (step 55c), which had been applied to chr11 only.
  76h reruns FUSION on all 22 chromosomes (8,125 models, 8,111 with a P; Bonferroni 0.05/8,111 = 6.2 × 10⁻⁶;
  FUSION.post_process per chromosome with a 100-kb window, MHC models included) twice: with the input as reported and
  with step 55c's correction applied genome-wide (1,056 records replaced: 489 strand-ambiguous or not a single A/C/G/T,
  567 allele-inconsistent; 15 had no allele-consistent record). The chr11 files reproduce step 55c's chr11 input and
  output exactly.
  - The correction changed the TWAS Z of 1,056 models (largest |ΔZ| 0.52, GTF2H5) but none of the jointly included
    genes: both runs include the same 24 genes in the same order, and NR1H3 is third (JOINT Z −7.7, P 1.2 × 10⁻¹⁴)
    behind LOXL1-AS1 (−13.7) and AFAP1 (−9.7). AFAP1's model (lasso) gives zero weight to its one corrected SNP.
  - The rerun reproduces the June 2025 ranks 1–18 (LOXL1-AS1 to ZBTB34) and its marginal Z-scores to the two
    significant figures reported. SF3A3 and NBPF2P swap places 19 and 20 (both JOINT P 4.4 × 10⁻⁶ at two significant
    figures). In the MHC, June included HLA-DQB1-AS1 (TWAS Z −5.6); in the rerun its Z is −3.9 under both inputs
    and HLA-DQB1 (Z −5.1) is included instead. The cause of this MHC difference is not known; it does not involve
    the allele correction or any gene ranked above 19.
- **Manuscript** (Final Draft v3): Code Availability names this release; the RNA-seq Methods sentence on a batch plus
  group model was removed.

## Changes for release v10 (27 October 2026)

- **Retinal colocalisation removed.** The manuscript no longer reports the colocalisation of NR1H3 retinal eQTLs
  (Regensburg retina eQTL mega-analysis, Strunz et al. 2020) with glaucoma, so step 72 is no longer run or bundled,
  and its input, the only one that had to be downloaded by hand, is no longer needed: step 00 now downloads every
  input. Supplementary Tables v2.15 drop the retina row of S9 Panel B (the other cells are identical to v2.14), and
  `reference/` holds v2.15. No other output changes.

## Changes for release v9 (24 October 2026)

- **New step 55c** (`55c_audit_gwas_inputs.py`, section 7b): input audit of the TWAS, conditional TWAS and cis-MR
  inputs of every GWAS (FECD, glaucoma, IOP, refractive error); see *Notes for reviewers*. It only reads the
  pipeline's inputs and outputs and writes `analysis/input_audit/` (and `results/input_audit_sensitivity/`), so no
  reported result changes. The master script runs it with `METAXCAN_PY` (`PY=... RUN`), as it uses numpy and scipy.
- **Supplementary Tables v2.14** (`reference/`): the S14 legend no longer states the DESeq2 subgroup model (batch +
  CTG18.1 group); every value is identical to v2.13, and step 140 finds 0 mismatches.
- `analysis_outputs/` now includes the input-audit tables of steps 55, 55b and 55c.
- **Project copies of the glaucoma and IOP S-PrediXcan inputs rebuilt (steps 46, 51).** Running the audit in the
  project and in the v8 fresh-folder copy showed that the project's two files dated from 9–10 June, when the
  rsID-to-model-variant map was built from a different set of model files. They differed from the current steps'
  output only at model variants that S-PrediXcan does not use (glaucoma: 10 rows written under a variant ID absent
  from the current models instead of an allele-inconsistent one; IOP: 3,674 rows at whole-blood and brain-cortex
  variants, which the IOP analysis does not use), and the S-PrediXcan outputs from the two versions are
  byte-identical (v8 fresh-folder test). The June files are kept as `*.bak_june_20261024`.
- **Check of step 55c.** Run in the project and in the v8 fresh-folder copy (all inputs regenerated from an empty
  folder), the seven audit tables and the sensitivity FUSION outputs are byte-identical. Step 140 with the v2.14
  workbook: 0 mismatches.

## Changes for Final Draft v8 (21 October 2026)

- **S1, S2, S11 and Figure 2 regenerated with this release.** Steps 17 → 95a → 95 → 19b → 95b were rerun with the
  released code. Only P2RX5 and TXN2 change (for example 3-tissue S-MultiXcan P 0.19 → 0.17 and 0.28 → 0.13; both genes
  gain a model SNP in some tissues), with the BH FDR values, ranks and per-tissue gene counts that depend on them.
  NR1H3 and ZNF408 values and ranks (17th and 21st of 17,285; 34th and 52nd of 22,222), every gene count and
  Bonferroni threshold in the manuscript, and Tables 1–2 are identical. Previous outputs:
  `results/bak_pre_s1s2s11_rerun_20261021/`.
- **S3 regenerated (step 82 rerun).** The fresh-folder test of this release (below) reproduced every table except
  S3. The project's European-GWAS FUSION input (`data/gwas/fecd_eur_fusion.tsv`, built 28 June) stored Z to
  6 significant figures and kept 54,023 duplicate rsID rows (multi-allelic sites; FUSION uses the first, which is
  identical in both builds), whereas step 82 writes full-precision Z with one row per rsID. FUSION TWAS Z-scores
  therefore differed in the 5th–6th significant figure (largest relative difference 0.5%, for a gene with Z ≈ 0.002),
  which moved the last printed digit of a few hundred S3 cells. Step 82 and the steps that use its outputs (69b, 24b,
  112a, 112, 112b, 113) were rerun in the project. The Bonferroni-significant genes in each tissue, NR1H3 (fibroblast
  Z −3.74, P 1.85 × 10⁻⁴; skin Z −3.53, P 4.21 × 10⁻⁴), ZNF408 (oesophagus Z 3.70), S5, S13, Figure 1 data and every
  number in the text are unchanged. Previous outputs: `results/bak_pre_fusion_eur_rerun_20261021/`.
- **New step 140** (`140_check_reported_tables.py`): writes S1, S2, S3, S10, S11 and the S4 power column exactly as
  reported (row order, rounding, derived Bonferroni/FDR columns) and checks every value in manuscript Tables 1–3 and
  Supplementary Tables S1–S14 against the run (exact for the exported tables; otherwise the reported value must be
  the output or a rounding of it to ≥ 2 significant figures). S3 now gives gene positions (P0, P1) exactly; the
  workbook had rounded them to 4 significant figures.
- **S4 power (step 84c).** Power is now computed for each gene and tissue from its own discovery Z, read from the
  S-PrediXcan outputs; earlier versions applied NR1H3's values to the ZNF408 rows. ZNF408: 24% (H18.4) and 51%
  (H18.5); NR1H3 unchanged (17%/21% and 36%/45%); NR1H3 oesophagus 6% and 11%. Written to
  `analysis/nr1h3_finngen_lookup_power.tsv`.
- **Table 1 (step 19c).** GWAS P-values below 1e-300 (stored at the double-precision floor, 1.1e-308) are written as
  "<1e-300" instead of a spurious number.
- **Step 110** writes its sample sheet (`analysis/nr1h3_endothelium_samples.tsv`, S14).
- **Supplementary Tables v2.12** built from these outputs (S1, S2, S10, S11, S4; S9 Panel A FECD fibroblast TWAS P
  1.8e-4, as FUSION prints 1.85e-04 and Figure 5 shows; S3, S4, S9 and S14 caption notes), and **v2.13** with S3
  from the step-140 export.
- **Reproducibility fixes found by the fresh-folder test:**
  - FOCUS needs Python 3.11 (pandas 1.5.3 has no wheels for later versions). `FOCUS_PY` is now documented in
    `config.sh` and `config.local.sh.example`, and `00_install_focus.sh` stops with a message for any other version
    (the first test stopped there because the default `python3` was 3.12).
  - `config.sh` no longer defines the SMR binary.
  - Fixed seeds for the ggrepel labels (69b, 49, 54, 59) and fixed TwoSampleMR IDs (83, 91); see *Notes for
    reviewers*. Before this, the PDF versions of the volcano plots and the ID columns of two MR tables differed
    between runs (label positions and random strings only).
  - No time limit for the ggrepel labels (`max.time = Inf`, same steps). The final test found that the
    refractive-error volcano plot (step 59, not in the manuscript) still placed two labels differently, because
    ggrepel stops after 0.5 s by default. Figures 1 and 2 and the other volcano plots are pixel-identical with and
    without the limit.

## Reproducibility check (21 October 2026, v8 release)

This release was run twice from an empty project folder on an Apple M1 laptop with 8 GB RAM, using only the
bundled scripts, `data/derived/` and the public inputs. To save disk and time, the large raw downloads (GWAS files,
MASHR models and covariances, LD reference, FUSION weights, eQTLGen, retina eQTLs, liftover chains) were linked from the project's
earlier download; every other file was regenerated by `run_full_pipeline.sh`. Step 140 then checked every reported
table, and every output file was compared with the file behind the manuscript.

- **First run** (bundle before the fixes above): FOCUS stopped under Python 3.12 (now `FOCUS_PY`, Python 3.11). After
  that was fixed, every table matched except S3 (the 28 June FUSION input, see *S3 regenerated*), and the MR
  ID columns and volcano-plot PDFs differed (random IDs and label positions).
- **Final run** (this release): completed with no failed step in 2 h 38 min of compute (longest steps: 11, GWAS
  preparation, 39 min; 82, FUSION, 35 min). Step 140: manuscript Tables 1–3 and Supplementary Tables S1–S14
  (25 tables or panels; S1, S2, S3, S10 and S11 cell by cell) match the outputs with **0 mismatches**.
- **Comparison with the project's outputs** (498 files in `analysis/` and `results/`): 479 are byte-identical,
  including every table, every MR, colocalisation, fine-mapping and FUSION result, and every PNG figure. The other
  19 differ only in metadata: 10 PDF figures (creation date; the rendered pages are pixel-identical after the
  `max.time` fix), 4 FUSION `.report` files (absolute paths; identical once the project folder is replaced), 1 gzip
  file (header timestamp; identical content) and 4 log files.
- Two inputs are not public (see *Data sources*); the run used the derived tables in `data/derived/` in their place,
  as a user of this bundle would.

## Changes for Final Draft v7 (21 October 2026)

- **Removed from the pipeline because the manuscript no longer reports them:** SMR/HEIDI (steps 90–93, with their
  manual downloads: the SMR binary, the BESD eQTL panels and the merged LD panel step 00 built for HEIDI) and the
  GTEx-instrumented MR for FECD (steps 31, 32) and glaucoma (step 78). The primary eQTLGen-instrumented cis-MR is
  unchanged. The scripts remain in the project but are not part of this release; earlier sections of this README
  that mention them describe earlier releases.
- **Supplementary Tables v2.11.** S9: glaucoma S-MultiXcan without brain cortex P = 8.4 × 10⁻⁶ (the output is
  8.45 × 10⁻⁶; the workbook had rounded it to 8.5 × 10⁻⁶); Panel D "Refractive error" label.

## Changes for Final Draft v6 (21 October 2026)

- **Figures redrawn from unchanged outputs** (no reported number changes):
  - Figure 2 (step 69b): ZNF408 highlighted with NR1H3; x-axis labelled "Mean per-tissue Z-score (S-PrediXcan, 3 tissues)".
  - Figure 5 (76e): signed FUSION TWAS Z-scores; the Z ± 1.96 bars (no inferential target) are removed; axis "FUSION TWAS
    Z-score"; "Refractive error" (positive = more hyperopic refraction) replaces "myopia".
  - Figure 6 (89b): "Refractive error" label and sign note; `analysis/nr1h3_cismr_standardised_v2.tsv` is
    byte-identical before and after (md5 377a1b88…).
  - Figure 7 (70b): the glaucoma cis-MR OR and P are now read from the step-89c Table 3 output; the hard-coded
    P = 1.5 × 10⁻⁵ predated the per-variant exposure-N fix and is now 2.0 × 10⁻⁵. The result text fits its box,
    "(protective)" is removed and the right-hand column is labelled as the glaucoma analysis it shows.
- **Supplementary Tables v2.10.** S9: I²GX defined as the MR-Egger regression-dilution (NOME) diagnostic;
  "Refractive error" label. S12: exact coordinates restored from `analysis/fecd_reported_loci.tsv` and
  `analysis/fecd_twas_gene_locus_class.tsv` (the earlier workbook export had kept four significant figures, e.g.
  NR1H3 TSS 47,250,000 instead of 47,248,300). S13: FOCUS gene sums labelled as descriptive sums of model-level PIPs.
  New S14: corneal endothelial RNA-seq samples (step 110 sample sheet and the GEO series matrices). The underlying
  analysis files are unchanged.
- **README.** Checksum wording (verified where recorded and supplied); P2RX5/TXN2 table references corrected; provenance
  section added.

## Changes since the 13 October 2026 (v3.10) bundle

- **Refractive-error GWAS merge (step 55) rebuilt on variant identity.** The original Hysi et al. METAL file
  (Z-scores; GRCh37 `MarkerName`) had been linked to the GWAS Catalog harmonised file (rsID, GRCh38) on rounded
  p-value and allele frequency, keeping the first harmonised record for a shared key. Step 55 now lifts each GRCh37
  `MarkerName` to GRCh38 (UCSC hg19ToHg38 chain, as step 50) and matches chromosome, position and the submitted
  alleles; p-value and frequency are checked for agreement only. Z is aligned explicitly to the harmonised effect
  allele, and missing, duplicated and allele-inconsistent records are rejected and counted
  (`analysis/myopia_gwas_merge_audit.tsv`): 8,675,333 of 8,754,054 records kept; every kept record's orientation agrees
  with the GWAS Catalog harmoniser. `FREQ` is now the A1 frequency and a per-variant `N` (METAL Weight) is carried.
- **What the old merge got wrong** (`analysis/myopia_gwas_merge_old_vs_new.tsv`): 714,794 variants were missing and
  249,068 rsIDs carried another variant's record (141,489 with a Z of the opposite sign). New step 55b checks the merge
  independently of liftOver and the GWAS Catalog (1000G EUR GRCh37 rsID/position → original record): 1,179,368 of
  1,179,368 comparable variants agree (old merge: 4,537 sign errors, 8,240 other mismatches, 3.9% missing).
- **Reruns.** 60 → 76d → 76e/76g (Figure 5, S8, S9 Panel A); 89c, 89b (Table 3, Figure 6, S9 Panel C);
  56 → 57 → 58 → 59 and 76 (refractive-error S-PrediXcan/S-MultiXcan, S9 Panel D). Glaucoma and IOP inputs and
  results are unchanged. 89b/89c take the refractive-error N from the merged file instead of four hard-coded values
  (identical for the four instruments). New step 76g builds the S8 table from the step-76d outputs.
- **Effect on reported refractive-error results.** Table 3 and Figure 6 are identical (same four instruments; the 221
  eQTLGen variants newly available as candidates do not change clumping). Whole-blood FUSION for NR1H3 (Figure 5, S8,
  S9 Panel A): TWAS Z 5.95 → 6.11 (P 1.0 × 10⁻⁹); conditional Z 1.51 → 3.54 (P 0.13 → 4.0 × 10⁻⁴), still not retained
  in the joint model (FUSION threshold |Z| > 3.87 for 463 chr11 models); at 11p11.2 the retained gene changes from
  SLC39A13 to C1QTNF4. In the conditional region (chr11:45.09–48.59 Mb, GRCh37; 17 models) the old input lacked
  77 of 1,117 reference SNPs, which FUSION imputed (`analysis/myopia_fusion_wb_11p11.2_input_audit.tsv`); one of
  NR1H3's three weighted SNPs (rs3824867) was among them. This file and `analysis/myopia_gwas_merge_old_vs_new*.tsv`
  were one-off comparisons with the pre-fix input (19 October) and are not produced by the pipeline; step 55c
  (release v9) produces the corresponding checks for the current inputs. NR1H3 S-MultiXcan P 0.94 → 0.93 (S9 Panel D).
- Backups of every replaced file: `results/bak_myopia_pvaleaf_merge_20261019/` and `*.bak_pvaleaf_20261019` in
  `data/gwas/`.
- **Exposure sample size in the cross-trait cis-MR (Table 3, Figure 6).** Steps 89b/89c converted the eQTLGen Z of
  each instrument to a standardised effect with a fixed N = 31,684, whereas the FECD cis-MR (83, 88, 91) uses each
  variant's `NrSamples`. The four instruments have NrSamples 31,676 (rs7395581), 31,684 (rs11038863), 31,346
  (rs3740694) and 21,502 (rs4647760; 21 of 37 cohorts), so the fixed N understated the exposure effect and SE of
  rs4647760 by a factor 1/1.214 (others ≤ 0.5%). 89b/89c now use the per-variant N, matching the FECD analysis;
  per-instrument values under both conventions are in `analysis/nr1h3_cismr_instruments_exposureN.tsv`. IVW estimates
  change by ≤ 0.6% (glaucoma log-OR −0.182 → −0.181, P 1.5 → 2.0 × 10⁻⁵; IOP −0.205 → −0.204, P 0.064 → 0.065;
  refractive error 0.0429 → 0.0431, P 1.1 × 10⁻⁶ → 1.3 × 10⁻⁷); I², F (597, already per-variant N) and all conclusions
  are unchanged. Previous outputs: `results/bak_exposureN_fixed_20261020/`.

## Changes since the 13 October 2026 (v3.9) bundle

- **Colocalisation on full cis summary statistics.** Steps 101c → 102c → 81c (Figure 3, S6) and 79b (glaucoma, S9
  Panel D) replace 83b, 89, 81b and 79. Gene lists are bundled in `data/tophits21/`.
- **Fine-mapping at 11p11.2 (new S13).** Steps 112a, 112, 112b and 113 (FOCUS, with `00_install_focus.sh`).
- **Figure 7** (step 70b) text updated for full-cis colocalisation and fine-mapping.
- 101c and 102c read the project folder and `tabix` from `config.sh` (they had hard-coded paths).
- **101c downloads are now checked.** A failed or interrupted `tabix` request used to be saved as an empty ("no data")
  or truncated file; in a clean-folder run on a congested connection this affected over half of the files. 101c now fetches
  each tissue's index once before the parallel downloads, retries each request five times, and stops with an error
  rather than saving a failed download. An empty file now always means the gene was not tested in that tissue
  (ENSG00000279672 in all three tissues; LINC01415 and SHISA9 in oesophagus mucosa).
- **PANO1 in 101c.** PANO1 is ENSG00000274897 in GENCODE v26 but ENSG00000288675 in the eQTL Catalogue; its S6 data had
  been fetched by a one-off step (101b) outside the master script, so a clean run lost PANO1 from S6. 101c now maps the ID.
- 112a creates `analysis/` if it does not yet exist. 112b downloads with a timeout and retries, via a temporary file.
- **Clean-folder test of the new and changed steps (13 October 2026).** Steps 00_install_focus, 100, 101c, 102c, 81c, 112a,
  112, 112b, 113, 78 and 79b were run from the unzipped bundle in an empty folder (raw GWAS files and MASHR models linked,
  as in the earlier test). Every output is byte-identical to the files behind the manuscript: S6 and Figure 3, all 81
  eQTL Catalogue downloads, the SuSiE/ABF, eQTL Catalogue SuSiE and FOCUS results (S13) and the glaucoma colocalisation;
  step 78 gives identical values (only TwoSampleMR's random ID columns differ). The test found and fixed the 101c, 112a,
  112b and PANO1 problems listed above. Compute time about 50 minutes, mostly the 101c downloads (20 min) and step 112 (15 min).

## Reproducibility check (13 October 2026, v3.9 bundle)

The bundle was run from an empty project folder on a MacBook with 8 GB RAM, using only the bundled
scripts, `data/derived/` and the public inputs. To save disk and time, the large raw downloads
(GWAS files, MASHR models, LD reference, BESD panels) were linked from an earlier download by step 00
rather than fetched again. Every other file was regenerated. The outputs were then compared with the files
behind the manuscript:

- **Identical** (byte-for-byte, or the same values for every key): Tables 1–3; S4 and S6–S9, S12; Figure 7 and the
  data behind Figures 3, 4 and 6; the FECD and cross-trait cis-MR (steps 83, 89b, 89c, 91, 94); the GTEx-instrumented
  MR (steps 31, 32 and 78); the colocalisation analyses (steps 72, 79, 83b, 89 and 108); every SMR result (steps
  90–93), including the reported FECD × eQTLGen estimate; and the harmonised GWAS and region files.
- **Equal within rounding:** the FUSION results: S3, S5, Figure 1 and the FUSION values in Figure 5. The largest
  difference in Z is 1 × 10⁻⁴, which is invisible at the reported precision.
- **Two genes differ:** P2RX5 and TXN2, in the genome-wide S-PrediXcan and S-MultiXcan tables S1, S2 and S11 and
  as two points in Figure 2 (S10 lists only NR1H3 and ZNF408 and is not affected; earlier versions of this README
  listed it in error). Step 17 no longer applies an rsID-to-varID remap, which mislabelled two multi-allelic
  indels in these genes' models and was not part of the intended analysis. Both genes stay far from significance
  (|Z| < 2 in every tissue; S-MultiXcan P > 0.12), and no number in the text changes.
- Compute time about 5.5 hours, plus downloads. The longest steps are 82 (FUSION, 1.6 h), 95b (49-tissue
  S-MultiXcan, 1.2 h), 92 (SMR, 50 min), 95 (35 min) and 83b (colocalisation, 30 min).

## Changes since the 12 October 2026 bundle

- **One settings file.** `scripts/config.sh` (+ optional `config.local.sh`) sets the project folder and every tool
  path; no script hard-codes `~/FECD_GWAS` any more.
- **Inputs and software.** `00_download_inputs.sh` downloads the public inputs and the pinned
  MetaXcan, FUSION and PLINK versions, verifies checksums where they are recorded and supplied to the download
  workflow, and lists the manual downloads. `00_install_R_packages.R` and
  `environment_metaxcan.yml` pin package versions.
- **Private inputs.** Steps 109, 109b and 76f fall back to `data/derived/` when the original file is absent.
- **Master script.** It stops at the first failed step. Step 00 was added, and step 100 is now called with its
  arguments. Step 109 now runs before 19c. Steps 24, 27, 40 and 61–64 were removed (no reported result).
- **Fixes found by the clean-folder run:**
  - 11: the run no longer fails out of memory on 8 GB.
  - 82: builds its own per-chromosome GWAS files.
  - 83a: writes the cis-MR region files (previously made by hand).
  - 83b: its gene list is bundled (`data/fecd_coloc_gene_set.tsv`), and it creates its output folder.
  - 84b: tolerates the pipe closing early when reading the header.
  - 24b: runs both minp settings itself.
  - 31/32: fixed bootstrap seed.
  - 78: weak instruments are flagged, not dropped.
  - 72: stops if the retina file is missing.
  - 49/54/59: dead code removed.
  - GTEx Portal requests are retried and fail loudly.
- **Table numbering.** Manuscript v3.9 swapped Tables 1 and 2 so they are cited in order; step 19c's labels follow the
  manuscript, and its output files are unchanged.
- **SMR.**
  - 90 keeps only single-rsID rows for glaucoma; SMR cannot parse FinnGen's `rsA,rsB` IDs, and the original
    file had been built this way.
  - 92 now exits with an error if any run fails; it previously warned and carried on.

## Changes since the 30 September 2026 bundle

- Pipeline re-targeted from manuscript v2.5 to v3.2. New scripts: 19c, 24b, 81b, 83a, 84b, 84c, 89c, 97b, 109a, 109b, 76f
  (most turn analyses previously run as one-off commands into scripts; each was checked against the reported values).
- Scripts added that v3.2 uses but the old bundle lacked: 19b, 95a, 95b, 31, 32, 75b, 76e, 89b, 108, 109, 110, SMR 90–93.
- Superseded and removed: 18/19 (un-remapped S-PrediXcan / S-MultiXcan), 69, 70, 76c, 81, 84, 85, 85b, 89 (non-per-SNP N),
  86, 88_nr1h3_cismr_forest_update, 71, 21, 23, 36, 65–68, 80, 90_coloc_eur, 92_coloc_nr1h3_3tissue_eur, 96, 97.
- 42 scripts that still pointed at the pre-rename directory `~/bk_analysis` (or at a FUSION copy that no longer exists)
  now use `~/FECD_GWAS` and `~/Downloads/fusion_twas-master`.
- Table 3 is now generated by a script (89c) with a fixed bootstrap seed and per-variant N for refractive error, so it
  agrees with Figure 6. The resulting small changes to Table 3 and the Results are tracked in manuscript v3.2.
