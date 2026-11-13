#!/bin/bash
# ============================================================
# NR1H3 / FECD manuscript pipeline — master script (release v14.9, 13 November 2026)
#
# Reproduces every result in:
#   "Integrative transcriptome-wide association, colocalisation, fine-mapping and Mendelian randomisation
#    analyses of Fuchs endothelial corneal dystrophy highlight 11p11.2 (ZNF408, NR1H3) as a suggestive locus"
#
# Maps to manuscript elements:
#   Fig 1   FUSION Manhattan (fibroblast, EUR GWAS)        [82 -> 69b]
#   Fig 2   S-MultiXcan volcano                            [95a, 95 -> 19b -> 69b]
#   Table 1 S-MultiXcan ranking; Table 2 per-tissue NR1H3   [19c]
#   Fig 3 / S6  colocalisation, full cis eQTL sumstats      [101c -> 102c -> 81c]
#   S6 S_MultiXcan_rank (European S-MultiXcan ranking)     [19d -> 19e -> 19f]
#   Fig 4 / S7  FECD cis-MR (EUR primary)                   [91, 83 -> 87 EUR; TSPAN10 rows of S7: 83a, 91, 83 _tspan10]
#   Fig 5 / S8  cross-trait FUSION whole-blood conditional  [76d -> 76e, 76g]
#   Glaucoma FUSION input matched to the LD-reference alleles  [60 -> 60b]
#   Genome-wide glaucoma rank of NR1H3 (whole blood)        [76f June 2025 run; 76h rerun, allele-matched + first-record]
#   Table 3 / S9C cross-trait cis-MR                        [83a eQTLGen allele frequencies -> 89c]
#   Fig 6   cross-trait cis-MR, each outcome on its own scale [89c -> 89b]
#   Fig 7   study flow diagram                              [70b]
#   S4 FinnGen look-up [84b, 84c]; S5 [24b]; S10-S11 [95, 95b]; S12 [109a -> 109]
#   S9 D coloc glaucoma [79b]; S13 fine-mapping 11p11.2 [112a -> 112, 112b, 113]; S14 RNA-seq samples [110]
#   Input audit of the TWAS, conditional TWAS and cis-MR inputs of every GWAS   [55, 55b, 55c]
#   Additional outputs requested in review (29 Oct 2026; analysis/additional_outputs/)   [139]
#   Figures 1-7 as journal files (figures/FigureN.pdf, .svg, .tif)   [139b]
#   Every value in Tables 1-3 and S1-S14 is checked against this run's outputs at the end   [140]
#
# All logs: logs/<script>.log  (+ master log logs/full_pipeline.log)
# ============================================================
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE etc.; machine-specific values in config.local.sh
BASE="$FECD_BASE"
cd "${BASE}"
mkdir -p "${BASE}/logs"
exec > >(tee -a "${BASE}/logs/full_pipeline.log") 2>&1

# RUN <script> "<note>" [args...] — interpreter from the extension, one log per script
RUN () {
  local s="$1" note="${2:-}"; shift 2 || shift $#
  local interp="bash"
  [[ "$s" == *.R  ]] && interp="Rscript"
  [[ "$s" == *.py ]] && interp="${PY:-python3}"     # PY=<interpreter> RUN ... for a step that needs numpy/scipy
  local log="${BASE}/logs/${s%.*}.log"
  echo ""; echo ">>> ${s}   ${note}"
  ${interp} "${SCRIPTS_DIR}/${s}" "$@" > "${log}" 2>&1
  # a step can exit 0 even though an embedded R or Python call failed: stop on either
  if grep -q -E "^Execution halted|^Traceback \(most recent call last\)" "${log}"; then
    echo "    FAILED: uncaught R/Python error in ${log}"; tail -5 "${log}"; exit 1
  fi
  echo "    done $(date)   [log: ${log}]"
}

echo "============================================================"
echo " NR1H3/FECD pipeline (release v14.9) started: $(date)"
echo "============================================================"

# ---- 0. INPUTS ----------------------------------------------
RUN 00_download_inputs.sh             "public inputs and pinned tools (checksums verified where recorded)"
[ -x "$FECD_BASE/software/focus_venv/bin/focus" ] || \
  RUN 00_install_focus.sh             "FOCUS 0.802 in its own venv, pinned (step 113)" "$FOCUS_PY"
RUN 11_prep_fecd_gwas.R               "harmonise FECD GWAS (Verma/MVP, GCST90480062)"
RUN 15_prep_fecd_fusion.sh            "multi-ancestry FECD -> FUSION format (step 24)"
RUN 17_prep_fecd_sprediXcan.sh        "FECD -> S-PrediXcan format"
RUN 95a_remap_fecd_gwas_alltissue.py  "remap GWAS to all 49 MASHR model varIDs"
RUN 39_download_eur_weights.sh        "GTEx v8 EUR FUSION weights"
RUN 45_download_glaucoma_gwas.sh      "FinnGen R12 H7_GLAUCOMA"
RUN 50_download_iop_gwas.sh           "IOP (GCST009413)"
RUN 55_download_myopia_gwas.sh        "refractive error (Hysi 2020)"
RUN 55b_verify_myopia_gwas.py         "refractive error: independent 1000G check + MR instrument records"
RUN 60_prep_fusion_gwas.sh            "cross-trait GWAS -> FUSION format"
RUN 60b_glaucoma_fusion_allele_match.py "glaucoma FUSION input matched to the LD-reference alleles (10 Nov 2026)"
RUN 83a_extract_eqtlgen_nr1h3.sh      "eQTLGen NR1H3 cis-eQTLs"
RUN 83a_extract_eqtlgen_tspan10.sh    "eQTLGen TSPAN10 cis-eQTLs (S7, release v14.8)"
RUN 100_fetch_eqtl_catalogue_region.sh "full cis GTEx sumstats, NR1H3 region (step 108)" ENSG00000025434 11 46248300 48248400
RUN 109a_fetch_gwascat_fecd.py        "GWAS Catalog FECD associations (step 109)"

# ---- 1. FECD TWAS -------------------------------------------
echo ">>> 95_alltissue_sprediXcan_fecd.sh   S-PrediXcan, 49 tissues, remapped GWAS (S1, S10)"
GWAS_FILE=data/gwas/fecd_alltissue_sprediXcan.txt.gz OUT_DIR=results/sprediXcan_alltissue_remap TAG=fecd \
  bash "${SCRIPTS_DIR}/95_alltissue_sprediXcan_fecd.sh" > logs/95_alltissue_sprediXcan_fecd.log 2>&1
RUN 19b_rerun_fecd_smultixcan_remapped.sh "3-tissue S-MultiXcan (S2)"
RUN 19d_prep_fecd_eur_sprediXcan.py  "European FECD GWAS -> S-PrediXcan format, matched by position + alleles (95a)"
RUN 19e_fecd_eur_smultixcan.sh       "European 3-tissue S-PrediXcan + S-MultiXcan (S6 ranks)"
RUN 19f_fecd_eur_smultixcan_rank.py  "European ranking -> S6 gene list ranks (data/tophits21/genes_coloc.tsv)"
RUN 95b_alltissue_smultixcan_fecd.sh  "49-tissue S-MultiXcan (S11)"
RUN 82_fecd_eur_matched_twas.sh       "FUSION, European GWAS (S3, Fig 1)"
RUN 109_fecd_reported_loci.py         "reported FECD loci + gene classification (S12)"   # before 19c: Table 2 uses the locus classification
RUN 19c_tables1_2_twas_summary.py     "Tables 1-2, per-tissue FDR"
RUN 69b_fig1_fig2_v2.R                "Figures 1-2"
RUN 24b_fusion_11p11_conditional.sh   "conditional within 11p11.2, European GWAS (S5)"
RUN 109b_insawang_11p11_lookup.py     "Insawang et al. 11p11.2 look-up"
RUN 97b_egene_status_nr1h3_znf408.py  "GTEx eGene status"
RUN 110_nr1h3_endothelium_fecd_vs_control.R "corneal endothelium FECD vs control"

# ---- 2. FINNGEN LOOK-UP (S4) --------------------------------
RUN 84b_finngen_replication_H184_H185.sh "NR1H3/ZNF408 in FinnGen H18.4 and H18.5"
RUN 84c_finngen_replication_power.py  "look-up power, per gene and tissue (S4)"

# ---- 3. COLOCALISATION (Fig 3, S6) --------------------------
# full GTEx v8 cis summary statistics (eQTL Catalogue r6); replaces the significant-pairs-only input of 83b/89/81b
export GENES=data/tophits21/genes_coloc.tsv
RUN 101c_fetch_coloc_allpairs.sh      "eQTL Catalogue all variant-gene pairs, S6 genes x 3 tissues"
unset GENES
RUN 102c_coloc_abf_allgenes.R         "coloc.abf, European GWAS (primary) + multi-ancestry"
RUN 81c_fecd_coloc_heatmap_fullcis.R  "Figure 3 + S6"

# ---- 3b. FINE-MAPPING AT 11p11.2 (S13) ----------------------
RUN 112a_ld_1000g_11p11.sh            "1000G GRCh38 European LD panel + European GWAS region"
RUN 112_finemap_11p11.R               "SuSiE-RSS + ABF credible sets"
RUN 112b_eqtlcat_susie_11p11.py       "eQTL Catalogue SuSiE eQTL credible sets"
RUN 113_focus_11p11.sh                "FOCUS gene-level fine-mapping (+ default-block sensitivity)"

# ---- 4. CAUSAL INFERENCE, FECD (Fig 4, S7) ------------------
RUN 91_mr_nr1h3_eur.R                 "primary cis-MR, European GWAS"
RUN 83_mr_nr1h3_eqtlgen_fecd.R        "cis-MR, multi-ancestry sensitivity"
RUN 91_mr_tspan10_eur.R               "TSPAN10 cis-MR, European GWAS (S7)"
RUN 83_mr_tspan10_eqtlgen_fecd.R      "TSPAN10 cis-MR, multi-ancestry sensitivity (S7)"
RUN 94_mr_ancestry_sidebyside.R       "ancestry comparison"
RUN 88_verify_nr1h3_fecd_mr.R         "independent verification"
RUN 87_nr1h3_fecd_mr_forest.R         "Figure 4" EUR

# ---- 5. CROSS-TRAIT TWAS (Fig 5, S8, S9 A-B, D) -------------
RUN 46_prep_glaucoma_sprediXcan.sh ; RUN 47_run_glaucoma_sprediXcan.sh
RUN 48_run_glaucoma_smultixcan.sh  ; RUN 49_analyse_glaucoma.R "glaucoma S-PrediXcan/S-MultiXcan"
RUN 51_prep_iop_sprediXcan.sh      ; RUN 52_run_iop_sprediXcan.sh
RUN 53_run_iop_smultixcan.sh       ; RUN 54_analyse_iop.R "IOP S-PrediXcan/S-MultiXcan"
RUN 56_prep_myopia_sprediXcan.sh   ; RUN 57_run_myopia_sprediXcan.sh
RUN 58_run_myopia_smultixcan.sh    ; RUN 59_analyse_myopia.R "refractive error S-PrediXcan/S-MultiXcan"
RUN 75_glaucoma_brain_panel.sh        "glaucoma: fibroblast + brain cortex + skin"
RUN 75b_glaucoma_smultixcan_nobrain.sh "glaucoma: without brain cortex"
RUN 77_glaucoma_wholeblood_panel.sh   "glaucoma: whole blood"
# genome-wide whole-blood FUSION + conditional/joint (scripts/fusion_genomewide.sh, the procedure step 82 uses for
# FECD): glaucoma, both inputs (helper scripts/76h_glaucoma_wb_genomewide_rank.py), then IOP and refractive error
# (helper scripts/76i_crosstrait_wb_genomewide_summary.py); step 76d takes their chromosome-11 results for S8
RUN 76h_glaucoma_wb_genomewide.sh     "genome-wide whole-blood FUSION + conditional, glaucoma (NR1H3 rank)"
RUN 76i_crosstrait_wb_genomewide.sh   "genome-wide whole-blood FUSION + conditional, IOP and refractive error"
RUN 76d_fusion_wb_conditional.sh      "chr11 whole-blood FUSION results from the genome-wide runs (S8)"
RUN 76g_s8_crosstrait_fusion_wb_table.R "S8 table"
RUN 76e_fig5_glaucoma_label.R         "Figure 5"
RUN 76f_glaucoma_wb_conditional_rank.py "genome-wide glaucoma rank (external 2025 results file)"

# ---- 6. CROSS-TRAIT CAUSAL INFERENCE (Table 3, Fig 6) -------
RUN 89c_table3_crosstrait_cismr.R     "Table 3 (+ S9 Panel C)"
RUN 89b_nr1h3_cismr_standardised_perSNP_N.R "Figure 6"
RUN 79b_coloc_nr1h3_glaucoma_fullcis.R "NR1H3 eQTL x glaucoma coloc, full cis (S9 Panel D)"
RUN 108_coloc_nr1h3_glaucoma_fullcis.R "FECD x glaucoma GWAS coloc at 11p11.2"

# ---- 7. FIGURE 7 --------------------------------------------
RUN 70b_flowchart_NR1H3_v2.R          "Figure 7"

# ---- 7b. INPUT AUDIT ----------------------------------------
# links every TWAS, conditional TWAS and cis-MR input (FECD, glaucoma, IOP, refractive error) back to its source GWAS
# records and to the variants it is matched against; re-derives model coverage, the 11p11.2 conditional regions and
# the cis-MR instruments (analysis/input_audit/). Uses numpy and scipy, so it runs in the MetaXcan Python.
PY="$METAXCAN_PY" RUN 55c_audit_gwas_inputs.py "input audit, all GWAS (analysis/input_audit/)"

# ---- 7c. ADDITIONAL OUTPUTS (peer review, 29 Oct 2026) -------
# collects the outputs requested in review with the tables derived for the request (analysis/additional_outputs/):
# glaucoma FUSION with the allele-matched input vs the first-record input of releases up to v13 (chr11 and 76h), the European S-MultiXcan
# behind S6 (19d-19f) and the basis for the RNA-seq batch assignment (110; fetches SRA run information once)
RUN 139_additional_outputs.py         "additional outputs requested in review"

# ---- 7d. FIGURE FILES -----------------------------------------
# Figures 1-7 re-plotted into figures/ as vector PDF (via headless Chrome), SVG and 600-dpi TIFF
RUN 139b_export_figures.sh            "Figures 1-7 as separate files (figures/)"

# ---- 8. CHECK THE REPORTED TABLES ---------------------------
# writes the S1, S2, S10, S11 and S4 tables (analysis/supplementary/) and compares every value in manuscript
# Tables 1-3 and Supplementary Tables S1-S14 (copies in reference/) with this run's outputs; fails on any mismatch
RUN 140_check_reported_tables.py      "Tables 1-3 and S1-S14 against this run"

echo ""
echo "============================================================"
echo " NR1H3/FECD pipeline (release v14.9) finished: $(date)"
echo "============================================================"
