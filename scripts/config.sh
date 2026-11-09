#!/bin/bash
# ============================================================
# Settings for the NR1H3/FECD pipeline. Every script sources this file.
# Edit the defaults below, export the variables before running, or put
# machine-specific values in scripts/config.local.sh (read first if present).
#
#   FECD_BASE    project folder holding data/, results/, analysis/, logs/, software/
#   METAXCAN_PY  Python interpreter with the MetaXcan dependencies
#   FOCUS_PY     Python 3.11 interpreter used to build the FOCUS venv (step 00_install_focus; pandas 1.5.3,
#                which FOCUS 0.802 needs, has no wheels for Python 3.12 or later)
#   TABIX        tabix binary (htslib), used for FinnGen region queries
#   H185_FILE    Insawang et al. 2026 meta-analysis summary statistics (not public)
#   GLAUCOMA_COND_XLSX  genome-wide FUSION whole-blood glaucoma conditional results (12 Jun 2025)
# External tools are expected under $FECD_BASE/software/ (see README):
#   MetaXcan/ (hakyimlab/MetaXcan v0.8.1), fusion_twas/ (gusevlab/fusion_twas e1ba5f7
#   + fusion_patch), plink/plink (v1.90b7.2).
# ============================================================
_cfg_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ -f "$_cfg_dir/config.local.sh" ] && source "$_cfg_dir/config.local.sh"

: "${FECD_BASE:=$HOME/FECD_GWAS}"
: "${METAXCAN_PY:=python3}"
: "${FOCUS_PY:=python3.11}"
: "${TABIX:=tabix}"
: "${H185_FILE:=$FECD_BASE/data/private/H185_MVP_FinnGen_UKB+AoU_SNPsinatleast2studies.gz}"
: "${GLAUCOMA_COND_XLSX:=$FECD_BASE/data/private/Eye Conditional results 12_6_25 FinnGEN GTEV8 whole blood.xlsx}"
FUSION_DIR="$FECD_BASE/software/fusion_twas"
METAXCAN_DIR="$FECD_BASE/software/MetaXcan/software"
PLINK="$FECD_BASE/software/plink/plink"
SCRIPTS_DIR="$_cfg_dir"
export FECD_BASE METAXCAN_PY FOCUS_PY TABIX H185_FILE GLAUCOMA_COND_XLSX FUSION_DIR METAXCAN_DIR PLINK SCRIPTS_DIR
