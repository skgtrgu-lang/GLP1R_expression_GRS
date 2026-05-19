#!/usr/bin/env bash
# =============================================================================
# GLP1R Expression Genetic Risk Score (GRS) — PLINK2 Scoring
# Exenatide-PD3 Pharmacogenomics
#
# Score:   GLP1R expression GRS
# Source:  Triozzi et al. eTable 2 / Supplementary Table S3
# Weights: 15 cis-eQTL variants from GTEx v8 whole blood (EUR)
# Build:   hg38 (GRCh38)
# Author:  Rowan Gurney
# Date:    2026-03-05
# =============================================================================

set -euo pipefail

# =============================================================================
# 01 — INPUT / OUTPUT PATHS
# Edit these to match your environment
# =============================================================================

VCF="chr6.vcf.gz"                    # Joint-called WGS chr6 VCF (hg38)
SCORE_FILE="glp1r_grs_score.txt"     # Score file (provided in this repo)
OUT_PREFIX="glp1r_grs"               # Output prefix for PLINK2 files

# =============================================================================
# 02 — SCORE
# Columns in score file: 1=variant_ID  2=effect_allele  3=beta
# Variant IDs are in chrN:POS:REF:ALT format matching the VCF
# No ID reformatting required
# =============================================================================

plink2 \
    --vcf "${VCF}" \
    --vcf-half-call missing \
    --score "${SCORE_FILE}" 1 2 3 header \
    --out "${OUT_PREFIX}"

echo "Scoring complete. Output: ${OUT_PREFIX}.sscore"
