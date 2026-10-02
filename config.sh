#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# All parameters in one place; any value can be overridden from the
# environment, e.g.  THREADS=8 REGION_END=500000 bash run_all.sh
# ---------------------------------------------------------------------------
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${THREADS:=2}"
: "${SEED:=7}"

# Reference: E. coli K-12 MG1655 (NC_000913.3), first REGION_END bp from NCBI
: "${REF_ACC:=NC_000913.3}"
: "${REGION_END:=200000}"
: "${REF_FASTA:=}"

# Cohort design: one row per simulated isolate (coverage, error rate,
# contamination, paralogs...). split=train -> model fitting with grouped CV,
# split=test -> never seen during training, used for the final report.
: "${SAMPLES_TSV:=${PROJECT_DIR}/config/samples.tsv}"
: "${N_SNPS:=150}"
: "${N_INDELS:=25}"

# Permissive calling: keep (almost) everything, let the model decide.
: "${CALL_PRIOR:=0.01}"      # bcftools call -P (default 1.1e-3); higher = more calls

# Baseline hard filter the model must beat (GATK-style, as in wgs-variant-pipeline)
: "${HF_MIN_QUAL:=30}"
: "${HF_MIN_DP:=10}"
: "${HF_MIN_MQ:=40}"
: "${HF_MIN_AF:=0.80}"

# CI acceptance: the chosen model must reach this F1 on the held-out samples
: "${MIN_TEST_F1:=0.90}"

DATA_DIR="${PROJECT_DIR}/data"
RESULTS_DIR="${PROJECT_DIR}/results"
REF="${DATA_DIR}/reference.fa"
SIM_DIR="${DATA_DIR}/cohort"
CALL_DIR="${RESULTS_DIR}/02_calls"
FEAT_DIR="${RESULTS_DIR}/03_features"
MODEL_DIR="${RESULTS_DIR}/04_models"
APPLY_DIR="${RESULTS_DIR}/05_filtered_vcfs"

export THREADS SEED
