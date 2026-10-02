#!/usr/bin/env bash
# Step 04 — train Logistic Regression + Random Forest, grouped CV on the
# train isolates, final evaluation on the held-out test isolates, compared
# with the GATK-style hard filter. Model + metrics + plots -> results/04_models
source "$(dirname "$0")/../lib/common.sh"
need_cmd python3
need_file "${FEAT_DIR}/features.tsv" "${FEAT_DIR}/truth_counts.tsv"
mkdir -p "${MODEL_DIR}"

python3 "${PROJECT_DIR}/scripts/train_models.py" \
    --features "${FEAT_DIR}/features.tsv" --truth-counts "${FEAT_DIR}/truth_counts.tsv" \
    --outdir "${MODEL_DIR}" --seed "${SEED}" \
    --hf-qual "${HF_MIN_QUAL}" --hf-dp "${HF_MIN_DP}" --hf-mq "${HF_MIN_MQ}" --hf-af "${HF_MIN_AF}" \
  | tee "${MODEL_DIR}/train.log"
