#!/usr/bin/env bash
# Step 05 — use the saved model as a real VCF filter on the held-out isolates:
# writes <sample>.ml_filtered.vcf.gz (INFO/ML_PROB + FILTER=ML_FAIL) and a
# per-isolate benchmark, then the run summary (results/summary.md).
source "$(dirname "$0")/../lib/common.sh"
need_cmd python3 bcftools
need_file "${MODEL_DIR}/variant_filter_model.joblib"
mkdir -p "${APPLY_DIR}"

while read -r s _ _ _ _ _ split; do
  [[ "${split}" == "test" ]] || continue
  python3 "${PROJECT_DIR}/scripts/apply_filter.py" --vcf "${CALL_DIR}/${s}.raw.vcf.gz" \
      --model "${MODEL_DIR}/variant_filter_model.joblib" \
      --truth "${CALL_DIR}/${s}.truth.norm.vcf.gz" --metrics "${APPLY_DIR}/${s}.metrics.json" \
    | bcftools view -Oz -o "${APPLY_DIR}/${s}.ml_filtered.vcf.gz"
  bcftools index -t "${APPLY_DIR}/${s}.ml_filtered.vcf.gz"
done < <(cohort)

python3 "${PROJECT_DIR}/scripts/make_report.py" --project "${PROJECT_DIR}" --out "${RESULTS_DIR}/summary.md"
log "summary -> ${RESULTS_DIR}/summary.md"
