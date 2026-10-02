#!/usr/bin/env bash
# Step 03 — VCF -> feature table, labelled against each isolate's truth set.
# Output: results/03_features/features.tsv (one row per call, all isolates)
# plus truth_counts.tsv (number of true variants per isolate, needed to
# compute recall: a classifier can only keep calls the caller produced).
source "$(dirname "$0")/../lib/common.sh"
need_cmd python3 bcftools
mkdir -p "${FEAT_DIR}"

printf 'sample\tsplit\tn_truth\n' > "${FEAT_DIR}/truth_counts.tsv"
tables=()
while read -r s _ _ _ _ _ split; do
  need_file "${CALL_DIR}/${s}.raw.vcf.gz"
  python3 "${PROJECT_DIR}/scripts/vcf_features.py" --vcf "${CALL_DIR}/${s}.raw.vcf.gz" \
      --sample "${s}" --split "${split}" --truth "${CALL_DIR}/${s}.truth.norm.vcf.gz" \
      --out "${FEAT_DIR}/${s}.features.tsv"
  tables+=("${FEAT_DIR}/${s}.features.tsv")
  printf '%s\t%s\t%s\n' "${s}" "${split}" \
      "$(bcftools view -H "${CALL_DIR}/${s}.truth.norm.vcf.gz" | wc -l)" >> "${FEAT_DIR}/truth_counts.tsv"
done < <(cohort)

# concatenate, keeping a single header
awk 'FNR==1 && NR!=1 {next} {print}' "${tables[@]}" > "${FEAT_DIR}/features.tsv"
log "feature table: $(( $(wc -l < "${FEAT_DIR}/features.tsv") - 1 )) labelled calls -> ${FEAT_DIR}/features.tsv"
