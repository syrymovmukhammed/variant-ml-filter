#!/usr/bin/env bash
# Step 01 — simulate a cohort of E. coli isolates with a known truth set each.
# Every row of config/samples.tsv becomes one isolate. Besides the planted
# SNPs/indels, isolates carry the confounders that create false positives in
# real bacterial WGS:
#   * low coverage + higher error rate   -> sequencing-error calls
#   * contamination with another strain  -> low allele-fraction calls
#   * plasmid paralogs (diverged copies of chromosomal genes that are absent
#     from the reference)                -> collapsed-repeat / mismapping calls
source "$(dirname "$0")/../lib/common.sh"
need_cmd python3
need_file "${REF}" "${SAMPLES_TSV}"
mkdir -p "${SIM_DIR}"

while read -r s cov err contam par seed _split; do
  python3 "${PROJECT_DIR}/scripts/wgs_sim.py" simulate \
    --reference "${REF}" --outdir "${SIM_DIR}" --prefix "${s}" --sample "${s}" \
    --n-snps "${N_SNPS}" --n-indels "${N_INDELS}" --coverage "${cov}" \
    --error-rate "${err}" --contam-frac "${contam}" --plasmid-paralogs "${par}" \
    --seed "${seed}"
done < <(cohort)
log "cohort ready: $(cohort | wc -l) isolates in ${SIM_DIR}"
