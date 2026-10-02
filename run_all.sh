#!/usr/bin/env bash
# Run the whole experiment in order. Any step can also be run on its own,
# e.g.  bash steps/05_call_variants.sh
source "$(dirname "$0")/lib/common.sh"
need_cmd bwa samtools bcftools python3 curl
cd "${PROJECT_DIR}" || exit 1
start=$(date +%s)
for step in steps/0*.sh; do
  log "==================== ${step} ===================="
  bash "${step}"
done
log "pipeline finished in $(( $(date +%s) - start )) s — see results/summary.md"
