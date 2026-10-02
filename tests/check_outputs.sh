#!/usr/bin/env bash
# Acceptance test used by CI. Fails (exit 1) if any artefact is missing, if
# the model is worse than the hard-filter baseline on the held-out isolates,
# or if its test F1 is below MIN_TEST_F1.
source "$(dirname "$0")/../lib/common.sh"
fail=0
check() { if eval "$2"; then echo "  ok    $1"; else echo "  FAIL  $1"; fail=1; fi; }

echo "== artefacts"
for f in "${FEAT_DIR}/features.tsv" "${MODEL_DIR}/variant_filter_model.joblib" \
         "${MODEL_DIR}/metrics.json" "${MODEL_DIR}/cv_folds.tsv" "${MODEL_DIR}/test_metrics.tsv" \
         "${MODEL_DIR}/feature_importance.tsv" "${MODEL_DIR}/roc_pr_test.png" \
         "${MODEL_DIR}/feature_importance.png" "${RESULTS_DIR}/summary.md"; do
  check "$(basename "$f")" "[[ -s '$f' ]]"
done

echo "== data sanity"
check "both classes present in training data" \
      "python3 -c \"import pandas as p; d=p.read_csv('${FEAT_DIR}/features.tsv',sep='\t'); exit(d.label.nunique()!=2)\""
while read -r s _ _ _ _ _ split; do
  [[ "${split}" == "test" ]] || continue
  vcf="${APPLY_DIR}/${s}.ml_filtered.vcf.gz"
  check "${s}: ML-filtered VCF parses and carries ML_PROB" \
        "bcftools view -h '${vcf}' | grep 'ID=ML_PROB' > /dev/null && bcftools view -H '${vcf}' | grep 'ML_PROB=' > /dev/null"
done < <(cohort)

echo "== model quality on held-out isolates"
read -r best ml_f1 hf_f1 < <(python3 - "${MODEL_DIR}/metrics.json" <<'PY'
import json, sys
m = json.load(open(sys.argv[1])); b = m["best_model_by_cv_PR_AUC"]
print(b, m["test"][b]["F1"], m["test"]["hard_filter_baseline"]["F1"])
PY
)
check "${best} F1 ${ml_f1} >= ${MIN_TEST_F1}" "awk 'BEGIN{exit !(${ml_f1} >= ${MIN_TEST_F1})}'"
check "${best} F1 ${ml_f1} >= hard-filter F1 ${hf_f1}" "awk 'BEGIN{exit !(${ml_f1} >= ${hf_f1})}'"

[[ ${fail} -eq 0 ]] && echo "ALL CHECKS PASSED" || { echo "SOME CHECKS FAILED"; exit 1; }
