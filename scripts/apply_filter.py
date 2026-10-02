#!/usr/bin/env python3
"""
apply_filter.py — use the trained model as a drop-in VCF filter.

Reads any bcftools VCF (split with `bcftools norm -m -any`), adds
INFO/ML_PROB (probability that the call is a real variant) and sets
FILTER=ML_FAIL where ML_PROB < threshold; other records get PASS.
Output is plain VCF on stdout -> pipe into `bcftools view -Oz`.

  python3 apply_filter.py --vcf calls.vcf.gz --model variant_filter_model.joblib \
      | bcftools view -Oz -o calls.ml.vcf.gz

With --truth, a one-line benchmark of the PASS set is written to stderr/--metrics.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys

import joblib

from vcf_features import parse_vcf, truth_keys


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--threshold", type=float, help="override the threshold stored with the model")
    ap.add_argument("--truth", help="normalised truth VCF for a quick benchmark")
    ap.add_argument("--metrics", help="write benchmark JSON here")
    a = ap.parse_args()

    bundle = joblib.load(a.model)
    thr = a.threshold if a.threshold is not None else bundle["threshold"]
    df = parse_vcf(a.vcf)
    prob = bundle["model"].predict_proba(df[bundle["features"]].to_numpy(float))[:, 1] if len(df) else []

    opener = gzip.open if a.vcf.endswith(".gz") else open
    out = sys.stdout
    i, kept = 0, set()
    with opener(a.vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                out.write('##INFO=<ID=ML_PROB,Number=1,Type=Float,Description='
                          f'"P(real variant) from {bundle["name"]} (variant-ml-filter)">\n')
                out.write(f'##FILTER=<ID=ML_FAIL,Description="ML_PROB < {thr}">\n')
                out.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            p = float(prob[i]); i += 1
            f[6] = "PASS" if p >= thr else "ML_FAIL"
            f[7] = f"{f[7]};ML_PROB={p:.4f}" if f[7] not in (".", "") else f"ML_PROB={p:.4f}"
            if p >= thr:
                kept.add((f[0], int(f[1]), f[3].upper(), f[4].upper()))
            out.write("\t".join(f) + "\n")

    msg = {"calls": i, "pass": len(kept), "ml_fail": i - len(kept), "model": bundle["name"], "threshold": thr}
    if a.truth:
        t = truth_keys(a.truth)
        tp = len(kept & t)
        p = tp / len(kept) if kept else 0.0
        r = tp / len(t) if t else 0.0
        msg.update(TP=tp, FP=len(kept) - tp, FN=len(t) - tp, precision=round(p, 4), recall=round(r, 4),
                   F1=round(2 * p * r / (p + r), 4) if p + r else 0.0)
    print(json.dumps(msg), file=sys.stderr)
    if a.metrics:
        with open(a.metrics, "w") as fh:
            json.dump(msg, fh, indent=2)


if __name__ == "__main__":
    main()
