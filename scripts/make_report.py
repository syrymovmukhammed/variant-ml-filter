#!/usr/bin/env python3
"""Collect the key numbers of a run into results/summary.md (also shown as
the GitHub Actions job summary)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = [f"{v:.4f}" if isinstance(v, float) else str(v) for v in r]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    proj = Path(a.project)
    res = proj / "results"
    src = (proj / "data" / "reference.source").read_text().strip()
    feats = pd.read_csv(res / "03_features" / "features.tsv", sep="\t")
    m = json.loads((res / "04_models" / "metrics.json").read_text())

    md = ["# ML variant filter — run summary", "", f"**Reference:** `{src}`"]
    if src.startswith("SYNTHETIC"):
        md.append("> ⚠️ NCBI was unreachable; this run used the synthetic offline fallback reference.")
    per = (feats.groupby(["sample", "split"])["label"].agg(calls="size", true_variants="sum")
           .reset_index())
    per["artefacts"] = per["calls"] - per["true_variants"]
    md += ["", "## Training data (permissive bcftools calls, labelled by the truth set)", "", md_table(per), ""]

    test = pd.read_csv(res / "04_models" / "test_metrics.tsv", sep="\t")
    test[["TP", "FP", "FN"]] = test[["TP", "FP", "FN"]].astype(int)
    keep = ["model", "TP", "FP", "FN", "precision", "recall", "F1"] + \
           [c for c in ("ROC_AUC", "PR_AUC") if c in test.columns]
    md += [f"## Held-out isolates ({', '.join(m['test_isolates'])})", "",
           md_table(test[keep].fillna("–")), "",
           f"Model saved for reuse: **{m['best_model_by_cv_PR_AUC']}** (best mean PR-AUC in grouped CV).", ""]

    imp = pd.read_csv(res / "04_models" / "feature_importance.tsv", sep="\t").head(8)
    md += ["## Most informative features", "", md_table(imp), ""]
    Path(a.out).write_text("\n".join(md) + "\n")


if __name__ == "__main__":
    main()
