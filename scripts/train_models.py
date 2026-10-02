#!/usr/bin/env python3
"""
train_models.py — learn to separate true variants from artefacts.

Models
  * Logistic Regression (standardised features, class-balanced) — linear,
    interpretable coefficients
  * Random Forest (class-balanced) — non-linear, captures interactions such
    as "high depth AND mid allele fraction" (collapsed paralogs)
  * Baseline: the GATK-style hard filter of the wgs-variant-pipeline repo
    (QUAL, DP, MQ, AF thresholds) — the thing the models have to beat

Validation
  * GroupKFold over TRAIN isolates (all calls of an isolate stay in the same
    fold) — a random row split would leak isolate-level properties
    (coverage, contamination) between train and validation
  * final, untouched evaluation on the TEST isolates

Recall is computed against the full truth set of each isolate (truth variants
the caller never produced count as misses), so numbers are comparable with
the wgs-variant-pipeline benchmark.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score, roc_curve
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from vcf_features import FEATURES


def make_models(seed: int) -> dict:
    return {
        "logistic_regression": make_pipeline(
            StandardScaler(), LogisticRegression(max_iter=2000, class_weight="balanced", C=1.0)),
        "random_forest": RandomForestClassifier(
            n_estimators=400, min_samples_leaf=3, class_weight="balanced_subsample",
            n_jobs=-1, random_state=seed),
    }


def hard_filter(df: pd.DataFrame, a) -> np.ndarray:
    return ((df["QUAL"] >= a.hf_qual) & (df["DP"] >= a.hf_dp) &
            (df["MQ"] >= a.hf_mq) & (df["AF"] >= a.hf_af)).astype(int).to_numpy()


def scores(y: np.ndarray, keep: np.ndarray, n_truth: int, prob: np.ndarray | None = None) -> dict:
    """precision over kept calls; recall over ALL truth variants of the isolates."""
    tp = int(((keep == 1) & (y == 1)).sum())
    fp = int(((keep == 1) & (y == 0)).sum())
    p = tp / (tp + fp) if tp + fp else 0.0
    r = tp / n_truth if n_truth else 0.0
    out = {"TP": tp, "FP": fp, "FN": n_truth - tp, "artefacts_removed": int(((keep == 0) & (y == 0)).sum()),
           "precision": p, "recall": r, "F1": 2 * p * r / (p + r) if p + r else 0.0}
    if prob is not None:
        out["ROC_AUC"] = roc_auc_score(y, prob)
        out["PR_AUC"] = average_precision_score(y, prob)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features", required=True)
    ap.add_argument("--truth-counts", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--hf-qual", type=float, default=30)
    ap.add_argument("--hf-dp", type=float, default=10)
    ap.add_argument("--hf-mq", type=float, default=40)
    ap.add_argument("--hf-af", type=float, default=0.8)
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(a.features, sep="\t")
    n_truth = pd.read_csv(a.truth_counts, sep="\t").set_index("sample")["n_truth"]
    train, test = df[df.split == "train"].reset_index(drop=True), df[df.split == "test"].reset_index(drop=True)
    X, y, groups = train[FEATURES].to_numpy(float), train["label"].to_numpy(), train["sample"].to_numpy()
    print(f"train: {len(train)} calls from {train['sample'].nunique()} isolates "
          f"({y.mean():.1%} true) | test: {len(test)} calls from {test['sample'].nunique()} isolates")

    # ---------------------------------------------------- grouped cross-validation
    n_splits = min(5, train["sample"].nunique())
    cv_rows = []
    for fold, (tr, va) in enumerate(GroupKFold(n_splits=n_splits).split(X, y, groups)):
        va_truth = int(n_truth[train.loc[va, "sample"].unique()].sum())
        cv_rows.append({"fold": fold, "model": "hard_filter_baseline",
                        **scores(y[va], hard_filter(train.loc[va], a), va_truth)})
        for name, model in make_models(a.seed).items():
            model.fit(X[tr], y[tr])
            prob = model.predict_proba(X[va])[:, 1]
            cv_rows.append({"fold": fold, "model": name,
                            **scores(y[va], (prob >= a.threshold).astype(int), va_truth, prob)})
    cv = pd.DataFrame(cv_rows)
    cv.to_csv(out / "cv_folds.tsv", sep="\t", index=False)
    cv_summary = cv.drop(columns="fold").groupby("model").agg(["mean", "std"]).round(4)
    cv_summary.to_csv(out / "cv_summary.tsv", sep="\t")

    # ---------------------------------------- final fit on all TRAIN, score TEST
    Xt, yt = test[FEATURES].to_numpy(float), test["label"].to_numpy()
    test_truth = int(n_truth[test["sample"].unique()].sum())
    raw_keep = np.ones_like(yt)
    results = {"no_filter (raw calls)": scores(yt, raw_keep, test_truth),
               "hard_filter_baseline": scores(yt, hard_filter(test, a), test_truth)}
    fitted, probs = {}, {}
    for name, model in make_models(a.seed).items():
        model.fit(X, y)
        fitted[name] = model
        probs[name] = model.predict_proba(Xt)[:, 1]
        results[name] = scores(yt, (probs[name] >= a.threshold).astype(int), test_truth, probs[name])
    test_tbl = pd.DataFrame(results).T
    test_tbl.index.name = "model"
    test_tbl.to_csv(out / "test_metrics.tsv", sep="\t", float_format="%.4f")

    best = max(fitted, key=lambda m: cv[cv.model == m]["PR_AUC"].mean())
    joblib.dump({"model": fitted[best], "features": FEATURES, "threshold": a.threshold, "name": best},
                out / "variant_filter_model.joblib")

    # ------------------------------------------------------------ interpretation
    rf = fitted["random_forest"]
    perm = permutation_importance(rf, Xt, yt, scoring="average_precision", n_repeats=10,
                                  random_state=a.seed, n_jobs=-1)
    imp = pd.DataFrame({"feature": FEATURES, "rf_gini_importance": rf.feature_importances_,
                        "permutation_importance_test_PR_AUC": perm.importances_mean,
                        "logreg_coef_standardised": fitted["logistic_regression"][-1].coef_[0]})
    imp = imp.sort_values("permutation_importance_test_PR_AUC", ascending=False)
    imp.to_csv(out / "feature_importance.tsv", sep="\t", index=False, float_format="%.4f")

    # ------------------------------------------------------------------- plots
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    colors = {"logistic_regression": "#1f77b4", "random_forest": "#d62728"}
    for name, p in probs.items():
        fpr, tpr, _ = roc_curve(yt, p)
        ax[0].plot(fpr, tpr, color=colors[name], label=f"{name} (AUC {results[name]['ROC_AUC']:.3f})")
        pr, rc, _ = precision_recall_curve(yt, p)
        ax[1].plot(rc, pr, color=colors[name], label=f"{name} (AP {results[name]['PR_AUC']:.3f})")
    hf = hard_filter(test, a)
    hf_tpr, hf_fpr = hf[yt == 1].mean(), hf[yt == 0].mean()
    hf_prec = yt[hf == 1].mean() if hf.sum() else 0
    ax[0].scatter([hf_fpr], [hf_tpr], marker="X", s=90, c="k", zorder=5, label="hard filter")
    ax[1].scatter([hf_tpr], [hf_prec], marker="X", s=90, c="k", zorder=5, label="hard filter")
    ax[0].plot([0, 1], [0, 1], ls=":", c="grey")
    ax[0].set(xlabel="False positive rate", ylabel="True positive rate", title="ROC — held-out isolates")
    ax[1].set(xlabel="Recall (of called variants)", ylabel="Precision", title="Precision-recall — held-out isolates")
    for x in ax:
        x.legend(loc="lower right", fontsize=8); x.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(out / "roc_pr_test.png", dpi=130); plt.close(fig)

    top = imp.head(12)[::-1]
    fig, axi = plt.subplots(figsize=(6, 4.5))
    axi.barh(top["feature"], top["permutation_importance_test_PR_AUC"], color="#d62728")
    axi.set(xlabel="drop in PR-AUC when shuffled", title="Random Forest — permutation importance (test)")
    axi.grid(axis="x", alpha=.3)
    fig.tight_layout(); fig.savefig(out / "feature_importance.png", dpi=130); plt.close(fig)

    summary = {"best_model_by_cv_PR_AUC": best, "threshold": a.threshold,
               "train_isolates": sorted(train["sample"].unique()), "test_isolates": sorted(test["sample"].unique()),
               "test": {k: {m: round(float(v), 4) for m, v in d.items()} for k, d in results.items()}}
    (out / "metrics.json").write_text(json.dumps(summary, indent=2))

    print("\nGrouped CV (train isolates), mean over folds:")
    print(cv.drop(columns="fold").groupby("model")[["precision", "recall", "F1", "PR_AUC"]].mean().round(4).to_string())
    print("\nHeld-out test isolates:")
    print(test_tbl[["TP", "FP", "FN", "precision", "recall", "F1"]].to_string(float_format=lambda v: f"{v:.4f}"))
    print(f"\nselected model: {best}  ->  {out / 'variant_filter_model.joblib'}")
    print("\ntop features (permutation importance on test):")
    print(imp.head(6)[["feature", "permutation_importance_test_PR_AUC"]].to_string(index=False))


if __name__ == "__main__":
    main()
