"""
Starter pipeline: blood (PBMC) gene expression -> HCC vs healthy.
Dataset: GEO GSE58208 (public, human).

Usage:
    pip install GEOparse pandas numpy scikit-learn
    python hcc_blood_classifier.py --inspect   # STEP 1: look at sample labels
    python hcc_blood_classifier.py             # STEP 2: run the model

Rules this script enforces:
  * Gene selection happens INSIDE cross-validation (no data leakage).
  * Labels come from sample text; any unclear sample stops the run.
  * A permutation test checks that the score is not luck.
"""
import argparse
import re
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.model_selection import (
    RepeatedStratifiedKFold,
    cross_val_score,
    permutation_test_score,
)
from sklearn.pipeline import Pipeline

GSE_ID = "GSE58208"
POSITIVE = "HCC"
NEGATIVE = "CHB"

# Keywords used to label samples. CHECK THESE with --inspect before trusting.
LABEL_RULES = {
    "HCC": ["hcc", "hepatocellular"],
    "healthy": ["healthy", "normal", "control"],
    "CHB": ["chb", "hepatitis b carrier", "chronic hepatitis b"],
}
SEED = 42


def load_series(gse_id=GSE_ID, dest="data"):
    """Download a GEO series. Returns (probes x samples matrix, metadata)."""
    import GEOparse  # imported here so the rest can be tested without it

    gse = GEOparse.get_GEO(geo=gse_id, destdir=dest, silent=True)
    expr = gse.pivot_samples("VALUE")
    rows = {}
    for name, gsm in gse.gsms.items():
        rows[name] = {
            "title": " ".join(gsm.metadata.get("title", [])),
            "characteristics": " | ".join(gsm.metadata.get("characteristics_ch1", [])),
        }
    meta = pd.DataFrame(rows).T
    return expr, meta


def label_samples(meta, rules=LABEL_RULES):
    """Give each sample exactly one label, or 'UNKNOWN' if 0 or 2+ match."""
    labels = {}
    for sample, row in meta.iterrows():
        text = f"{row['title']} {row['characteristics']}".lower()
        hits = [
            cls
            for cls, words in rules.items()
            if any(re.search(r"\b" + re.escape(w) + r"\b", text) for w in words)
        ]
        labels[sample] = hits[0] if len(hits) == 1 else "UNKNOWN"
    return pd.Series(labels, name="label")


def prepare(expr, meta):
    """Build X (samples x probes) and y (1 = HCC, 0 = healthy)."""
    labels = label_samples(meta)
    print("\nLabel counts:\n", labels.value_counts().to_string())

    keep = labels[labels.isin([POSITIVE, NEGATIVE])].index
    unknown = labels[labels == "UNKNOWN"]
    if len(unknown) > 0:
        print(f"\nWARNING: {len(unknown)} samples could not be labelled:")
        print(meta.loc[unknown.index].to_string())
        sys.exit("Fix LABEL_RULES, then run again. Nothing was modelled.")

    y = (labels.loc[keep] == POSITIVE).astype(int)
    X = expr[keep].T.astype(float)
    X = X.dropna(axis=1)  # drop probes with any missing value

    counts = y.value_counts()
    if len(counts) < 2 or counts.min() < 5:
        sys.exit(f"Too few samples in one class: {counts.to_dict()}")

    # Arrays should be on a log2 scale. Very large values suggest raw intensity.
    if X.values.max() > 100:
        print("Values look un-logged (max > 100). Applying log2(x + 1).")
        X = np.log2(X.clip(lower=0) + 1)
    return X, y


def build_model(k=50):
    # Every step below is fitted on TRAINING folds only.
    return Pipeline(
        [
            ("var", VarianceThreshold(threshold=0.0)),
            ("select", SelectKBest(f_classif, k=k)),
            (
                "rf",
                RandomForestClassifier(
                    n_estimators=300, random_state=SEED, n_jobs=-1
                ),
            ),
        ]
    )


def evaluate(X, y, k=50, n_perm=100):
    n_splits = int(min(5, y.value_counts().min()))
    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=10, random_state=SEED)
    k = min(k, X.shape[1])
    model = build_model(k)

    scores = cross_val_score(model, X, y, cv=cv, scoring="roc_auc")
    print(f"\nCV ROC-AUC: {scores.mean():.3f} +/- {scores.std():.3f}"
          f"  (samples={len(y)}, probes={X.shape[1]}, folds={n_splits}x10)")

    perm_cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=1, random_state=SEED)
    _, _, p = permutation_test_score(
        model, X, y, cv=perm_cv, scoring="roc_auc",
        n_permutations=n_perm, random_state=SEED, n_jobs=1,
    )
    print(f"Permutation p-value: {p:.4f}  (small = not luck)")
    return scores, p


def top_probes(X, y, k=50, n=15):
    """Exploration only: probes ranked on ALL data. Do not report as a result."""
    model = build_model(min(k, X.shape[1])).fit(X, y)
    kept = X.columns[model.named_steps["var"].get_support()]
    chosen = kept[model.named_steps["select"].get_support()]
    imp = pd.Series(model.named_steps["rf"].feature_importances_, index=chosen)
    return imp.sort_values(ascending=False).head(n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inspect", action="store_true", help="print sample metadata and stop")
    ap.add_argument("--k", type=int, default=50, help="number of probes to keep")
    args = ap.parse_args()

    expr, meta = load_series()
    print(f"Loaded {expr.shape[0]} probes x {expr.shape[1]} samples")

    if args.inspect:
        pd.set_option("display.max_colwidth", 120)
        print(meta.to_string())
        print("\nAuto-labels:\n", label_samples(meta).value_counts().to_string())
        return

    X, y = prepare(expr, meta)
    evaluate(X, y, k=args.k)
    print("\nTop probes (exploration only):\n", top_probes(X, y, k=args.k).to_string())


if __name__ == "__main__":
    main()
