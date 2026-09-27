"""
Cross-dataset validation for HCC vs healthy, blood (PBMC) gene expression.

Train on GSE49515 (larger, more balanced: 10 HCC + 10 healthy).
Test on GSE58208 (10 HCC + 5 healthy) -- a dataset the model
never sees during training or gene selection.

Both series use platform GPL570, so probe IDs line up directly.

Usage:
    python cross_dataset_validation.py
"""
import re

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from hcc_blood_classifier import build_model, load_series

TRAIN_GSE = "GSE49515"
TEST_GSE = "GSE58208"
K = 50
SEED = 42

TRAIN_RULES = {"HCC": ["hepatocellular"], "healthy": ["healthy"]}
TEST_RULES = {"HCC": ["hepatocellular"], "healthy": ["healthy control"]}


def label_from_rules(meta, rules):
    labels = {}
    for sample, row in meta.iterrows():
        text = f"{row['title']} {row['characteristics']}".lower()
        hits = [c for c, words in rules.items()
                if any(re.search(r"\b" + re.escape(w) + r"\b", text) for w in words)]
        labels[sample] = hits[0] if len(hits) == 1 else "OTHER"
    return pd.Series(labels)


def build_xy(expr, meta, rules):
    labels = label_from_rules(meta, rules)
    keep = labels[labels.isin(["HCC", "healthy"])].index
    y = (labels.loc[keep] == "HCC").astype(int)
    X = expr[keep].T.astype(float)
    return X, y


def zscore(X):
    """Standardize each gene (column) within its own dataset.
    Simple per-study centering, a common light batch correction."""
    return (X - X.mean(axis=0)) / X.std(axis=0).replace(0, 1)


def main():
    print(f"Loading TRAIN set: {TRAIN_GSE}")
    train_expr, train_meta = load_series(TRAIN_GSE)
    X_train, y_train = build_xy(train_expr, train_meta, TRAIN_RULES)
    print(f"Train: {X_train.shape[0]} samples, "
          f"{y_train.value_counts().to_dict()}")

    print(f"\nLoading TEST set: {TEST_GSE}")
    test_expr, test_meta = load_series(TEST_GSE)
    X_test, y_test = build_xy(test_expr, test_meta, TEST_RULES)
    print(f"Test: {X_test.shape[0]} samples, "
          f"{y_test.value_counts().to_dict()}")

    common = X_train.columns.intersection(X_test.columns)
    print(f"\nShared probes between platforms: {len(common)}")
    X_train, X_test = X_train[common], X_test[common]

    X_train = zscore(X_train)
    X_test = zscore(X_test)

    model = build_model(k=min(K, X_train.shape[1]))
    model.fit(X_train, y_train)

    proba = model.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, proba)
    print(f"\nExternal test ROC-AUC on {TEST_GSE}: {auc:.3f}")
    print("(This model was trained ONLY on", TRAIN_GSE, "-- it never saw",
          TEST_GSE, "until this line.)")

    preds = (proba >= 0.5).astype(int)
    correct = (preds == y_test.values).sum()
    print(f"Correct calls: {correct}/{len(y_test)}")


if __name__ == "__main__":
    main()
