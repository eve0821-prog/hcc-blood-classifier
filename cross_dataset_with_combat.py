"""
Cross-dataset validation, WITH ComBat batch correction.

Same train/test setup as cross_dataset_validation.py (train on GSE58208,
test on GSE49515), but first runs ComBat to remove batch (dataset)
differences before training.

IMPORTANT, read this:
  * ComBat uses the BATCH label (which dataset a sample is from) --
    never the disease label (HCC/healthy) -- to estimate the correction.
  * It is fit on train + test expression VALUES together. This is
    standard practice for combining cohorts in genomics, but it does
    mean the test set's expression data (not its labels) is used
    during preprocessing. Say this plainly in your write-up.
  * The disease labels are still never touched until scoring, so
    this is not label leakage -- just be transparent about the method.

Usage:
    pip install combat
    python cross_dataset_with_combat.py
"""
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from hcc_blood_classifier import build_model, load_series
from cross_dataset_validation import TRAIN_GSE, TEST_GSE, TRAIN_RULES, TEST_RULES, build_xy, K


def log2_if_needed(X):
    if X.values.max() > 100:
        return np.log2(X.clip(lower=0) + 1)
    return X


def main():
    from combat.pycombat import pycombat  # imported here: fails loudly if not installed

    print(f"Loading TRAIN set: {TRAIN_GSE}")
    train_expr, train_meta = load_series(TRAIN_GSE)
    X_train, y_train = build_xy(train_expr, train_meta, TRAIN_RULES)
    X_train = log2_if_needed(X_train)
    print(f"Train: {X_train.shape[0]} samples, {y_train.value_counts().to_dict()}")

    print(f"\nLoading TEST set: {TEST_GSE}")
    test_expr, test_meta = load_series(TEST_GSE)
    X_test, y_test = build_xy(test_expr, test_meta, TEST_RULES)
    X_test = log2_if_needed(X_test)
    print(f"Test: {X_test.shape[0]} samples, {y_test.value_counts().to_dict()}")

    common = X_train.columns.intersection(X_test.columns)
    print(f"\nShared probes: {len(common)}")
    X_train, X_test = X_train[common], X_test[common]

    # pycombat expects genes as ROWS, samples as COLUMNS.
    combined = pd.concat([X_train.T, X_test.T], axis=1)
    batch = ["train"] * X_train.shape[0] + ["test"] * X_test.shape[0]

    print("\nRunning ComBat (this can take a minute)...")
    corrected = pycombat(combined, batch)

    X_train_corr = corrected[X_train.index.tolist()].T
    X_test_corr = corrected[X_test.index.tolist()].T

    model = build_model(k=min(K, X_train_corr.shape[1]))
    model.fit(X_train_corr, y_train)

    proba = model.predict_proba(X_test_corr)[:, 1]
    auc = roc_auc_score(y_test, proba)
    preds = (proba >= 0.5).astype(int)
    correct = (preds == y_test.values).sum()

    print(f"\nExternal test ROC-AUC on {TEST_GSE} (after ComBat): {auc:.3f}")
    print(f"Correct calls: {correct}/{len(y_test)}")
    print("\nCompare this to the pre-ComBat AUC from cross_dataset_validation.py.")
    print("If this AUC is now above 0.5 and closer to (or above) 0.7,")
    print("ComBat likely fixed a real batch effect. If it is still far")
    print("from 0.5 in either direction, the two studies may not be")
    print("comparable enough for a simple pooled model.")


if __name__ == "__main__":
    main()
