"""
Small classifier using ONLY FLNA + CLU -- the 2-gene combination that
Puttipanyalears et al. 2021 (Sci Rep) found gave 86.2% accuracy in a
large clinical cohort (153 HCC, 76 healthy).

Train on one dataset, test on the other (both directions), using just
these 2 genes instead of a 50-gene black-box model.

Note: label rules are defined HERE, per real dataset, and do not
depend on cross_dataset_validation.py's current state (that file may
have been edited/swapped earlier and its meaning can change).

Usage:
    python two_gene_classifier.py
"""
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from hcc_blood_classifier import load_series
from cross_dataset_validation import build_xy
from check_published_markers import gene_to_probes

GENES = ["FLNA", "CLU"]

# Real label wording for each dataset, fixed here regardless of any
# other script's current state.
RULES_BY_DATASET = {
    "GSE58208": {"HCC": ["hepatocellular"], "healthy": ["healthy control"]},
    "GSE49515": {"HCC": ["hepatocellular"], "healthy": ["healthy"]},
}


def gene_features(expr, meta, rules, gse_id):
    _, y = build_xy(expr, meta, rules)
    probe_map = gene_to_probes(gse_id)
    cols = {}
    for gene in GENES:
        probes = probe_map.get(gene, [])
        if not probes:
            raise RuntimeError(f"{gene} not found on platform for {gse_id}")
        cols[gene] = expr.loc[probes, y.index].mean(axis=0)
    return pd.DataFrame(cols), y


def zscore(X):
    return (X - X.mean(axis=0)) / X.std(axis=0).replace(0, 1)


def run_direction(train_gse, test_gse, label):
    print(f"\n--- {label}: train on {train_gse}, test on {test_gse} ---")
    train_expr, train_meta = load_series(train_gse)
    X_train, y_train = gene_features(
        train_expr, train_meta, RULES_BY_DATASET[train_gse], train_gse
    )
    print(f"Train: {X_train.shape[0]} samples, {y_train.value_counts().to_dict()}")

    test_expr, test_meta = load_series(test_gse)
    X_test, y_test = gene_features(
        test_expr, test_meta, RULES_BY_DATASET[test_gse], test_gse
    )
    print(f"Test: {X_test.shape[0]} samples, {y_test.value_counts().to_dict()}")

    if y_test.nunique() < 2:
        print("ERROR: test set has only one class after labeling. "
              "Check RULES_BY_DATASET wording against --inspect output.")
        return None

    X_train_z = zscore(X_train)
    X_test_z = zscore(X_test)

    model = LogisticRegression()
    model.fit(X_train_z, y_train)

    proba = model.predict_proba(X_test_z)[:, 1]
    auc = roc_auc_score(y_test, proba)
    preds = (proba >= 0.5).astype(int)
    correct = (preds == y_test.values).sum()

    print(f"External test ROC-AUC: {auc:.3f}")
    print(f"Correct calls: {correct}/{len(y_test)}")
    return auc


def main():
    auc1 = run_direction("GSE58208", "GSE49515", "Direction A")
    auc2 = run_direction("GSE49515", "GSE58208", "Direction B")
    print(f"\nSummary: Direction A AUC={auc1}, Direction B AUC={auc2}")
    print("Compare both to the 50-gene model's 0.165 and 0.180.")
    print("Closer to 0.5 or above 0.7 in BOTH directions = the 2-gene")
    print("model generalizes better than the full classifier did.")


if __name__ == "__main__":
    main()
