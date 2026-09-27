"""
Check three published HCC blood markers (FLNA, CLU, CAP1) in our two
GEO datasets. These genes were found by Puttipanyalears et al. 2021
(Scientific Reports, DOI: 10.1038/s41598-021-94330-1) using the same
two public datasets we are using, then validated in a large clinical
cohort (153 HCC, 76 healthy) with lab-based qRT-PCR.

This is a much safer test than a full trained classifier: it just asks
"does this specific, already-published gene move the same direction
in HCC vs healthy, in our data too?"

Usage:
    python check_published_markers.py
"""
import pandas as pd

from hcc_blood_classifier import load_series
from cross_dataset_validation import (
    TRAIN_GSE, TEST_GSE, TRAIN_RULES, TEST_RULES, build_xy,
)

GENES_OF_INTEREST = ["FLNA", "CLU", "CAP1"]


def gene_to_probes(gse_id, dest="data"):
    """Map gene symbols -> probe IDs using this series' own platform table."""
    import GEOparse

    gse = GEOparse.get_GEO(geo=gse_id, destdir=dest, silent=True)
    gpl_id = list(gse.gpls.keys())[0]
    table = gse.gpls[gpl_id].table
    symbol_col = next(
        (c for c in ["Gene Symbol", "gene_assignment", "Symbol"] if c in table.columns),
        None,
    )
    if symbol_col is None:
        raise RuntimeError(f"No gene symbol column found for {gse_id}")
    out = {}
    for gene in GENES_OF_INTEREST:
        matches = table[table[symbol_col] == gene]["ID"].tolist()
        out[gene] = matches
    return out


def summarize(name, expr, meta, rules, probe_map):
    _, y = build_xy(expr, meta, rules)
    print(f"\n=== {name} ===")
    for gene, probes in probe_map.items():
        if not probes:
            print(f"{gene}: no matching probe found on this platform")
            continue
        vals = expr.loc[probes, y.index].mean(axis=0)
        hcc_mean = vals[y == 1].mean()
        healthy_mean = vals[y == 0].mean()
        direction = "UP in HCC" if hcc_mean > healthy_mean else "DOWN in HCC"
        print(f"{gene} (probes: {probes}): "
              f"HCC mean={hcc_mean:.2f}, healthy mean={healthy_mean:.2f} "
              f"-> {direction}")


def main():
    print("Loading TRAIN set:", TRAIN_GSE)
    train_expr, train_meta = load_series(TRAIN_GSE)
    train_probes = gene_to_probes(TRAIN_GSE)
    summarize(TRAIN_GSE, train_expr, train_meta, TRAIN_RULES, train_probes)

    print("\nLoading TEST set:", TEST_GSE)
    test_expr, test_meta = load_series(TEST_GSE)
    test_probes = gene_to_probes(TEST_GSE)
    summarize(TEST_GSE, test_expr, test_meta, TEST_RULES, test_probes)

    print("\nPublished result (Puttipanyalears et al. 2021, Sci Rep):")
    print("FLNA and CLU were UP in HCC vs healthy, in a large validated cohort.")
    print("Check above: do our two datasets agree with 'UP in HCC' for these genes?")


if __name__ == "__main__":
    main()
