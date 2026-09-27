"""
Map Affymetrix probe IDs (e.g. 234993_at) to gene symbols.
Run this AFTER hcc_blood_classifier.py, using its --inspect run to
confirm the platform, or just run this directly (it re-downloads
the series, which GEOparse caches locally so it is fast the 2nd time).

Usage:
    python map_probes_to_genes.py
"""
import GEOparse
import pandas as pd

from hcc_blood_classifier import GSE_ID

# Paste probe IDs from your "Top probes" output here:
PROBE_IDS = [
    "234993_at", "233506_at", "235675_at", "201565_s_at", "232883_at",
    "230192_at", "203406_at", "231005_at", "1559895_x_at", "223598_at",
    "202887_s_at", "225985_at", "205239_at", "1556060_a_at", "205094_at",
]


def main():
    gse = GEOparse.get_GEO(geo=GSE_ID, destdir="data", silent=True)

    # A GEO series can use more than one platform (GPL). Usually there is one.
    gpl_ids = list(gse.gpls.keys())
    print("Platform(s) in this series:", gpl_ids)

    gpl_table = gse.gpls[gpl_ids[0]].table
    print("\nAnnotation columns available:", list(gpl_table.columns))

    # Try common column names for the gene symbol field.
    symbol_col = next(
        (c for c in ["Gene Symbol", "gene_assignment", "Symbol"] if c in gpl_table.columns),
        None,
    )
    if symbol_col is None:
        print("\nCould not auto-find a gene symbol column.")
        print("Look at the columns above and set symbol_col manually.")
        return

    lookup = gpl_table.set_index("ID")[symbol_col]

    result = pd.DataFrame({
        "probe_id": PROBE_IDS,
        "gene_symbol": [lookup.get(p, "NOT_FOUND") for p in PROBE_IDS],
    })
    print("\nProbe -> gene mapping:")
    print(result.to_string(index=False))
    result.to_csv("probe_to_gene.csv", index=False)
    print("\nSaved to probe_to_gene.csv")


if __name__ == "__main__":
    main()
