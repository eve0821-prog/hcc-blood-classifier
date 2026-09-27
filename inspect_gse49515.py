"""
STEP 1 for the second dataset: look at GSE49515's samples before
trusting any label. Same idea as hcc_blood_classifier.py --inspect.

Usage:
    python inspect_gse49515.py
"""
import GEOparse
import pandas as pd

GSE_ID = "GSE49515"


def main():
    gse = GEOparse.get_GEO(geo=GSE_ID, destdir="data", silent=True)
    rows = {}
    for name, gsm in gse.gsms.items():
        rows[name] = {
            "title": " ".join(gsm.metadata.get("title", [])),
            "characteristics": " | ".join(gsm.metadata.get("characteristics_ch1", [])),
        }
    meta = pd.DataFrame(rows).T
    pd.set_option("display.max_colwidth", 120)
    print(f"Loaded {len(meta)} samples\n")
    print(meta.to_string())
    print("\nPlatform(s):", list(gse.gpls.keys()))


if __name__ == "__main__":
    main()
