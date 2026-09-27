cd ~/hcc-project
cat > README.md << 'EOF'
# Blood Gene Classifier for Hepatocellular Carcinoma (HCC)

A student passion project exploring whether blood gene expression
(peripheral blood mononuclear cells, PBMCs) can detect liver cancer
(HCC), using only public data.

**Why HCC:** liver cancer is a leading cause of death in Mongolia —
Mongolia has the world's highest age-standardized liver cancer death
rate (Global Burden of Disease study, 2019).

## What this repo does

1. Downloads two public gene expression datasets from NCBI GEO
   (`GSE58208`, `GSE49515`) — both PBMC samples from HCC patients and
   healthy people, on the same microarray platform (GPL570).
2. Trains a machine learning classifier (random forest) to tell HCC
   apart from healthy blood samples.
3. Honestly tests whether that classifier generalizes to a *second*,
   independent dataset it never saw during training.
4. Compares that to a small classifier built from **two published,
   independently validated genes** (FLNA, CLU).

## Key finding

A classifier trained and tested *within one dataset* looked perfect
(ROC-AUC = 1.0). But when tested on a **second, independent dataset**
it had never seen, it performed *worse than random guessing*
(ROC-AUC = 0.165–0.180 depending on direction). Standard batch
correction (ComBat) did not fix this.

This is a well-documented failure mode in blood gene expression
research: a model can learn technical differences between studies
(batch effects) instead of real biology (Leek et al., 2010).

Switching to two genes already validated in a large, independent
clinical cohort — **FLNA** and **CLU**, from Puttipanyalears et al.
(2021, *Scientific Reports*) — fixed this. That 2-gene model
generalized well in both directions:

| Train → Test | ROC-AUC | Correct calls |
|---|---|---|
| GSE58208 → GSE49515 | 0.940 | 14/20 |
| GSE49515 → GSE58208 | 1.000 | 12/15 |

**Takeaway:** for small, public blood gene expression datasets, a
handful of biologically-motivated, independently-validated genes can
generalize far better than an automatically-selected multi-gene model.

## Repository contents

| File | Purpose |
|---|---|
| `hcc_blood_classifier.py` | Loads GSE58208, labels samples, trains/cross-validates a 50-gene random forest classifier within one dataset |
| `map_probes_to_genes.py` | Converts Affymetrix probe IDs to gene symbols |
| `inspect_gse49515.py` | Inspects sample labels in the second dataset |
| `cross_dataset_validation.py` | Trains on one dataset, tests on the other (the honest, out-of-study test) |
| `cross_dataset_with_combat.py` | Same test, with ComBat batch correction applied first |
| `check_published_markers.py` | Checks whether FLNA, CLU, CAP1 move the same direction (HCC vs healthy) in both datasets |
| `two_gene_classifier.py` | A small classifier using only FLNA + CLU, tested in both cross-dataset directions |

## How to run

```bash
python3 -m venv venv
source venv/bin/activate
pip install GEOparse pandas numpy scikit-learn

python3 hcc_blood_classifier.py --inspect   # look at sample labels first
python3 hcc_blood_classifier.py             # train + cross-validate
python3 cross_dataset_validation.py         # honest cross-dataset test
python3 check_published_markers.py          # check FLNA/CLU/CAP1 direction
python3 two_gene_classifier.py              # final 2-gene model

# Optional (exploratory -- did not fix the cross-dataset problem here):
pip install combat
python3 cross_dataset_with_combat.py
```

Downloads are cached in a local `data/` folder after the first run.

## Data sources

- GSE58208 — Hui, K.M. *Gene expression profiling of PBMC from normal
  individuals, chronic hepatitis B carriers and hepatocellular
  carcinoma patients.* NCBI GEO (2014).
- GSE49515 — Hui, K.M. *Expression profiling of PBMC from patients
  with hepatocellular carcinoma.* NCBI GEO (2013).

## Key references

- Puttipanyalears, C. *et al.* The expression of FLNA and CLU in
  PBMCs as a novel screening marker for hepatocellular carcinoma.
  *Sci Rep* **11**, 14969 (2021).
  https://doi.org/10.1038/s41598-021-94330-1
- Leek, J.T. *et al.* Tackling the widespread and critical impact of
  batch effects in high-throughput data. *Nat Rev Genet* **11**,
  733–739 (2010).
- Johnson, W.E., Li, C. & Rabinovic, A. Adjusting batch effects in
  microarray expression data using empirical Bayes methods.
  *Biostatistics* **8**, 118–127 (2007). (the ComBat method)

## Limitations

- Sample sizes are small (15–20 samples per cross-dataset test).
  Results should be treated as preliminary, not clinical-grade.
- Both public datasets are from populations outside Mongolia. Results
  may not directly transfer; local validation would be needed before
  any clinical use.
- This is exploratory/educational work, not a diagnostic tool.

## Status / next steps

- [x] Data pipeline (download, label, clean)
- [x] Within-dataset classifier
- [x] Honest cross-dataset validation
- [x] Batch-effect investigation (ComBat)
- [x] Validated 2-gene model (FLNA + CLU)
- [ ] Test on a third, independent dataset
- [ ] Write up as a short research paper
- [ ] Explore relevance to liver cancer / hepatitis B burden in Mongolia
EOF
## Data sources

- GSE58208 — Hui, K.M. *Gene expression profiling of PBMC from normal
  individuals, chronic hepatitis B carriers and hepatocellular
  carcinoma patients.* NCBI GEO (2014).
- GSE49515 — Hui, K.M. *Expression profiling of PBMC from patients
  with hepatocellular carcinoma.* NCBI GEO (2013).

## Key references

- Puttipanyalears, C. *et al.* The expression of FLNA and CLU in
  PBMCs as a novel screening marker for hepatocellular carcinoma.
  *Sci Rep* **11**, 14969 (2021).
  https://doi.org/10.1038/s41598-021-94330-1
- Leek, J.T. *et al.* Tackling the widespread and critical impact of
  batch effects in high-throughput data. *Nat Rev Genet* **11**,
  733–739 (2010).
- Johnson, W.E., Li, C. & Rabinovic, A. Adjusting batch effects in
  microarray expression data using empirical Bayes methods.
  *Biostatistics* **8**, 118–127 (2007). (the ComBat method)

## Limitations

- Sample sizes are small (15–20 samples per cross-dataset test).
  Results should be treated as preliminary, not clinical-grade.
- Both public datasets are from populations outside Mongolia. Results
  may not directly transfer; local validation would be needed before
  any clinical use.
- This is exploratory/educational work, not a diagnostic tool.

## Status / next steps

- [x] Data pipeline (download, label, clean)
- [x] Within-dataset classifier
- [x] Honest cross-dataset validation
- [x] Batch-effect investigation (ComBat)
- [x] Validated 2-gene model (FLNA + CLU)
- [ ] Test on a third, independent dataset
- [ ] Write up as a short research paper
- [ ] Explore relevance to liver cancer / hepatitis B burden in Mongolia
head -1 README.md
tail -1 README.md
git add README.md
git commit -m "Clarify optional ComBat step in README"
git push
