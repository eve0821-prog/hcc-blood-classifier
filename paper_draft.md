# Cross-Study Generalization of Blood Gene Expression Classifiers for Hepatocellular Carcinoma: A Cautionary Case Study

**Otgontugs Sugar**

*Draft — student project*

## Abstract

Hepatocellular carcinoma (HCC) is a leading cause of cancer death
worldwide and is especially common in Mongolia, where it is linked to
a high burden of hepatitis B infection. Blood-based gene expression
from peripheral blood mononuclear cells (PBMCs) has been proposed as
a low-cost, non-invasive screening approach. We built a machine
learning classifier using two public PBMC gene expression datasets
(GSE58208, GSE49515) and tested whether it generalized across
studies. A classifier using 50 automatically selected genes achieved
perfect discrimination within a single dataset (ROC-AUC = 1.00) but
performed worse than chance when tested on an independent dataset
(ROC-AUC = 0.165–0.180 depending on training direction), and standard
batch correction (ComBat) did not resolve this. In contrast, a
2-gene model built from FLNA and CLU — markers independently
validated in a large clinical cohort by Puttipanyalears et al. (2021)
— generalized well across both datasets (ROC-AUC = 0.940 and 1.000).
These results illustrate a well-known but easily overlooked risk in
small-sample blood transcriptomics: automatically selected multi-gene
signatures can fail to generalize even when they appear highly
accurate, while biologically-motivated, independently validated
markers may be more robust.

## 1. Introduction

HCC accounts for a disproportionate share of cancer mortality in
Mongolia; it has been reported to have the world's highest
age-standardized liver cancer mortality rate, driven substantially by
hepatitis B and C infection (Global Burden of Disease Study, 2019).
Current screening tools (alpha-fetoprotein blood tests, ultrasound)
have known sensitivity and specificity limitations. Blood-based gene
expression from PBMCs has been explored as an inexpensive,
minimally-invasive alternative, motivated by the observation that
tumors can influence circulating immune cell gene expression.

A recurring challenge in this line of research is that gene
expression datasets are small (tens of samples) and collected in
different labs, at different times, using different protocols —
producing "batch effects" that can be mistaken for real biological
signal (Leek et al., 2010). This project set out to build a basic HCC
blood classifier, then explicitly test how well it holds up outside
the dataset it was trained on.

## 2. Data

Two public datasets from NCBI GEO were used, both PBMC samples on the
Affymetrix GPL570 platform:

- **GSE58208**: 10 HCC, 12 chronic hepatitis B carriers, 5 healthy
  controls.
- **GSE49515**: 10 HCC, 10 healthy controls, plus 3 gastric and 3
  pancreatic cancer samples (not used here).

Only HCC and healthy-control samples were used for the classifiers
described below.

## 3. Methods

### 3.1 Within-dataset classifier

A random forest classifier (300 trees) was trained on GSE58208 (10
HCC, 5 healthy). Gene selection (ANOVA F-test, top 50 probes) was
performed *inside* a repeated stratified 5-fold cross-validation loop
(10 repeats) to avoid data leakage. Performance was measured by
ROC-AUC, and a permutation test (100 label permutations) was used to
estimate whether the score exceeded chance.

### 3.2 Cross-dataset validation

The same modeling approach (variance filtering, top-50-gene ANOVA
selection, random forest) was trained on the full HCC/healthy subset
of one dataset and evaluated, without refitting, on the full
HCC/healthy subset of the other dataset. Each dataset was
standardized (z-scored) independently before combination, to reduce
gross scale differences. Both directions were tested (GSE58208 →
GSE49515 and GSE49515 → GSE58208).

### 3.3 Batch correction

ComBat (Johnson, Li & Rabinovic, 2007), a standard empirical Bayes
batch-correction method, was applied to the combined expression
matrix (genes as rows, samples as columns) using dataset identity as
the batch variable, before repeating the cross-dataset validation.

### 3.4 Published-marker replication

Puttipanyalears et al. (2021) identified 187 genes significantly
upregulated in HCC across both GSE58208 and GSE49515, and validated
three (FLNA, CLU, CAP1) by qRT-PCR in a large independent cohort (153
HCC, 76 healthy), reporting 86.2% accuracy for the FLNA + CLU
combination. We checked whether FLNA, CLU, and CAP1 showed the same
direction of change (up in HCC) in our own re-processing of both
public datasets, then trained a simple logistic regression classifier
using only FLNA and CLU expression, evaluated with the same
cross-dataset procedure as above.

## 4. Results

**Within-dataset performance.** The 50-gene random forest achieved
ROC-AUC = 1.00 (cross-validated) within GSE58208, with a permutation
test p-value of 0.030, suggesting the score was unlikely to be due to
chance *within this dataset*.

**Cross-dataset performance.** The same model, trained on GSE58208
and tested on GSE49515, achieved ROC-AUC = 0.165. Training in the
reverse direction (GSE49515 → GSE58208) gave ROC-AUC = 0.180. Both
values are well below 0.5, indicating the model's predictions were
systematically *inverted* relative to the true labels in the unseen
dataset — a pattern consistent with the classifier having learned
dataset-specific (batch) structure rather than a disease signal that
holds across studies.

**After ComBat correction.** Applying ComBat batch correction before
retraining did not improve cross-dataset performance (ROC-AUC =
0.110), suggesting the mismatch was not simply an additive or
multiplicative scaling difference between batches that ComBat's
model could remove.

**Published-marker replication.** FLNA, CLU, and CAP1 were all
upregulated in HCC relative to healthy controls in both GSE58208 and
GSE49515 independently, consistent with the direction reported in the
published validation cohort. A logistic regression classifier using
only FLNA and CLU achieved ROC-AUC = 0.940 (GSE58208 → GSE49515) and
1.000 (GSE49515 → GSE58208) under the same cross-dataset evaluation
that the 50-gene model failed.

## 5. Discussion

The central finding of this project is methodological as much as
biological: a classifier that appears excellent by within-dataset
cross-validation (ROC-AUC = 1.00) can fail completely — including
in a directionally inverted way — when evaluated on an independent
dataset. This is consistent with prior reports of batch effects
undermining blood-based gene expression classifiers in other
diseases (e.g., lupus; not part of this project's original data, but
noted in the broader literature as a similar cautionary example).

By contrast, restricting the model to two genes with independent,
large-cohort clinical validation (FLNA, CLU) produced a classifier
that generalized well across both public datasets. This suggests
that, at least for these small public datasets, biologically
plausible, externally-validated markers were a more robust foundation
than an automatically selected, larger gene panel.

## 6. Limitations

- Sample sizes are small (15–20 samples per cross-dataset test),
  so confidence intervals on these AUC estimates are wide; single
  additional or removed samples could shift results noticeably.
- Both datasets are from non-Mongolian cohorts. No Mongolian patient
  data was used in this project. Any future clinical relevance to
  Mongolia's HCC burden would require local validation.
- This project used only two public datasets; the batch-effect vs.
  generalizable-signal distinction would benefit from a third,
  independent dataset (in progress).
- This is exploratory computational work, not a validated diagnostic
  method, and should not be interpreted as clinical guidance.

## 7. Data and code availability

All code is available at [GitHub link — add once published]. Data
are publicly available from NCBI GEO under accessions GSE58208 and
GSE49515.

## References

1. Puttipanyalears, C. *et al.* The expression of FLNA and CLU in
   PBMCs as a novel screening marker for hepatocellular carcinoma.
   *Sci Rep* **11**, 14969 (2021).
   https://doi.org/10.1038/s41598-021-94330-1
2. Leek, J.T. *et al.* Tackling the widespread and critical impact of
   batch effects in high-throughput data. *Nat Rev Genet* **11**,
   733–739 (2010).
3. Johnson, W.E., Li, C. & Rabinovic, A. Adjusting batch effects in
   microarray expression data using empirical Bayes methods.
   *Biostatistics* **8**, 118–127 (2007).
4. Hui, K.M. Gene expression profiling of PBMC from normal
   individuals, chronic hepatitis B carriers and hepatocellular
   carcinoma patients. NCBI GEO accession GSE58208 (2014).
5. Hui, K.M. Expression profiling of PBMC from patients with
   hepatocellular carcinoma. NCBI GEO accession GSE49515 (2013).
6. GBD 2019 Cancer Collaboration and related Global Burden of Disease
   Study 2019 results for Mongolia (liver cancer mortality and
   burden). Institute for Health Metrics and Evaluation.
