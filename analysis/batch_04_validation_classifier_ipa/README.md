# Batch 04: Subtype validation, LASSO classifier, and IPA preparation

## Purpose

This batch documents validation and follow-up analyses for the sporadic ALS subtype framework. It includes subtype sensitivity analyses, sparse LASSO subtype classifier development, and generation of IPA-ready differential-expression tables.

## Scripts used

```text
scripts/documented/validation/SensitivityAnalysis.R
scripts/refactored/validation/01_sensitivity_analysis.R
scripts/documented/validation/LassoClassifier.R
scripts/refactored/validation/02_lasso_classifier.R
scripts/documented/validation/IPAPrepare.R
scripts/refactored/validation/03_ipa_prepare.R
```

## Scientific questions

```text
1. Are subtype or axis findings driven by differentiation quality, glial contamination, proliferation, or outliers?
2. Can sparse Hallmark-axis gene panels classify sporadic ALS molecular subtypes?
3. Which genes contribute most to one-vs-rest subtype classification?
4. Are key signals such as CCT2 robust across samples rather than driven by a few outliers?
5. What input table should be passed to Ingenuity Pathway Analysis for S-ALS-3 vs S-ALS-2?
```

## Inputs

Primary upstream objects from Batch 03:

```text
expr_clean_mat
ax_clean_dedup
axis_scores
df_adj
df_nfl_uniq
sporadic_final_check
sporadic_subtype_assignments
res_final_12
res_final_13
res_final_23
```

Additional local-only metadata may include:

```text
clinical metadata with ALSFRS-R progression
NfL measurements
PRS scores
CCT2 expression matrix or expression object
subtype DESeq2 result objects
```

## Outputs

Sensitivity-analysis outputs:

```text
diff_quality_by_subtype.pdf
diff_quality_sensitivity_results.csv
sensitivity plots and subtype confounder summaries
CCT2/AXIS_METABOLIC/PRS/NfL/ALSFRS-R association summaries
```

LASSO classifier outputs:

```text
one-vs-rest subtype classifier performance tables
cross-validation metrics
sparse gene panels
classifier coefficient tables
CCT2 follow-up plots/tables
```

IPA preparation outputs:

```text
IPA-ready S-ALS-3 vs S-ALS-2 DESeq2 table
CCT2 outlier/uniformity checks
ranked gene lists for pathway upload
```

## Upstream dependencies

This batch should run after Batch 03. It depends on axis-scoring and sporadic-subtyping outputs.

## How to run

Sensitivity analysis:

```bash
Rscript scripts/refactored/validation/01_sensitivity_analysis.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_04_validation_classifier_ipa/results/sensitivity_analysis \
  --objects outputs/objects/03_axis_scoring_objects.rds,outputs/objects/01_sporadic_subtyping_objects.rds
```

LASSO classifier:

```bash
Rscript scripts/refactored/validation/02_lasso_classifier.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_04_validation_classifier_ipa/results/lasso_classifier \
  --objects outputs/objects/03_axis_scoring_objects.rds,outputs/objects/01_sporadic_subtyping_objects.rds
```

IPA preparation:

```bash
Rscript scripts/refactored/validation/03_ipa_prepare.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_04_validation_classifier_ipa/results/ipa_prepare \
  --objects outputs/objects/03_axis_scoring_objects.rds,outputs/objects/01_sporadic_subtyping_objects.rds
```

## Run record

```text
Date run:
Machine/HPC node:
Git branch:
Git commit:
R version:
Run completed: yes/no
```

## Quality-control notes

Record:

```text
Confounders tested:
Subtype sample counts after filtering:
Classifier folds/repeats:
Nested CV used: yes/no
Model family and alpha/lambda settings:
Outlier criteria:
IPA upload table dimensions:
```

## Main findings

Add validation and classifier summary after local execution.

## Limitations

Potential limitations:

```text
small subtype-specific sample sizes
risk of classifier overfitting
feature correlation within Hallmark gene sets
outlier sensitivity for individual marker genes
IPA database version and settings
```
