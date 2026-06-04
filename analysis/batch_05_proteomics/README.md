# Batch 05: Proteomics normalization, DE, batch correction, and sensitivity analysis

## Purpose

This batch documents proteomics workflows, including normalization-method comparison, protein differential abundance, batch-correction comparison, and sensitivity analysis comparing C9orf72 vs control results before and after batch correction.

## Scripts used

```text
scripts/documented/proteomics/CompareProteomicNormalization.R
scripts/refactored/proteomics/01_compare_normalization.R
scripts/documented/proteomics/proteomicsDE.R
scripts/refactored/proteomics/02_proteomics_de.R
scripts/documented/proteomics/CompareBatchCorrectionProteomics.R
scripts/refactored/proteomics/03_compare_batch_correction.R
scripts/documented/proteomics/sensitivity_analysis_prebatch.py
scripts/refactored/proteomics/04_sensitivity_analysis_prebatch.py
```

## Scientific questions

```text
1. Which proteomic normalization method gives the most interpretable sample structure?
2. What proteins differ by ALS mutation/subtype or case-control group?
3. How sensitive are C9orf72 vs control protein findings to batch-correction strategy?
4. Which proteins remain significant before and after batch correction?
5. Do proteomic signatures align with transcriptomic or pathway-level signatures?
```

## Inputs

Local-only inputs:

```text
AnswerALS proteomics protein matrices
pre-batch-corrected protein matrix
post-batch-corrected/imputed protein matrix
protein-to-gene annotation files
mutation subtype table
clinical/sample metadata
```

Common matrix paths used in the original project included:

```text
proteomics/4_matrix/AnswerALS-447-P_proteomics-protein-matrix_prebatchcorrected.txt
proteomics/4_matrix/AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt
ALS_Mutation_Subtypes_with_ALS_Control.csv
```

## Outputs

Normalization comparison outputs:

```text
PCA plots for log2, quantile, z-score, and median-centered normalization
normalization comparison notes/tables
```

Proteomics DE outputs:

```text
participant-level collapsed protein matrix
limma contrast tables
Hallmark fgsea tables
PCA plots
targeted protein plots
protein-to-gene mapped result tables
```

Batch correction outputs:

```text
pre-batch vs post-batch comparisons
shipment-adjusted model results
ComBat strategy results
balanced-shipment sensitivity outputs
batch correction comparison plots/tables
```

Pre/post batch sensitivity outputs:

```text
DE_results_C9orf72_vs_Control_PREBATCH.tsv
DE_results_C9orf72_vs_Control_POSTBATCH.tsv
comparison of significant proteins at FDR < 0.05 and FDR < 0.01
lists of proteins lost/gained after correction
```

## Upstream dependencies

This batch is mostly independent of transcriptomic object handoffs, but depends on consistent mutation subtype/sample metadata. Cross-omics summaries may use outputs from Batch 02 or Batch 03.

## How to run

Normalization comparison:

```bash
Rscript scripts/refactored/proteomics/01_compare_normalization.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_05_proteomics/results/normalization
```

Proteomics DE:

```bash
Rscript scripts/refactored/proteomics/02_proteomics_de.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_05_proteomics/results/proteomics_de
```

Batch-correction comparison:

```bash
Rscript scripts/refactored/proteomics/03_compare_batch_correction.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_05_proteomics/results/batch_correction
```

Pre/post batch Python sensitivity analysis:

```bash
python scripts/refactored/proteomics/04_sensitivity_analysis_prebatch.py \
  --prebatch /path/to/private/ALS/proteomics/4_matrix/AnswerALS-447-P_proteomics-protein-matrix_prebatchcorrected.txt \
  --postbatch /path/to/private/ALS/proteomics/4_matrix/AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt \
  --subtypes /path/to/private/ALS/ALS_Mutation_Subtypes_with_ALS_Control.csv \
  --outdir analysis/batch_05_proteomics/results/prebatch_sensitivity
```

## Run record

```text
Date run:
Machine/HPC node:
Git branch:
Git commit:
R/Python version:
Package versions:
Run completed: yes/no
```

## Quality-control notes

Record:

```text
Number of proteins in each matrix:
Number of samples in each matrix:
C9orf72 sample count:
Control sample count:
Missingness filtering:
Normalization method selected:
Batch variables tested:
Shipment/batch balance:
```

## Main findings

Add proteomics and batch-sensitivity findings after local execution.

## Limitations

Potential limitations:

```text
pre- and post-batch matrices may contain different sample/protein sets
protein IDs may require careful gene-symbol mapping
batch-correction can remove biological signal if confounded with group
simple t-tests are sensitivity analyses and not replacements for full model-based DE
```
