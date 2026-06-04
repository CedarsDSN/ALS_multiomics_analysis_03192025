# Batch 01: Transcriptomics progression-class RNA-seq analysis

## Purpose

This batch documents the initial RNA-seq quality-control and differential-expression workflow comparing ALS progression classes. The analysis evaluates fast, intermediate, and slow progressors, generates exploratory QC plots, and creates pairwise DESeq2 result tables.

## Scripts used

```text
scripts/documented/transcriptomics/run_RNAseq_DE.R
scripts/refactored/transcriptomics/01_progression_rnaseq_de.R
```

## Scientific questions

```text
1. Do ALS progression groups show detectable transcriptomic differences?
2. Are those differences robust after adjustment for race or other clinical covariates?
3. Are there PCA, sample-distance, or lower-PC patterns suggesting batch or sample-quality structure?
```

## Inputs

These files should be available locally or on the HPC/Box data system, but should not be committed to GitHub.

```text
transcriptomics/4_matrix/AnswerALS-651-T-v1-release6_raw-counts.csv
transcriptomics/4_matrix/Sample_Mapping_Information/Sample_Mapping_File_Feb_2024.csv
ALS_Progression_Rates_Full.csv
metadata/aals_dataportal_datatable.csv
```

## Outputs

Expected output files include:

```text
DGE_Fast_vs_Slow_with_race.csv
DGE_Intermediate_vs_Slow_with_race.csv
DGE_Fast_vs_Intermediate_with_race.csv
DGE_Fast_vs_Slow_dispersion_model.csv
DGE_Intermediate_vs_Slow_dispersion_model.csv
DGE_Fast_vs_Intermediate_dispersion_model.csv
PCA_ProgressorClass.pdf
Sample_Distance_Heatmap.pdf
Top20_Loadings_<PC>.csv
Lower_PC_Plots_Top500.pdf
Batch_visualization.pdf
```

Generated outputs should be written under:

```text
analysis/batch_01_transcriptomics_progression/results/
```

This directory is ignored by Git by default. Only non-sensitive summary notes should be committed.

## Upstream dependencies

None. This is an entry-point transcriptomics analysis, assuming the raw count matrix and metadata files are present.

## How to run

Documented script:

```bash
Rscript scripts/documented/transcriptomics/run_RNAseq_DE.R
```

Refactored script:

```bash
Rscript scripts/refactored/transcriptomics/01_progression_rnaseq_de.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_01_transcriptomics_progression/results/progression_rnaseq_de
```

## Run record

```text
Date run:
Machine/HPC node:
Git branch:
Git commit:
R version:
Key package versions:
Run completed: yes/no
```

## Quality-control notes

Record the following after running:

```text
Number of samples retained:
Number of genes retained after filtering:
Progression classes represented:
PCA findings:
Sample-distance outliers:
Batch or covariate concerns:
```

## Main findings

Add final summarized findings here after local execution. Do not paste patient-level data.

## Limitations

Potential limitations to track:

```text
progression-class definitions
sample-size imbalance across groups
race/sex/age covariate completeness
batch effects or hidden structure
multiple-testing burden
```
