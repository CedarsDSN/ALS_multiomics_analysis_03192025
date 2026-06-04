# Batch 03: Axis scoring and sporadic ALS subtype discovery

## Purpose

This batch documents the transcriptomic axis-scoring and sporadic ALS subtyping workflow. Axis scoring computes Hallmark-derived metabolic, inflammation, and cell-cycle signatures. Sporadic subtyping uses those axis scores to identify stable sporadic ALS molecular subtypes, followed by subtype DE, pathway enrichment, visualization, clinical association, and sensitivity checks.

## Scripts used

```text
scripts/documented/transcriptomics/AxisScoring.R
scripts/refactored/transcriptomics/03_axis_scoring.R
scripts/documented/subtyping/sporadic_subtyping.R
scripts/refactored/subtyping/01_sporadic_subtyping.R
```

## Scientific questions

```text
1. Can ALS samples be represented by interpretable Hallmark-derived axes?
2. Are metabolic, inflammatory, or cell-cycle axes associated with mutation subtype, NfL, disease duration, or ALSFRS-R progression?
3. Do sporadic ALS samples cluster into reproducible molecular subtypes based on these axes?
4. What genes and pathways distinguish the sporadic subtypes?
5. Are subtype signals robust to outlier, differentiation-quality, glial, EMT, or marker-driven confounding?
```

## Inputs

Local-only inputs for axis scoring:

```text
transcriptomics/4_matrix/AnswerALS-651-T-v1-release6_raw-counts.csv
transcriptomics/4_matrix/Sample_Mapping_Information/Sample_Mapping_File_Feb_2024.csv
metadata/clinical/aals_dataportal_datatable.csv
ALS_Mutation_Subtypes_with_ALS_Control.csv
answerals/NfL/NfL_Seimens_Atellica_Immunoassay_Proteomics_2022_07.csv
MSigDB Hallmark gene sets through msigdbr
```

Primary objects consumed by subtyping:

```text
df_adj
axis_scores
df_nfl_uniq
expr_symbol
counts_filtered
annot_als_filtered
clean_svs
axis_sets or axis_sets_clean
```

## Outputs

Axis scoring outputs:

```text
axis_boxplots_adjusted.pdf
metabolic_axis_vs_progression.pdf
metabolic_quartiles_vs_progression.pdf
metabolic_partial_effect_vs_progression_annotated.pdf
axis_metabolic_sex.pdf
axis_covariate_screen_continuous_ALS_only.csv
axis_covariate_screen_clinical_ALS_only.csv
axis_covariate_screen_marker_ALS_only.csv
axis_covariate_screen_top_hits_table.csv
axis_covariate_screen_interactions_ALS_only.csv
```

Sporadic subtyping outputs:

```text
sporadic_subtype_heatmap.pdf
sporadic_umap.pdf
sporadic_subtype_heatmap_FINAL.pdf
pca_ssgsea_subtypes.pdf
pca_sals2_internal.pdf
degs_sporadic_subtype1_vs_2.csv
degs_sporadic_subtype1_vs_3.csv
degs_sporadic_subtype2_vs_3.csv
nfl_by_sporadic_subtype.pdf
alsfrs_progression_by_subtype.pdf
nfl_by_subtype_FINAL.pdf
nfl_by_subtype_FINAL_v2.pdf
sporadic_subtype_pathway_enrichment.csv
sporadic_subtype_pathway_enrichment_clean.csv
sporadic_subtype_pathway_enrichment_FINAL.csv
sporadic_subtype_pathway_enrichment_EMTcorrected.csv
sporadic_subtype_pathway_enrichment_FINAL_scaled.csv
deseq2_downsample_sensitivity.csv
```

Recommended handoff objects:

```text
axis_scores.rds
df_adj.rds
df_nfl_uniq.rds
expr_symbol.rds
expr_clean_mat.rds
counts_filtered.rds
annot_als_filtered.rds
svobj.rds
clean_svs.rds
axis_sets.rds
ax_clean_df.rds
ax_clean_dedup.rds
sporadic_subtype_assignments.csv
axis_scores_final.rds
```

## Upstream dependencies

`AxisScoring.R` is an upstream producer for many downstream analyses. `sporadic_subtyping.R` should run after `AxisScoring.R` or after loading its saved object bundle.

## How to run

Axis scoring:

```bash
Rscript scripts/refactored/transcriptomics/03_axis_scoring.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_03_axis_scoring_sporadic_subtyping/results/axis_scoring
```

Sporadic subtyping:

```bash
Rscript scripts/refactored/subtyping/01_sporadic_subtyping.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_03_axis_scoring_sporadic_subtyping/results/sporadic_subtyping \
  --objects analysis/batch_03_axis_scoring_sporadic_subtyping/objects/03_axis_scoring_objects.rds
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
Number of ALS samples:
Number of sporadic ALS samples:
Number of controls:
Genes used in each axis:
Subtype cluster counts:
Consensus clustering settings:
PAC/ICL or stability metrics:
UMAP/PCA separation:
Outliers excluded:
```

## Main findings

Add final subtype summary here after local execution.

Suggested summary format:

```text
Subtype label:
Dominant axis pattern:
Top DE genes/pathways:
Clinical/biomarker associations:
Robustness notes:
```

## Limitations

Potential limitations:

```text
axis definitions depend on Hallmark gene-set mapping
subtype labels depend on consensus clustering choices
sample-level clinical metadata completeness
hidden differentiation or cell-composition effects
object handoff required for downstream scripts
```
