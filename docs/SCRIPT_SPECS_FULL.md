# Script specifications: purpose, inputs, outputs, dependencies, run commands

This is the primary GitHub documentation file for the ALS project scripts.

## `run_RNAseq_DE.R`

**Category:** Transcriptomics / QC / DE

**Purpose**

Initial RNA-seq differential-expression and QC workflow. Compares ALS progression classes using DESeq2, including a race-adjusted model and a simpler progression-only model, then produces PCA, sample-distance, lower-PC, and batch-visualization plots.

**Inputs**

- `transcriptomics/4_matrix/AnswerALS-651-T-v1-release6_raw-counts.csv` — raw gene counts with `Geneid` column.
- `transcriptomics/4_matrix/Sample_Mapping_Information/Sample_Mapping_File_Feb_2024.csv` — sample-to-participant map.
- `ALS_Progression_Rates_Full.csv` — progression class table.
- `metadata/aals_dataportal_datatable.csv` — clinical metadata.

**Outputs**

- `DGE_Fast_vs_Slow_with_race.csv`
- `DGE_Intermediate_vs_Slow_with_race.csv`
- `DGE_Fast_vs_Intermediate_with_race.csv`
- `DGE_Fast_vs_Slow_dispersion_model.csv`
- `DGE_Intermediate_vs_Slow_dispersion_model.csv`
- `DGE_Fast_vs_Intermediate_dispersion_model.csv`
- `PCA_ProgressorClass.pdf`
- `Sample_Distance_Heatmap.pdf`
- `Top20_Loadings_<PC>.csv`
- `Lower_PC_Plots_Top500.pdf`
- `Batch_visualization.pdf`

**Upstream dependencies**

- None. This is an entry-point transcriptomics script, assuming all input files are present.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/transcriptomics/run_RNAseq_DE.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/transcriptomics/01_progression_rnaseq_de.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/transcriptomics/01_progression_rnaseq_de
```

**Refactor notes**

- Replace `setwd(mainDir)` with config-driven paths.
- Write outputs under `results/transcriptomics/progression_de/` and `figures/transcriptomics/qc/`.
- Save the DESeq2 dataset and VST object if downstream QC is reused.

## `C9SOD1Convergence.R`

**Category:** Transcriptomics / familial ALS convergence

**Purpose**

Familial ALS transcriptomic analysis focused on C9orf72 and SOD1. Performs count filtering, VST/PCA QC, surrogate-variable analysis, model comparison, C9orf72-vs-control and SOD1-vs-control DESeq2 contrasts, convergence/divergence summaries, pathway enrichment, ALS gene-set guardrail scoring, and RNA/protein signature comparisons.

**Inputs**

- `transcriptomics/4_matrix/AnswerALS-651-T-v1-release6_raw-counts.csv`.
- `transcriptomics/4_matrix/Sample_Mapping_Information/Sample_Mapping_File_Feb_2024.csv`.
- `metadata/clinical/aals_dataportal_datatable.csv`.
- `ALS_Mutation_Subtypes_with_ALS_Control.csv`.
- `sod1_sample_mutations.csv` for SOD1 mutation labels.
- Optional external paper/supplement table read through `readxl::read_excel()` in later cross-reference sections.

**Outputs**

- SVA/model summaries: `sva_sensitivity_all_svs.csv`, `sva_filtered_summary.csv`, `model_comparison_simple_vs_filtered.csv`, `complete_three_model_comparison.csv`.
- Pathway outputs: `pathway_GO_BP_results.csv`, `pathway_KEGG_results.csv`, `pathway_Reactome_results.csv`, `pathway_shared_C9_SOD1.csv`.
- Familial signatures: `familial_signature_broad.csv`, `familial_signature_strict.csv`, `familial_signature_lenient.csv`, `familial_signature_strict_sva_markers.csv`.
- Guardrail and sporadic follow-up outputs: `sporadic_guardrail_High_vs_Low_DE.csv`, `GSEA_GO_BP_sporadic_guardrail.csv`, `sporadic_samples_vcf.txt`, `sporadic_samples_present.txt`.
- Cross-omics summaries: `rna_signature_summary.csv`, `protein_signature_summary.csv`, `cross_omics_feature.csv`, `cross_omics_pathway.csv`, `integrated_signature_results_matrix.csv`.
- Figures including dispersion plots, PCA plots, GSEA plots, mutation/sex/age/progression distributions, and C9/SOD1 divergence plots.

**Upstream dependencies**

- None for the first RNA-seq sections; later sections depend on objects created earlier in the same script.
- Downstream consumers may use `counts_filtered`, `annot_als_filtered`, `svobj`, `clean_svs`, and familial signature tables.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/transcriptomics/C9SOD1Convergence.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/transcriptomics/02_c9_sod1_convergence.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/transcriptomics/02_c9_sod1_convergence
```

**Refactor notes**

- Split into `01_load_qc.R`, `02_sva_model_selection.R`, `03_c9_sod1_de.R`, `04_convergence_metrics.R`, `05_pathway_analysis.R`, and `06_cross_omics_summary.R`.
- Save reusable intermediates as RDS: `counts_filtered.rds`, `annot_als_filtered.rds`, `svobj.rds`, `clean_svs.rds`, `expr_vst.rds`, `familial_signatures.rds`.
- Move all plots/tables into task-specific `results/` and `figures/` subfolders.

## `AxisScoring.R`

**Category:** Transcriptomics / axis scoring

**Purpose**

Computes Hallmark-derived metabolic, inflammation, and cell-cycle axis scores using transformed RNA-seq expression, removes mutation-confounded surrogate variables, and tests axis associations with mutation subtype, sporadic ALS, NfL, disease duration, ALSFRS-R progression, sex, and other clinical/marker covariates.

**Inputs**

- `transcriptomics/4_matrix/AnswerALS-651-T-v1-release6_raw-counts.csv`.
- `transcriptomics/4_matrix/Sample_Mapping_Information/Sample_Mapping_File_Feb_2024.csv`.
- `metadata/clinical/aals_dataportal_datatable.csv`.
- `ALS_Mutation_Subtypes_with_ALS_Control.csv`.
- `answerals/NfL/NfL_Seimens_Atellica_Immunoassay_Proteomics_2022_07.csv`.
- MSigDB Hallmark sets via `msigdbr`.

**Outputs**

- Figures: `axis_boxplots_adjusted.pdf`, `metabolic_axis_vs_progression.pdf`, `metabolic_quartiles_vs_progression.pdf`, `metabolic_partial_effect_vs_progression_annotated.pdf`, `axis_metabolic_sex.pdf`, plus one covariate-specific plot created through `outfile`.
- Tables: `axis_covariate_screen_continuous_ALS_only.csv`, `axis_covariate_screen_clinical_ALS_only.csv`, `axis_covariate_screen_marker_ALS_only.csv`, `axis_covariate_screen_top_hits_table.csv`, `axis_covariate_screen_interactions_ALS_only.csv`.
- Recommended new intermediates: `axis_scores.rds`, `df_adj.rds`, `df_nfl_uniq.rds`, `expr_symbol.rds`, `expr_clean_mat.rds`, `counts_filtered.rds`, `annot_als_filtered.rds`, `svobj.rds`, `clean_svs.rds`, `axis_sets.rds`.

**Upstream dependencies**

- None beyond data files. This should become the main upstream producer for sporadic subtyping, PRSet, eQTL, LASSO, IPA, and sensitivity scripts.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/transcriptomics/AxisScoring.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/transcriptomics/03_axis_scoring.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/transcriptomics/03_axis_scoring
```

**Refactor notes**

- Add explicit `saveRDS()` calls for all shared objects.
- Move ssGSEA axis gene definitions into `R/axis_definitions.R` so LASSO/genomics/PRS use exactly the same gene sets.
- Create one table `results/axis_scores/axis_scores_with_metadata.csv` for collaborators and one RDS for downstream code.

## `sporadic_subtyping.R`

**Category:** Sporadic ALS subtyping

**Purpose**

Defines transcriptomic subtypes among sporadic ALS samples from metabolic/inflammation/cell-cycle axis scores. Uses consensus clustering, assigns S-ALS subtype labels, generates heatmap/UMAP visualizations, runs pairwise subtype DESeq2 contrasts, performs pathway enrichment, evaluates NfL and ALSFRS-R progression by subtype, and runs EMT/S100B/downsampling sensitivity analyses.

**Inputs**

- Workspace objects from `AxisScoring.R`: `df_adj`, `axis_scores`, `df_nfl_uniq`, `expr_symbol`, `counts_filtered`, `annot_als_filtered`, `clean_svs`.
- Gene-set objects such as `axis_sets` or `axis_sets_clean` for marker-filtered re-scoring.
- Optional genomic burden object `burden_all` for exploratory genomic subtyping section.

**Outputs**

- Subtype figures: `sporadic_subtype_heatmap.pdf`, `sporadic_umap.pdf`, `sporadic_subtype_heatmap_FINAL.pdf`, `pca_ssgsea_subtypes.pdf`, `pca_sals2_internal.pdf`.
- DE tables: `degs_sporadic_subtype1_vs_2.csv`, `degs_sporadic_subtype1_vs_3.csv`, `degs_sporadic_subtype2_vs_3.csv`.
- Clinical plots: `nfl_by_sporadic_subtype.pdf`, `alsfrs_progression_by_subtype.pdf`, `nfl_by_subtype_FINAL.pdf`, `nfl_by_subtype_FINAL_v2.pdf`.
- Pathway tables: `sporadic_subtype_pathway_enrichment.csv`, `sporadic_subtype_pathway_enrichment_clean.csv`, `sporadic_subtype_pathway_enrichment_FINAL.csv`, `sporadic_subtype_pathway_enrichment_EMTcorrected.csv`, `sporadic_subtype_pathway_enrichment_FINAL_scaled.csv`.
- Sensitivity table: `deseq2_downsample_sensitivity.csv`.
- Recommended new intermediates: `ax_clean_df.rds`, `ax_clean_dedup.rds`, `sporadic_subtype_assignments.csv`, `expr_clean_mat.rds`, `axis_scores_final.rds`.

**Upstream dependencies**

- Must run after `AxisScoring.R` or load its saved intermediate RDS files.
- For DE sections, also needs the filtered count matrix and clean SV list from the same transcriptomics pipeline.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/subtyping/sporadic_subtyping.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/subtyping/01_sporadic_subtyping.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/subtyping/01_sporadic_subtyping --objects outputs/objects/axis_scoring_objects.rds,outputs/objects/sporadic_subtyping_objects.rds
```

**Refactor notes**

- Split into clustering, subtype DE, pathway enrichment, clinical associations, and sensitivity scripts.
- Make `k = 3`, PAC/ICL settings, and subtype label names explicit parameters.
- Save final subtype assignments before any sensitivity analyses so downstream scripts do not depend on a live workspace.

## `SensitivityAnalysis.R`

**Category:** Subtype validation / sensitivity

**Purpose**

Tests whether transcriptomic subtype and axis findings are artifacts of differentiation quality, glial contamination, proliferation, marker outliers, or related confounders. Also examines CCT2, PRS, NfL, AXIS_METABOLIC, and ALSFRS-R relationships.

**Inputs**

- Required base objects: `expr_clean_mat` and `ax_clean_dedup`.
- Later PRS/CCT2 sections require `analysis_df`, `expr_sym` or `expr_clean_mat`, and subtype/clinical metadata with `cluster_final`, axis scores, sex, NfL, ALSFRS-R slope, and progression variables.

**Outputs**

- `diff_quality_by_subtype.pdf`.
- `diff_quality_sensitivity_results.csv`.
- `sensitivity_summary.csv`.
- Additional scatter/forest/model outputs may be produced by later sections depending on available objects.

**Upstream dependencies**

- Must run after `AxisScoring.R` and `sporadic_subtyping.R`.
- PRS-specific checks require `PRSetTesting.R` to have produced `analysis_df`.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/validation/SensitivityAnalysis.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/validation/01_sensitivity_analysis.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/validation/01_sensitivity_analysis --objects outputs/objects/axis_scoring_objects.rds,outputs/objects/sporadic_subtyping_objects.rds
```

**Refactor notes**

- Split into `01_differentiation_quality.R`, `02_glial_residualization.R`, `03_cct2_checks.R`, `04_prs_cct2_nfl.R`, and `05_axis_clinical_models.R`.
- Save each model summary as a tidy CSV instead of only printing to console.
- Define a single metadata contract for NfL, ALSFRS slope, age, sex, and subtype columns.

## `LassoClassifier.R`

**Category:** Subtype classifier

**Purpose**

Builds one-vs-rest LASSO classifiers for S-ALS-1, S-ALS-2, and S-ALS-3 using Hallmark axis genes. Runs nested cross-validation, extracts sparse gene panels, performs S-ALS-3 stability selection, and checks CCT2 subtype/NfL relationships.

**Inputs**

- `expr_clean_mat` — symbol-keyed expression matrix, genes x samples.
- `ax_clean_dedup` — deduplicated subtype table with `Participant_ID` and `cluster_final`.
- MSigDB Hallmark gene sets via `msigdbr`; should use the exact same axis sets as `AxisScoring.R`.

**Outputs**

- `CCT2_by_subtype.pdf`.
- `CCT2_vs_NfL.pdf`.
- `S-ALS-3_stability_selection.csv`.
- `panel_axis_composition.csv`.
- `lasso_results.rds`.

**Upstream dependencies**

- Must run after `AxisScoring.R` and `sporadic_subtyping.R`, or after loading `expr_clean_mat.rds` and `ax_clean_dedup.rds`.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/validation/LassoClassifier.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/validation/02_lasso_classifier.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/validation/02_lasso_classifier --objects outputs/objects/axis_scoring_objects.rds,outputs/objects/sporadic_subtyping_objects.rds
```

**Refactor notes**

- Move axis gene definitions into one shared file and import them here.
- Set and record random seeds, fold IDs, lambda settings, and training sample IDs.
- Write model performance metrics into a CSV, not only an RDS.

## `IPAPrepare.R`

**Category:** IPA export / subtype DE

**Purpose**

Creates an IPA-ready DESeq2 result table for the S-ALS-3 vs S-ALS-2 contrast using clean SVs and sex adjustment. Also checks whether CCT2 elevation is uniform across S-ALS-3 or driven by outlier participants.

**Inputs**

- Workspace objects: `counts_filtered`, `annot_als_filtered`, `svobj`, `clean_svs`.
- Subtype object: one of `ax_clean_dedup`, `ax_clean_df`, `sporadic_meta`, or `sporadic_final_check`, with a cluster column such as `cluster_final`.

**Outputs**

- `ipa_inputs/DEG_S3_vs_S2_IPA_upload.tsv`.
- `ipa_inputs/DEG_S3_vs_S2_full.tsv`.
- `ipa_inputs/DEG_S3_vs_S2_full.rds`.
- `cct2_distribution_check.pdf`.

**Upstream dependencies**

- Must run after `AxisScoring.R` and `sporadic_subtyping.R`.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/validation/IPAPrepare.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/validation/03_ipa_prepare.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/validation/03_ipa_prepare --objects outputs/objects/axis_scoring_objects.rds,outputs/objects/sporadic_subtyping_objects.rds
```

**Refactor notes**

- Convert workspace object discovery into explicit RDS arguments.
- Keep `ipa_inputs/` under `results/ipa_inputs/`.
- Add an option to choose contrast, e.g. `--case S-ALS-3 --reference S-ALS-2`.

## `CompareProteomicNormalization.R`

**Category:** Proteomics / normalization QC

**Purpose**

Compares proteomics normalization strategies: log2(x+1), quantile normalization, row-wise z-score normalization, and median-centering. Produces normalized matrices and a PCA comparison PDF.

**Inputs**

- `proteomics/4_matrix/AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt`.

**Outputs**

- `Proteomics_Log2Plus1.csv`.
- `Proteomics_QuantileNormalized.csv`.
- `Proteomics_ZscoreNormalized.csv`.
- `Proteomics_MedianCentered.csv`.
- `PCA_Normalization_Comparison.pdf`.

**Upstream dependencies**

- None. This is a QC/preprocessing comparison script.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/proteomics/CompareProteomicNormalization.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/proteomics/01_compare_normalization.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/proteomics/01_compare_normalization
```

**Refactor notes**

- Parameterize the input matrix and output directory.
- Use a consistent sample/protein orientation check and write orientation to log.
- Move outputs to `results/proteomics/normalization_qc/`.

## `proteomicsDE.R`

**Category:** Proteomics / differential abundance

**Purpose**

End-to-end proteomics differential abundance workflow. Loads corrected/imputed protein matrix, collapses technical/repeated proteomics runs to participant-level medians, builds mutation subtype metadata, normalizes abundance values, maps proteins to gene symbols, computes proteomics guardrail score, runs limma contrasts, runs Hallmark fgsea, and generates PCA/diagnostic plots.

**Inputs**

- `AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt`.
- `Sample_Mapping_Information/Sample_Mapping_File_Feb_2024.csv` in the proteomics matrix directory.
- `metadata/clinical/aals_dataportal_datatable.csv`.
- `ALS_Mutation_Subtypes_with_ALS_Control.csv`.
- Transcriptomics sample mapping and raw counts for sample filtering.
- `familial_signature_strict.csv` with columns `ensembl_id` and `direction`.

**Outputs**

- Under `outputs/proteomics_guardrail/`: `proteomics_guardrail_scores.tsv`, `guardrail_score_boxplot.png`, `limma_C9orf72_vs_Control.tsv`, `limma_SOD1_vs_Control.tsv`, `limma_Sporadic_vs_Control.tsv`, `fgsea_hallmark_all.tsv`, `fgsea_hallmark_top20_per_comparison.tsv`, `proteomics_pcs.tsv`, `proteomics_PCA_PC1_PC2.png`, `proteomics_PCA_PC1_PC2_colored_by_batch.png`, `proteomics_PCA_PC1_PC2_after_reg_sex.png`, `CBX1_protein_abundance_boxplot_C9_vs_Control.pdf`.

**Upstream dependencies**

- Uses familial guardrail signature from `C9SOD1Convergence.R`.
- Uses mutation subtype and clinical metadata shared with transcriptomics scripts.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/proteomics/proteomicsDE.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/proteomics/02_proteomics_de.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/proteomics/02_proteomics_de
```

**Refactor notes**

- Remove or separate the initial exploratory block before the `#!/usr/bin/env Rscript` section.
- Replace `MulticoreParam(20)` with a config parameter.
- Cache UniProt/Ensembl mapping tables so biomaRt failures do not break reruns.

## `CompareBatchCorrectionProteomics.R`

**Category:** Proteomics / batch correction sensitivity

**Purpose**

Compares C9orf72-vs-control proteomics results across pre-batch and post-batch correction strategies. Maps protein IDs to genes, collapses to gene-level values, runs limma variants with/without shipment adjustment or ComBat, compares t-statistic rank correlations, and runs Hallmark camera pathway analysis.

**Inputs**

- `AnswerALS-447-P_proteomics-protein-matrix_prebatchcorrected.txt`.
- `AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt`.
- `metadata/clinical/aals_dataportal_datatable.csv`.
- `ALS_Mutation_Subtypes_with_ALS_Control.csv`.
- MSigDB Hallmark sets via `msigdbr`.

**Outputs**

- Folder `C9vsCtrl_proteoDA_genelevel_outputs/`: `protein_to_gene_mapping_prebatch.csv`, `protein_to_gene_mapping_postbatch.csv`, `limma_genelevel_<variant>.csv`, `summary_genelevel_across_variants.csv`, `tstat_spearman_corr_genelevel.csv`, `tstat_spearman_corr_genelevel_heatmap.png`, `camera_hallmark_genelevel_across_variants.csv`, `camera_hallmark_genelevel_top20_per_variant.csv`.

**Upstream dependencies**

- No required upstream script if all data files exist. It reuses the same metadata conventions as proteomicsDE.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/proteomics/CompareBatchCorrectionProteomics.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/proteomics/03_compare_batch_correction.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/proteomics/03_compare_batch_correction
```

**Refactor notes**

- Make shipment variable and C9/control group definitions explicit.
- Remove dead/commented RNA count code unless it is truly needed.
- Save the aligned participant metadata used for each variant.

## `sensitivity_analysis_prebatch.py`

**Category:** Proteomics / batch correction sensitivity

**Purpose**

Runs a simple protein-level C9orf72-vs-control sensitivity analysis on pre-batch-corrected and post-batch-corrected proteomics matrices, using independent t-tests per protein, FDR correction, and overlap summaries.

**Inputs**

- Pre-batch proteomics matrix: `AnswerALS-447-P_proteomics-protein-matrix_prebatchcorrected.txt`.
- Post-batch proteomics matrix: `AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt`.
- `ALS_Mutation_Subtypes_with_ALS_Control.csv`.

**Outputs**

- `DE_results_C9orf72_vs_Control_PREBATCH.tsv`.
- `DE_results_C9orf72_vs_Control_POSTBATCH.tsv`.
- Console summary of sample sizes, significant protein counts, top proteins, and proteins gained/lost after correction.

**Upstream dependencies**

- No upstream script required beyond proteomics matrices and mutation metadata.

**How to run**

Legacy documented version:

```bash
python scripts/documented/proteomics/sensitivity_analysis_prebatch.py
```

Refactored command-line version:

```bash
python scripts/refactored/proteomics/04_sensitivity_analysis_prebatch.py \
  --prebatch /path/to/prebatch.txt \
  --postbatch /path/to/postbatch.txt \
  --subtypes /path/to/ALS_Mutation_Subtypes_with_ALS_Control.csv \
  --outdir results/proteomics/prebatch_sensitivity
```

**Refactor notes**

- Use `argparse` instead of hard-coded paths.
- Write a comparison summary TSV in addition to console output.
- Optionally replace t-tests with limma-style modeling if covariate adjustment is required.

## `MAF_allele_frequency.R`

**Category:** Genomics / QC

**Purpose**

Summarizes minor allele frequency distributions from a GATK joint VCF-derived allele-frequency file. Computes MAF from AF, reports allele-number and frequency-tier summaries, and produces site-frequency spectrum plots.

**Inputs**

- `all_chr_AF.txt` with columns `CHROM`, `POS`, `REF`, `ALT`, `AF`, and `AN`.

**Outputs**

- `MAF_full_SFS.pdf`.
- `MAF_rare_zoom.pdf`.
- `MAF_tier_composition.pdf`.

**Upstream dependencies**

- Requires prior extraction of AF/AN fields from the joint VCF. No R upstream dependency.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/genomics/MAF_allele_frequency.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/genomics/01_maf_allele_frequency.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/genomics/01_maf_allele_frequency
```

**Refactor notes**

- Parameterize the AF file and MAF cutoff.
- Do not hard-code plot colors if journal/slide style will be controlled centrally.
- Write frequency-tier counts to CSV for reproducible reporting.

## `genomic_analysis.R`

**Category:** Genomics / rare variant burden

**Purpose**

Parses ANNOVAR-style exon annotation files, creates sample-level rare damaging/LoF burden metrics, compares burden across sporadic ALS subtypes, associates LoF burden with NfL, builds gene-level burden matrices, and explores ALS-panel burden/clustering.

**Inputs**

- Exon annotation files such as `ALS.chr2.anno_report_exon.txt` and corresponding files for all chromosomes.
- Workspace object `sporadic_meta` with sample/participant/CGND IDs and subtype labels.
- NfL and ALSFRS fields in merged subtype/clinical metadata for downstream associations.

**Outputs**

- `genomic_burden_per_sample_all_chrom.csv`.
- `genomic_burden_by_subtype.pdf`.
- `nfl_vs_lof_burden_scatter.pdf`.
- Important in-memory objects currently not saved consistently: `rare_burden_all`, `rare_burden_sub`, `gene_burden_combined`, `gene_burden_filt`, `col_to_id`, `use_cols`, `crosswalk`.

**Upstream dependencies**

- Run after `sporadic_subtyping.R` has produced sporadic subtype metadata.
- Requires WGS/ANNOVAR preprocessing outside R.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/genomics/genomic_analysis.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/genomics/02_genomic_analysis.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/genomics/02_genomic_analysis
```

**Refactor notes**

- Separate sample-level burden, gene-level burden, ALS-panel analysis, and clinical association into separate scripts.
- Save burden matrices and crosswalk objects as RDS for TRiC/MYC/DEG burden scripts.
- Make variant class names and MAF column names configurable because ANNOVAR headers vary.

## `genomic_stratification.R`

**Category:** Genomics / network stratification

**Purpose**

Builds a genotype-derived stratification model from rare damaging variants in transcriptomic-axis genes and ALS-panel genes. Constructs patient-by-gene burden matrices, builds PPI/coexpression networks, scores edge conservation/perturbation, and clusters patients using SNF.

**Inputs**

- ANNOVAR exon files for all chromosomes.
- MSigDB Hallmark gene sets and curated ALS panel genes.
- Expression matrix from the transcriptomics axis pipeline for motor-neuron coexpression.
- Subtype/crosswalk metadata from sporadic transcriptomic workflow.
- STRINGdb data access for PPI network.

**Outputs**

- `gene_tier_membership.csv`.
- `patient_gene_matrix_sporadic.rds` and `.csv`.
- `genes_with_variants.rds`.
- `mapped_df.rds`, `interactions_genes.rds`, `coexpr_edges_df.rds`, `validated_edges.rds`.
- `als_mn_network_final.rds`, `patient_gene_matrix_network.rds`, `als_mn_network_edges.csv`, `als_mn_network_nodes.csv`.
- `edge_conservation_matrix.rds`.

**Upstream dependencies**

- Requires WGS annotations, axis gene definitions, expression matrix, and sporadic subtype metadata.
- Can consume burden/crosswalk objects from `genomic_analysis.R` after refactoring.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/genomics/genomic_stratification.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/genomics/03_genomic_stratification.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/genomics/03_genomic_stratification
```

**Refactor notes**

- Split into `01_define_gene_universe.R`, `02_extract_target_variants.R`, `03_build_network.R`, `04_score_edges.R`, `05_snf_cluster.R`.
- Cache STRING mappings and coexpression edges.
- Add deterministic seeds for SNF/spectral clustering.

## `CreateFilesForPRSet.R`

**Category:** PRS / PRSet prep

**Purpose**

Creates PLINK/PRSet input files for pathway-level PRS testing in sporadic ALS: keep file, phenotype file containing axis scores, and covariate file containing sex. Also deduplicates longitudinal samples to the earliest transcriptomic timepoint.

**Inputs**

- `sporadic_final_check` with `NYGC_CGND_ID` and clinical covariates.
- `ax_clean_df` with subtype assignments and axis scores.
- `vcf_sample_ids.txt`.

**Outputs**

- `sporadic_keep.txt`.
- `prset_phenotype.txt`.
- `prset_covariates.txt`.
- A no-FID phenotype/covariate variant written later in the script.

**Upstream dependencies**

- Must run after `sporadic_subtyping.R` has generated `ax_clean_df` and `sporadic_final_check`.
- Requires external genotype/VCF sample IDs.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/prs/CreateFilesForPRSet.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/prs/01_create_files_for_prset.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/prs/01_create_files_for_prset --objects outputs/objects/axis_scoring_objects.rds,outputs/objects/sporadic_subtyping_objects.rds
```

**Refactor notes**

- Convert `vcf_sample_ids.txt` path to a config entry.
- Write a manifest describing exact sample counts and excluded samples.
- Keep PRSet inputs under `inputs_for_external_tools/prset/` or `results/prs/prset_inputs/`.

## `PRSetTesting.R`

**Category:** PRS / PRSet results

**Purpose**

Tests whether pathway-specific PRS values predict matching transcriptomic axis scores. Reads PRSet `.best` files, merges PRS with sporadic axis/subtype metadata, fits all 3 x 3 PRS-to-axis regressions, and visualizes specificity.

**Inputs**

- `prset_axis_metabolic.best`.
- `prset_axis_cell_cycle.best`.
- `prset_axis_inflammation.best`.
- `vcf_sample_ids.txt`.
- Workspace objects `sporadic_final_check` and `ax_clean_df`.

**Outputs**

- `prset_specificity_matrix.pdf` and `.png`.
- `prs_cell_cycle_vs_transcriptomic_cell_cycle.pdf`.
- `prset_specificity_matrix_slides.png` and `.pdf`.
- `prset_diagonal_vs_offdiagonal_barplot.png`.
- `prset_pairwise_r2_barplot.png`.
- `prset_diagonal_regression_panel.png` and `.pdf`.
- `prset_all_9_regressions.png`.
- `prset_summary_table.csv`.
- Recommended new output: `analysis_df.rds` for sensitivity scripts.

**Upstream dependencies**

- Run after external PRSet/PRSice finishes using files from `CreateFilesForPRSet.R`.
- Requires sporadic subtype/axis metadata from `sporadic_subtyping.R`.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/prs/PRSetTesting.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/prs/02_prset_testing.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/prs/02_prset_testing --objects outputs/objects/axis_scoring_objects.rds,outputs/objects/sporadic_subtyping_objects.rds
```

**Refactor notes**

- Add covariate-adjusted models if required by analysis plan.
- Save `r2_mat`, `p_mat`, `b_mat`, and `analysis_df` as RDS/CSV.
- Separate publication figures from exploratory slide figures.

## `ciseQTLTesting.R`

**Category:** eQTL / CCT2

**Purpose**

Tests whether common variants in the CCT2 cis-window predict CCT2 expression and checks whether top cis-eQTL variants overlap ALS GWAS association signals.

**Inputs**

- `cct2_eqtl/cct2_cis_dosage.raw` from PLINK `--recode A`.
- `expr_clean_mat` containing row `CCT2`.
- `sporadic_final_check` containing `NYGC_CGND_ID`, `Participant_ID`, `SEX`, and `SV1`–`SV6`.
- `van_rheenen2021_als.tsv` GWAS summary statistics.

**Outputs**

- `cct2_cis_eqtl_results.csv`.
- `cct2_eqtl_gwas_lookup.csv`.
- `cct2_cis_eqtl_manhattan.pdf`.

**Upstream dependencies**

- Run after `AxisScoring.R`/`sporadic_subtyping.R` produces expression and covariate objects.
- Requires PLINK dosage generation outside R.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/eqtl/ciseQTLTesting.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/eqtl/01_cct2_cis_eqtl.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/eqtl/01_cct2_cis_eqtl --objects outputs/objects/axis_scoring_objects.rds,outputs/objects/sporadic_subtyping_objects.rds
```

**Refactor notes**

- Use explicit file inputs for expression and covariates instead of workspace objects.
- Make cis-window, MAF filter, missingness filter, and SV count configurable.
- Write a model specification block to the output directory for reproducibility.

## `DEGBurdenTesting.R`

**Category:** Genomics / DEG-informed burden

**Purpose**

Tests whether rare damaging variant burden in subtype-defining DEG sets is enriched in the matching subtype and decomposes LoF-to-NfL associations by axis gene-set burden.

**Inputs**

- DESeq2 result objects for subtype contrasts: `res_1v2`, `res_1v3`, `res_2v3`.
- Gene burden objects from `genomic_analysis.R`: `gene_burden_combined` and/or `gene_burden_filt`.
- Subtype/clinical metadata including `rare_burden_sub`, subtype labels, NfL, age, and sex.
- Annotation mapping through `org.Hs.eg.db`.

**Outputs**

- `targeted_deg_burden.pdf`.
- `axis_burden_nfl_forest.pdf`.
- `axis_burden_nfl_scatter.pdf`.
- `axis_decomposed_lof_nfl.csv`.

**Upstream dependencies**

- Requires `sporadic_subtyping.R` subtype DE objects and `genomic_analysis.R` burden matrices.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/genomics/DEGBurdenTesting.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/genomics/04_deg_burden_testing.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/genomics/04_deg_burden_testing --objects outputs/objects/axis_scoring_objects.rds,outputs/objects/sporadic_subtyping_objects.rds
```

**Refactor notes**

- Save subtype DESeq2 results from `sporadic_subtyping.R` as RDS/CSV and load them explicitly.
- Unify gene ID conversion before burden filtering.
- Return all model results as tidy tables, not only plots.

## `tric_burden.R`

**Category:** Genomics / targeted TRiC burden

**Purpose**

Performs targeted rare damaging/LoF variant burden tests for CCT2 and the TRiC chaperonin complex, including subtype association, S-ALS-3 membership models, and CCT2-expression association.

**Inputs**

- Exon annotation files `ALS.chrN.anno_report_exon.txt` in `EXON_DIR`.
- Workspace objects `use_cols`, `col_to_id`, `rare_burden_sub`, and `cct2_df`.

**Outputs**

- Folder `tric_burden_results/`: `tric_burden_full.csv`, `tric_burden_clean.csv`, `tric_per_gene_fisher.csv`, `tric_models.rds`.

**Upstream dependencies**

- Requires `genomic_analysis.R` to define sample columns, ID crosswalk, and rare burden metadata.
- Requires CCT2 expression table from sensitivity/eQTL/expression workflow.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/genomics/tric_burden.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/genomics/05_tric_burden.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/genomics/05_tric_burden --objects outputs/objects/axis_scoring_objects.rds,outputs/objects/sporadic_subtyping_objects.rds
```

**Refactor notes**

- Turn TRiC gene list, MAF threshold, variant classes, and chromosome list into parameters.
- Load `use_cols`, `col_to_id`, and `rare_burden_sub` from RDS files.
- Write per-test tidy model summaries alongside the RDS.

## `TRICBurdenMYCBurden.R`

**Category:** Genomics / targeted MYC burden

**Purpose**

Tests whether rare damaging burden across HALLMARK_MYC_TARGETS_V1 could explain subtype-associated CCT2 expression variation as a trans-regulatory mechanism.

**Inputs**

- Exon annotation files for all autosomes and chrX.
- MYC target genes from `axis_sets_clean`, `axis_sets`, or `msigdbr`.
- Workspace objects `use_cols`, `col_to_id`, `rare_burden_sub`, and `cct2_df`.

**Outputs**

- Folder `myc_burden_results/`: `myc_burden_full.csv`, `myc_burden_clean.csv`, `myc_per_gene_dmg.csv`, optional `myc_per_gene_lof.csv`, `myc_models.rds`.

**Upstream dependencies**

- Requires `genomic_analysis.R` and CCT2 expression objects.
- Usually run after `tric_burden.R` as a trans-regulatory follow-up.

**How to run**

Legacy documented version:

```bash
Rscript scripts/documented/genomics/TRICBurdenMYCBurden.R
```

Refactored command-line version:

```bash
Rscript scripts/refactored/genomics/06_myc_burden.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/genomics/06_myc_burden --objects outputs/objects/axis_scoring_objects.rds,outputs/objects/sporadic_subtyping_objects.rds
```

**Refactor notes**

- Share burden-matrix helper functions with `tric_burden.R` instead of duplicating them.
- Save gene universe used for MYC testing to the output folder.
- Add a single `targeted_pathway_burden.R` template to support TRiC, MYC, and future pathways.
