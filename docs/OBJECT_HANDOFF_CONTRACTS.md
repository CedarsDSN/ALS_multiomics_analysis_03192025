# Object handoff contracts

The refactored R scripts load object bundles with `--objects` and save selected objects into `<outdir>/<script_id>_objects.rds` if those objects exist after successful execution.

## `01_progression_rnaseq_de_objects.rds`

Produced by: `scripts/refactored/transcriptomics/01_progression_rnaseq_de.R`

Saved objects, if present:

- `dds`
- `dds_no_race`
- `dds_disp`
- `vsd`
- `vst_data`
- `pca_data`
- `pca_result`
- `merged`
- `top_genes`

## `02_c9_sod1_convergence_objects.rds`

Produced by: `scripts/refactored/transcriptomics/02_c9_sod1_convergence.R`

Saved objects, if present:

- `counts_filtered`
- `annot_als_filtered`
- `svobj`
- `clean_svs`
- `dds`
- `dds_filtered`
- `res_C9_filtered`
- `res_SOD1_filtered`
- `res_Sporadic_filtered`
- `res_simple_c9`
- `res_simple_sod1`
- `convergent_high_conf`
- `sig_strict`
- `sig_lenient`
- `sig_details`
- `rna_signature_summary`
- `protein_signature_summary`

## `03_axis_scoring_objects.rds`

Produced by: `scripts/refactored/transcriptomics/03_axis_scoring.R`

Saved objects, if present:

- `counts_filtered`
- `annot_als_filtered`
- `annot_sva_filtered`
- `svobj`
- `clean_svs`
- `dds_filtered`
- `res_C9_filtered`
- `res_SOD1_filtered`
- `res_Sporadic_filtered`
- `expr`
- `expr_symbol`
- `expr_clean_mat`
- `axis_sets`
- `axis_sets_clean`
- `axis_scores`
- `df_adj`
- `df_nfl_uniq`
- `marker_data`
- `cont_res`
- `clinical_res`
- `marker_res`
- `top_hits`
- `interaction_res`

## `01_sporadic_subtyping_objects.rds`

Produced by: `scripts/refactored/subtyping/01_sporadic_subtyping.R`

Saved objects, if present:

- `axis_scores`
- `df_adj`
- `ax_clean_df`
- `ax_clean_dedup`
- `sporadic_final_check`
- `cluster_final`
- `res_final_12`
- `res_final_13`
- `res_final_23`
- `degs_12`
- `degs_13`
- `degs_23`
- `all_ora`
- `all_ora_clean`
- `all_ora_final`
- `ds_df`

## `01_sensitivity_analysis_objects.rds`

Produced by: `scripts/refactored/validation/01_sensitivity_analysis.R`

Saved objects, if present:

- `meta`
- `sensitivity_summary`
- `quality_results`
- `analysis_df`
- `cct2_df`
- `prs_scores`
- `cct2_meta`
- `cct2_models`

## `02_lasso_classifier_objects.rds`

Produced by: `scripts/refactored/validation/02_lasso_classifier.R`

Saved objects, if present:

- `results`
- `stability_df`
- `panel_composition`
- `gene_to_axis`
- `cct2_df`

## `03_ipa_prepare_objects.rds`

Produced by: `scripts/refactored/validation/03_ipa_prepare.R`

Saved objects, if present:

- `dds`
- `res`
- `res_ipa`
- `cct2_dist`
- `subtype_frame`

## `01_compare_normalization_objects.rds`

Produced by: `scripts/refactored/proteomics/01_compare_normalization.R`

Saved objects, if present:

- `raw_data`
- `datExpr_log`
- `datExpr_quantile`
- `datExpr_zscore`
- `datExpr_median_norm`

## `02_proteomics_de_objects.rds`

Produced by: `scripts/refactored/proteomics/02_proteomics_de.R`

Saved objects, if present:

- `prot_matrix`
- `meta_prot2`
- `score_df`
- `res_c9`
- `res_sod1`
- `res_spor`
- `fg_all`
- `top20`
- `pcs`
- `fit`
- `design`

## `03_compare_batch_correction_objects.rds`

Produced by: `scripts/refactored/proteomics/03_compare_batch_correction.R`

Saved objects, if present:

- `pre_mat`
- `post_mat`
- `map_pre`
- `map_post`
- `summary_df`
- `corr_mat`
- `camera_df`
- `top_cam`

## `01_maf_allele_frequency_objects.rds`

Produced by: `scripts/refactored/genomics/01_maf_allele_frequency.R`

Saved objects, if present:

- `af_df`
- `maf_summary`
- `tier_summary`

## `02_genomic_analysis_objects.rds`

Produced by: `scripts/refactored/genomics/02_genomic_analysis.R`

Saved objects, if present:

- `burden_all`
- `burden_long`
- `rare_burden_sub`
- `gene_burden_combined`
- `gene_burden_filt`
- `use_cols`
- `col_to_id`

## `03_genomic_stratification_objects.rds`

Produced by: `scripts/refactored/genomics/03_genomic_stratification.R`

Saved objects, if present:

- `gene_tier_df`
- `patient_gene_matrix_sporadic`
- `genes_with_variants_final`
- `mapped_df`
- `interactions_genes`
- `coexpr_edges_df`
- `validated_edges`
- `als_mn_network_final`
- `patient_gene_matrix_network`
- `edge_conservation_matrix`
- `snf_clusters`

## `01_create_files_for_prset_objects.rds`

Produced by: `scripts/refactored/prs/01_create_files_for_prset.R`

Saved objects, if present:

- `keep_df`
- `pheno_out`
- `cov_out`
- `pheno_out_nofid`
- `cov_out_nofid`
- `ax_clean_df`
- `sporadic_final_check`

## `02_prset_testing_objects.rds`

Produced by: `scripts/refactored/prs/02_prset_testing.R`

Saved objects, if present:

- `analysis_df`
- `results_df`
- `summary_table`
- `prs_met`
- `prs_cc`
- `prs_inf`

## `01_cct2_cis_eqtl_objects.rds`

Produced by: `scripts/refactored/eqtl/01_cct2_cis_eqtl.R`

Saved objects, if present:

- `dosage`
- `cct2_df`
- `crosswalk_v2`
- `geno_pheno`
- `eqtl_results`
- `top_hit`
- `gwas_hits`

## `04_deg_burden_testing_objects.rds`

Produced by: `scripts/refactored/genomics/04_deg_burden_testing.R`

Saved objects, if present:

- `targeted_results`
- `axis_results`
- `p_forest`
- `p_scatter`

## `05_tric_burden_objects.rds`

Produced by: `scripts/refactored/genomics/05_tric_burden.R`

Saved objects, if present:

- `tric_df`
- `tric_clean`
- `per_gene`
- `t4_unadj`
- `t4_adj`
- `t4_lof`
- `t5`
- `cct2_model`

## `06_myc_burden_objects.rds`

Produced by: `scripts/refactored/genomics/06_myc_burden.R`

Saved objects, if present:

- `myc_df`
- `myc_clean`
- `per_gene_dmg`
- `per_gene_lof`
- `t4a`
- `t4b`
- `t4c`
- `t5`
- `cct2_model`
