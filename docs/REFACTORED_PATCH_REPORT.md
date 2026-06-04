# Refactored script patch report

The `scripts/refactored/` files are best-effort command-line wrappers around the uploaded scripts. They preserve the analysis body, add config/CLI/bootstrap support, redirect relative outputs to `--outdir`, and save selected handoff objects when present.

They are not scientifically rewritten. Model formulas, thresholds, contrasts, and exploratory analysis logic are retained from the originals.

## `run_RNAseq_DE.R` -> `scripts/refactored/transcriptomics/01_progression_rnaseq_de.R`

- Script ID: `01_progression_rnaseq_de`
- Output handoff attempted: `01_progression_rnaseq_de_objects.rds`
- Objects saved if present: dds, dds_no_race, dds_disp, vsd, vst_data, pca_data, pca_result, merged, top_genes.

## `C9SOD1Convergence.R` -> `scripts/refactored/transcriptomics/02_c9_sod1_convergence.R`

- Script ID: `02_c9_sod1_convergence`
- Output handoff attempted: `02_c9_sod1_convergence_objects.rds`
- Objects saved if present: counts_filtered, annot_als_filtered, svobj, clean_svs, dds, dds_filtered, res_C9_filtered, res_SOD1_filtered, res_Sporadic_filtered, res_simple_c9, res_simple_sod1, convergent_high_conf, sig_strict, sig_lenient, sig_details, rna_signature_summary, protein_signature_summary.

## `AxisScoring.R` -> `scripts/refactored/transcriptomics/03_axis_scoring.R`

- Script ID: `03_axis_scoring`
- Output handoff attempted: `03_axis_scoring_objects.rds`
- Objects saved if present: counts_filtered, annot_als_filtered, annot_sva_filtered, svobj, clean_svs, dds_filtered, res_C9_filtered, res_SOD1_filtered, res_Sporadic_filtered, expr, expr_symbol, expr_clean_mat, axis_sets, axis_sets_clean, axis_scores, df_adj, df_nfl_uniq, marker_data, cont_res, clinical_res, marker_res, top_hits, interaction_res.

## `sporadic_subtyping.R` -> `scripts/refactored/subtyping/01_sporadic_subtyping.R`

- Script ID: `01_sporadic_subtyping`
- Output handoff attempted: `01_sporadic_subtyping_objects.rds`
- Objects saved if present: axis_scores, df_adj, ax_clean_df, ax_clean_dedup, sporadic_final_check, cluster_final, res_final_12, res_final_13, res_final_23, degs_12, degs_13, degs_23, all_ora, all_ora_clean, all_ora_final, ds_df.

## `SensitivityAnalysis.R` -> `scripts/refactored/validation/01_sensitivity_analysis.R`

- Script ID: `01_sensitivity_analysis`
- Output handoff attempted: `01_sensitivity_analysis_objects.rds`
- Objects saved if present: meta, sensitivity_summary, quality_results, analysis_df, cct2_df, prs_scores, cct2_meta, cct2_models.

## `LassoClassifier.R` -> `scripts/refactored/validation/02_lasso_classifier.R`

- Script ID: `02_lasso_classifier`
- Output handoff attempted: `02_lasso_classifier_objects.rds`
- Objects saved if present: results, stability_df, panel_composition, gene_to_axis, cct2_df.

## `IPAPrepare.R` -> `scripts/refactored/validation/03_ipa_prepare.R`

- Script ID: `03_ipa_prepare`
- Output handoff attempted: `03_ipa_prepare_objects.rds`
- Objects saved if present: dds, res, res_ipa, cct2_dist, subtype_frame.

## `CompareProteomicNormalization.R` -> `scripts/refactored/proteomics/01_compare_normalization.R`

- Script ID: `01_compare_normalization`
- Output handoff attempted: `01_compare_normalization_objects.rds`
- Objects saved if present: raw_data, datExpr_log, datExpr_quantile, datExpr_zscore, datExpr_median_norm.

## `proteomicsDE.R` -> `scripts/refactored/proteomics/02_proteomics_de.R`

- Script ID: `02_proteomics_de`
- Output handoff attempted: `02_proteomics_de_objects.rds`
- Objects saved if present: prot_matrix, meta_prot2, score_df, res_c9, res_sod1, res_spor, fg_all, top20, pcs, fit, design.

## `CompareBatchCorrectionProteomics.R` -> `scripts/refactored/proteomics/03_compare_batch_correction.R`

- Script ID: `03_compare_batch_correction`
- Output handoff attempted: `03_compare_batch_correction_objects.rds`
- Objects saved if present: pre_mat, post_mat, map_pre, map_post, summary_df, corr_mat, camera_df, top_cam.

## `sensitivity_analysis_prebatch.py` -> `scripts/refactored/proteomics/04_sensitivity_analysis_prebatch.py`

- Script ID: `04_sensitivity_analysis_prebatch`
- Output handoff attempted: `04_sensitivity_analysis_prebatch_objects.rds`
- Objects saved if present: None configured.

## `MAF_allele_frequency.R` -> `scripts/refactored/genomics/01_maf_allele_frequency.R`

- Script ID: `01_maf_allele_frequency`
- Output handoff attempted: `01_maf_allele_frequency_objects.rds`
- Objects saved if present: af_df, maf_summary, tier_summary.

## `genomic_analysis.R` -> `scripts/refactored/genomics/02_genomic_analysis.R`

- Script ID: `02_genomic_analysis`
- Output handoff attempted: `02_genomic_analysis_objects.rds`
- Objects saved if present: burden_all, burden_long, rare_burden_sub, gene_burden_combined, gene_burden_filt, use_cols, col_to_id.

## `genomic_stratification.R` -> `scripts/refactored/genomics/03_genomic_stratification.R`

- Script ID: `03_genomic_stratification`
- Output handoff attempted: `03_genomic_stratification_objects.rds`
- Objects saved if present: gene_tier_df, patient_gene_matrix_sporadic, genes_with_variants_final, mapped_df, interactions_genes, coexpr_edges_df, validated_edges, als_mn_network_final, patient_gene_matrix_network, edge_conservation_matrix, snf_clusters.

## `CreateFilesForPRSet.R` -> `scripts/refactored/prs/01_create_files_for_prset.R`

- Script ID: `01_create_files_for_prset`
- Output handoff attempted: `01_create_files_for_prset_objects.rds`
- Objects saved if present: keep_df, pheno_out, cov_out, pheno_out_nofid, cov_out_nofid, ax_clean_df, sporadic_final_check.

## `PRSetTesting.R` -> `scripts/refactored/prs/02_prset_testing.R`

- Script ID: `02_prset_testing`
- Output handoff attempted: `02_prset_testing_objects.rds`
- Objects saved if present: analysis_df, results_df, summary_table, prs_met, prs_cc, prs_inf.

## `ciseQTLTesting.R` -> `scripts/refactored/eqtl/01_cct2_cis_eqtl.R`

- Script ID: `01_cct2_cis_eqtl`
- Output handoff attempted: `01_cct2_cis_eqtl_objects.rds`
- Objects saved if present: dosage, cct2_df, crosswalk_v2, geno_pheno, eqtl_results, top_hit, gwas_hits.

## `DEGBurdenTesting.R` -> `scripts/refactored/genomics/04_deg_burden_testing.R`

- Script ID: `04_deg_burden_testing`
- Output handoff attempted: `04_deg_burden_testing_objects.rds`
- Objects saved if present: targeted_results, axis_results, p_forest, p_scatter.

## `tric_burden.R` -> `scripts/refactored/genomics/05_tric_burden.R`

- Script ID: `05_tric_burden`
- Output handoff attempted: `05_tric_burden_objects.rds`
- Objects saved if present: tric_df, tric_clean, per_gene, t4_unadj, t4_adj, t4_lof, t5, cct2_model.

## `TRICBurdenMYCBurden.R` -> `scripts/refactored/genomics/06_myc_burden.R`

- Script ID: `06_myc_burden`
- Output handoff attempted: `06_myc_burden_objects.rds`
- Objects saved if present: myc_df, myc_clean, per_gene_dmg, per_gene_lof, t4a, t4b, t4c, t5, cct2_model.

## Additional output patch

- Relative `outdir <- "..."` and `out_dir <- "..."` assignments in refactored R scripts are redirected under `ALSPROJ_OUTDIR`.
- `R/project_io.R` now warns instead of creating a missing private `base_dir`.
