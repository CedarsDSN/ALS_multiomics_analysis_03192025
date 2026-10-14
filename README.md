# Molecular subtyping of sporadic ALS

Analysis code for a study defining three transcriptomic subtypes of
sporadic ALS (S-ALS-1, S-ALS-2, S-ALS-3) from AnswerALS iPSC-derived
motor neuron multi-omics. The subtypes come from metabolic, cell-cycle
and inflammatory axes scored by ssGSEA on Hallmark gene sets and
partitioned with ConsensusClusterPlus. CCT2 is the single-gene
classifier for S-ALS-3. Supporting layers cover genotype differential
expression anchored on C9orf72 and SOD1, DIA proteomics, germline
burden and cis-eQTL, clinical endpoints, pathway-database
independence, IPA mechanism, and Target ALS post-mortem spatial
validation.

Cohort after removing the 98 replicated reference-line samples
(`CTRL-NEUEU392AE8`): 24,253 genes by 644 samples, 185 Control / 46
C9orf72 / 13 SOD1 / 400 Sporadic, giving 391 sporadic patients at
patient level and subtypes of 85 / 204 / 102.

## Requirements

R 4.3 or later with DESeq2, sva, limma, GSVA, ConsensusClusterPlus,
fgsea, msigdbr, glmnet, survival, org.Hs.eg.db, AnnotationDbi, Seurat,
data.table, dplyr, ggplot2, patchwork.

External tools for the germline arm: ANNOVAR, bcftools, PLINK.
Mechanism figures need IPA exports.

## How the scripts are built

Three rules hold throughout.

1. Inputs are pinned. Each script starts with `readRDS()` of a named
   file. No script depends on an object that happens to be in the
   session.
2. Outputs are saved. Each script ends with `saveRDS()` to a named
   file.
3. Reruns are gated. Before overwriting a canonical `_FINAL` file, a
   script recomputes its results and compares them to what is on disk.
   On a mismatch it writes `*_AUDITCHECK.rds` and leaves the canonical
   file untouched.

All `_FINAL` files are treated as immutable.

## Pinned inputs

Not in this repository. Expected alongside the project data directory.

| File | Contents |
|---|---|
| `cohort_inputs_FINAL.rds` | corrected counts, annotation with SV1 to SV6, 742-sample SVA annotation |
| `backbone_FINAL.rds` | axis scores, `axis_sets`, round-1 consensus, final labels |
| `patient_level_FINAL.rds` | 391 patients, earliest timepoint, clinical and NfL and CCT2 |
| `genotype_de_corrected_FINAL.rds` | genotype DE on clean SVs 1, 2, 6, 7, 8, 9 |
| `proteomics_inputs_FINAL.rds` | `res_c9`, `res_sod1`, `res_spor`, `fg_c9`, `fg_sod1`, `fg_spor` |
| `tals_de_results.rds` | `res_C9_tals`, `res_SOD1_tals`, `res_Spor_tals` and metadata |
| `genomic_burden_FINAL.rds` | burden matrices, blind clustering, `landmarks` |

## Run order

### Core pipeline

Run 01 and 02 first. Scripts 03 through 08 are independent of each
other. Script 16 needs 18, 14 and 16a.

| # | Script | Reads | Writes |
|---|---|---|---|
| 01 | `01_rebuild_backbone_AUDITED.R` | `cohort_inputs_FINAL.rds` | `backbone_FINAL.rds` |
| 02 | `02_patient_level_AUDITED.R` | backbone, clinical, NfL | `patient_level_FINAL.rds` |
| 03 | `03_differentiation_sensitivity.R` | backbone, patient level | `diff_quality_sensitivity_FINAL.rds` |
| 03b | `03b_marker_rescue_rerun.R` | as 03, plus pinned counts | marker-rescue adjudication |
| 04 | `04_lasso_AUDITED.R` | backbone, patient level | `lasso_FINAL.rds` |
| 05 | `05_downstream_ladders_AUDITED.R` | patient level, clinical, PRS | `downstream_analyses_FINAL.rds` |
| 06 | `06_clinical_screen_AUDITED.R` | patient level, portal CSVs, wave monolith | `clinical_screen_combined_FINAL.rds` |
| 07 | `07_ipa_prep_AUDITED.R` | backbone, pinned counts | `ipa_inputs_FINAL/`, `*_noXY.tsv` |
| 08 | `08_cct2_protein_test.R` | proteomics matrix, patient level | `cct2_protein_by_subtype.csv` |
| 14 | `14_tals_rebuild.R` | Target ALS bulk counts, metadata | `tals_de_results.rds` |
| 16a | `16a_pin_proteomics_inputs.R` | proteomics TSVs or session | `proteomics_inputs_FINAL.rds` |
| 16 | `16_convergence_corrected.R` | 18, 16a, 14 | `convergence_results_corrected/`, `proteomics_cct2_FINAL.rds` |
| 18 | `18_full_recovery_and_de_AUDITED.R` | `cohort_inputs_FINAL.rds` | `genotype_de_corrected_FINAL.rds` |

Script 06 sources `AllClinicalVariableTesting_newlabels.R`, which holds
the five clinical screen waves.

### Germline arm

Run from the directory holding the ANNOVAR
`ALS.chrN.anno_report_exon.txt` files, because `genomic_analysis.R`
reads them with `list.files(".")`.

```r
source("genomic_preamble.R")     # supplies sporadic_meta, ax_clean_df,
                                 # df_nfl_uniq, cct2_df
source("genomic_analysis.R")     # genome-wide burden
source("tric_burden.R")          # CCT2 and TRiC rare damaging burden
source("TRICBurdenMYCBurden.R")  # MYC burden
source("ciseQTLTesting.R")       # CCT2 cis-eQTL and GWAS lookup
```

| Script | Writes |
|---|---|
| `genomic_analysis.R` | `genomic_burden_FINAL.rds` |
| `tric_burden.R` | `tric_burden_results/tric_models.rds` |
| `TRICBurdenMYCBurden.R` | `myc_burden_results/myc_models.rds` |
| `ciseQTLTesting.R` | `cct2_cis_eqtl_results.csv`, `cct2_eqtl_gwas_lookup.csv` |

Rare variant definition: MAF below 0.01 on the ANNOVAR "Max Freq of 3"
column. Damaging classes are `nonsynonymous_SNV`, `stopgain`,
`stoploss`, `frameshift_deletion`, `frameshift_insertion`.

The GWAS lookup needs the harmonised GRCh38 release of the van Rheenen
summary statistics (GCST90027164, the `.h.tsv.gz` file) and its `hm_*`
columns. The GRCh37 release will appear to run and match almost
nothing.

### Pathway robustness

```r
source("30_pathway_preamble.R")
source("31_pathway_endpoint_screen.R")
source("32_database_independence.R")
```

30 supplies every object 31 and 32 need, from pinned files only. 31
writes the `phase2_*`, `phase6_*` and `phase7_*` CSVs. 32 writes
`database_independence_FINAL.rds` and `step3_ARI_matrix.csv`.

`RUN_CLEANED` in 32 stays `FALSE`. It guards an abandoned
residualisation that produces a degenerate solution.

### Spatial

```r
source("spatialdata.R")   # writes spatial_results_FINAL.rds
```

Frontal cortex only. Motor cortex contributes four samples, too few
for a group comparison.

### Figures

Source `00_figure_common.R` first, then any figure script. It provides
`p_()`, `load_pin()`, `get_backbone()`, `get_patients()`,
`get_group()`, `fmt_p()`, `save_fig()` and `figure_readiness()`, and
sets the subtype palette.

| Script | Figure |
|---|---|
| `20_figure1_cohort.R` | design, cohort sizes, sex, age at onset |
| `21_figure2_genotype_de.R` | familial to sporadic convergence |
| `22_figure3_axes.R` | axis activity, correlation, dimensionality |
| `23_figure4_consensus.R` | consensus matrix, PAC, axis heatmap |
| `24_figure5_cct2.R` | CCT2 classifier, ROC, bootstrap AUC, ECDF |
| `25_figure6_mechanism.R` | IPA mechanism, mTORC1 to LARP1, MYC |
| `26_figure7_genomic.R` | germline burden, cis-eQTL QQ, genomic clustering |
| `27_figure8_clinical.R` | clinical screen and CCT2 effects |
| `28_figure9_spatial.R` | Target ALS frontal cortex |

Figure 6 needs the IPA export CSVs from the noSexChr analysis. Figure 7
needs the germline arm to have run.

### Tables

```r
source("make_tables.R")   # Table 1, Table 2, Supplementary S1 to S4
```

Where a value has no saved object behind it, the console prints
`TYPED:` so it can be checked by hand.

## Verification

1. Fresh R session, no workspace restore.
2. `source()` the script top to bottom. No manual steps and no console
   detours between the seed and any stochastic step.
3. Compare the printed landmarks to the values in the script header.
