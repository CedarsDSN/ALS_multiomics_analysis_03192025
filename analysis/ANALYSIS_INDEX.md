# Analysis index

This file maps the major ALS project analysis batches to scripts, inputs, outputs, and downstream dependencies.

| Batch | Folder | Purpose | Scripts used | Upstream dependencies | Primary downstream use |
|---|---|---|---|---|---|
| 01 | `batch_01_transcriptomics_progression/` | RNA-seq QC and ALS progression-class DE | `run_RNAseq_DE.R` | Raw counts, sample map, progression class table, clinical metadata | Baseline transcriptomics QC and progression comparisons |
| 02 | `batch_02_familial_c9_sod1_convergence/` | Familial ALS C9orf72/SOD1 DE, convergence, pathways, cross-omics summaries | `C9SOD1Convergence.R` | Raw counts, sample map, clinical metadata, mutation subtype table | Familial ALS signatures and SVA/model objects |
| 03 | `batch_03_axis_scoring_sporadic_subtyping/` | Hallmark axis scoring and sporadic ALS subtype discovery | `AxisScoring.R`, `sporadic_subtyping.R` | Transcriptomics objects, Hallmark gene sets, clinical metadata, NfL | Subtype assignments, subtype DE, axis scores |
| 04 | `batch_04_validation_classifier_ipa/` | Sensitivity analyses, LASSO subtype classifier, IPA-ready outputs | `SensitivityAnalysis.R`, `LassoClassifier.R`, `IPAPrepare.R` | Axis scoring and sporadic subtype objects | Validation evidence, sparse gene panels, IPA tables |
| 05 | `batch_05_proteomics/` | Proteomic normalization, differential abundance, batch correction, pre/post batch sensitivity | `CompareProteomicNormalization.R`, `proteomicsDE.R`, `CompareBatchCorrectionProteomics.R`, `sensitivity_analysis_prebatch.py` | Proteomics matrices, mutation subtype table, clinical/sample metadata | Proteomic DE and batch-robust protein results |
| 06 | `batch_06_genomics_prs_eqtl/` | MAF, rare variant burden, genomic stratification, PRS, cis-eQTL | `MAF_allele_frequency.R`, `genomic_analysis.R`, `genomic_stratification.R`, `CreateFilesForPRSet.R`, `PRSetTesting.R`, `ciseQTLTesting.R`, `DEGBurdenTesting.R`, `tric_burden.R`, `TRICBurdenMYCBurden.R` | Variant annotations, subtype assignments, axis scores, genotype/PRS files | Genetic burden, PRS specificity, CCT2 eQTL evidence |

## Recommended execution order

```text
1. batch_01_transcriptomics_progression
2. batch_02_familial_c9_sod1_convergence
3. batch_03_axis_scoring_sporadic_subtyping
4. batch_04_validation_classifier_ipa
5. batch_05_proteomics
6. batch_06_genomics_prs_eqtl
```

The strongest dependency chain is:

```text
AxisScoring.R -> sporadic_subtyping.R -> SensitivityAnalysis.R / LassoClassifier.R / IPAPrepare.R / PRSetTesting.R / ciseQTLTesting.R / DEGBurdenTesting.R
```

For full object-level dependencies, see `docs/OBJECT_HANDOFF_CONTRACTS.md`.
