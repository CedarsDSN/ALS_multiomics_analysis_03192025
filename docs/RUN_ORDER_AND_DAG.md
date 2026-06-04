# Workflow DAG and run order

```mermaid
flowchart TD
  raw_tx[Raw RNA-seq counts + metadata] --> rnaseq[run_RNAseq_DE]
  raw_tx --> fam[C9SOD1Convergence]
  raw_tx --> axis[AxisScoring]
  axis --> subtype[sporadic_subtyping]
  subtype --> sensitivity[SensitivityAnalysis]
  subtype --> lasso[LassoClassifier]
  subtype --> ipa[IPAPrepare]
  raw_prot[Proteomics matrices] --> norm[CompareProteomicNormalization]
  raw_prot --> protde[proteomicsDE]
  raw_prot --> batch[CompareBatchCorrectionProteomics]
  raw_prot --> pybatch[sensitivity_analysis_prebatch]
  raw_wgs[ANNOVAR + VCF files] --> maf[MAF_allele_frequency]
  raw_wgs --> genomics[genomic_analysis]
  genomics --> strat[genomic_stratification]
  subtype --> prprep[CreateFilesForPRSet]
  prprep --> external[External PLINK/PRSice PRSet]
  external --> prstest[PRSetTesting]
  axis --> eqtl[ciseQTLTesting]
  genomics --> deg[DEGBurdenTesting]
  genomics --> tric[tric_burden]
  tric --> myc[TRICBurdenMYCBurden]
```

## Recommended execution order

1. `scripts/refactored/transcriptomics/01_progression_rnaseq_de.R`
2. `scripts/refactored/transcriptomics/02_c9_sod1_convergence.R`
3. `scripts/refactored/transcriptomics/03_axis_scoring.R`
4. `scripts/refactored/subtyping/01_sporadic_subtyping.R`
5. `scripts/refactored/validation/03_ipa_prepare.R`
6. `scripts/refactored/validation/02_lasso_classifier.R`
7. `scripts/refactored/validation/01_sensitivity_analysis.R`
8. Proteomics scripts: normalization, differential abundance, batch correction, Python pre/post sensitivity.
9. Genomics scripts: MAF, burden, network stratification, DEG burden, TRiC, MYC.
10. PRS/eQTL scripts after external PLINK/PRSice and dosage generation.

Use `--objects` for downstream scripts that previously required objects to exist in the R session.
