# Batch 02: Familial C9orf72/SOD1 convergence transcriptomic analysis

## Purpose

This batch documents the familial ALS transcriptomic workflow focused on C9orf72 and SOD1 mutation subtypes. The analysis runs RNA-seq QC, surrogate-variable analysis, DESeq2 contrasts, convergence/divergence summaries, pathway enrichment, ALS gene-set guardrail scoring, and cross-omics signature summaries.

## Scripts used

```text
scripts/documented/transcriptomics/C9SOD1Convergence.R
scripts/refactored/transcriptomics/02_c9_sod1_convergence.R
```

## Scientific questions

```text
1. What genes are differentially expressed in C9orf72 ALS vs controls?
2. What genes are differentially expressed in SOD1 ALS vs controls?
3. Which transcriptomic signals converge between C9orf72 and SOD1?
4. Which signals diverge between mutation groups?
5. What pathways are enriched in the shared or mutation-specific signatures?
6. Do familial ALS signatures connect to later sporadic subtype or cross-omics findings?
```

## Inputs

Controlled or local-only inputs:

```text
transcriptomics/4_matrix/AnswerALS-651-T-v1-release6_raw-counts.csv
transcriptomics/4_matrix/Sample_Mapping_Information/Sample_Mapping_File_Feb_2024.csv
metadata/clinical/aals_dataportal_datatable.csv
ALS_Mutation_Subtypes_with_ALS_Control.csv
sod1_sample_mutations.csv
optional external paper/supplement tables used for cross-reference sections
```

## Outputs

Expected table outputs include:

```text
sva_sensitivity_all_svs.csv
sva_filtered_summary.csv
model_comparison_simple_vs_filtered.csv
complete_three_model_comparison.csv
pathway_GO_BP_results.csv
pathway_KEGG_results.csv
pathway_Reactome_results.csv
pathway_shared_C9_SOD1.csv
familial_signature_broad.csv
familial_signature_strict.csv
familial_signature_lenient.csv
familial_signature_strict_sva_markers.csv
sporadic_guardrail_High_vs_Low_DE.csv
GSEA_GO_BP_sporadic_guardrail.csv
sporadic_samples_vcf.txt
sporadic_samples_present.txt
rna_signature_summary.csv
protein_signature_summary.csv
cross_omics_feature.csv
cross_omics_pathway.csv
integrated_signature_results_matrix.csv
```

Expected plot outputs include dispersion plots, PCA plots, GSEA plots, mutation/sex/age/progression distributions, and C9/SOD1 divergence summaries.

Recommended handoff objects:

```text
counts_filtered.rds
annot_als_filtered.rds
svobj.rds
clean_svs.rds
expr_vst.rds
familial_signatures.rds
```

## Upstream dependencies

None for the first RNA-seq sections. Later sections depend on objects generated earlier in the same script.

## How to run

Documented script:

```bash
Rscript scripts/documented/transcriptomics/C9SOD1Convergence.R
```

Refactored script:

```bash
Rscript scripts/refactored/transcriptomics/02_c9_sod1_convergence.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_02_familial_c9_sod1_convergence/results/c9_sod1_convergence
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
C9orf72 sample count:
SOD1 sample count:
Control sample count:
Genes retained:
Surrogate variables estimated:
Surrogate variables removed or retained:
Batch/covariate concerns:
PCA findings:
```

## Main findings

Add final summarized results here after local execution.

## Limitations

Potential limitations:

```text
small mutation-subtype sample sizes
family-specific and mutation-specific heterogeneity
SVA interpretation
cross-omics sample overlap
multiple exploratory sections in one script
```
