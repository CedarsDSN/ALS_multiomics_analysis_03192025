# Batch 06: Genomics, rare variant burden, PRS, and CCT2 cis-eQTL

## Purpose

This batch documents the genomics side of the ALS project. It includes minor allele frequency summaries, rare damaging/LoF burden analyses, patient-gene matrices, PPI/coexpression network stratification, PRSet/PRS testing, CCT2 cis-eQTL testing, DEG-set burden testing, TRiC burden, and MYC-target burden follow-up.

## Scripts used

```text
scripts/documented/genomics/MAF_allele_frequency.R
scripts/refactored/genomics/01_maf_allele_frequency.R
scripts/documented/genomics/genomic_analysis.R
scripts/refactored/genomics/02_genomic_analysis.R
scripts/documented/genomics/genomic_stratification.R
scripts/refactored/genomics/03_genomic_stratification.R
scripts/documented/prs/CreateFilesForPRSet.R
scripts/refactored/prs/01_create_files_for_prset.R
scripts/documented/prs/PRSetTesting.R
scripts/refactored/prs/02_prset_testing.R
scripts/documented/eqtl/ciseQTLTesting.R
scripts/refactored/eqtl/01_cct2_cis_eqtl.R
scripts/documented/genomics/DEGBurdenTesting.R
scripts/refactored/genomics/04_deg_burden_testing.R
scripts/documented/genomics/tric_burden.R
scripts/refactored/genomics/05_tric_burden.R
scripts/documented/genomics/TRICBurdenMYCBurden.R
scripts/refactored/genomics/06_myc_burden.R
```

## Scientific questions

```text
1. What is the allele-frequency distribution of variants used in the genomic analyses?
2. Is rare damaging or LoF burden associated with ALS subtype or NfL?
3. Do patient-level rare variant patterns cluster into genomic strata?
4. Are axis-specific PRS scores associated with matching transcriptomic axes?
5. Are CCT2 expression differences supported by local cis-eQTL evidence?
6. Is rare burden enriched in subtype-defining DEGs or specific pathways such as TRiC/MYC targets?
```

## Inputs

Local-only inputs:

```text
GATK allele-frequency output
exon annotation files
rare variant annotation tables
patient/sample ID mapping files
sporadic subtype assignments
axis scores
NfL and clinical metadata
PLINK genotype files
PRS/PRSet score files
CCT2 expression values
dosage files for cis-eQTL testing
ALS GWAS summary statistics for overlap checks
```

Do not commit genotype files, dosage files, PLINK files, VCF/BCF files, or patient-level burden matrices.

## Outputs

MAF outputs:

```text
site-frequency spectrum plots
MAF distribution summaries
```

Burden/genomic analysis outputs:

```text
damaging/LoF burden tables
burden-by-subtype tests
NfL burden association summaries
```

Genomic stratification outputs:

```text
rare-variant patient-gene matrices
PPI/coexpression network summaries
edge perturbation/conservation metrics
SNF clustering outputs
```

PRS/PRSet outputs:

```text
PLINK/PRSet keep files
phenotype files
covariate files
axis-specific PRS association tables
3x3 PRS specificity matrices
PRS plots
```

cis-eQTL outputs:

```text
CCT2 cis-eQTL association table
variant dosage/expression plots
ALS GWAS overlap summary
```

DEG/TRiC/MYC burden outputs:

```text
subtype-DEG burden enrichment tables
axis burden decomposition for NfL
TRiC/CCT2 burden tables
HALLMARK_MYC_TARGETS_V1 burden results
```

## Upstream dependencies

This batch depends heavily on transcriptomic and subtype outputs from Batch 03.

Primary upstream objects/files:

```text
axis_scores
ax_clean_dedup
sporadic_subtype_assignments.csv
expr_clean_mat
subtype DEG result tables
clinical/NfL metadata
```

## How to run

MAF summary:

```bash
Rscript scripts/refactored/genomics/01_maf_allele_frequency.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_06_genomics_prs_eqtl/results/maf
```

Rare burden analysis:

```bash
Rscript scripts/refactored/genomics/02_genomic_analysis.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_06_genomics_prs_eqtl/results/genomic_analysis \
  --objects outputs/objects/03_axis_scoring_objects.rds,outputs/objects/01_sporadic_subtyping_objects.rds
```

Genomic stratification:

```bash
Rscript scripts/refactored/genomics/03_genomic_stratification.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_06_genomics_prs_eqtl/results/genomic_stratification
```

Create PRSet files:

```bash
Rscript scripts/refactored/prs/01_create_files_for_prset.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_06_genomics_prs_eqtl/results/prset_files \
  --objects outputs/objects/03_axis_scoring_objects.rds
```

Test PRSet/PRS associations:

```bash
Rscript scripts/refactored/prs/02_prset_testing.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_06_genomics_prs_eqtl/results/prset_testing \
  --objects outputs/objects/03_axis_scoring_objects.rds
```

CCT2 cis-eQTL:

```bash
Rscript scripts/refactored/eqtl/01_cct2_cis_eqtl.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_06_genomics_prs_eqtl/results/cct2_cis_eqtl \
  --objects outputs/objects/03_axis_scoring_objects.rds,outputs/objects/01_sporadic_subtyping_objects.rds
```

DEG burden:

```bash
Rscript scripts/refactored/genomics/04_deg_burden_testing.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_06_genomics_prs_eqtl/results/deg_burden \
  --objects outputs/objects/03_axis_scoring_objects.rds,outputs/objects/01_sporadic_subtyping_objects.rds
```

TRiC burden:

```bash
Rscript scripts/refactored/genomics/05_tric_burden.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_06_genomics_prs_eqtl/results/tric_burden \
  --objects outputs/objects/03_axis_scoring_objects.rds,outputs/objects/01_sporadic_subtyping_objects.rds
```

MYC-target burden:

```bash
Rscript scripts/refactored/genomics/06_myc_burden.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_06_genomics_prs_eqtl/results/myc_burden \
  --objects outputs/objects/03_axis_scoring_objects.rds,outputs/objects/01_sporadic_subtyping_objects.rds
```

## Run record

```text
Date run:
Machine/HPC node:
Git branch:
Git commit:
R version:
External tools and versions:
Run completed: yes/no
```

## Quality-control notes

Record:

```text
Variant filtering criteria:
MAF thresholds:
Annotation source/version:
Damaging/LoF definitions:
Number of patients with genomic data:
Number of patients overlapping subtype data:
Genotype/sample ID matching issues:
PRS model/covariates:
eQTL covariates and SVs:
```

## Main findings

Add summarized genomics/PRS/eQTL findings after local execution.

## Limitations

Potential limitations:

```text
small sample sizes after omics overlap
rare variant burden depends on annotation/filtering choices
population structure and ancestry adjustment
genotype-expression sample matching
eQTL power limitations
multiple testing across genes/pathways/subtypes
```
