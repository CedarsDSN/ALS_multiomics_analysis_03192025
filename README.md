# ALS Multiomics Project

## Project description

This repository documents the ALS multi-omics analysis project using the AnswerALS dataset. The project aims to identify and characterize molecular subtypes of ALS by integrating transcriptomics, proteomics, genomics, pathway enrichment, and clinical metadata. The analyses include differential expression, subtype discovery, subtype validation, proteomic normalization and batch sensitivity, rare variant burden analyses, PRS testing, and CCT2 cis-eQTL follow-up.

The repository is organized to preserve the original project folders requested for the GitHub handoff:

```text
analysis/      Detailed analysis records, organized by analysis batch.
scripts/       All scripts used in the project.
large_data/    Empty placeholder only; controlled-access data are not committed.
```

Additional support folders were added for reproducibility:

```text
config/        Example configuration file. The real config/config.yml stays local.
R/             Shared helper functions for refactored R scripts.
docs/          Cross-script documentation, run order, object handoff contracts, and script specs.
outputs/       Placeholder for generated object handoffs; ignored by Git.
results/       Placeholder for generated results; ignored by Git.
.github/       Basic repository smoke-test workflow.
```

## Data

AnswerALS contains clinical and omics data, including genomics, transcriptomics, epigenomics, proteomics, and clinical metadata. Access requires a research proposal and data user agreement.

Project data access timeline:

```text
Data requested from AnswerALS:       March 27, 2025
Data user agreement received:        April 9, 2025
Access granted:                      April 28, 2025
```

The Box folder used for data information was:

```text
Projects_Internal/ALS_multiomics_data/Data/
```

Within that folder, clinical data were under `metadata/` and transcriptomics data were under `transcriptomic/`.

The HPC project folder used for analysis was:

```text
/common/compbiomed-dsn/projects/ALS_multiomics_analysis_03192025/analysis
```

Raw AnswerALS data and patient-level files should not be committed to this repository. The `large_data/` folder is retained only as an empty placeholder.

## Repository layout

```text
ALS_multiomics_project/
├── README.md
├── analysis/
│   ├── README.md
│   ├── ANALYSIS_INDEX.md
│   ├── _template/
│   ├── batch_01_transcriptomics_progression/
│   ├── batch_02_familial_c9_sod1_convergence/
│   ├── batch_03_axis_scoring_sporadic_subtyping/
│   ├── batch_04_validation_classifier_ipa/
│   ├── batch_05_proteomics/
│   └── batch_06_genomics_prs_eqtl/
├── scripts/
│   ├── README.md
│   ├── ALL_SCRIPTS_USED.csv
│   ├── documented/
│   └── refactored/
├── large_data/
│   ├── README.md
│   └── .gitkeep
├── config/
│   └── config.example.yml
├── R/
│   └── project_io.R
├── docs/
├── outputs/
└── results/
```

## Analysis batches

The `analysis/` folder is the project notebook/logbook for the repository. Each analysis batch has its own subfolder with a README that records the goal, scripts used, inputs, outputs, run commands, status, and interpretation notes.

| Batch folder | Main purpose | Primary scripts |
|---|---|---|
| `analysis/batch_01_transcriptomics_progression/` | RNA-seq QC and progression-class DE | `run_RNAseq_DE.R` |
| `analysis/batch_02_familial_c9_sod1_convergence/` | C9orf72/SOD1 familial ALS convergence | `C9SOD1Convergence.R` |
| `analysis/batch_03_axis_scoring_sporadic_subtyping/` | Hallmark axis scoring and sporadic subtype discovery | `AxisScoring.R`, `sporadic_subtyping.R` |
| `analysis/batch_04_validation_classifier_ipa/` | Subtype sensitivity checks, LASSO classifier, IPA preparation | `SensitivityAnalysis.R`, `LassoClassifier.R`, `IPAPrepare.R` |
| `analysis/batch_05_proteomics/` | Proteomic normalization, DE, batch correction, pre/post batch sensitivity | `CompareProteomicNormalization.R`, `proteomicsDE.R`, `CompareBatchCorrectionProteomics.R`, `sensitivity_analysis_prebatch.py` |
| `analysis/batch_06_genomics_prs_eqtl/` | MAF, burden, genomic stratification, PRS, cis-eQTL | `MAF_allele_frequency.R`, `genomic_analysis.R`, `genomic_stratification.R`, `CreateFilesForPRSet.R`, `PRSetTesting.R`, `ciseQTLTesting.R`, `DEGBurdenTesting.R`, `tric_burden.R`, `TRICBurdenMYCBurden.R` |

## Scripts

All scripts used in the project are under `scripts/`.

```text
scripts/documented/    Script files with standardized headers explaining purpose, inputs, outputs, dependencies, and run commands.
scripts/refactored/    Best-effort command-line refactors using --config, --base-dir, --outdir, and optional --objects.
```

Use `scripts/ALL_SCRIPTS_USED.csv` for a compact index of the full script set.

The original scientific analysis logic was preserved in the documented/refactored versions. Refactoring focused on reproducibility structure: explicit paths, config-based execution, output directories, and object handoff files.

## How to run locally

Create a private local config file from the public template:

```bash
cp config/config.example.yml config/config.yml
```

Edit `config/config.yml` locally so paths point to the private ALS data location. Do not commit `config/config.yml`.

Example command:

```bash
Rscript scripts/refactored/transcriptomics/03_axis_scoring.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_03_axis_scoring_sporadic_subtyping/results/axis_scoring
```

Downstream scripts can receive saved object bundles through `--objects`:

```bash
Rscript scripts/refactored/subtyping/01_sporadic_subtyping.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_03_axis_scoring_sporadic_subtyping/results/sporadic_subtyping \
  --objects outputs/objects/03_axis_scoring_objects.rds
```

## Data safety

Do not commit:

```text
raw AnswerALS data
clinical metadata
NfL tables
RNA-seq count matrices
proteomics matrices
genotype files
VCF/BCF/BAM/CRAM/FASTQ files
PLINK files
patient-level .rds/.RData objects
real config/config.yml files
generated results containing patient-level data
```

The `.gitignore` is configured to keep `large_data/` as an empty placeholder and to block common raw data, genotype, intermediate-object, and generated-output file types.

## Main documentation files

```text
analysis/ANALYSIS_INDEX.md             Batch-level analysis index.
scripts/ALL_SCRIPTS_USED.csv           Script inventory.
docs/SCRIPT_SPECS_FULL.md              Purpose, inputs, outputs, dependencies, and run command for every script.
docs/RUN_ORDER_AND_DAG.md              Suggested execution order.
docs/OBJECT_HANDOFF_CONTRACTS.md       Shared intermediate object contracts.
docs/REFACTORED_PATCH_REPORT.md        Summary of refactoring changes.
SECURITY_AND_DATA_POLICY.md            Data safety policy.
PUBLISH_TO_GITHUB.md                   GitHub publishing instructions.
```
