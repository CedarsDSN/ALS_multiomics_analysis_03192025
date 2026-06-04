# Analysis folder

This folder is the analysis logbook for the ALS multi-omics project.

Per the repository convention, all analysis details should live here, not only in commit messages or scattered notes. Each subfolder represents a distinct batch of analysis. A batch may be a single analysis script, a related group of scripts, or a complete analysis theme such as proteomics or genomics.

## Folder convention

```text
analysis/
├── README.md
├── ANALYSIS_INDEX.md
├── _template/
│   └── analysis_record_template.md
├── batch_01_transcriptomics_progression/
├── batch_02_familial_c9_sod1_convergence/
├── batch_03_axis_scoring_sporadic_subtyping/
├── batch_04_validation_classifier_ipa/
├── batch_05_proteomics/
└── batch_06_genomics_prs_eqtl/
```

## What each batch README should contain

Each batch folder should document:

```text
Purpose
Question being asked
Scripts used
Inputs
Outputs
Upstream dependencies
How to run
Run status
Main findings
Quality-control notes
Known limitations
Next steps
```

Use `analysis/_template/analysis_record_template.md` when creating new analysis batches.

## What belongs in analysis/

Commit these when they do not expose controlled data:

```text
README files
run notes
analysis decisions
non-sensitive summary tables
workflow diagrams
small aggregate plots cleared for sharing
batch-level interpretation notes
```

Do not commit these unless explicitly cleared:

```text
raw AnswerALS data
patient-level metadata
sample-level clinical tables
count matrices
protein matrices
genotype files
large generated outputs
patient-level RDS objects
```

Generated results should usually be written to `analysis/<batch>/results/` or `results/`, but those paths are ignored by Git by default. If a small summary output is safe to include, move it into the batch README or a clearly named summary file after checking that it contains no protected patient-level information.

## Batch index

See `analysis/ANALYSIS_INDEX.md` for the current batch list and script mapping.
