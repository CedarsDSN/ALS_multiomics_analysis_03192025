# Scripts folder

This folder contains all scripts used in the ALS multi-omics project.

Two script forms are included:

```text
scripts/documented/    Analysis scripts with standardized documentation headers.
scripts/refactored/    Command-line versions using --config, --base-dir, --outdir, and optional --objects.
```

The documented scripts preserve the original analysis style and file names. The refactored scripts provide cleaner entrypoints for reproducible local or HPC execution.

Use `scripts/ALL_SCRIPTS_USED.csv` as the script inventory.

## Script domains

```text
transcriptomics/    RNA-seq DE, familial ALS convergence, and axis scoring.
subtyping/          Sporadic ALS subtype discovery.
validation/         Sensitivity analysis, LASSO classifier, and IPA preparation.
proteomics/         Normalization, proteomic DE, batch correction, and pre/post batch sensitivity.
genomics/           MAF, rare burden, genomic stratification, DEG/TRiC/MYC burden analyses.
prs/                PRSet file creation and PRS-axis testing.
eqtl/               CCT2 cis-eQTL testing.
```

## Rule for adding new scripts

When adding a new script:

1. Put the script under the appropriate domain in `scripts/documented/` or `scripts/refactored/`.
2. Add or update its entry in `scripts/ALL_SCRIPTS_USED.csv`.
3. Document the associated run in the matching `analysis/batch_*` folder.
4. Avoid hard-coded private paths when possible; use `config/config.yml` or command-line arguments.
5. Do not commit raw data or patient-level outputs.

## Minimal documented header

Every new script should start with:

```text
Purpose
Inputs
Outputs
Upstream dependencies
How to run
```

Full script-level details are also maintained in `scripts/README.md`.
