# Analysis batch: <short descriptive name>

## Purpose

Describe the scientific question and why this analysis was run.

## Date / owner

```text
Date started:
Date completed:
Analyst:
Git branch:
Git commit:
```

## Scripts used

```text
scripts/documented/<domain>/<script_name>
scripts/refactored/<domain>/<script_name>
```

## Inputs

List input files, their private location, and whether they are controlled-access.

```text
Input file:
Source:
Description:
Controlled data: yes/no
```

## Outputs

List tables, figures, RDS objects, and logs.

```text
Output:
Location:
Description:
Committed to GitHub: yes/no
Reason if not committed:
```

## Upstream dependencies

List scripts or objects that must be created first.

## How to run

```bash
Rscript scripts/refactored/<domain>/<script>.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/<batch>/results/<analysis_name>
```

## Run notes

Record package issues, file path corrections, warnings, or deviations from the planned workflow.

## Quality control

Summarize QC plots, sample exclusions, outlier checks, batch effects, and sensitivity checks.

## Main findings

Summarize high-level findings without exposing patient-level data.

## Limitations

Document caveats, unresolved issues, and assumptions.

## Next steps

List follow-up analyses or validation tasks.
