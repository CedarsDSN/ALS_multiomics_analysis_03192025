# Migration guide: keep existing `analysis/`, `scripts/`, and `large_data/` folders

This guide updates the existing GitHub repository while preserving the folder structure requested for the project handoff.

## Target repository structure

```text
analysis/      Keep all analysis details here, with subfolders for batches.
scripts/       Include all scripts used.
large_data/    Keep as an empty placeholder only.
```

Support folders are also added:

```text
config/        Public example config only.
R/             Shared helper functions.
docs/          Script specs, run order, object contracts, and refactor notes.
outputs/       Ignored placeholder for generated object handoffs.
results/       Ignored placeholder for generated results.
.github/       Basic smoke-test workflow.
```

## Step 1: clone the existing repository

```bash
git clone https://github.com/priyanka8590/YOUR_REPO_NAME.git
cd YOUR_REPO_NAME
```

Replace `YOUR_REPO_NAME` with the actual repository name.

## Step 2: create an update branch

```bash
git checkout -b keep-analysis-scripts-layout
```

## Step 3: confirm `large_data/` is only a placeholder

```bash
find large_data -maxdepth 3 -type f -print
```

Expected output should only include placeholder files such as:

```text
large_data/README.md
large_data/.gitkeep
```

If other files appear, do not commit them until you verify they are safe.

## Step 4: copy this package into the repo

From inside the cloned repository, copy the contents of this package into the repo root:

```bash
rsync -av --exclude='.git' /path/to/als_existing_repo_keep_folders_v6/ ./
```

## Step 5: check what changed

```bash
git status
git diff --stat
```

Check that raw data are not staged:

```bash
git status --short
```

You should not see raw AnswerALS files, genotype files, count matrices, protein matrices, clinical metadata, `.rds`, `.RData`, or real `config/config.yml`.

## Step 6: stage the GitHub-safe files

```bash
git add README.md \
  MIGRATION_KEEP_EXISTING_FOLDERS.md \
  PUBLISH_TO_GITHUB.md \
  SECURITY_AND_DATA_POLICY.md \
  LICENSE_NOTE.md \
  .gitignore \
  requirements.txt \
  analysis \
  scripts \
  large_data \
  config \
  R \
  docs \
  outputs \
  results \
  .github
```

Review staged files:

```bash
git diff --cached --stat
git status
```

## Step 7: commit and push

```bash
git commit -m "Document ALS analyses while preserving analysis and scripts folders"
git push -u origin keep-analysis-scripts-layout
```

## Step 8: open a pull request

Open a pull request from:

```text
keep-analysis-scripts-layout -> main
```

Suggested PR title:

```text
Document ALS multi-omics analyses and add script inventory
```

Suggested PR summary:

```text
This update preserves the existing analysis/, scripts/, and large_data/ repository folders. It adds batch-level analysis documentation under analysis/, an inventory of all scripts used under scripts/, refactored command-line script versions, script specifications, run-order documentation, object handoff notes, and a data-safety policy. large_data/ remains an empty placeholder and controlled-access AnswerALS files remain excluded through .gitignore.
```
