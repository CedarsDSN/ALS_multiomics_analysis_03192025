# Publish/update the ALS repository on GitHub

This guide assumes the repository already exists and should preserve these top-level folders:

```text
analysis/
scripts/
large_data/
```

Use a private GitHub repository unless the PI/institution explicitly approves public release.

## 1. Clone the existing repository

```bash
git clone https://github.com/priyanka8590/YOUR_REPO_NAME.git
cd YOUR_REPO_NAME
```

## 2. Create a branch

```bash
git checkout -b keep-analysis-scripts-layout
```

## 3. Confirm `large_data/` is empty or placeholder-only

```bash
find large_data -maxdepth 3 -type f -print
```

The folder should not contain controlled-access AnswerALS files.

## 4. Copy the update package into the repository root

```bash
rsync -av --exclude='.git' /path/to/als_existing_repo_keep_folders_v6/ ./
```

## 5. Review before committing

```bash
git status
git diff --stat
```

Check that none of the following are staged:

```text
raw AnswerALS files
clinical metadata
NfL tables
RNA-seq count matrices
proteomics matrices
VCF/BCF/BAM/CRAM/FASTQ files
PLINK files
patient-level .rds/.RData objects
config/config.yml
```

## 6. Stage GitHub-safe files

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

## 7. Commit and push

```bash
git commit -m "Document ALS analyses while preserving analysis and scripts folders"
git push -u origin keep-analysis-scripts-layout
```

## 8. Open a pull request

Open a pull request on GitHub:

```text
keep-analysis-scripts-layout -> main
```

Suggested PR title:

```text
Document ALS multi-omics analyses and add script inventory
```

Suggested PR summary:

```text
This update preserves the existing analysis/, scripts/, and large_data/ repository folders. It adds batch-level analysis documentation under analysis/, an inventory of all scripts used under scripts/, documented and refactored script versions, run-order documentation, object handoff notes, and a data-safety policy. The large_data/ folder remains an empty placeholder and controlled-access AnswerALS data remain excluded by .gitignore.
```

## 9. Local execution after cloning

```bash
cp config/config.example.yml config/config.yml
```

Edit `config/config.yml` locally to point to the private ALS data directory. Do not commit it.

Example:

```bash
Rscript scripts/refactored/transcriptomics/03_axis_scoring.R \
  --config config/config.yml \
  --base-dir /path/to/private/ALS \
  --outdir analysis/batch_03_axis_scoring_sporadic_subtyping/results/axis_scoring
```
