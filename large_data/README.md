# large_data placeholder

This folder is intentionally kept as an empty placeholder because the original GitHub repository included `large_data/`.

Do not commit raw or controlled-access AnswerALS data here.

Examples of files that must stay local/private:

```text
raw RNA-seq count matrices
proteomics matrices
clinical metadata
NfL tables
mutation subtype tables if restricted
VCF/BCF files
BAM/CRAM/FASTQ files
PLINK genotype files
patient-level .rds/.RData objects
large generated result files
```

Use `config/config.yml` locally to point scripts to the private data location. The public repository should only contain `large_data/README.md` and `large_data/.gitkeep`.
