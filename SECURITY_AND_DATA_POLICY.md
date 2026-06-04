# Security and data policy

This repository should contain code and documentation only.

Do not commit:

- AnswerALS raw data or downloaded matrices.
- Patient-level clinical metadata.
- NfL tables, sample-level phenotype tables, or mutation subtype tables if they are restricted.
- VCF/BCF/BAM/CRAM/FASTQ files or PLINK genotype files.
- Generated `.rds`, `.RData`, `.h5`, `.parquet`, or large intermediate objects.
- Real `config/config.yml` files containing local private data paths.
- Local absolute paths, usernames, tokens, or credentials.

Use `config/config.example.yml` as the public template. Keep the real config file local and ignored by Git.
