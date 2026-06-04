# ============================================================
# Refactored entrypoint for original script: CreateFilesForPRSet.R
# Use --config, --base-dir, --outdir, and --objects. See scripts/README.md.
# ============================================================
# ---- ALS project bootstrap ----
.get_script_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- "--file="
  hit <- grep(file_arg, args, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub(file_arg, "", hit[[1]]), mustWork = FALSE)))
  getwd()
}
.find_project_root <- function(start) {
  path <- normalizePath(start, mustWork = FALSE)
  for (i in seq_len(8)) {
    if (file.exists(file.path(path, "R", "project_io.R"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }
  stop("Cannot find R/project_io.R. Run this script from the repository root.", call. = FALSE)
}
ALSPROJ_ROOT <- .find_project_root(.get_script_dir())
source(file.path(ALSPROJ_ROOT, "R", "project_io.R"))
setwd(ALSPROJ_ROOT)
ALS_ARGS <- als_parse_args()
if (isTRUE(ALS_ARGS$help)) { als_print_help("01_create_files_for_prset"); quit(save = "no", status = 0) }
CFG <- als_load_config(ALS_ARGS$config)
ALSPROJ_BASE_DIR <- als_choose(ALS_ARGS$base_dir, CFG, "project.base_dir", getwd())
ALSPROJ_OUTDIR <- als_choose(ALS_ARGS$outdir, CFG, "outputs.root", file.path(getwd(), "results", "01_create_files_for_prset"))
ALSPROJ_SCRIPT_ID <- "01_create_files_for_prset"
als_init(ALSPROJ_SCRIPT_ID, ALSPROJ_BASE_DIR, ALSPROJ_OUTDIR, ALS_ARGS$redirect_outputs)
als_register_output_wrappers()
als_load_objects(ALS_ARGS$objects)
on.exit(als_write_session_info(ALSPROJ_SCRIPT_ID), add = TRUE)
# ---- End bootstrap ----

# Build keep list with -b38 suffix to match VCF IDs
# PLINK format: FID IID (we'll use IID=IID, FID=0 since no family structure)
sporadic_cgnd <- sporadic_final_check %>%
  filter(!is.na(NYGC_CGND_ID)) %>%
  mutate(VCF_ID = paste0(NYGC_CGND_ID, "-b38"))

# Check — does sporadic_final_check have cluster_final?
colnames(sporadic_final_check)

# Write keep file (PLINK format: FID IID)
keep_df <- data.frame(
  FID = sporadic_cgnd$VCF_ID,
  IID = sporadic_cgnd$VCF_ID
)

write.table(keep_df, 
            file = "sporadic_keep.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

cat("Keep file written:", nrow(keep_df), "samples\n")

# Check ax_clean_df join key
colnames(ax_clean_df)
head(ax_clean_df)  # see what ID columns exist

# Attempt join
sporadic_with_clusters <- sporadic_final_check %>%
  left_join(
    ax_clean_df %>% select(sample_id, cluster_final),
    by = c("Sample_ID" = "sample_id")
  )

table(sporadic_with_clusters$cluster_final, useNA = "always")


# Check what's in prs_meta for these duplicated IDs
dup_ids <- c("CGND-HDA-01134-b38", "CGND-HDA-01160-b38", "CGND-HDA-03444-b38",
             "CGND-HDA-03445-b38", "CGND-HDA-01374-b38", "CGND-HDA-01372-b38",
             "CGND-HDA-01061-b38", "CGND-HDA-00187-b38", "CGND-HDA-01338-b38")

# prs_meta should have one row per patient — check
prs_meta[prs_meta$VCF_ID %in% dup_ids, c("VCF_ID", "Participant_ID", "cluster_final", 
                                         "AXIS_METABOLIC", "AXIS_CELL_CYCLE", "AXIS_INFLAMMATION")]


# How many unique Participant_IDs vs rows in prs_meta?
cat("Rows in prs_meta:", nrow(prs_meta), "\n")
cat("Unique Participant_IDs:", length(unique(prs_meta$Participant_ID)), "\n")
cat("Unique VCF_IDs:", length(unique(prs_meta$VCF_ID)), "\n")

ax_clean_df <- ax_clean_df %>%
  left_join(
    sporadic_meta %>% select(Sample_ID, Participant_ID),
    by = c("sample_id" = "Sample_ID")
  )
# For the subtype-discordant pairs specifically — which sample was used in final clustering?
# Check ax_clean_df which has the definitive cluster assignments
ax_clean_df[ax_clean_df$Participant_ID %in% 
              c("CASE-NEUAW157NMJ", "CASE-NEUBW008RJ5", "CASE-NEUEB422WW0",
                "CASE-NEUGE540TC4", "CASE-NEUHG644RYB", "CASE-NEUJK720JCL",
                "CASE-NEULZ548ZXV", "CASE-NEUMN012EVP", "CASE-NEUZY975XKL"),
            c("Participant_ID", "cluster_final")]


# What columns does ax_clean_df have?
names(ax_clean_df)
head(ax_clean_df, 3)

# How many total rows vs unique patients?
cat("Total rows:", nrow(ax_clean_df), "\n")
cat("Unique Participant_IDs:", length(unique(ax_clean_df$Participant_ID)), "\n")

# For the 9 duplicated patients — what distinguishes the two rows?
# (sample ID, timepoint, library batch?)
ax_clean_df[ax_clean_df$Participant_ID %in% 
              c("CASE-NEUAW157NMJ", "CASE-NEUBW008RJ5", "CASE-NEUEB422WW0",
                "CASE-NEUGE540TC4", "CASE-NEUHG644RYB", "CASE-NEUJK720JCL",
                "CASE-NEULZ548ZXV", "CASE-NEUMN012EVP", "CASE-NEUZY975XKL"), ]


# Extract timepoint number from sample_id suffix
ax_clean_df$timepoint <- as.integer(sub(".*-(\\d+)-T$", "\\1", ax_clean_df$sample_id))

# For each Participant_ID, keep the row with the lowest timepoint number
library(dplyr)
ax_clean_dedup <- ax_clean_df %>%
  group_by(Participant_ID) %>%
  slice_min(timepoint, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Rows before dedup:", nrow(ax_clean_df), "\n")
cat("Rows after dedup:", nrow(ax_clean_dedup), "\n")
cat("Unique patients:", length(unique(ax_clean_dedup$Participant_ID)), "\n")

# Verify the 9 duplicated patients now have one row each
ax_clean_dedup[ax_clean_dedup$Participant_ID %in%
                 c("CASE-NEUAW157NMJ", "CASE-NEUBW008RJ5", "CASE-NEUEB422WW0",
                   "CASE-NEUGE540TC4", "CASE-NEUHG644RYB", "CASE-NEUJK720JCL",
                   "CASE-NEULZ548ZXV", "CASE-NEUMN012EVP", "CASE-NEUZY975XKL"),
               c("Participant_ID", "sample_id", "cluster_final", "AXIS_METABOLIC")]


# Rebuild vcf_df from the saved VCF sample IDs file
vcf_ids_raw <- read.table(
  "vcf_sample_ids.txt",
  header = FALSE, stringsAsFactors = FALSE
)$V1

vcf_df <- data.frame(
  VCF_ID    = vcf_ids_raw,
  CGND_bare = sub("-b38$", "", vcf_ids_raw),
  stringsAsFactors = FALSE
)

# Deduplicate sporadic_final_check to one row per Participant_ID
# (just need covariates — AGE and SEX don't change between timepoints)
sporadic_meta_dedup <- sporadic_final_check %>%
  select(Participant_ID, NYGC_CGND_ID, AGE_AT_SYMPTOM_ONSET, SEX) %>%
  distinct(Participant_ID, .keep_all = TRUE)

cat("sporadic_meta_dedup rows:", nrow(sporadic_meta_dedup), "\n")
cat("Duplicate Participant_IDs:", sum(duplicated(sporadic_meta_dedup$Participant_ID)), "\n")

# Rebuild prs_meta_clean with deduped metadata
prs_meta_clean <- ax_clean_dedup %>%
  left_join(sporadic_meta_dedup, by = "Participant_ID") %>%
  mutate(CGND_bare = sub("-b38$", "", NYGC_CGND_ID)) %>%
  left_join(vcf_df, by = "CGND_bare") %>%
  filter(!is.na(VCF_ID))

cat("prs_meta_clean rows:", nrow(prs_meta_clean), "\n")
cat("Duplicate VCF_IDs:", sum(duplicated(prs_meta_clean$VCF_ID)), "\n")

# Phenotype file
pheno_out <- data.frame(
  FID = 0,
  IID = prs_meta_clean$VCF_ID,
  AXIS_METABOLIC    = prs_meta_clean$AXIS_METABOLIC,
  AXIS_CELL_CYCLE   = prs_meta_clean$AXIS_CELL_CYCLE,
  AXIS_INFLAMMATION = prs_meta_clean$AXIS_INFLAMMATION
)
cat("Duplicate IIDs in phenotype:", sum(duplicated(pheno_out$IID)), "\n")

write.table(pheno_out,
            "prset_phenotype.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# Rewrite covariate file without FID column
cov_out_nofid <- data.frame(
  IID = prs_meta_clean$VCF_ID,
  SEX = prs_meta_clean$SEX
)

write.table(cov_out_nofid,
            "prset_covariates.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# Same for phenotype file
pheno_out_nofid <- data.frame(
  IID = prs_meta_clean$VCF_ID,
  AXIS_METABOLIC    = prs_meta_clean$AXIS_METABOLIC,
  AXIS_CELL_CYCLE   = prs_meta_clean$AXIS_CELL_CYCLE,
  AXIS_INFLAMMATION = prs_meta_clean$AXIS_INFLAMMATION
)

write.table(pheno_out_nofid,
            "prset_phenotype.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

cat("Done\n")
head(cov_out_nofid, 3)
