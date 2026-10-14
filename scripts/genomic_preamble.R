################################################################################
## genomic_preamble.R
##
## Supplies the three objects genomic_analysis.R expects, rebuilt on the
## FINAL subtype labels. Run this from the directory holding the ANNOVAR
## ALS.chrN.anno_report_exon.txt files, then:
##
##   source("genomic_analysis.R")      # builds burden_all, rare_burden_sub
##   source("tric_burden.R")           # CCT2 / TRiC rare damaging burden
##   source("TRICBurdenMYCBurden.R")
##
## Provides:
##   sporadic_meta   Participant_ID + NYGC_CGND_ID + covariates
##   ax_clean_df     sample_id + cluster_clean (new labels)
##   df_nfl_uniq     sample_id + logNfL + clinical + covariates
##   cct2_df         Participant_ID + CCT2 (for tric_burden.R step 5)
##
## Note: genomic_analysis.R reads the ANNOVAR files from the working
## directory via list.files("."), so run it from where they live.
################################################################################

suppressPackageStartupMessages({ library(dplyr); library(data.table) })
set.seed(42)

BACKBONE <- "backbone_FINAL.rds"
PATIENTS <- "patient_level_FINAL.rds"
for (f in c(BACKBONE, PATIENTS))
  if (!file.exists(f)) stop("not found in working directory: ", f,
                            "\nCopy it here or set an absolute path.")

b    <- readRDS(BACKBONE)
pl   <- readRDS(PATIENTS)
meta <- as.data.frame(b$meta)
pat  <- as.data.frame(pl$pat)

cat("ANNOVAR files present:",
    length(list.files(".", pattern = "^ALS\\.chr.*anno_report_exon\\.txt$")), "\n")

## ---------------------------------------------------------------------------
## ax_clean_df  -- sample level, cluster_clean carries the new labels
## ---------------------------------------------------------------------------
ax_clean_df <- as.data.frame(b$ax_round1)
if (is.null(ax_clean_df$Participant_ID))
  ax_clean_df$Participant_ID <- meta$Participant_ID[
    match(ax_clean_df$sample_id, meta$sample_id)]
ax_clean_df$cluster_clean <- ax_clean_df$cluster_r1
ax_clean_df$cluster_final <- ax_clean_df$cluster_r1
cat("\nax_clean_df:", nrow(ax_clean_df), "sporadic samples\n")
print(table(ax_clean_df$cluster_clean))

## ---------------------------------------------------------------------------
## sporadic_meta  -- one row per sporadic participant, with the CGND ID
## ---------------------------------------------------------------------------
cg <- grep("CGND", colnames(meta), value = TRUE, ignore.case = TRUE)[1]
if (is.na(cg)) stop("no CGND column in b$meta")
cat("\nusing CGND column:", cg, "\n")

sp_pid <- unique(ax_clean_df$Participant_ID)
keep <- c("Participant_ID", cg, "SEX", "AGE_AT_SYMPTOM_ONSET",
          "ALSFRS_R_PROGRESSION_SLOPE", "AGE_AT_DEATH",
          paste0("SV", 1:6))
keep <- keep[keep %in% colnames(meta)]

sporadic_meta <- meta[meta$Participant_ID %in% sp_pid, keep, drop = FALSE]
sporadic_meta <- sporadic_meta[!duplicated(sporadic_meta$Participant_ID), ]
names(sporadic_meta)[names(sporadic_meta) == cg] <- "NYGC_CGND_ID"
sporadic_meta$NYGC_CGND_ID <- sub("-b[0-9]+$", "",
                                  as.character(sporadic_meta$NYGC_CGND_ID))

cat("sporadic_meta:", nrow(sporadic_meta), "participants,",
    sum(!is.na(sporadic_meta$NYGC_CGND_ID) &
          sporadic_meta$NYGC_CGND_ID != ""), "with a CGND ID\n")
cat("  Participant_ID example:", sporadic_meta$Participant_ID[1], "\n")
cat("  NYGC_CGND_ID example  :",
    sporadic_meta$NYGC_CGND_ID[!is.na(sporadic_meta$NYGC_CGND_ID)][1], "\n")

## ---------------------------------------------------------------------------
## df_nfl_uniq  -- sample-keyed clinical frame
## genomic_analysis.R strips -N-T from sample_id to get Participant_ID,
## so sample_id must be the full form.
## ---------------------------------------------------------------------------
clin <- c("SEX", "AGE_AT_SYMPTOM_ONSET", "ALSFRS_R_PROGRESSION_SLOPE",
          "COGNITIVE_PROGRESSION_SLOPE", paste0("SV", 1:6))
clin <- clin[clin %in% colnames(meta)]

df_nfl_uniq <- pat[, c("Participant_ID", "sample_id",
                       intersect(c("logNfL","CCT2"), colnames(pat))), drop = FALSE]
df_nfl_uniq <- merge(df_nfl_uniq,
                     unique(meta[, c("Participant_ID", clin)]),
                     by = "Participant_ID", all.x = TRUE)
df_nfl_uniq$cluster_clean <- pat$subtype[match(df_nfl_uniq$Participant_ID,
                                               pat$Participant_ID)]
df_nfl_uniq$cluster_final <- df_nfl_uniq$cluster_clean
df_nfl_uniq$group <- "Sporadic"

cat("\ndf_nfl_uniq:", nrow(df_nfl_uniq), "patients\n")
for (v in c("logNfL", "ALSFRS_R_PROGRESSION_SLOPE"))
  if (v %in% colnames(df_nfl_uniq))
    cat("  ", sprintf("%-30s %d\n", v, sum(!is.na(df_nfl_uniq[[v]]))))

## ---------------------------------------------------------------------------
## cct2_df  -- tric_burden.R step 5
## ---------------------------------------------------------------------------
cct2_df <- data.frame(Participant_ID = pat$Participant_ID,
                      CCT2 = pat$CCT2, stringsAsFactors = FALSE)
cct2_df <- cct2_df %>% distinct(Participant_ID, .keep_all = TRUE)
cat("\ncct2_df:", nrow(cct2_df), "patients\n")

## ---------------------------------------------------------------------------
## sanity: will the crosswalk actually match the ANNOVAR sample columns?
## ---------------------------------------------------------------------------
f1 <- list.files(".", pattern = "^ALS\\.chr1\\.anno_report_exon\\.txt$")[1]
if (!is.na(f1)) {
  hdr <- names(fread(f1, nrows = 0))
  scols <- grep("^(CASE|CTRL)-", hdr, value = TRUE)
  ov <- intersect(scols, sporadic_meta$Participant_ID)
  cat("\nANNOVAR sample columns:", length(scols),
      "| matching sporadic participants:", length(ov), "\n")
  if (length(ov) < 100)
    warning("low overlap: check that Participant_ID format matches the ",
            "annotation column names")
} else cat("\nALS.chr1.anno_report_exon.txt not in working directory\n")

cat("\n--------------------------------------------------------\n")
cat("Ready. Now source genomic_analysis.R, then tric_burden.R.\n")
cat("genomic_analysis.R may contain its own object-building block;\n")
cat("if it errors on a stopifnot, tell me which object it wants.\n")
cat("--------------------------------------------------------\n")