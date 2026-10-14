# =====================================================================
# 30_pathway_preamble.R   (rebuilt 2026-08-26)
#
# Supplies every object 31_ and 32_ need, from PINNED files only.
# Nothing is taken from the session.
#
#   source("30_pathway_preamble.R")
#   source("31_pathway_endpoint_screen.R")
#   source("32_database_independence.R")
#
# WHAT THE OLD pathway_go_preamble.R GOT WRONG
#   1. `universe` was never created, but mohamadFeedback1.R lists it as
#      required -- the script died at its own stopifnot on line 34.
#   2. `counts_filtered` likewise: only a warning was printed, then
#      feedback2 Step 2 failed on it.
#   3. clin_frame carried NO `group` column, so Phase 1's
#      filter(group == "Sporadic") errored.
#   4. df_als_only was set to a PATIENT-level frame (391 rows) while
#      Phase 2 joins sample-level ssGSEA scores by sample_id. Mixing
#      levels silently pseudoreplicates. This version is patient-level
#      END TO END, which is the convention the rest of the pipeline uses
#      after the eQTL n=369 vs 360 finding.
#   5. axis_sets was assumed to exist; it comes from backbone_FINAL.rds.
#
# SV NOTE: expr_resid is residualised on SEX + SV1-6, which is the
# CLUSTERING backbone specification. The genotype DE layer uses clean
# SVs 1,2,6,7,8,9 instead. The two are different by design; do not
# "reconcile" them.
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(msigdbr)
})
set.seed(42)

## ---- pinned inputs ----------------------------------------------------
find_first <- function(cands, label) {
  cands <- path.expand(cands)
  hit <- cands[file.exists(cands)]
  if (!length(hit)) stop(label, " not found. Tried:\n  ",
                         paste(cands, collapse = "\n  "))
  hit[1]
}
BACKBONE <- find_first(c("backbone_FINAL.rds",
                         "~/Documents/DSN/ALS/backbone_FINAL.rds"), "backbone_FINAL.rds")
PATIENTS <- find_first(c("patient_level_FINAL.rds",
                         "~/Documents/DSN/ALS/patient_level_FINAL.rds"), "patient_level_FINAL.rds")
DE_FILE  <- find_first(c("~/Documents/convergence_results/genotype_de_corrected_FINAL.rds",
                         "genotype_de_corrected_FINAL.rds"), "genotype_de_corrected_FINAL.rds")

b    <- readRDS(BACKBONE)
pl   <- readRDS(PATIENTS)
de   <- readRDS(DE_FILE)
meta <- as.data.frame(b$meta)
AXV  <- c("AXIS_METABOLIC", "AXIS_INFLAMMATION", "AXIS_CELL_CYCLE")
axis_vars <- AXV
cat("backbone:", BACKBONE, "\npatients:", PATIENTS, "\nDE      :", DE_FILE, "\n\n")

## ---- expression -------------------------------------------------------
expr_symbol <- b$expr_symbol
cat("expr_symbol:", paste(dim(expr_symbol), collapse = " x "), "\n")

sv_vars <- intersect(paste0("SV", 1:6), colnames(meta))
if (length(sv_vars) < 6)
  warning("only ", length(sv_vars), " SV columns found in meta: ",
          paste(sv_vars, collapse = ", "))

m2 <- meta[match(colnames(expr_symbol), meta$sample_id), , drop = FALSE]
stopifnot(nrow(m2) == ncol(expr_symbol), !any(is.na(m2$sample_id)))
cov1 <- model.matrix(as.formula(paste("~ SEX +", paste(sv_vars, collapse = " + "))),
                     data = m2)
P <- cov1 %*% solve(t(cov1) %*% cov1) %*% t(cov1)
expr_resid <- expr_symbol - expr_symbol %*% P
dimnames(expr_resid) <- dimnames(expr_symbol)
cat("expr_resid (SEX +", paste(sv_vars, collapse = "+"), "):",
    paste(dim(expr_resid), collapse = " x "), "\n")

expr_mat <- as.matrix(expr_symbol); storage.mode(expr_mat) <- "double"
expr_resid_mat <- as.matrix(expr_resid); storage.mode(expr_resid_mat) <- "double"

## ---- gene sets ---------------------------------------------------------
.msig <- function(...) {
  a <- list(...)
  tryCatch(do.call(msigdbr, c(list(species = "Homo sapiens"), a)),
           error = function(e) {
             names(a)[names(a) == "category"]    <- "collection"
             names(a)[names(a) == "subcategory"] <- "subcollection"
             do.call(msigdbr, c(list(species = "Homo sapiens"), a))
           })
}
h <- .msig(category = "H") %>% dplyr::select(gs_name, gene_symbol) %>% distinct()
hallmark_list <- split(h$gene_symbol, h$gs_name)
cat("hallmark_list:", length(hallmark_list), "sets\n")

if (is.null(b$axis_sets)) stop("backbone_FINAL.rds has no axis_sets")
axis_sets <- b$axis_sets
cat("axis_sets:", paste(names(axis_sets), lengths(axis_sets), sep = "=",
                        collapse = "  "), "\n")

## ---- universe: the DE gene space, version-stripped ---------------------
## feedback1 and feedback2 both need this and neither preamble made it.
res_Sporadic <- as.data.frame(de$res_Sporadic)
res_C9       <- as.data.frame(de$res_C9)
res_SOD1     <- as.data.frame(de$res_SOD1)
universe     <- rownames(res_Sporadic)
universe_ens <- sub("\\..*$", "", universe)
cat("universe:", length(universe), "genes (DE space)\n")

## aliases the feedback scripts refer to
res_C9_filtered <- res_C9_corrected <- res_C9
res_SOD1_filtered <- res_SOD1_corrected <- res_SOD1
res_Sporadic_filtered <- res_Sporadic_corrected <- res_Sporadic

## counts_filtered is only ever used for its ROWNAMES (the universe).
## Supply a zero-row stand-in with the right rownames rather than
## loading a 500 MB matrix nothing else needs.
counts_filtered <- matrix(numeric(0), nrow = length(universe), ncol = 0,
                          dimnames = list(universe, NULL))

## ---- subtype labels ----------------------------------------------------
ax_clean_df <- as.data.frame(b$ax_round1)
if (is.null(ax_clean_df$Participant_ID))
  ax_clean_df$Participant_ID <- meta$Participant_ID[
    match(ax_clean_df$sample_id, meta$sample_id)]
ax_clean_df$cluster_clean <- ax_clean_df$cluster_r1
ax_clean_df$cluster_final <- ax_clean_df$cluster_r1
ax_clean_df$timepoint <- as.integer(sub(".*-(\\d+)-T$", "\\1", ax_clean_df$sample_id))
cat("\nax_clean_df:", nrow(ax_clean_df), "sporadic samples\n")
print(table(ax_clean_df$cluster_clean))

ax_clean_dedup <- ax_clean_df %>%
  group_by(Participant_ID) %>% slice_min(timepoint, n = 1, with_ties = FALSE) %>%
  ungroup() %>% as.data.frame()
cat("ax_clean_dedup:", nrow(ax_clean_dedup), "patients\n")

## ---- differentiation markers -------------------------------------------
mk <- intersect(c("PCT_SMI32","PCT_ISL1","PCT_NKX61","PCT_TUJ1",
                  "PCT_S100b","PCT_Nestin"), colnames(meta))
marker_data <- meta[, c("sample_id", mk)]
for (m in mk) marker_data[[m]] <- suppressWarnings(as.numeric(marker_data[[m]]))
cat("\nmarker_data:", nrow(marker_data), "rows,", length(mk), "markers,",
    sum(complete.cases(marker_data[, mk])), "complete\n")

## ---- clinical frame, PATIENT LEVEL -------------------------------------
## One row per participant, earliest timepoint. `group` is set explicitly
## because this table is sporadic-only by construction; the old preamble
## omitted it and Phase 1 died on filter(group == "Sporadic").
pat <- as.data.frame(pl$pat)
pat$Participant_ID <- as.character(pat$Participant_ID)

extra <- c("COGNITIVE_PROGRESSION_SLOPE", sv_vars)
extra <- extra[extra %in% colnames(meta)]
addl <- meta[!duplicated(meta$Participant_ID), c("Participant_ID", extra), drop = FALSE]
clin_frame <- merge(pat, addl, by = "Participant_ID", all.x = TRUE)
clin_frame$group <- "Sporadic"
if (!"cluster_final" %in% names(clin_frame))
  clin_frame$cluster_final <- clin_frame$subtype
clin_frame$cluster_clean <- clin_frame$cluster_final

df_nfl_uniq <- clin_frame
df_als_only <- clin_frame

endpoints <- intersect(c("logNfL", "ALSFRS_R_PROGRESSION_SLOPE",
                         "COGNITIVE_PROGRESSION_SLOPE"), colnames(clin_frame))
cat("\nclin_frame:", nrow(clin_frame), "patients (PATIENT level, earliest timepoint)\n")
for (e in endpoints)
  cat(sprintf("  %-30s %d non-missing\n", e, sum(!is.na(clin_frame[[e]]))))
cat("  covariates: SEX, AGE_AT_SYMPTOM_ONSET,", paste(sv_vars, collapse = ", "), "\n")
if (!length(sv_vars))
  warning("no SV columns in the clinical frame; the endpoint models will be ",
          "unadjusted for surrogate variables")

## ---- all-group axis scores, for the familial-proximity phase -----------
## backbone$scores_r1 is 3 x 644 and spans Control/C9/SOD1/Sporadic, so
## axis_merged can be built here rather than left missing.
sc <- as.data.frame(scale(t(b$scores_r1)))
sc$sample_id <- rownames(sc)
gcols <- grep("^SUBJECT_GROUP|Mutation_Subtype|^group$", names(meta), value = TRUE)
if (length(gcols)) {
  gv  <- as.character(meta[[gcols[1]]][match(sc$sample_id, meta$sample_id)])
  ssg <- grep("^SUBJECT_SUBGROUP", names(meta), value = TRUE)[1]
  sv  <- if (!is.na(ssg)) as.character(meta[[ssg]][match(sc$sample_id, meta$sample_id)])
         else rep("", nrow(sc))
  sc$group <- ifelse(grepl("control", gv, ignore.case = TRUE), "Control",
              ifelse(grepl("c9", sv, ignore.case = TRUE), "C9orf72",
              ifelse(grepl("sod1", sv, ignore.case = TRUE), "SOD1",
              ifelse(grepl("als|mnd", gv, ignore.case = TRUE), "Sporadic", NA))))
  axis_merged <- sc
  cat("\naxis_merged:", nrow(axis_merged), "samples\n")
  print(table(axis_merged$group, useNA = "ifany"))
} else {
  axis_merged <- NULL
  cat("\naxis_merged not built: no genotype column in meta\n")
}

cat("\n---------------------------------------------------------------\n")
cat("Ready. Provided: expr_symbol, expr_resid, expr_mat, expr_resid_mat,\n")
cat("hallmark_list, axis_sets, axis_vars, ax_clean_df, ax_clean_dedup,\n")
cat("marker_data, universe, universe_ens, counts_filtered, clin_frame,\n")
cat("df_nfl_uniq, df_als_only, sv_vars, endpoints, axis_merged,\n")
cat("res_C9/_SOD1/_Sporadic (+ _filtered and _corrected aliases).\n")
cat("---------------------------------------------------------------\n")
