# 14_tals_rebuild.R   (CORRECTED 2026-08-26)
# Rebuild and SAVE res_C9_tals / res_SOD1_tals / res_Spor_tals.
#
# Why this exists: C9ConvergenceTargetALS.R contains no saveRDS anywhere
# in its 1604 lines, so the three TALS results objects only ever lived
# in the old session. This is the minimal deterministic path back to
# them, ending in a saved rds so this never has to run again.
#
# Cohort definition made explicit (the original file is ambiguous: the
# meta_de2 built at line 264 was filtered to Control + C9 only, so the
# literal top-to-bottom file cannot produce the SOD1 and Sporadic
# contrasts; the session must have used the later all-group rebuild).
# Encoded here, matching the dds_tals block (lines 1105-1153):
#   - Motor cortex (Lateral + Medial)
#   - group from c9_status plus conservative parsing of
#     reported_genomic_mutations for SOD1/C9 positives
#   - design ~ sex2 + group
#   - gene filter: >= 10 counts in >= 10% of samples
#
# FIXES IN THIS VERSION
#   1. ORDERING BUG: the original built tals_gene_cor_lfc BEFORE the
#      rna_tbl() calls that create tals_c9 / tals_sod1 / tals_sp, so it
#      died with "object 'tals_c9' not found" after the DE had already
#      run. The three tables are now built first.
#   2. The DE section is skipped if tals_de_results.rds already exists,
#      so a rerun after the failure does not repeat ~10 minutes of DESeq.
#   3. setwd() removed; paths are explicit so the script cannot move the
#      session out from under whatever called it.
#   4. is.finite(tstat) filter applied consistently before correlating.

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(DESeq2)
  library(BiocParallel)
})
register(MulticoreParam(8))
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::mutate,
                               dplyr::slice, base::intersect, base::union,
                               base::setdiff, base::unname, .quiet = TRUE)
}
set.seed(42)

## ---- 0. Paths (explicit; no setwd) ----
mainDir <- path.expand("~/Documents")
out_dir <- file.path(mainDir, "convergence_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

COUNTS_PARQUET <- file.path(mainDir, "counts_matrix.parquet")
SAMPLE_SHEET   <- file.path(mainDir, "sample_sheet.csv")
METADATA_CSV   <- file.path(mainDir,
                            "collections.postmortem_tissue_core.rnaseq_metadata-export-2026-01-27.csv")

RESULT_FILE <- file.path(out_dir, "tals_de_results.rds")

## =====================================================================
## PART A: differential expression (skipped if already saved)
## =====================================================================
if (file.exists(RESULT_FILE)) {
  cat("Found", RESULT_FILE, "-- loading instead of refitting.\n")
  tals <- readRDS(RESULT_FILE)
  res_C9_tals   <- tals$res_C9_tals
  res_SOD1_tals <- tals$res_SOD1_tals
  res_Spor_tals <- tals$res_Spor_tals
  cat("  design:", tals$design, "| cohort:", tals$cohort, "\n")
} else {
  
  ## ---- 1. Counts ----
  cat("Loading counts ...\n")
  mat_df <- as.data.frame(read_parquet(COUNTS_PARQUET))
  stopifnot("gene" %in% names(mat_df))
  count_mat <- as.matrix(mat_df[, setdiff(names(mat_df), "gene")])
  rownames(count_mat) <- mat_df$gene
  storage.mode(count_mat) <- "double"
  count_mat <- round(count_mat, 0)
  storage.mode(count_mat) <- "integer"
  
  ## ---- 2. Metadata ----
  samples <- fread(SAMPLE_SHEET)
  meta    <- fread(METADATA_CSV)
  setnames(meta, old = intersect("externalsampleid", names(meta)),
           new = "data_file_id")
  smeta <- as_tibble(merge(samples, meta, by = "data_file_id", all.x = TRUE))
  
  tissue_col <- if ("tissue.x" %in% names(smeta)) "tissue.x" else "tissue"
  
  smeta <- smeta %>%
    mutate(
      tissue_use = .data[[tissue_col]],
      condition = case_when(
        subject_group %in% c("ALS Spectrum MND",
                             "ALS Spectrum MND, Other Neurological Disorders") ~ "ALS",
        subject_group == "Non-Neurological Control" ~ "Control",
        TRUE ~ "Other"
      ),
      sex2 = case_when(
        tolower(as.character(sex)) == "male" ~ "Male",
        tolower(as.character(sex)) == "female" ~ "Female",
        TRUE ~ NA_character_
      ),
      c9x = trimws(tolower(c9orf72_repeat_expansion)),
      c9_status = case_when(
        c9x %in% c("yes", "abnormal repeat expansion") ~ "C9",
        c9x %in% c("no", "negative") ~ "NonC9",
        TRUE ~ NA_character_
      ),
      mut_raw = str_to_lower(str_trim(as.character(reported_genomic_mutations))),
      sod1_negative = str_detect(mut_raw, "negative for sod1"),
      sod1_positive = str_detect(mut_raw,
                                 "positive for sod1|sod1\\s*a4v|sod1\\s*\\(|sod1\\s*[a-z0-9]+") &
        !sod1_negative,
      c9_positive = (!is.na(c9_status) & c9_status == "C9") |
        (str_detect(mut_raw, "positive for c9orf72|c9orf72\\b|co9rf72") &
           !str_detect(mut_raw, "negative for c9orf72")),
      group = case_when(
        condition == "Control" ~ "Control",
        condition == "ALS" & c9_positive ~ "C9",
        condition == "ALS" & !c9_positive & sod1_positive ~ "SOD1",
        condition == "ALS" & !c9_positive & !sod1_positive ~ "Sporadic",
        TRUE ~ NA_character_
      )
    )
  
  ## ---- 3. Motor cortex cohort ----
  meta_tals <- smeta %>%
    filter(tissue_use %in% c("Motor Cortex Lateral", "Motor Cortex Medial")) %>%
    filter(!is.na(group), !is.na(sex2)) %>%
    mutate(group = factor(group, levels = c("Control", "C9", "SOD1", "Sporadic")),
           sex2 = factor(sex2))
  
  cat("\nGroup counts (motor cortex, samples):\n")
  print(table(meta_tals$group))
  cat("\nGroup x sex:\n")
  print(table(meta_tals$group, meta_tals$sex2))
  
  common <- intersect(colnames(count_mat), meta_tals$data_file_id)
  count_use <- count_mat[, common, drop = FALSE]
  meta_use  <- meta_tals[match(common, meta_tals$data_file_id), , drop = FALSE]
  stopifnot(identical(colnames(count_use), meta_use$data_file_id))
  meta_use <- droplevels(meta_use)
  cat("\nSamples with counts:", ncol(count_use), "\n")
  print(table(meta_use$group))
  stopifnot(all(table(meta_use$group) > 0))
  
  ## Gene filter: >= 10 counts in >= 10% of samples
  keep <- rowSums(count_use >= 10) >= round(0.1 * ncol(count_use))
  count_use <- count_use[keep, , drop = FALSE]
  cat("Genes kept:", nrow(count_use), "\n")
  
  ## ---- 4. DE ----
  dds_tals <- DESeqDataSetFromMatrix(countData = count_use,
                                     colData = as.data.frame(meta_use),
                                     design = ~ sex2 + group)
  dds_tals <- DESeq(dds_tals, parallel = TRUE)
  
  res_C9_tals   <- results(dds_tals, contrast = c("group", "C9", "Control"))
  res_SOD1_tals <- results(dds_tals, contrast = c("group", "SOD1", "Control"))
  res_Spor_tals <- results(dds_tals, contrast = c("group", "Sporadic", "Control"))
  
  saveRDS(list(res_C9_tals = res_C9_tals,
               res_SOD1_tals = res_SOD1_tals,
               res_Spor_tals = res_Spor_tals,
               meta_use = as.data.frame(meta_use),
               design = "~ sex2 + group",
               cohort = "Motor Cortex Lateral+Medial",
               gene_filter = ">=10 counts in >=10% of samples"),
          RESULT_FILE)
  cat("\nSaved", RESULT_FILE, "\n")
}

## =====================================================================
## PART B: gene-level correlations
## ORDER MATTERS: build the three tables BEFORE any correlation block.
## =====================================================================
strip_ens_version <- function(x) sub("\\..*$", "", x)

rna_tbl <- function(res, label) {
  df <- as.data.frame(res)
  df$gene <- strip_ens_version(rownames(df))
  df %>%
    transmute(gene, logFC = log2FoldChange, tstat = stat,
              pval = pvalue, FDR = padj, contrast = label) %>%
    filter(!is.na(gene) & gene != "" & is.finite(tstat))
}

tals_c9   <- rna_tbl(res_C9_tals,   "C9_vs_CTRL")
tals_sod1 <- rna_tbl(res_SOD1_tals, "SOD1_vs_CTRL")
tals_sp   <- rna_tbl(res_Spor_tals, "SPORADIC_vs_CTRL")
cat("\ngenes after filtering:",
    paste(sapply(list(tals_c9, tals_sod1, tals_sp), nrow), collapse = " / "), "\n")

pair_cor <- function(df1, df2, col = "tstat") {
  m <- inner_join(df1 %>% dplyr::select(gene, s1 = all_of(col)),
                  df2 %>% dplyr::select(gene, s2 = all_of(col)), by = "gene")
  cor(m$s1, m$s2, method = "spearman", use = "pairwise.complete.obs")
}
topK_dir_agree <- function(df1, df2, K = 500) {
  d1 <- df1 %>% arrange(desc(abs(tstat))) %>% slice(1:K) %>%
    transmute(gene, s1 = sign(tstat))
  d2 <- df2 %>% transmute(gene, s2 = sign(tstat))
  m <- inner_join(d1, d2, by = "gene")
  mean(m$s1 == m$s2)
}

## PRIMARY metric: logFC, matching the manuscript convention.
tals_gene_cor_lfc <- tibble(
  pair = c("C9 vs SOD1", "C9 vs Sporadic", "SOD1 vs Sporadic"),
  spearman_logFC = c(pair_cor(tals_c9, tals_sod1, "logFC"),
                     pair_cor(tals_c9, tals_sp,   "logFC"),
                     pair_cor(tals_sod1, tals_sp, "logFC"))
)
cat("\nlogFC-based (primary):\n")
print(tals_gene_cor_lfc)
write.csv(tals_gene_cor_lfc,
          file.path(out_dir, "tals_gene_level_correlations_logFC.csv"),
          row.names = FALSE)

## SECONDARY: Wald stat. Kept for continuity, not for the paper.
tals_gene_cor <- tibble(
  pair = c("C9 vs SOD1", "C9 vs Sporadic", "SOD1 vs Sporadic"),
  spearman = c(pair_cor(tals_c9, tals_sod1),
               pair_cor(tals_c9, tals_sp),
               pair_cor(tals_sod1, tals_sp)),
  top500_dir_agree = c(topK_dir_agree(tals_c9, tals_sod1),
                       topK_dir_agree(tals_c9, tals_sp),
                       topK_dir_agree(tals_sod1, tals_sp))
)
cat("\nWald-stat based (secondary) -- handover targets",
    "C9~SOD1 0.013, C9~Sporadic 0.289, SOD1~Sporadic 0.260:\n")
print(tals_gene_cor)
write.csv(tals_gene_cor,
          file.path(out_dir, "tals_gene_level_correlations.csv"),
          row.names = FALSE)

cat("\nNOTE: the end-stage TALS values (0.487 / 0.751 / 0.751 in the\n")
cat("handover) are the ones that motivate the cell-type-matched design.\n")
cat("Check which metric and ordering those refer to before quoting them.\n")
cat("\nres_C9_tals / res_SOD1_tals / res_Spor_tals are now in the session,\n")
cat("so Section 11 of 16_convergence_corrected.R will run.\n")