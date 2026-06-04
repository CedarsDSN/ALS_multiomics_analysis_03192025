# ============================================================
# Script: IPAPrepare.R
# Project: ALS multi-omics analysis
# Source status: documented/sanitized for GitHub handoff
# ============================================================
# Documentation summary from scripts/README.md
# 
# ## `IPAPrepare.R`
# 
# **Category:** IPA export / subtype DE
# 
# **Purpose**
# 
# Creates an IPA-ready DESeq2 result table for the S-ALS-3 vs S-ALS-2 contrast using clean SVs and sex adjustment. Also checks whether CCT2 elevation is uniform across S-ALS-3 or driven by outlier participants.
# 
# **Inputs**
# 
# - Workspace objects: `counts_filtered`, `annot_als_filtered`, `svobj`, `clean_svs`.
# - Subtype object: one of `ax_clean_dedup`, `ax_clean_df`, `sporadic_meta`, or `sporadic_final_check`, with a cluster column such as `cluster_final`.
# 
# **Outputs**
# 
# - `ipa_inputs/DEG_S3_vs_S2_IPA_upload.tsv`.
# - `ipa_inputs/DEG_S3_vs_S2_full.tsv`.
# - `ipa_inputs/DEG_S3_vs_S2_full.rds`.
# - `cct2_distribution_check.pdf`.
# 
# **Upstream dependencies**
# 
# - Must run after `AxisScoring.R` and `sporadic_subtyping.R`.
# 
# **How to run**
# 
# Legacy documented version:
# 
# ```bash
# Rscript scripts/documented/validation/IPAPrepare.R
# ```
# 
# Refactored command-line version:
# 
# ```bash
# Rscript scripts/refactored/validation/03_ipa_prepare.R --config config/config.yml --base-dir /path/to/private/ALS --outdir analysis/batch_04_validation_classifier_ipa/results/ipa_prepare --objects analysis/batch_03_axis_scoring_sporadic_subtyping/objects/03_axis_scoring_objects.rds,analysis/batch_03_axis_scoring_sporadic_subtyping/objects/01_sporadic_subtyping_objects.rds
# ```
# 
# **Refactor notes**
# 
# - Convert workspace object discovery into explicit RDS arguments.
# - Keep `ipa_inputs/` under `results/ipa_inputs/`.
# - Add an option to choose contrast, e.g. `--case S-ALS-3 --reference S-ALS-2`.
# ============================================================

# ============================================================
# DESeq2 contrast: S-ALS-3 vs S-ALS-2 for IPA Upstream Regulator
#
# Uses the same pipeline structure as Nick's existing axis-scoring
# workflow: clean SVs as covariates, SEX adjustment, sporadic-only.
#
# Inputs from workspace (variables from existing pipeline):
#   - counts_filtered   : gene x sample filtered raw counts
#   - annot_als_filtered: sample metadata (Sample_ID, Mutation_Subtype, SEX, etc.)
#   - svobj             : sva output (svobj$sv is samples x 10)
#   - clean_svs         : index vector of "clean" SVs (not correlated w/ Mutation_Subtype)
#   - ax_clean_dedup    : dedup frame with cluster_final assignments
#     (or whatever your subtype-assignment frame is called)
#
# Output:
#   ipa_inputs/DEG_S3_vs_S2_IPA_upload.tsv  -- upload to IPA
#   ipa_inputs/DEG_S3_vs_S2_full.tsv        -- full table for record
#   ipa_inputs/DEG_S3_vs_S2_full.rds        -- DESeq2 object + results
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

setwd(file.path(ALSPROJ_BASE_DIR))   # adjust if needed
OUT_DIR <- "ipa_inputs"
dir.create(OUT_DIR, showWarnings = FALSE)

# ---------------- Sanity checks -------------------------------
stopifnot(exists("counts_filtered"))
stopifnot(exists("annot_als_filtered"))
stopifnot(exists("svobj"))
stopifnot(exists("clean_svs"))

# Find the subtype assignment frame -- adjust this if yours is named differently
subtype_frame_candidates <- c("ax_clean_dedup", "ax_clean_df",
                              "sporadic_meta", "sporadic_final_check")
subtype_frame_name <- subtype_frame_candidates[
  sapply(subtype_frame_candidates, exists)
][1]
if (is.na(subtype_frame_name))
  stop("Cannot find subtype frame. Adjust subtype_frame_candidates ",
       "or load your subtype assignment object first.")
cat("Using subtype frame:", subtype_frame_name, "\n")
subtype_df <- get(subtype_frame_name)

# Find the cluster column
cluster_col <- intersect(c("cluster_final", "cluster_clean", "cluster"),
                         names(subtype_df))[1]
if (is.na(cluster_col))
  stop("No cluster column found in ", subtype_frame_name)
cat("Using cluster column:", cluster_col, "\n")

# Find the sample ID column in subtype_df
sid_col <- intersect(c("sample_id", "Sample_ID", "Participant_ID"),
                     names(subtype_df))[1]
if (is.na(sid_col))
  stop("No sample ID column in ", subtype_frame_name)
cat("Subtype frame keyed on:", sid_col, "\n\n")

# ---------------- Step 1: Build contrast metadata -----------
# Restrict to sporadic patients with S-ALS-2 or S-ALS-3 assignment.

subtype_clean <- subtype_df %>%
  as.data.frame() %>%
  dplyr::select(all_of(c(sid_col, cluster_col))) %>%
  rename(join_id = !!sid_col, cluster = !!cluster_col) %>%
  filter(cluster %in% c("S-ALS-2", "S-ALS-3")) %>%
  filter(!is.na(cluster))

cat("Subtype frame rows in S-ALS-2/3:", nrow(subtype_clean), "\n")

# Join to annot_als_filtered. The join key needs to match Sample_ID in annot.
# If subtype_df is keyed on Sample_ID, direct join. If keyed on Participant_ID,
# join through annot_als_filtered$Participant_ID.
if (sid_col == "Sample_ID") {
  contrast_meta <- annot_als_filtered %>%
    inner_join(subtype_clean, by = c("Sample_ID" = "join_id"))
} else if (sid_col == "Participant_ID") {
  contrast_meta <- annot_als_filtered %>%
    inner_join(subtype_clean, by = c("Participant_ID" = "join_id"))
} else if (sid_col == "sample_id") {
  contrast_meta <- annot_als_filtered %>%
    inner_join(subtype_clean, by = c("Sample_ID" = "join_id"))
}

# Sporadic only (subtype was defined on sporadic; this is a safety filter)
contrast_meta <- contrast_meta %>%
  filter(Mutation_Subtype == "Sporadic") %>%
  filter(!is.na(SEX))

# Deduplicate to one row per patient (earliest timepoint if longitudinal)
contrast_meta <- contrast_meta %>%
  distinct(Participant_ID, .keep_all = TRUE)

contrast_meta$cluster <- factor(contrast_meta$cluster,
                                levels = c("S-ALS-2", "S-ALS-3"))
contrast_meta$SEX <- factor(contrast_meta$SEX)

cat("\nFinal contrast metadata:\n")
print(table(contrast_meta$cluster))

# ---------------- Step 2: Subset count matrix ---------------
keep_samples <- intersect(contrast_meta$Sample_ID, colnames(counts_filtered))
cat("Samples with counts available:", length(keep_samples), "\n")

contrast_meta <- contrast_meta %>%
  filter(Sample_ID %in% keep_samples) %>%
  arrange(match(Sample_ID, keep_samples))

counts_sub <- counts_filtered[, contrast_meta$Sample_ID, drop = FALSE]
stopifnot(identical(colnames(counts_sub), contrast_meta$Sample_ID))

cat("Final n:", ncol(counts_sub), "samples,", nrow(counts_sub), "genes\n\n")
cat("Subtype distribution:\n")
print(table(contrast_meta$cluster))

# ---------------- Step 3: Attach clean SVs ------------------
# Pull the same clean SVs that the original full-cohort SVA gave us,
# subset to the patients in this contrast. Critically: we do NOT re-run
# SVA on the subset -- that would give different SVs and break consistency
# with the prior analyses.

# svobj$sv is samples x 10, in the order of annot_als_filtered (full cohort)
# We need to subset by row to match our contrast samples.

# Find the row indices in the original SVA matrix corresponding to our samples
orig_sample_order <- annot_als_filtered$Sample_ID
sv_row_idx <- match(contrast_meta$Sample_ID, orig_sample_order)
if (any(is.na(sv_row_idx)))
  stop("Some contrast samples not found in original SVA frame.")

sv_matrix <- svobj$sv[sv_row_idx, clean_svs, drop = FALSE]
colnames(sv_matrix) <- paste0("SV", seq_len(ncol(sv_matrix)))

# Attach to contrast_meta
for (i in seq_len(ncol(sv_matrix))) {
  contrast_meta[[colnames(sv_matrix)[i]]] <- sv_matrix[, i]
}

n_clean_sv <- ncol(sv_matrix)
cat("\nClean SVs attached:", n_clean_sv, "(SV1..SV", n_clean_sv, ")\n", sep = "")

# ---------------- Step 4: Build DESeq2 model ----------------
# Design: ~ SV1 + ... + SVk + SEX + cluster
# cluster LAST so the contrast extracts cleanly via name.

sv_terms <- paste0("SV", seq_len(n_clean_sv), collapse = " + ")
design_formula <- as.formula(
  paste("~", sv_terms, "+ SEX + cluster")
)
cat("Design formula:", deparse(design_formula), "\n\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData   = contrast_meta,
  design    = design_formula
)

# Light gene filter (same as your original pipeline: >=10 counts in >=10%)
min_count <- 10
min_samples <- round(0.1 * ncol(dds))
keep <- rowSums(counts(dds) >= min_count) >= min_samples
cat("Genes passing 10/10% pre-filter:", sum(keep), "of", nrow(dds), "\n")
dds <- dds[keep, ]

# ---------------- Step 5: Run DESeq2 ------------------------
cat("\nRunning DESeq2...\n")
dds <- DESeq(dds)
dds <- nbinomWaldTest(dds, maxit = 500)
n_unconverged <- sum(!mcols(dds)$betaConv)
if (n_unconverged > 0)
  cat("WARNING:", n_unconverged, "genes failed beta convergence\n")

# Extract the S-ALS-3 vs S-ALS-2 contrast
res <- results(
  dds,
  contrast = c("cluster", "S-ALS-3", "S-ALS-2"),
  independentFiltering = TRUE,
  alpha = 0.05
)

cat("\nDESeq2 summary:\n")
summary(res)

# Shrunken LFC for better ranking (apeglm preferred; ashr fallback)
coef_name <- "cluster_S.ALS.3_vs_S.ALS.2"
if (!coef_name %in% resultsNames(dds)) {
  # Some R versions use different separator
  coef_name <- grep("cluster.*S-ALS-3.*S-ALS-2|cluster.*S\\.ALS\\.3",
                    resultsNames(dds), value = TRUE)[1]
  cat("Using coef name:", coef_name, "\n")
}

res_shrunk <- tryCatch({
  lfcShrink(dds, coef = coef_name, type = "apeglm")
}, error = function(e) {
  message("apeglm failed: ", e$message, " -- falling back to ashr")
  lfcShrink(dds, contrast = c("cluster", "S-ALS-3", "S-ALS-2"), type = "ashr")
})

# ---------------- Step 6: ENSEMBL -> SYMBOL mapping ---------
# Your count matrix uses ENSEMBL IDs (per the rownames in your pipeline).
# IPA accepts symbols directly; using symbol avoids any ID-conversion
# issues inside IPA.

ens_ids <- rownames(res)
ens_clean <- sub("\\..*$", "", ens_ids)  # strip version suffix

symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = ens_clean,
  keytype  = "ENSEMBL",
  columns  = c("SYMBOL", "ENTREZID")
)
symbol_map <- symbol_map[!is.na(symbol_map$SYMBOL), ]
symbol_map <- symbol_map[!duplicated(symbol_map$ENSEMBL), ]

# Build results table
res_df <- data.frame(
  ensembl_id   = ens_clean,
  log2FC_raw   = res$log2FoldChange,
  log2FC_shrunk= res_shrunk$log2FoldChange,
  baseMean     = res$baseMean,
  lfcSE        = res$lfcSE,
  stat         = res$stat,
  pvalue       = res$pvalue,
  padj         = res$padj,
  stringsAsFactors = FALSE
)

# Attach symbols
res_df <- res_df %>%
  left_join(symbol_map %>%
              rename(ensembl_id = ENSEMBL,
                     gene_symbol = SYMBOL,
                     entrez_id = ENTREZID),
            by = "ensembl_id")

# Drop genes with no symbol or no padj
res_df <- res_df %>%
  filter(!is.na(padj), !is.na(gene_symbol), gene_symbol != "")

# If multiple ENSEMBL map to same symbol, keep most significant
res_df <- res_df %>%
  group_by(gene_symbol) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(pvalue)

cat("\nFinal IPA table:", nrow(res_df), "unique-symbol genes\n")
cat("Significant at padj < 0.05:", sum(res_df$padj < 0.05), "\n")
cat("  Up in S-ALS-3:  ",
    sum(res_df$log2FC_shrunk > 0 & res_df$padj < 0.05), "\n")
cat("  Down in S-ALS-3:",
    sum(res_df$log2FC_shrunk < 0 & res_df$padj < 0.05), "\n")

# ---------------- Step 7: Write IPA upload file -------------
# IPA's preferred format: tab-delimited, gene ID + Expr Log Ratio +
# Expr p-value + (optional) Expr FDR.

ipa_upload <- res_df %>%
  transmute(
    `ID`                                  = gene_symbol,
    `Expr Log Ratio`                      = round(log2FC_shrunk, 4),
    `Expr p-value`                        = signif(pvalue, 4),
    `Expr False Discovery Rate (q-value)` = signif(padj, 4)
  )

fwrite(ipa_upload,
       file.path(OUT_DIR, "DEG_S3_vs_S2_IPA_upload.tsv"),
       sep = "\t", quote = FALSE)

fwrite(res_df,
       file.path(OUT_DIR, "DEG_S3_vs_S2_full.tsv"),
       sep = "\t")

saveRDS(list(dds = dds,
             res = res,
             res_shrunk = res_shrunk,
             contrast_meta = contrast_meta,
             clean_svs_used = clean_svs),
        file.path(OUT_DIR, "DEG_S3_vs_S2_full.rds"))

# ---------------- Step 8: Sanity checks ---------------------
cat("\n========================================\n")
cat("SANITY CHECKS\n")
cat("========================================\n\n")

# CCT2 should be UP in S-ALS-3 with strong significance
cct2 <- res_df %>% filter(gene_symbol == "CCT2")
if (nrow(cct2) > 0) {
  cat("CCT2 (expected: UP in S-ALS-3, p < 1e-20):\n")
  cat("  log2FC (raw):   ", round(cct2$log2FC_raw, 3), "\n")
  cat("  log2FC (shrunk):", round(cct2$log2FC_shrunk, 3), "\n")
  cat("  p-value:        ", signif(cct2$pvalue, 3), "\n")
  cat("  padj:           ", signif(cct2$padj, 3), "\n\n")
  if (cct2$log2FC_shrunk < 0)
    cat("  WARNING: CCT2 is DOWN in S-ALS-3 -- check contrast direction!\n")
  if (cct2$padj > 0.01)
    cat("  WARNING: CCT2 padj weaker than expected -- check covariates\n")
} else {
  cat("WARNING: CCT2 not in results -- check ENSEMBL mapping\n\n")
}

# Top 15 up + down
cat("Top 15 UP in S-ALS-3 (by p-value):\n")
top_up <- res_df %>%
  filter(log2FC_shrunk > 0) %>%
  arrange(pvalue) %>%
  head(15) %>%
  dplyr::select(gene_symbol, log2FC_shrunk, pvalue, padj)
print(as.data.frame(top_up))

cat("\nTop 15 DOWN in S-ALS-3 (by p-value):\n")
top_dn <- res_df %>%
  filter(log2FC_shrunk < 0) %>%
  arrange(pvalue) %>%
  head(15) %>%
  dplyr::select(gene_symbol, log2FC_shrunk, pvalue, padj)
print(as.data.frame(top_dn))

cat("\n\nOutput files written to:", OUT_DIR, "\n")
cat("  DEG_S3_vs_S2_IPA_upload.tsv  <-- upload to IPA\n")
cat("  DEG_S3_vs_S2_full.tsv        <-- record for paper\n")
cat("  DEG_S3_vs_S2_full.rds        <-- full DESeq2 + metadata\n")






# ============================================================
# Is CCT2 expression in S-ALS-3 driven by a few outlier patients
# or uniformly elevated across the subtype?
#
# Two ways to fail "uniform elevation":
#   1. A small number of S-ALS-3 patients have extreme CCT2 while
#      the rest sit at S-ALS-2 levels (driver-patient pattern)
#   2. CCT2 has a bimodal distribution within S-ALS-3 itself
#      (sub-subtype possibility)
#
# Inputs:
#   - expr_symbol  : gene x sample VST-normalized matrix (symbol-keyed)
#   - df_adj OR ax_clean_dedup : has cluster_final + sample_id mapping
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# ---------------- Pull CCT2 expression with subtype ----------------
# Use whichever frame in your workspace has the subtype assignments
subtype_frame <- if (exists("ax_clean_dedup")) ax_clean_dedup else
  if (exists("df_adj")) df_adj else
    stop("Need a subtype assignment frame in workspace")

cluster_col <- intersect(c("cluster_final", "cluster_clean"),
                         names(subtype_frame))[1]
sid_col     <- intersect(c("Sample_ID", "sample_id"),
                         names(subtype_frame))[1]

cct2_by_pt <- data.frame(
  sample_id = colnames(expr_symbol),
  CCT2      = expr_symbol["CCT2", ],
  stringsAsFactors = FALSE
) %>%
  inner_join(
    subtype_frame %>%
      dplyr::select(all_of(sid_col), all_of(cluster_col)) %>%
      rename(sample_id = !!sid_col, subtype = !!cluster_col),
    by = "sample_id"
  ) %>%
  filter(subtype %in% c("S-ALS-2", "S-ALS-3")) %>%
  mutate(subtype = factor(subtype, levels = c("S-ALS-2", "S-ALS-3")))

cat("Patients per subtype:\n")
print(table(cct2_by_pt$subtype))

# ---------------- Distribution stats ----------------
stats <- cct2_by_pt %>%
  group_by(subtype) %>%
  summarise(
    n        = n(),
    mean     = mean(CCT2),
    median   = median(CCT2),
    sd       = sd(CCT2),
    min      = min(CCT2),
    q25      = quantile(CCT2, 0.25),
    q75      = quantile(CCT2, 0.75),
    max      = max(CCT2),
    iqr      = q75 - q25
  )
cat("\nCCT2 distribution by subtype:\n")
print(stats)

# Key check: what's the overlap?
# If S-ALS-3 is uniformly elevated, almost no S-ALS-3 patients should
# fall within the S-ALS-2 range.
s2 <- cct2_by_pt %>% filter(subtype == "S-ALS-2") %>% pull(CCT2)
s3 <- cct2_by_pt %>% filter(subtype == "S-ALS-3") %>% pull(CCT2)

s2_q90 <- quantile(s2, 0.90)
s2_q95 <- quantile(s2, 0.95)
s2_max <- max(s2)

s3_below_s2_q90 <- sum(s3 < s2_q90)
s3_below_s2_q95 <- sum(s3 < s2_q95)
s3_below_s2_max <- sum(s3 < s2_max)
s3_above_s2_max <- sum(s3 > s2_max)

cat("\n========================================\n")
cat("CRITICAL CHECK: Distributional overlap\n")
cat("========================================\n")
cat("S-ALS-2 90th percentile CCT2:", round(s2_q90, 3), "\n")
cat("S-ALS-2 max CCT2:            ", round(s2_max, 3), "\n\n")
cat("S-ALS-3 patients below S-ALS-2 90th pctile:",
    s3_below_s2_q90, "/", length(s3),
    sprintf("(%.1f%%)", 100*s3_below_s2_q90/length(s3)), "\n")
cat("S-ALS-3 patients below S-ALS-2 max:        ",
    s3_below_s2_max, "/", length(s3),
    sprintf("(%.1f%%)", 100*s3_below_s2_max/length(s3)), "\n")
cat("S-ALS-3 patients ABOVE S-ALS-2 max:        ",
    s3_above_s2_max, "/", length(s3),
    sprintf("(%.1f%%)", 100*s3_above_s2_max/length(s3)), "\n\n")

# If most S-ALS-3 patients sit above S-ALS-2's range, the elevation
# is uniform. If only a few S-ALS-3 patients are extreme and the rest
# overlap with S-ALS-2, the elevation is driver-patient.

# ---------------- Bimodality test within S-ALS-3 ----------------
# If CCT2 within S-ALS-3 is bimodal, that's evidence of a sub-subtype
if (requireNamespace("diptest", quietly = TRUE)) {
  dip <- diptest::dip.test(s3)
  cat("Hartigan's dip test for bimodality WITHIN S-ALS-3:\n")
  cat("  D =", round(dip$statistic, 4), "  p =", signif(dip$p.value, 3), "\n")
  cat("  (p < 0.05 = significantly non-unimodal)\n\n")
}

# ---------------- Driver-patient check ----------------
# If we drop the top N most extreme S-ALS-3 patients, does the contrast
# survive? This is the formal driver-patient test.
cat("Driver-patient test: drop top-N CCT2 patients, re-test mean diff\n")
cct2_s3_sorted <- sort(s3, decreasing = TRUE)
for (n_drop in c(0, 3, 5, 10, 15)) {
  if (n_drop >= length(s3)) next
  s3_trimmed <- cct2_s3_sorted[(n_drop+1):length(s3)]
  diff <- mean(s3_trimmed) - mean(s2)
  tt <- t.test(s3_trimmed, s2)
  cat(sprintf("  Drop top %2d S-ALS-3 patients: mean diff = %.3f, t = %.2f, p = %.2g, n_S3_remaining = %d\n",
              n_drop, diff, tt$statistic, tt$p.value, length(s3_trimmed)))
}

# ---------------- Distribution plots ----------------
pdf("cct2_distribution_check.pdf", width = 10, height = 8)

# Strip plot: every patient is a dot
p1 <- ggplot(cct2_by_pt, aes(x = subtype, y = CCT2, color = subtype)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_boxplot(fill = NA, outlier.shape = NA, color = "black", width = 0.3) +
  geom_hline(yintercept = s2_max, linetype = "dashed", color = "red") +
  annotate("text", x = 0.6, y = s2_max + 0.02,
           label = "S-ALS-2 max", color = "red", size = 3, hjust = 0) +
  labs(title = "Every patient as a point",
       subtitle = paste0("S-ALS-2 n=", length(s2),
                         ", S-ALS-3 n=", length(s3)),
       y = "CCT2 (VST normalized)") +
  theme_bw() + theme(legend.position = "none")

# Density overlay
p2 <- ggplot(cct2_by_pt, aes(x = CCT2, fill = subtype)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density distribution overlap",
       subtitle = "Uniform elevation = shifted distributions; driver pattern = long S-ALS-3 tail",
       x = "CCT2 (VST)") +
  theme_bw()

# Ranked CCT2 within S-ALS-3 (bimodality check)
p3 <- cct2_by_pt %>%
  filter(subtype == "S-ALS-3") %>%
  arrange(CCT2) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(x = rank, y = CCT2)) +
  geom_point(size = 2, color = "darkred") +
  geom_smooth(method = "loess", se = FALSE) +
  labs(title = "S-ALS-3 patients ranked by CCT2 (bimodality check)",
       subtitle = "Smooth curve = uniform; step/break = sub-subtype",
       x = "Patient rank within S-ALS-3", y = "CCT2 (VST)") +
  theme_bw()

# ECDF comparison
p4 <- ggplot(cct2_by_pt, aes(x = CCT2, color = subtype)) +
  stat_ecdf(size = 1.2) +
  labs(title = "Cumulative distribution",
       subtitle = "Parallel curves = uniform shift; crossing curves = mixture",
       x = "CCT2 (VST)", y = "Cumulative proportion") +
  theme_bw()

# Print all four
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()
cat("\nPlots saved to: cct2_distribution_check.pdf\n")

# ---------------- Summary verdict ----------------
cat("\n========================================\n")
cat("HOW TO INTERPRET\n")
cat("========================================\n")
cat("Verdict logic:\n")
cat("  - If >75% of S-ALS-3 patients sit above S-ALS-2's median  -> UNIFORM ELEVATION\n")
cat("  - If <50% of S-ALS-3 patients sit above S-ALS-2's max     -> DRIVER PATTERN\n")
cat("  - If dip test p < 0.05                                     -> SUB-SUBTYPE candidate\n")
cat("  - If contrast survives dropping top 10 S-ALS-3 patients    -> UNIFORM (not driver)\n")
