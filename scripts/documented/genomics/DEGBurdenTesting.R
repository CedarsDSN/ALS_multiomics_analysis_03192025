# ============================================================
# Script: DEGBurdenTesting.R
# Project: ALS multi-omics analysis
# Source status: documented/sanitized for GitHub handoff
# ============================================================
# Documentation summary from docs/SCRIPT_SPECS_FULL.md
# 
# ## `DEGBurdenTesting.R`
# 
# **Category:** Genomics / DEG-informed burden
# 
# **Purpose**
# 
# Tests whether rare damaging variant burden in subtype-defining DEG sets is enriched in the matching subtype and decomposes LoF-to-NfL associations by axis gene-set burden.
# 
# **Inputs**
# 
# - DESeq2 result objects for subtype contrasts: `res_1v2`, `res_1v3`, `res_2v3`.
# - Gene burden objects from `genomic_analysis.R`: `gene_burden_combined` and/or `gene_burden_filt`.
# - Subtype/clinical metadata including `rare_burden_sub`, subtype labels, NfL, age, and sex.
# - Annotation mapping through `org.Hs.eg.db`.
# 
# **Outputs**
# 
# - `targeted_deg_burden.pdf`.
# - `axis_burden_nfl_forest.pdf`.
# - `axis_burden_nfl_scatter.pdf`.
# - `axis_decomposed_lof_nfl.csv`.
# 
# **Upstream dependencies**
# 
# - Requires `sporadic_subtyping.R` subtype DE objects and `genomic_analysis.R` burden matrices.
# 
# **How to run**
# 
# Legacy documented version:
# 
# ```bash
# Rscript scripts/documented/genomics/DEGBurdenTesting.R
# ```
# 
# Refactored command-line version:
# 
# ```bash
# Rscript scripts/refactored/genomics/04_deg_burden_testing.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/genomics/04_deg_burden_testing --objects outputs/objects/axis_scoring_objects.rds,outputs/objects/sporadic_subtyping_objects.rds
# ```
# 
# **Refactor notes**
# 
# - Save subtype DESeq2 results from `sporadic_subtyping.R` as RDS/CSV and load them explicitly.
# - Unify gene ID conversion before burden filtering.
# - Return all model results as tidy tables, not only plots.
# ============================================================

# ============================================================
# TARGETED GENOMIC BURDEN IN SUBTYPE DEG SETS
# Teammate Comment 1 — transcriptomically-informed variant test
# ============================================================

library(DESeq2)
library(dplyr)
library(ggplot2)

# ---- STEP 1: Get subtype DEGs from DESeq2 results ----------
# Assumes you have DESeq2 results objects already
# res_1v2, res_1v3, res_2v3 — standard DESeq2 result objects
# If not, load them:
# res_1v2 <- results(dds, contrast = c("cluster_final", "S-ALS-1", "S-ALS-2"))

# Extract significant DEGs per comparison (FDR < 0.05, |LFC| > 0.5)
get_degs <- function(res, name) {
  as.data.frame(res) |>
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5) |>
    mutate(comparison = name,
           gene = rownames(res)[!is.na(res$padj) & 
                                  res$padj < 0.05 & 
                                  abs(res$log2FoldChange) > 0.5])
}

# These are the objects to plug in as res_1v2, res_1v3, res_2v3
res_1v2 <- res_final_12   # S-ALS-1 vs S-ALS-2
res_1v3 <- res_final_13   # S-ALS-1 vs S-ALS-3
res_2v3 <- res_final_23   # S-ALS-2 vs S-ALS-3

degs_1v2 <- get_degs(res_1v2, "S-ALS-1_vs_S-ALS-2")
degs_1v3 <- get_degs(res_1v3, "S-ALS-1_vs_S-ALS-3")
degs_2v3 <- get_degs(res_2v3, "S-ALS-2_vs_S-ALS-3")

# Subtype-specific gene sets
# Genes that distinguish each subtype from BOTH others
sals1_genes <- intersect(
  rownames(degs_1v2[degs_1v2$log2FoldChange > 0, ]),  # up in S-ALS-1 vs 2
  rownames(degs_1v3[degs_1v3$log2FoldChange > 0, ])   # up in S-ALS-1 vs 3
)

sals2_genes <- intersect(
  rownames(degs_1v2[degs_1v2$log2FoldChange < 0, ]),  # up in S-ALS-2 vs 1
  rownames(degs_2v3[degs_2v3$log2FoldChange > 0, ])   # up in S-ALS-2 vs 3
)

sals3_genes <- intersect(
  rownames(degs_1v3[degs_1v3$log2FoldChange < 0, ]),  # up in S-ALS-3 vs 1
  rownames(degs_2v3[degs_2v3$log2FoldChange < 0, ])   # up in S-ALS-3 vs 2
)

cat("S-ALS-1 specific DEGs:", length(sals1_genes), "\n")
cat("S-ALS-2 specific DEGs:", length(sals2_genes), "\n")
cat("S-ALS-3 specific DEGs:", length(sals3_genes), "\n")

# ---- STEP 2: Filter genomic burden matrix to DEG sets -------
# rare_burden_all: rows = samples, columns include gene-level counts
# Assumes gene-level damaging/LoF counts are available
# If your burden matrix is sample x gene, filter columns

# Load gene-level burden matrix (samples x genes)
# This should be your per-gene rare damaging count matrix
# gene_burden_matrix: rows = Participant_ID, cols = gene symbols

# ============================================================
# USE EXISTING gene_burden_combined FROM YOUR PIPELINE
# ============================================================

# gene_burden_combined is genes x samples (from parse_gene_burden_chr)
# Need: samples x genes, with Participant_IDs as rownames

# Step 1: Transpose to samples x genes
gene_burden_t <- as.data.frame(t(gene_burden_combined))
# rownames are CASE-NEU col_ids at this point
cat("Dimensions:", nrow(gene_burden_t), "samples x",
    ncol(gene_burden_t), "genes\n")

# Step 2: Map col_id → Participant_ID using your crosswalk
# col_to_id maps col_id → cgnd_id
# crosswalk maps cgnd_id → expr_sample_id (Participant_ID)
cgnd_to_participant <- setNames(
  crosswalk$expr_sample_id,   # Participant_ID (CASE-NEU short form)
  crosswalk$cgnd_id
)

# Current rownames are col_ids (CASE-NEU format from annotation file)
# These ARE the Participant_IDs already — check:
cat("First 3 rownames of gene_burden_t:\n")
print(head(rownames(gene_burden_t), 3))
cat("First 3 sporadic_final_check Participant_IDs:\n")
print(head(sporadic_final_check$Participant_ID, 3))

# Direct match check
direct_match <- intersect(rownames(gene_burden_t),
                          sporadic_final_check$Participant_ID)
cat("Direct matches:", length(direct_match), "\n")

# If direct match works, filter to sporadic samples
gene_burden_matrix <- gene_burden_t[
  rownames(gene_burden_t) %in% sporadic_final_check$Participant_ID, ]

# Add zero rows for any missing sporadic patients
missing_pts <- setdiff(sporadic_final_check$Participant_ID,
                       rownames(gene_burden_matrix))
if (length(missing_pts) > 0) {
  cat("Adding", length(missing_pts), "zero-variant patients\n")
  zero_mat <- matrix(0L,
                     nrow = length(missing_pts),
                     ncol = ncol(gene_burden_matrix),
                     dimnames = list(missing_pts, colnames(gene_burden_matrix))
  )
  gene_burden_matrix <- rbind(gene_burden_matrix, zero_mat)
}

cat("gene_burden_matrix ready:",
    nrow(gene_burden_matrix), "patients x",
    ncol(gene_burden_matrix), "genes\n")

# ============================================================
# CONVERT DEG ENSEMBL IDs TO GENE SYMBOLS
# DESeq2 returns Ensembl IDs; ANNOVAR uses gene symbols
# ============================================================

library(AnnotationDbi)
library(org.Hs.eg.db)

convert_to_symbols <- function(ensembl_ids) {
  ids_clean <- sub("\\..*$", "", ensembl_ids)
  mapped <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = unique(ids_clean),
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  )
  mapped <- mapped[!is.na(mapped$SYMBOL), ]
  mapped <- mapped[!duplicated(mapped$ENSEMBL), ]
  unique(mapped$SYMBOL)
}

cat("Converting DEG Ensembl IDs to gene symbols...\n")
sals1_symbols <- convert_to_symbols(sals1_genes)
sals2_symbols <- convert_to_symbols(sals2_genes)
sals3_symbols <- convert_to_symbols(sals3_genes)

cat("S-ALS-1:", length(sals1_genes), "Ensembl →",
    length(sals1_symbols), "symbols\n")
cat("S-ALS-2:", length(sals2_genes), "Ensembl →",
    length(sals2_symbols), "symbols\n")
cat("S-ALS-3:", length(sals3_genes), "Ensembl →",
    length(sals3_symbols), "symbols\n")

# How many of these symbols are in your burden matrix?
cat("\nDEG symbols present in gene_burden_matrix:\n")
cat("S-ALS-1:", sum(sals1_symbols %in% colnames(gene_burden_matrix)),
    "/", length(sals1_symbols), "\n")
cat("S-ALS-2:", sum(sals2_symbols %in% colnames(gene_burden_matrix)),
    "/", length(sals2_symbols), "\n")
cat("S-ALS-3:", sum(sals3_symbols %in% colnames(gene_burden_matrix)),
    "/", length(sals3_symbols), "\n")

# ============================================================
# NOW RUN filter_burden_to_genes WITH SYMBOL GENE SETS
# ============================================================

burden_sals1_genes <- filter_burden_to_genes(
  gene_burden_matrix, sals1_symbols, "SALS1_DEGs")
burden_sals2_genes <- filter_burden_to_genes(
  gene_burden_matrix, sals2_symbols, "SALS2_DEGs")
burden_sals3_genes <- filter_burden_to_genes(
  gene_burden_matrix, sals3_symbols, "SALS3_DEGs")

# ---- STEP 3: Merge with subtype labels ----------------------
targeted_burden <- rare_burden_sub |>  # your existing df with cluster_final
  left_join(burden_sals1_genes, by = "Participant_ID") |>
  left_join(burden_sals2_genes, by = "Participant_ID") |>
  left_join(burden_sals3_genes, by = "Participant_ID")

# ---- STEP 4: Test — does each subtype carry more burden 
#              in its own DEG genes? --------------------------
# This is the key question: does S-ALS-3 carry more rare variants
# specifically in S-ALS-3-defining genes?

run_targeted_test <- function(df, burden_col, subtype_label) {
  df$is_target <- as.integer(df$cluster_final == subtype_label)
  df$burden <- df[[burden_col]]
  
  # ANOVA across all three subtypes
  aov_result <- aov(burden ~ cluster_final, data = df)
  aov_p <- summary(aov_result)[[1]][["Pr(>F)"]][1]
  
  # Wilcoxon: target subtype vs all others
  target <- df$burden[df$is_target == 1]
  others <- df$burden[df$is_target == 0]
  wt <- wilcox.test(target, others)
  
  cat("\n---", subtype_label, "DEG-targeted burden ---\n")
  cat("ANOVA p:", round(aov_p, 4), "\n")
  cat("Wilcoxon (target vs rest) p:", round(wt$p.value, 4), "\n")
  cat("Median burden in target subtype:", median(target, na.rm=TRUE), "\n")
  cat("Median burden in other subtypes:", median(others, na.rm=TRUE), "\n")
}

run_targeted_test(targeted_burden, "burden_SALS1_DEGs", "S-ALS-1")
run_targeted_test(targeted_burden, "burden_SALS2_DEGs", "S-ALS-2")
run_targeted_test(targeted_burden, "burden_SALS3_DEGs", "S-ALS-3")

# ---- STEP 5: Visualise ------------------------------------
targeted_long <- targeted_burden |>
  select(Participant_ID, cluster_final, 
         burden_SALS1_DEGs, burden_SALS2_DEGs, burden_SALS3_DEGs) |>
  tidyr::pivot_longer(
    cols = starts_with("burden_"),
    names_to = "gene_set",
    values_to = "burden"
  ) |>
  mutate(gene_set = gsub("burden_", "", gene_set))

ggplot(targeted_long, aes(x = cluster_final, y = log1p(burden), 
                          fill = cluster_final)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~gene_set, scales = "free_y") +
  scale_fill_manual(values = c("S-ALS-1"="#E41A1C",
                               "S-ALS-2"="#377EB8",
                               "S-ALS-3"="#4DAF4A")) +
  labs(title = "Targeted genomic burden in subtype-specific DEG sets",
       x = "Subtype", y = "log(1 + rare damaging burden)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("targeted_deg_burden.pdf", width = 10, height = 4)







# ══════════════════════════════════════════════════════════════════════════════
# PATHWAY-DECOMPOSED LoF → NfL ANALYSIS
# Goal: Does the pan-sporadic LoF → NfL survival bias signal (β=-0.40, p=0.0004)
#       concentrate in a specific biological axis, or is it evenly distributed?
# ══════════════════════════════════════════════════════════════════════════════

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ── 1. Check whether a LoF-specific gene-level matrix exists ─────────────────
# gene_burden_matrix is confirmed as rare damaging (391 × 18808)
# Check if a parallel LoF matrix was built

cat("Checking for LoF-specific gene burden matrix...\n")
lof_matrix_candidates <- c("gene_lof_matrix", "lof_burden_matrix", 
                           "gene_burden_lof", "lof_mat")
for (obj in lof_matrix_candidates) {
  if (exists(obj)) cat("  FOUND:", obj, "—", dim(get(obj)), "\n")
}

# Also check gene_burden_filt and gene_burden_combined
if (exists("gene_burden_filt")) {
  cat("gene_burden_filt dims:", dim(gene_burden_filt), "\n")
  cat("  rownames sample:", head(rownames(gene_burden_filt)), "\n")
}
if (exists("gene_burden_combined")) {
  cat("gene_burden_combined dims:", dim(gene_burden_combined), "\n")
  cat("  colnames:", head(colnames(gene_burden_combined)), "\n")
}

# ── 2. Extract axis gene sets and check overlap with gene_burden_matrix ───────
metabolic_genes  <- axis_sets_clean[["AXIS_METABOLIC"]]
cellcycle_genes  <- axis_sets_clean[["AXIS_CELL_CYCLE"]]
inflam_genes     <- axis_sets_clean[["AXIS_INFLAMMATION"]]

# How many of each axis are in gene_burden_matrix?
gb_genes <- colnames(gene_burden_matrix)

met_avail  <- intersect(metabolic_genes, gb_genes)
cc_avail   <- intersect(cellcycle_genes, gb_genes)
inf_avail  <- intersect(inflam_genes,    gb_genes)

cat(sprintf("\nAxis gene availability in gene_burden_matrix:\n"))
cat(sprintf("  AXIS_METABOLIC:    %d / %d genes (%.1f%%)\n", 
            length(met_avail),  length(metabolic_genes),
            100 * length(met_avail) / length(metabolic_genes)))
cat(sprintf("  AXIS_CELL_CYCLE:   %d / %d genes (%.1f%%)\n",
            length(cc_avail),   length(cellcycle_genes),
            100 * length(cc_avail) / length(cellcycle_genes)))
cat(sprintf("  AXIS_INFLAMMATION: %d / %d genes (%.1f%%)\n",
            length(inf_avail),  length(inflam_genes),
            100 * length(inf_avail) / length(inflam_genes)))

# ── 3. Compute per-patient per-axis rare damaging burden ──────────────────────
# gene_burden_matrix rows = patients (confirmed: CASE-NEU... IDs)

cat("\nComputing per-axis burden from gene_burden_matrix...\n")

axis_burden_df <- data.frame(
  Participant_ID     = rownames(gene_burden_matrix),
  burden_met_dam     = rowSums(gene_burden_matrix[, met_avail, drop = FALSE]),
  burden_cc_dam      = rowSums(gene_burden_matrix[, cc_avail,  drop = FALSE]),
  burden_inf_dam     = rowSums(gene_burden_matrix[, inf_avail, drop = FALSE])
)

cat("Axis burden summary (rare damaging):\n")
print(summary(axis_burden_df[, -1]))

# ── 4. Merge with targeted_burden (has logNfL, rare_lof, cluster, covariates) ─
# targeted_burden already has Participant_ID confirmed

analysis_df <- targeted_burden |>
  select(Participant_ID, cluster_final, logNfL,
         rare_lof, rare_damaging,
         AGE_AT_SYMPTOM_ONSET, SEX) |>
  inner_join(axis_burden_df, by = "Participant_ID")

cat(sprintf("\nMerged dataset: %d patients total\n", nrow(analysis_df)))
cat(sprintf("With non-NA NfL: %d\n", sum(!is.na(analysis_df$logNfL))))
cat(sprintf("Subtype distribution:\n"))
print(table(analysis_df$cluster_final))

# Filter to NfL-complete for modelling
nfl_df <- analysis_df |> filter(!is.na(logNfL))
cat(sprintf("\nModelling on n = %d with complete NfL\n", nrow(nfl_df)))

# ── 5. Replicate pan-sporadic signal first ────────────────────────────────────
# Confirm β ≈ -0.40, p ≈ 0.0004 before decomposing

m_lof_total <- lm(logNfL ~ rare_lof + AGE_AT_SYMPTOM_ONSET + SEX,
                  data = nfl_df)
m_dam_total <- lm(logNfL ~ rare_damaging + AGE_AT_SYMPTOM_ONSET + SEX,
                  data = nfl_df)

cat("\n─── REPLICATION CHECK ───────────────────────────────────────\n")
cat(sprintf("Total rare LoF     → NfL: β = %+.3f, p = %.4f\n",
            coef(m_lof_total)["rare_lof"],
            summary(m_lof_total)$coefficients["rare_lof", "Pr(>|t|)"]))
cat(sprintf("Total rare damaging → NfL: β = %+.3f, p = %.4f\n",
            coef(m_dam_total)["rare_damaging"],
            summary(m_dam_total)$coefficients["rare_damaging", "Pr(>|t|)"]))
cat("─────────────────────────────────────────────────────────────\n")

# ── 6. Per-axis models (rare damaging) → NfL ─────────────────────────────────

run_axis_model <- function(burden_col, axis_label, df) {
  f <- as.formula(paste("logNfL ~", burden_col, 
                        "+ AGE_AT_SYMPTOM_ONSET + SEX"))
  m <- lm(f, data = df)
  s <- summary(m)$coefficients
  data.frame(
    axis      = axis_label,
    n_genes   = c(met = length(met_avail), 
                  cc  = length(cc_avail), 
                  inf = length(inf_avail))[
                    c(met = "met", cc = "cc", inf = "inf")[
                      sub("burden_(.*)_dam", "\\1", burden_col)]],
    beta      = round(s[burden_col, "Estimate"],   4),
    se        = round(s[burden_col, "Std. Error"], 4),
    t_stat    = round(s[burden_col, "t value"],    3),
    p_value   = round(s[burden_col, "Pr(>|t|)"],  4),
    r2_adj    = round(summary(lm(f, data = df))$adj.r.squared, 4)
  )
}

axis_results <- bind_rows(
  run_axis_model("burden_met_dam", "AXIS_METABOLIC",    nfl_df),
  run_axis_model("burden_inf_dam", "AXIS_INFLAMMATION", nfl_df),
  run_axis_model("burden_cc_dam",  "AXIS_CELL_CYCLE",   nfl_df)
)

# Fix n_genes column (above indexing is a bit awkward — do it cleanly)
axis_results$n_genes <- c(length(met_avail), length(inf_avail), length(cc_avail))

cat("\n══════════════════════════════════════════════════════════════\n")
cat("PER-AXIS RARE DAMAGING BURDEN → NfL\n")
cat("══════════════════════════════════════════════════════════════\n")
print(axis_results[, c("axis", "n_genes", "beta", "se", "p_value")])

# ── 7. Mutual adjustment — independent effects when all axes included ─────────
m_mutual <- lm(logNfL ~ burden_met_dam + burden_inf_dam + burden_cc_dam +
                 AGE_AT_SYMPTOM_ONSET + SEX,
               data = nfl_df)

cat("\n─── MUTUALLY ADJUSTED (all axes simultaneously) ──────────────\n")
mut_coefs <- summary(m_mutual)$coefficients
print(round(mut_coefs[c("burden_met_dam", "burden_inf_dam", "burden_cc_dam"), ], 4))

# ── 8. Subtype interaction test ───────────────────────────────────────────────
# Is the signal pan-sporadic (consistent) or concentrated in one subtype?

cat("\n─── SUBTYPE INTERACTION TESTS ────────────────────────────────\n")
for (burden_col in c("burden_met_dam", "burden_inf_dam", "burden_cc_dam")) {
  
  axis_label <- gsub("burden_(.*)_dam", "AXIS_\\1", 
                     toupper(burden_col)) |>
    gsub("AXIS_MET", "AXIS_METABOLIC", x = _) |>
    gsub("AXIS_INF", "AXIS_INFLAMMATION", x = _) |>
    gsub("AXIS_CC",  "AXIS_CELL_CYCLE", x = _)
  
  f_main <- as.formula(paste("logNfL ~", burden_col, 
                             "+ cluster_final + AGE_AT_SYMPTOM_ONSET + SEX"))
  f_int  <- as.formula(paste("logNfL ~", burden_col, 
                             "* cluster_final + AGE_AT_SYMPTOM_ONSET + SEX"))
  
  m_main <- lm(f_main, data = nfl_df)
  m_int  <- lm(f_int,  data = nfl_df)
  int_p  <- anova(m_main, m_int)[["Pr(>F)"]][2]
  
  cat(sprintf("%s × subtype interaction: p = %.4f  %s\n",
              axis_label, int_p,
              ifelse(int_p < 0.05, "← SIGNIFICANT", "(pan-sporadic)")))
}

# ── 9. Figures ────────────────────────────────────────────────────────────────

# 9a. Forest plot of per-axis effect sizes
forest_df <- axis_results |>
  mutate(
    ci_lo      = beta - 1.96 * se,
    ci_hi      = beta + 1.96 * se,
    sig        = ifelse(p_value < 0.05, "p < 0.05", "n.s."),
    axis_label = factor(axis,
                        levels = c("AXIS_METABOLIC", 
                                   "AXIS_INFLAMMATION",
                                   "AXIS_CELL_CYCLE"),
                        labels = c("Metabolic\n(AXIS_METABOLIC)",
                                   "Inflammatory\n(AXIS_INFLAMMATION)",
                                   "Cell Cycle\n(AXIS_CELL_CYCLE)"))
  )

# Replace scale_y_continuous with this — x is discrete, don't constrain it
p_forest <- ggplot(forest_df, 
                   aes(x = axis_label, y = beta, colour = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), 
                width = 0.15, linewidth = 1.1) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("b=%.4f\np=%.3f\n(%d genes)", 
                                beta, p_value, n_genes)),
            vjust = -0.8, size = 3, lineheight = 0.9) +
  scale_colour_manual(values = c("p < 0.05" = "#E74C3C", 
                                 "n.s."     = "#95A5A6")) +
  scale_x_discrete() +                          # <- fix: explicitly discrete
  scale_y_continuous(limits = c(-0.008, 0.004)) + # <- fix: correct y range
  labs(
    title    = "Per-axis rare damaging burden -> log NfL",
    subtitle = sprintf("Adjusted for age at onset and sex  |  n = %d with NfL", 
                       nrow(nfl_df)),
    x        = NULL,
    y        = "Beta coefficient (95% CI)",
    colour   = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        plot.title      = element_text(face = "bold"))

print(p_forest)
ggsave("axis_burden_nfl_forest.pdf", p_forest, width = 7, height = 5)

# Fix scatter y-axis label encoding issue
make_scatter <- function(burden_col, axis_title, df) {
  ggplot(df, aes(x = .data[[burden_col]], y = logNfL)) +
    geom_point(aes(colour = cluster_final), alpha = 0.45, size = 1.8) +
    geom_smooth(method = "lm", colour = "black", se = TRUE, linewidth = 1) +
    scale_colour_manual(values = subtype_cols) +
    labs(title  = axis_title,
         x      = "Rare damaging burden (gene count)",
         y      = "log NfL",          # plain ASCII avoids encoding warning
         colour = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 11))
}
ggsave("axis_burden_nfl_forest.pdf", p_forest, width = 7, height = 5)

# 9b. Scatter: per-axis burden vs NfL coloured by subtype
subtype_cols <- c("S-ALS-1" = "#E74C3C", 
                  "S-ALS-2" = "#3498DB", 
                  "S-ALS-3" = "#2ECC71")

make_scatter <- function(burden_col, axis_title, df) {
  ggplot(df, aes_string(x = burden_col, y = "logNfL")) +
    geom_point(aes(colour = cluster_final), alpha = 0.45, size = 1.8) +
    geom_smooth(method = "lm", colour = "black", se = TRUE, linewidth = 1) +
    scale_colour_manual(values = subtype_cols) +
    labs(title  = axis_title,
         x      = "Rare damaging burden (gene count)",
         y      = "log₁₀ NfL",
         colour = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 11))
}

p_scatter <- 
  make_scatter("burden_met_dam", "AXIS_METABOLIC",    nfl_df) +
  make_scatter("burden_inf_dam", "AXIS_INFLAMMATION", nfl_df) +
  make_scatter("burden_cc_dam",  "AXIS_CELL_CYCLE",   nfl_df) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

p_scatter <- p_scatter + 
  plot_annotation(
    title    = "Per-axis rare damaging burden vs plasma NfL",
    subtitle = "Black line: pan-sporadic linear fit  |  Coloured by transcriptomic subtype"
  )

print(p_scatter)
ggsave("axis_burden_nfl_scatter.pdf", p_scatter, width = 13, height = 5)

# ── 10. Paper-ready summary ───────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PAPER-READY SUMMARY\n")
cat("══════════════════════════════════════════════════════════════\n")
cat(sprintf("Pan-sporadic replication — rare LoF → NfL: β = %.3f, p = %.4f\n",
            coef(m_lof_total)["rare_lof"],
            summary(m_lof_total)$coefficients["rare_lof", "Pr(>|t|)"]))
cat(sprintf("Pan-sporadic rare damaging → NfL:          β = %.3f, p = %.4f\n",
            coef(m_dam_total)["rare_damaging"],
            summary(m_dam_total)$coefficients["rare_damaging", "Pr(>|t|)"]))
cat("\nPer-axis decomposition (rare damaging):\n")
print(axis_results[, c("axis", "n_genes", "beta", "se", "p_value")])

# Save
write.csv(axis_results, "axis_decomposed_lof_nfl.csv", row.names = FALSE)
cat("\nSaved: axis_decomposed_lof_nfl.csv\n")

