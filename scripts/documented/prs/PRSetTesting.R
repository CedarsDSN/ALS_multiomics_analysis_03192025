# ============================================================
# Script: PRSetTesting.R
# Project: ALS multi-omics analysis
# Source status: documented/sanitized for GitHub handoff
# ============================================================
# Documentation summary from scripts/README.md
# 
# ## `PRSetTesting.R`
# 
# **Category:** PRS / PRSet results
# 
# **Purpose**
# 
# Tests whether pathway-specific PRS values predict matching transcriptomic axis scores. Reads PRSet `.best` files, merges PRS with sporadic axis/subtype metadata, fits all 3 x 3 PRS-to-axis regressions, and visualizes specificity.
# 
# **Inputs**
# 
# - `prset_axis_metabolic.best`.
# - `prset_axis_cell_cycle.best`.
# - `prset_axis_inflammation.best`.
# - `vcf_sample_ids.txt`.
# - Workspace objects `sporadic_final_check` and `ax_clean_df`.
# 
# **Outputs**
# 
# - `prset_specificity_matrix.pdf` and `.png`.
# - `prs_cell_cycle_vs_transcriptomic_cell_cycle.pdf`.
# - `prset_specificity_matrix_slides.png` and `.pdf`.
# - `prset_diagonal_vs_offdiagonal_barplot.png`.
# - `prset_pairwise_r2_barplot.png`.
# - `prset_diagonal_regression_panel.png` and `.pdf`.
# - `prset_all_9_regressions.png`.
# - `prset_summary_table.csv`.
# - Recommended new output: `analysis_df.rds` for sensitivity scripts.
# 
# **Upstream dependencies**
# 
# - Run after external PRSet/PRSice finishes using files from `CreateFilesForPRSet.R`.
# - Requires sporadic subtype/axis metadata from `sporadic_subtyping.R`.
# 
# **How to run**
# 
# Legacy documented version:
# 
# ```bash
# Rscript scripts/documented/prs/PRSetTesting.R
# ```
# 
# Refactored command-line version:
# 
# ```bash
# Rscript scripts/refactored/prs/02_prset_testing.R --config config/config.yml --base-dir /path/to/private/ALS --outdir analysis/batch_06_genomics_prs_eqtl/results/prset_testing --objects analysis/batch_03_axis_scoring_sporadic_subtyping/objects/03_axis_scoring_objects.rds,analysis/batch_03_axis_scoring_sporadic_subtyping/objects/01_sporadic_subtyping_objects.rds
# ```
# 
# **Refactor notes**
# 
# - Add covariate-adjusted models if required by analysis plan.
# - Save `r2_mat`, `p_mat`, `b_mat`, and `analysis_df` as RDS/CSV.
# - Separate publication figures from exploratory slide figures.
# ============================================================

library(dplyr)

mainDir <- file.path(ALSPROJ_BASE_DIR, "prset")
setwd(mainDir)

# Each .best file's own-axis PRS is at the best threshold for that phenotype
# Extract only the relevant PRS column from each file

prs_met <- read.table(file.path("prset_axis_metabolic.best"),
                      header=TRUE, stringsAsFactors=FALSE) %>%
  filter(In_Regression == "Yes") %>%
  select(IID, PRS_METABOLIC = AXIS_METABOLIC)

prs_cc <- read.table(file.path("prset_axis_cell_cycle.best"),
                     header=TRUE, stringsAsFactors=FALSE) %>%
  filter(In_Regression == "Yes") %>%
  select(IID, PRS_CELL_CYCLE = AXIS_CELL_CYCLE)

prs_inf <- read.table(file.path("prset_axis_inflammation.best"),
                      header=TRUE, stringsAsFactors=FALSE) %>%
  filter(In_Regression == "Yes") %>%
  select(IID, PRS_INFLAMMATION = AXIS_INFLAMMATION)

cat("Metabolic PRS samples:", nrow(prs_met), "\n")
cat("Cell cycle PRS samples:", nrow(prs_cc), "\n")
cat("Inflammation PRS samples:", nrow(prs_inf), "\n")

# Merge all three PRS into one dataframe
prs_all <- prs_met %>%
  inner_join(prs_cc,  by = "IID") %>%
  inner_join(prs_inf, by = "IID")

cat("Merged PRS samples:", nrow(prs_all), "\n")


# Rebuild vcf_df if needed
vcf_ids_raw <- read.table(
  file.path(ALSPROJ_BASE_DIR, "vcf_sample_ids.txt"),
  header = FALSE, stringsAsFactors = FALSE
)$V1

vcf_df <- data.frame(
  VCF_ID    = vcf_ids_raw,
  CGND_bare = sub("-b38$", "", vcf_ids_raw),
  stringsAsFactors = FALSE
)

# Deduplicate sporadic_final_check to one row per Participant_ID
sporadic_meta_dedup <- sporadic_final_check %>%
  select(Participant_ID, NYGC_CGND_ID, AGE_AT_SYMPTOM_ONSET, SEX) %>%
  distinct(Participant_ID, .keep_all = TRUE)


# ax_clean_df must be in memory — check first
cat("ax_clean_df rows:", nrow(ax_clean_df), "\n")
cat("Columns:", names(ax_clean_df), "\n")

# Extract timepoint number from sample_id suffix (e.g. CASE-NEUAW157NMJ-5550-T → 5550)
ax_clean_df$timepoint <- as.integer(sub(".*-(\\d+)-T$", "\\1", ax_clean_df$sample_id))

# Keep earliest timepoint per patient
library(dplyr)
ax_clean_dedup <- ax_clean_df %>%
  group_by(Participant_ID) %>%
  slice_min(timepoint, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Rows before dedup:", nrow(ax_clean_df), "\n")   # should be 400
cat("Rows after dedup:", nrow(ax_clean_dedup), "\n")  # should be 391
cat("Unique patients:", length(unique(ax_clean_dedup$Participant_ID)), "\n")  # should be 391
cat("Duplicate Participant_IDs:", sum(duplicated(ax_clean_dedup$Participant_ID)), "\n")  # should be 0

# Build prs_meta_clean from ax_clean_dedup (earliest timepoint, 391 patients)
prs_meta_clean <- ax_clean_dedup %>%
  left_join(sporadic_meta_dedup, by = "Participant_ID") %>%
  mutate(CGND_bare = sub("-b38$", "", NYGC_CGND_ID)) %>%
  left_join(vcf_df, by = "CGND_bare") %>%
  filter(!is.na(VCF_ID))

cat("prs_meta_clean rows:", nrow(prs_meta_clean), "\n")        # should be 360
cat("Duplicate VCF_IDs:", sum(duplicated(prs_meta_clean$VCF_ID)), "\n")  # should be 0
# Join to ssGSEA axis scores
analysis_df <- prs_meta_clean %>%
  select(VCF_ID, AXIS_METABOLIC, AXIS_CELL_CYCLE, AXIS_INFLAMMATION, cluster_final) %>%
  inner_join(prs_all, by = c("VCF_ID" = "IID"))

cat("Analysis dataframe:", nrow(analysis_df), "samples\n")

# ── 3x3 specificity matrix ────────────────────────────────────────────────
prs_vars  <- c("PRS_METABOLIC", "PRS_CELL_CYCLE", "PRS_INFLAMMATION")
axis_vars <- c("AXIS_METABOLIC", "AXIS_CELL_CYCLE", "AXIS_INFLAMMATION")

r2_mat <- matrix(NA, 3, 3, dimnames = list(prs_vars, axis_vars))
p_mat  <- matrix(NA, 3, 3, dimnames = list(prs_vars, axis_vars))
b_mat  <- matrix(NA, 3, 3, dimnames = list(prs_vars, axis_vars))

for (prs in prs_vars) {
  for (axis in axis_vars) {
    fit <- lm(as.formula(paste(axis, "~", prs)), data = analysis_df)
    s <- summary(fit)
    r2_mat[prs, axis] <- s$r.squared
    p_mat[prs, axis]  <- coef(s)[prs, "Pr(>|t|)"]
    b_mat[prs, axis]  <- coef(s)[prs, "Estimate"]
  }
}

cat("\n=== R² matrix ===\n");       print(round(r2_mat, 5))
cat("\n=== P-value matrix ===\n");  print(round(p_mat, 4))
cat("\n=== Beta matrix ===\n");     print(round(b_mat, 3))

cat("\nDiagonal (on-target) results:\n")
for (i in 1:3) {
  cat(sprintf("  %s -> %s: R²=%.4f, p=%.4f\n",
              prs_vars[i], axis_vars[i], r2_mat[i,i], p_mat[i,i]))
}
cat(sprintf("\nMean off-diagonal R²: %.5f\n",
            mean(r2_mat[row(r2_mat) != col(r2_mat)])))


library(ggplot2)
library(reshape2)

# Build long-format dataframe for plotting
r2_df <- melt(r2_mat, varnames = c("PRS", "ssGSEA_Axis"), value.name = "R2")
p_df  <- melt(p_mat,  varnames = c("PRS", "ssGSEA_Axis"), value.name = "P")
plot_df <- left_join(r2_df, p_df, by = c("PRS", "ssGSEA_Axis")) %>%
  mutate(
    sig_label = case_when(
      P < 0.001 ~ "***",
      P < 0.01  ~ "**",
      P < 0.05  ~ "*",
      P < 0.1   ~ "†",
      TRUE      ~ ""
    ),
    is_diagonal = (PRS == "PRS_METABOLIC"    & ssGSEA_Axis == "AXIS_METABOLIC") |
      (PRS == "PRS_CELL_CYCLE"   & ssGSEA_Axis == "AXIS_CELL_CYCLE") |
      (PRS == "PRS_INFLAMMATION" & ssGSEA_Axis == "AXIS_INFLAMMATION"),
    # Clean labels
    PRS = recode(PRS,
                 PRS_METABOLIC    = "Metabolic PRS",
                 PRS_CELL_CYCLE   = "Cell cycle PRS",
                 PRS_INFLAMMATION = "Inflammation PRS"),
    ssGSEA_Axis = recode(ssGSEA_Axis,
                         AXIS_METABOLIC    = "Metabolic axis",
                         AXIS_CELL_CYCLE   = "Cell cycle axis",
                         AXIS_INFLAMMATION = "Inflammation axis"),
    # Order factors
    PRS = factor(PRS, levels = c("Metabolic PRS", "Cell cycle PRS", "Inflammation PRS")),
    ssGSEA_Axis = factor(ssGSEA_Axis, levels = c("Metabolic axis", "Cell cycle axis", "Inflammation axis"))
  )

# Plot
p <- ggplot(plot_df, aes(x = ssGSEA_Axis, y = PRS, fill = R2)) +
  geom_tile(aes(color = is_diagonal), linewidth = 1.2) +
  scale_color_manual(values = c("TRUE" = "#1D9E75", "FALSE" = "transparent"), guide = "none") +
  geom_text(aes(label = paste0("R²=", round(R2, 4), "\n", sig_label)),
            size = 3.5, lineheight = 1.2) +
  scale_fill_gradient(low = "white", high = "#3B8BD4",
                      name = "R²", limits = c(0, 0.02)) +
  labs(
    title = "PRSet cross-pathway specificity matrix",
    subtitle = "ALS GWAS pathway PRS (rows) vs transcriptomic axis ssGSEA scores (cols)",
    x = "Transcriptomic axis (ssGSEA score)",
    y = "Pathway-specific PRS",
    caption = "† p<0.1  * p<0.05  ** p<0.01 | van Rheenen 2021 GWAS · n=360 sporadic ALS · unadjusted"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )

print(p)

ggsave("prset_specificity_matrix.pdf", p, width = 7, height = 5)
ggsave("prset_specificity_matrix.png", p, width = 7, height = 5, dpi = 300)
cat("Saved\n")


# Does cell cycle PRS differ across subtypes?
summary(aov(PRS_CELL_CYCLE ~ cluster_final, data = analysis_df))
TukeyHSD(aov(PRS_CELL_CYCLE ~ cluster_final, data = analysis_df))

pdf("prs_cell_cycle_vs_transcriptomic_cell_cycle.pdf", width = 6, height = 5)
# Also useful — scatter plot of PRS_CELL_CYCLE vs AXIS_CELL_CYCLE coloured by subtype
ggplot(analysis_df, aes(x = PRS_CELL_CYCLE, y = AXIS_CELL_CYCLE, color = cluster_final)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_manual(values = c("S-ALS-1" = "#D85A30", "S-ALS-2" = "#888780", "S-ALS-3" = "#BA7517")) +
  labs(x = "Cell cycle pathway PRS", y = "Cell cycle ssGSEA score",
       title = "Cell cycle PRS vs transcriptomic cell cycle activation",
       color = "Subtype") +
  theme_minimal(base_size = 12)
dev.off()



library(dplyr)
library(ggplot2)
library(reshape2)

plot_df <- melt(r2_mat, varnames = c("PRS", "Axis"), value.name = "R2") %>%
  left_join(
    melt(p_mat, varnames = c("PRS", "Axis"), value.name = "P"),
    by = c("PRS", "Axis")
  ) %>%
  left_join(
    melt(b_mat, varnames = c("PRS", "Axis"), value.name = "Beta"),
    by = c("PRS", "Axis")
  ) %>%
  mutate(
    sig_class = case_when(
      P < 0.05 ~ "p < 0.05",
      P < 0.10 ~ "p < 0.10",
      TRUE ~ "NS"
    ),
    is_diagonal = (PRS == "PRS_METABOLIC"    & Axis == "AXIS_METABOLIC") |
      (PRS == "PRS_CELL_CYCLE"   & Axis == "AXIS_CELL_CYCLE") |
      (PRS == "PRS_INFLAMMATION" & Axis == "AXIS_INFLAMMATION"),
    PRS = recode(PRS,
                 PRS_METABOLIC    = "Metabolic PRS",
                 PRS_CELL_CYCLE   = "Cell cycle PRS",
                 PRS_INFLAMMATION = "Inflammation PRS"
    ),
    Axis = recode(Axis,
                  AXIS_METABOLIC    = "Metabolic axis",
                  AXIS_CELL_CYCLE   = "Cell cycle axis",
                  AXIS_INFLAMMATION = "Inflammation axis"
    ),
    PRS = factor(PRS, levels = c("Metabolic PRS", "Cell cycle PRS", "Inflammation PRS")),
    Axis = factor(Axis, levels = c("Metabolic axis", "Cell cycle axis", "Inflammation axis")),
    label = paste0("R² = ", sprintf("%.3f", R2), "\nβ = ", sprintf("%.3f", Beta))
  )

p_heatmap2 <- ggplot(plot_df, aes(x = Axis, y = PRS, fill = R2)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_point(
    data = subset(plot_df, sig_class == "p < 0.05"),
    shape = 21, size = 4, stroke = 1, fill = "white", color = "black"
  ) +
  geom_point(
    data = subset(plot_df, sig_class == "p < 0.10"),
    shape = 4, size = 3, stroke = 1, color = "black"
  ) +
  geom_text(aes(label = label), size = 3.2, lineheight = 1.1) +
  scale_fill_gradient(low = "white", high = "#3B8BD4", limits = c(0, max(plot_df$R2, na.rm = TRUE))) +
  labs(
    title = "PRSet cross-pathway specificity matrix",
    subtitle = "White circle: p < 0.05   |   Cross: p < 0.10",
    x = "Transcriptomic axis",
    y = "Pathway-specific PRS",
    fill = expression(R^2)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)
  )

ggsave("prset_specificity_matrix_slides.png", p_heatmap2, width = 7.2, height = 5.2, dpi = 300)
ggsave("prset_specificity_matrix_slides.pdf", p_heatmap2, width = 7.2, height = 5.2)
p_heatmap2


diag_summary <- plot_df %>%
  mutate(pair_type = if_else(is_diagonal, "Diagonal", "Off-diagonal")) %>%
  group_by(pair_type) %>%
  summarise(
    mean_R2 = mean(R2, na.rm = TRUE),
    median_R2 = median(R2, na.rm = TRUE),
    n_p05 = sum(P < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

p_diag_bar <- ggplot(diag_summary, aes(x = pair_type, y = mean_R2)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("Mean R² = %.3f", mean_R2)), vjust = -0.5, size = 4) +
  labs(
    title = "On-target PRS–axis pairs show stronger signal",
    x = NULL,
    y = "Mean R²"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("prset_diagonal_vs_offdiagonal_barplot.png", p_diag_bar, width = 6, height = 4.5, dpi = 300)
p_diag_bar

pair_bar_df <- plot_df %>%
  mutate(
    pair = paste(PRS, "→", Axis),
    pair = factor(pair, levels = pair)
  )

p_pair_bar <- ggplot(pair_bar_df, aes(x = reorder(pair, R2), y = R2, fill = is_diagonal)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("p=%.3f", P)), hjust = -0.1, size = 3.3) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#1D9E75", "FALSE" = "#BFC7D5"), guide = "none") +
  labs(
    title = "Pairwise PRS–axis association strength",
    x = NULL,
    y = expression(R^2)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("prset_pairwise_r2_barplot.png", p_pair_bar, width = 7.5, height = 5.5, dpi = 300)
p_pair_bar

library(patchwork)

plot_diag_reg <- function(df, prs_col, axis_col, color_var = "cluster_final", title_text = NULL, xlab = NULL, ylab = NULL) {
  ggplot(df, aes_string(x = prs_col, y = axis_col, color = color_var)) +
    geom_point(alpha = 0.55, size = 1.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    scale_color_manual(values = c(
      "S-ALS-1" = "#D85A30",
      "S-ALS-2" = "#888780",
      "S-ALS-3" = "#BA7517"
    )) +
    labs(
      title = title_text,
      x = xlab,
      y = ylab,
      color = "Subtype"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank()
    )
}

p1 <- plot_diag_reg(
  analysis_df,
  "PRS_METABOLIC", "AXIS_METABOLIC",
  title_text = "Metabolic PRS vs metabolic axis",
  xlab = "Metabolic pathway PRS",
  ylab = "Metabolic ssGSEA score"
)

p2 <- plot_diag_reg(
  analysis_df,
  "PRS_CELL_CYCLE", "AXIS_CELL_CYCLE",
  title_text = "Cell cycle PRS vs cell cycle axis",
  xlab = "Cell cycle pathway PRS",
  ylab = "Cell cycle ssGSEA score"
)

p3 <- plot_diag_reg(
  analysis_df,
  "PRS_INFLAMMATION", "AXIS_INFLAMMATION",
  title_text = "Inflammation PRS vs inflammation axis",
  xlab = "Inflammation pathway PRS",
  ylab = "Inflammation ssGSEA score"
)

p_reg_panel <- (p1 | p2 | p3) +
  plot_annotation(title = "Diagonal PRS–axis associations")

ggsave("prset_diagonal_regression_panel.png", p_reg_panel, width = 13, height = 4.5, dpi = 300)
ggsave("prset_diagonal_regression_panel.pdf", p_reg_panel, width = 13, height = 4.5)
p_reg_panel

all_pairs <- expand.grid(
  prs = c("PRS_METABOLIC", "PRS_CELL_CYCLE", "PRS_INFLAMMATION"),
  axis = c("AXIS_METABOLIC", "AXIS_CELL_CYCLE", "AXIS_INFLAMMATION"),
  stringsAsFactors = FALSE
)

all_pair_plots <- lapply(seq_len(nrow(all_pairs)), function(i) {
  prs_col <- all_pairs$prs[i]
  axis_col <- all_pairs$axis[i]
  
  ggplot(analysis_df, aes_string(x = prs_col, y = axis_col, color = "cluster_final")) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
    scale_color_manual(values = c(
      "S-ALS-1" = "#D85A30",
      "S-ALS-2" = "#888780",
      "S-ALS-3" = "#BA7517"
    )) +
    labs(
      title = paste(gsub("PRS_", "", prs_col), "→", gsub("AXIS_", "", axis_col)),
      x = prs_col,
      y = axis_col
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    )
})

p_all9 <- wrap_plots(all_pair_plots, ncol = 3)
ggsave("prset_all_9_regressions.png", p_all9, width = 12, height = 10, dpi = 300)
p_all9


summary_table <- plot_df %>%
  mutate(
    pair_type = if_else(is_diagonal, "Diagonal", "Off-diagonal"),
    R2 = round(R2, 4),
    Beta = round(Beta, 3),
    P = signif(P, 3)
  ) %>%
  arrange(desc(is_diagonal), desc(R2)) %>%
  select(PRS, Axis, pair_type, Beta, R2, P)

print(summary_table)

write.csv(summary_table, "prset_summary_table.csv", row.names = FALSE)
