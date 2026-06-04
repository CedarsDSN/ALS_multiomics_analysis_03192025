# ============================================================
# Script: SensitivityAnalysis.R
# Project: ALS multi-omics analysis
# Source status: documented/sanitized for GitHub handoff
# ============================================================
# Documentation summary from scripts/README.md
# 
# ## `SensitivityAnalysis.R`
# 
# **Category:** Subtype validation / sensitivity
# 
# **Purpose**
# 
# Tests whether transcriptomic subtype and axis findings are artifacts of differentiation quality, glial contamination, proliferation, marker outliers, or related confounders. Also examines CCT2, PRS, NfL, AXIS_METABOLIC, and ALSFRS-R relationships.
# 
# **Inputs**
# 
# - Required base objects: `expr_clean_mat` and `ax_clean_dedup`.
# - Later PRS/CCT2 sections require `analysis_df`, `expr_sym` or `expr_clean_mat`, and subtype/clinical metadata with `cluster_final`, axis scores, sex, NfL, ALSFRS-R slope, and progression variables.
# 
# **Outputs**
# 
# - `diff_quality_by_subtype.pdf`.
# - `diff_quality_sensitivity_results.csv`.
# - `sensitivity_summary.csv`.
# - Additional scatter/forest/model outputs may be produced by later sections depending on available objects.
# 
# **Upstream dependencies**
# 
# - Must run after `AxisScoring.R` and `sporadic_subtyping.R`.
# - PRS-specific checks require `PRSetTesting.R` to have produced `analysis_df`.
# 
# **How to run**
# 
# Legacy documented version:
# 
# ```bash
# Rscript scripts/documented/validation/SensitivityAnalysis.R
# ```
# 
# Refactored command-line version:
# 
# ```bash
# Rscript scripts/refactored/validation/01_sensitivity_analysis.R --config config/config.yml --base-dir /path/to/private/ALS --outdir analysis/batch_04_validation_classifier_ipa/results/sensitivity_analysis --objects analysis/batch_03_axis_scoring_sporadic_subtyping/objects/03_axis_scoring_objects.rds,analysis/batch_03_axis_scoring_sporadic_subtyping/objects/01_sporadic_subtyping_objects.rds
# ```
# 
# **Refactor notes**
# 
# - Split into `01_differentiation_quality.R`, `02_glial_residualization.R`, `03_cct2_checks.R`, `04_prs_cct2_nfl.R`, and `05_axis_clinical_models.R`.
# - Save each model summary as a tidy CSV instead of only printing to console.
# - Define a single metadata contract for NfL, ALSFRS slope, age, sex, and subtype columns.
# ============================================================

# ============================================================
# DIFFERENTIATION QUALITY SENSITIVITY ANALYSIS
# Purpose: confirm transcriptomic subtypes are not artifacts of
# differentiation efficiency or glial contamination
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(mclust)              # for adjustedRandIndex
  library(ConsensusClusterPlus)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ------------------------------------------------------------
# SETUP — confirm inputs and build a symbol-keyed expression mat
# ------------------------------------------------------------

# Sanity checks on the objects we rely on
stopifnot(exists("expr_clean_mat"))
stopifnot(exists("ax_clean_dedup"))

cat("expr_clean_mat dims:", dim(expr_clean_mat), "\n")
cat("ax_clean_dedup rows:", nrow(ax_clean_dedup), "\n")
cat("First 3 expr rownames:", head(rownames(expr_clean_mat), 3), "\n")
cat("First 3 expr colnames:", head(colnames(expr_clean_mat), 3), "\n")

# Convert Ensembl -> Symbol if rownames look like ENSG...
# If your rownames are ALREADY gene symbols, skip this block and
# set expr_sym <- expr_clean_mat
if (any(grepl("^ENSG", rownames(expr_clean_mat)))) {
  ens_clean <- sub("\\..*$", "", rownames(expr_clean_mat))  # strip version suffix
  sym_map   <- mapIds(org.Hs.eg.db,
                      keys      = ens_clean,
                      column    = "SYMBOL",
                      keytype   = "ENSEMBL",
                      multiVals = "first")
  
  keep_idx  <- !is.na(sym_map) & !duplicated(sym_map)
  expr_sym  <- expr_clean_mat[keep_idx, ]
  rownames(expr_sym) <- sym_map[keep_idx]
  cat("Mapped to symbols:", nrow(expr_sym), "genes\n")
} else {
  expr_sym <- expr_clean_mat
  cat("Rownames already appear to be symbols; skipping conversion\n")
}

# Build a sample -> Participant_ID -> subtype table
# This assumes colnames(expr_sym) are sample_ids and ax_clean_dedup
# has a matching key. Adjust the join key if your setup differs.
if ("sample_id" %in% colnames(ax_clean_dedup)) {
  sample_key <- "sample_id"
} else if ("Sample_ID" %in% colnames(ax_clean_dedup)) {
  sample_key <- "Sample_ID"
} else {
  stop("Need a sample_id column in ax_clean_dedup to match expr colnames")
}

meta <- ax_clean_dedup %>%
  filter(.data[[sample_key]] %in% colnames(expr_sym)) %>%
  arrange(match(.data[[sample_key]], colnames(expr_sym)))

expr_sym <- expr_sym[, meta[[sample_key]]]
stopifnot(all(meta[[sample_key]] == colnames(expr_sym)))
cat("Aligned matrix:", nrow(expr_sym), "genes x", ncol(expr_sym), "samples\n")
cat("Subtype counts:\n"); print(table(meta$cluster_final))


# ============================================================
# ANALYSIS 1 — DIFFERENTIATION QUALITY SCORE BY SUBTYPE
# ============================================================

# Maturity markers: mature neuronal + spinal MN identity
maturity_markers   <- c("MAP2","SYN1","NEFL","NEFM","SNAP25","RBFOX3",
                        "STMN2","CHAT","ISL1","MNX1","SLC18A3",
                        "TUBB3","ELAVL3","ELAVL4")

# Progenitor / proliferation markers: if high, cells are still dividing
progenitor_markers <- c("PAX6","SOX2","NES","MKI67","TOP2A","PCNA",
                        "VIM","SOX1","MSI1")

present <- function(g) intersect(g, rownames(expr_sym))
mat_g <- present(maturity_markers)
pro_g <- present(progenitor_markers)

cat("Maturity markers found:", length(mat_g), "/", length(maturity_markers), "\n")
cat("Progenitor markers found:", length(pro_g), "/", length(progenitor_markers), "\n")
cat("Missing maturity:", setdiff(maturity_markers, mat_g), "\n")
cat("Missing progenitor:", setdiff(progenitor_markers, pro_g), "\n")

# Z-score each gene across samples before averaging — prevents any single
# highly-expressed gene from dominating the score
zrow <- function(x) (x - mean(x)) / sd(x)
expr_z <- t(apply(expr_sym, 1, zrow))

maturity_score   <- colMeans(expr_z[mat_g, , drop = FALSE])
progenitor_score <- colMeans(expr_z[pro_g, , drop = FALSE])
diff_quality     <- maturity_score - progenitor_score

meta$maturity_score   <- maturity_score
meta$progenitor_score <- progenitor_score
meta$diff_quality     <- diff_quality

# ANOVA: does differentiation quality vary by subtype?
cat("\n=== ANOVA: differentiation quality ~ subtype ===\n")
aov1 <- aov(diff_quality ~ cluster_final, data = meta)
print(summary(aov1))

cat("\n=== Tukey HSD ===\n")
print(TukeyHSD(aov1))

cat("\n=== Group means ± SD ===\n")
print(meta %>%
        group_by(cluster_final) %>%
        summarise(n        = n(),
                  mat_mean = mean(maturity_score),
                  pro_mean = mean(progenitor_score),
                  dq_mean  = mean(diff_quality),
                  dq_sd    = sd(diff_quality)))

# Effect size — eta^2
ss     <- summary(aov1)[[1]]$`Sum Sq`
eta_sq <- ss[1] / sum(ss)
cat("\nEta-squared (subtype on diff quality):", round(eta_sq, 3), "\n")
cat("Interpretation: <0.01 negligible, 0.06 medium, 0.14 large\n")

# Plot — save for supplementary figure
p1 <- ggplot(meta, aes(cluster_final, diff_quality, fill = cluster_final)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  labs(x = NULL, y = "Maturity − Progenitor (z-score mean)",
       title = "Differentiation quality by subtype") +
  theme_bw() + theme(legend.position = "none")
ggsave("diff_quality_by_subtype.pdf", p1, width = 5, height = 4)


# ============================================================
# ANALYSIS 2 — BROADER GLIAL RESIDUALIZATION
# ============================================================

astrocyte_markers    <- c("GFAP","AQP4","ALDH1L1","S100B","SLC1A3",
                          "SLC1A2","GJA1","AGT")
microglia_markers    <- c("AIF1","CX3CR1","P2RY12","TMEM119","CD68",
                          "ITGAM","TYROBP")
oligodendrocyte_mk   <- c("OLIG2","MBP","MOG","PLP1","SOX10","GALC")

glia_all <- c(astrocyte_markers, microglia_markers, oligodendrocyte_mk)
glia_g   <- present(glia_all)
cat("\nGlial markers found:", length(glia_g), "/", length(glia_all), "\n")
cat("Missing:", setdiff(glia_all, glia_g), "\n")

# Per-sample composite scores (separate, so we can report each coefficient)
meta$astro_score  <- colMeans(expr_z[present(astrocyte_markers), , drop = FALSE])
meta$micro_score  <- colMeans(expr_z[present(microglia_markers),  , drop = FALSE])
meta$oligo_score  <- colMeans(expr_z[present(oligodendrocyte_mk), , drop = FALSE])

# Residualize each axis on all three glial scores + differentiation quality
# (diff_quality in the regression makes this a joint test)
resid_axis <- function(y) {
  fit <- lm(y ~ astro_score + micro_score + oligo_score + diff_quality,
            data = meta)
  cat("  Adjusted R^2:", round(summary(fit)$adj.r.squared, 3), "\n")
  residuals(fit)
}

cat("\n=== Residualizing AXIS_METABOLIC ===\n")
meta$AXIS_METABOLIC_resid   <- resid_axis(meta$AXIS_METABOLIC)
cat("=== Residualizing AXIS_CELL_CYCLE ===\n")
meta$AXIS_CELL_CYCLE_resid  <- resid_axis(meta$AXIS_CELL_CYCLE)
cat("=== Residualizing AXIS_INFLAMMATION ===\n")
meta$AXIS_INFLAMMATION_resid <- resid_axis(meta$AXIS_INFLAMMATION)

# Re-cluster on residualized axes — same pipeline as original
axis_mat_resid <- meta %>%
  select(AXIS_METABOLIC_resid, AXIS_CELL_CYCLE_resid, AXIS_INFLAMMATION_resid) %>%
  as.matrix() %>% scale() %>% t()
colnames(axis_mat_resid) <- meta[[sample_key]]

dir.create("ccp_resid_broadglial", showWarnings = FALSE)
ccp_res <- ConsensusClusterPlus(axis_mat_resid,
                                maxK        = 6,
                                reps        = 1000,
                                pItem       = 0.8,
                                pFeature    = 1,
                                clusterAlg  = "km",
                                distance    = "euclidean",
                                seed        = 1234,
                                title       = "ccp_resid_broadglial",
                                plot        = "pdf")

meta$cluster_resid <- ccp_res[[3]]$consensusClass

# Align labels via majority vote so ARI isn't punished by label permutation
contingency <- table(meta$cluster_final, meta$cluster_resid)
cat("\n=== Contingency: original (rows) vs residualized (cols) ===\n")
print(contingency)

ari <- adjustedRandIndex(meta$cluster_final, meta$cluster_resid)
cat("\n>>> ARI (original vs broad-glial-residualized):", round(ari, 3), "\n")
cat("Prior S100b-only ARI was 0.791; target is similar or better\n")
cat("ARI > 0.70 = subtypes robust; 0.50-0.70 = partial; <0.50 = problem\n")


# ============================================================
# ANALYSIS 3 — CCT2 ORTHOGONALITY TO PROLIFERATION / MATURITY
# ============================================================

stopifnot("CCT2" %in% rownames(expr_sym))
meta$CCT2_expr <- expr_sym["CCT2", ]

prolif_markers   <- c("MKI67","PCNA","TOP2A","MCM2","CCNB1","CDK1","BIRC5")
prolif_g         <- present(prolif_markers)
cat("\nProliferation markers found:", length(prolif_g), "/", length(prolif_markers), "\n")

# Per-gene correlations with CCT2
cat("\n=== CCT2 correlations with proliferation markers ===\n")
prolif_cors <- sapply(prolif_g, function(g) {
  ct <- cor.test(meta$CCT2_expr, expr_sym[g, ], method = "pearson")
  c(r = ct$estimate, p = ct$p.value)
})
print(t(prolif_cors))

cat("\n=== CCT2 correlations with maturity markers ===\n")
mat_cors <- sapply(mat_g, function(g) {
  ct <- cor.test(meta$CCT2_expr, expr_sym[g, ], method = "pearson")
  c(r = ct$estimate, p = ct$p.value)
})
print(t(mat_cors))

# Composite test — CCT2 ~ proliferation composite + maturity composite + diff quality
meta$prolif_score <- colMeans(expr_z[prolif_g, , drop = FALSE])

cat("\n=== Joint regression: CCT2 ~ proliferation + maturity + diff quality ===\n")
cct2_fit <- lm(CCT2_expr ~ prolif_score + maturity_score + diff_quality,
               data = meta)
print(summary(cct2_fit))

cat("\n=== Does CCT2 still separate S-ALS-3 after adjusting for these? ===\n")
cct2_adj_fit <- lm(CCT2_expr ~ cluster_final + prolif_score + maturity_score +
                     diff_quality, data = meta)
print(summary(cct2_adj_fit))

# Compare unadjusted vs adjusted effect size for S-ALS-3
cct2_unadj <- lm(CCT2_expr ~ cluster_final, data = meta)
cat("\n=== Unadjusted S-ALS-3 coefficient (reference = S-ALS-1) ===\n")
print(summary(cct2_unadj)$coefficients)

cat("\n=== Adjusted S-ALS-3 coefficient ===\n")
print(summary(cct2_adj_fit)$coefficients)


# ============================================================
# SAVE CONSOLIDATED OUTPUTS
# ============================================================

write.csv(meta, "diff_quality_sensitivity_results.csv", row.names = FALSE)

sensitivity_summary <- data.frame(
  test = c("Diff quality ANOVA p",
           "Diff quality eta^2",
           "ARI (original vs broad-glial-resid)",
           "CCT2 ~ prolif_score p",
           "CCT2 ~ maturity_score p",
           "S-ALS-3 CCT2 effect (unadjusted)",
           "S-ALS-3 CCT2 effect (adjusted)"),
  value = c(summary(aov1)[[1]]$`Pr(>F)`[1],
            eta_sq,
            ari,
            summary(cct2_fit)$coefficients["prolif_score", 4],
            summary(cct2_fit)$coefficients["maturity_score", 4],
            coef(cct2_unadj)["cluster_finalS-ALS-3"],      # label may differ
            coef(cct2_adj_fit)["cluster_finalS-ALS-3"])    # label may differ
)
print(sensitivity_summary)
write.csv(sensitivity_summary, "sensitivity_summary.csv", row.names = FALSE)

cat("\n=== DONE — share console output for interpretation ===\n")


# S-ALS-3 vs S-ALS-2 CCT2 test, dropping the confounded S-ALS-1
meta_s23 <- meta %>% filter(cluster_final %in% c("S-ALS-2", "S-ALS-3")) %>%
  mutate(cluster_final = droplevels(cluster_final))

cat("=== CCT2 ~ subtype (S-ALS-2 vs S-ALS-3 only, unadjusted) ===\n")
print(summary(lm(CCT2_expr ~ cluster_final, data = meta_s23))$coefficients)

cat("\n=== CCT2 ~ subtype + prolif + maturity + diff_quality ===\n")
print(summary(lm(CCT2_expr ~ cluster_final + prolif_score + maturity_score + 
                   diff_quality, data = meta_s23))$coefficients)

cat("\n=== And with full glial panel ===\n")
print(summary(lm(CCT2_expr ~ cluster_final + prolif_score + maturity_score +
                   diff_quality + astro_score + micro_score + oligo_score, 
                 data = meta_s23))$coefficients)





# ============================================================
# PRS_CELL_CYCLE → CCT2 expression
# Closing the gene-level link in the genome → transcriptome chain
# ============================================================

# Confirm the objects are still in env
stopifnot(exists("analysis_df"))   # 360 patients with PRS + axis scores
stopifnot(exists("expr_clean_mat") | exists("expr_sym"))

# Build the regression dataframe
# We need: PRS_CELL_CYCLE, CCT2 expression, SEX, cluster_final, ancestry PCs if available
if (exists("expr_sym")) {
  cct2_vec <- expr_sym["CCT2", ]
} else {
  # If you only have ENSG IDs, look up CCT2's ENSG (ENSG00000166226) and pull
  cct2_ensg <- "ENSG00000166226"
  match_idx <- which(sub("\\..*$", "", rownames(expr_clean_mat)) == cct2_ensg)
  cct2_vec  <- expr_clean_mat[match_idx, ]
}

# Make a sample-keyed dataframe for joining
cct2_df <- data.frame(
  sample_id = names(cct2_vec),
  CCT2_expr = as.numeric(cct2_vec)
)

# Adjust the join key to match analysis_df — likely Participant_ID or sample_id
# Check first what keys exist
cat("analysis_df columns:\n"); print(colnames(analysis_df))
cat("\nFirst 3 cct2_df sample_ids:", head(cct2_df$sample_id, 3), "\n")
cat("First 3 analysis_df sample IDs:\n")
print(head(analysis_df[, intersect(c("sample_id","Sample_ID","Participant_ID","VCF_ID"), 
                                   colnames(analysis_df))], 3))

map_df <- sporadic_final_check %>%
  select(Participant_ID, NYGC_CGND_ID, SEX) %>%
  distinct() %>%
  group_by(NYGC_CGND_ID) %>%
  slice(1) %>%
  ungroup()
map_df <- unique(map_df)
analysis_df$NYGC_CGND_ID <- sub("-b38$", "", analysis_df$VCF_ID)

analysis_df2 <- merge(analysis_df, map_df, by = "NYGC_CGND_ID")
# Assuming Participant_ID is the join key (adjust if needed):
# If sample_id has the long format with -N-T suffix, strip it
cct2_df$Participant_ID <- sub("-[0-9]+-T$", "", cct2_df$sample_id)

# Join — should produce 360 rows (patients with PRS)
df <- merge(
  analysis_df2,
  cct2_df[, c("Participant_ID", "CCT2_expr")],
  by = "Participant_ID"
)
cat("\nMerged rows:", nrow(df), "(expecting 360)\n")
cat("CCT2 missingness:", sum(is.na(df$CCT2_expr)), "\n")

# ============================================================
# THE PRIMARY TEST
# ============================================================

cat("\n=== Test 1: PRS_CELL_CYCLE → CCT2 (univariate) ===\n")
fit1 <- lm(CCT2_expr ~ PRS_CELL_CYCLE, data = df)
print(summary(fit1)$coefficients)
cat("R^2:", round(summary(fit1)$r.squared, 4), "\n")

cat("\n=== Test 2: PRS_CELL_CYCLE → CCT2 + SEX ===\n")
fit2 <- lm(CCT2_expr ~ PRS_CELL_CYCLE + SEX, data = df)
print(summary(fit2)$coefficients)

cat("\n=== Test 3: PRS_CELL_CYCLE → CCT2 + SEX + subtype ===\n")
# Conditional test: does PRS predict CCT2 *beyond* what subtype assignment captures?
fit3 <- lm(CCT2_expr ~ PRS_CELL_CYCLE + SEX + cluster_final, data = df)
print(summary(fit3)$coefficients)

# ============================================================
# SPECIFICITY CHECK — does this work for the right PRS only?
# ============================================================

cat("\n=== Specificity: each PRS → CCT2 ===\n")
for (prs in c("PRS_METABOLIC","PRS_CELL_CYCLE","PRS_INFLAMMATION")) {
  cat("\n--", prs, "--\n")
  f <- lm(as.formula(paste("CCT2_expr ~", prs, "+ SEX")), data = df)
  print(summary(f)$coefficients[2, ])
}

# ============================================================
# ALTERNATIVE FRAMING — predict S-ALS-3 membership directly
# ============================================================

cat("\n=== Logistic: PRS_CELL_CYCLE → P(S-ALS-3) ===\n")
df$is_sals3 <- as.integer(df$cluster_final == "S-ALS-3")
glm1 <- glm(is_sals3 ~ PRS_CELL_CYCLE + SEX, data = df, family = "binomial")
print(summary(glm1)$coefficients)

# Restricted to S-ALS-2 vs S-ALS-3 (clean half) — matches Test 3 framing
cat("\n=== Logistic restricted to S-ALS-2 + S-ALS-3 ===\n")
df_clean <- df[df$cluster_final %in% c("S-ALS-2","S-ALS-3"), ]
df_clean$is_sals3 <- as.integer(df_clean$cluster_final == "S-ALS-3")
glm2 <- glm(is_sals3 ~ PRS_CELL_CYCLE + SEX, data = df_clean, family = "binomial")
print(summary(glm2)$coefficients)



# CCT2 -> NFl connection?
# ============================================================
# CCT2 → log NfL in the clean subset (S-ALS-2 + S-ALS-3)
# Closes the second link in the chain:
#   PRS_CELL_CYCLE → CCT2 → NfL → disease duration
# ============================================================

# --- Build the analysis frame ---
# Pull CCT2 expression from the VST matrix (genes × samples)
cct2_vec <- expr_clean_mat["CCT2", ]
# Sample_ID → Participant_ID (strip the -N-T suffix)
pid_from_sample <- sub("-[0-9]+-T$", "", names(cct2_vec))
cct2_df <- data.frame(
  Participant_ID = pid_from_sample,
  CCT2_expr      = as.numeric(cct2_vec),
  stringsAsFactors = FALSE
)

# Attach subtype labels and NfL
df_nfl_cct2 <- merge(cct2_df,
                     ax_clean_dedup[, c("Participant_ID", "cluster_final")],
                     by = "Participant_ID")
df_nfl_cct2 <- merge(df_nfl_cct2, df_nfl_uniq, by = "Participant_ID")
# (If AGE_AT_SYMPTOM_ONSET / SEX aren't on df_nfl_uniq, merge them in here from
#  whichever metadata frame holds them — clinical_df, meta, etc.)

cat("Total with CCT2 + NfL:", nrow(df_nfl_cct2), "\n")
print(table(df_nfl_cct2$cluster_final))

# Restrict to the clean half — S-ALS-1 excluded (differentiation confound)
df_clean <- subset(df_nfl_cct2, cluster_final %in% c("S-ALS-2", "S-ALS-3"))
df_clean$cluster_final <- droplevels(factor(df_clean$cluster_final))
cat("\nClean subset n:", nrow(df_clean), "\n")

# ----------------------------------------------------------------
# Test 1: Univariate
# Headline: does CCT2 expression alone predict log NfL?
# ----------------------------------------------------------------
cat("\n=== Test 1: log NfL ~ CCT2 (univariate) ===\n")
fit1 <- lm(logNfL ~ CCT2_expr, data = df_clean)
print(summary(fit1)$coefficients)
cat("R^2:", round(summary(fit1)$r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 2: Standard biomarker adjustment (age + sex)
# NfL rises with age and differs by sex — these always go in.
# ----------------------------------------------------------------
cat("\n=== Test 2: log NfL ~ CCT2 + age + sex ===\n")
fit2 <- lm(logNfL ~ CCT2_expr + AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
print(summary(fit2)$coefficients)
cat("Adj R^2:", round(summary(fit2)$adj.r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 3: + subtype — the load-bearing test
# Does CCT2 carry biomarker signal BEYOND the discrete subtype label?
# Expect collinearity (CCT2 alone classifies S-ALS-3 at AUC 0.918 in
# the clean half), so attenuation is expected. The question is
# whether CCT2 retains an independent direct effect.
# ----------------------------------------------------------------
cat("\n=== Test 3: log NfL ~ CCT2 + age + sex + subtype ===\n")
fit3 <- lm(logNfL ~ CCT2_expr + AGE_AT_SYMPTOM_ONSET + SEX + cluster_final,
           data = df_clean)
print(summary(fit3)$coefficients)
cat("Adj R^2:", round(summary(fit3)$adj.r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 4: Within-subtype regression
# Forces CCT2 to operate as a continuous trait at fixed subtype.
# If significant in either group, CCT2 carries graded signal that
# the categorical label misses.
# ----------------------------------------------------------------
cat("\n=== Test 4: Within-subtype CCT2 → log NfL ===\n")
for (sub in levels(df_clean$cluster_final)) {
  d_sub <- df_clean[df_clean$cluster_final == sub, ]
  cat("\n--", sub, "(n =", nrow(d_sub), ") --\n")
  f <- lm(logNfL ~ CCT2_expr + AGE_AT_SYMPTOM_ONSET + SEX, data = d_sub)
  print(summary(f)$coefficients["CCT2_expr", ])
}

# ----------------------------------------------------------------
# Test 5: Model comparison (BIC)
# Which encoding of S-ALS-3 biology best explains NfL — the 
# discrete label, the continuous gene, or both? BIC penalizes
# extra parameters, so 'both' wins only if CCT2 adds real signal.
# ----------------------------------------------------------------
cat("\n=== Test 5: Model comparison (BIC) ===\n")
m_sub  <- lm(logNfL ~ cluster_final + AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
m_cct2 <- lm(logNfL ~ CCT2_expr     + AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
m_both <- lm(logNfL ~ CCT2_expr + cluster_final + AGE_AT_SYMPTOM_ONSET + SEX,
             data = df_clean)
bic_tab <- data.frame(
  model  = c("subtype only", "CCT2 only", "subtype + CCT2"),
  BIC    = c(BIC(m_sub), BIC(m_cct2), BIC(m_both)),
  adj_R2 = c(summary(m_sub)$adj.r.squared,
             summary(m_cct2)$adj.r.squared,
             summary(m_both)$adj.r.squared))
print(bic_tab)

# ----------------------------------------------------------------
# Test 6: Sensitivity — full cohort (include S-ALS-1)
# Confirms the relationship isn't an artifact of restricting to the
# clean half. Direction and magnitude should be similar.
# ----------------------------------------------------------------
cat("\n=== Test 6: Full cohort (all three subtypes) ===\n")
fit_full <- lm(logNfL ~ CCT2_expr + AGE_AT_SYMPTOM_ONSET + SEX,
               data = df_nfl_cct2)
print(summary(fit_full)$coefficients)






# ============================================================
# AXIS_METABOLIC → log NfL on the clean subset (S-ALS-2 + S-ALS-3)
# Companion to CCT2 → NfL. Tests whether pathway-level signal is
# graded where gene-level signal was categorical.
#
# Predicted two-level story:
#   - Gene resolution (CCT2):     categorical classifier of subtype
#   - Pathway resolution (AXIS):  continuous biomarker correlate
# ============================================================

# --- Build the analysis frame ---
df_axis <- merge(
  ax_clean_dedup[, c("Participant_ID",
                     "cluster_final")],
  df_nfl_uniq, by = "Participant_ID")

# Bring CCT2 in for the head-to-head test
cct2_vec <- expr_clean_mat["CCT2", ]
pid_from_sample <- sub("-[0-9]+-T$", "", names(cct2_vec))
cct2_df <- data.frame(Participant_ID = pid_from_sample,
                      CCT2_expr = as.numeric(cct2_vec))
df_axis <- merge(df_axis, cct2_df, by = "Participant_ID")

# (Merge AGE_AT_SYMPTOM_ONSET / SEX from clinical_df here if not already on df_nfl_uniq)

# Clean subset
df_clean <- subset(df_axis, cluster_final %in% c("S-ALS-2","S-ALS-3"))
df_clean$cluster_final <- droplevels(factor(df_clean$cluster_final))
cat("Clean subset n:", nrow(df_clean), "\n")
print(table(df_clean$cluster_final))

# ----------------------------------------------------------------
# Test 1: Univariate AXIS_METABOLIC → log NfL
# ----------------------------------------------------------------
cat("\n=== Test 1: log NfL ~ AXIS_METABOLIC (univariate) ===\n")
fit1 <- lm(logNfL ~ AXIS_METABOLIC, data = df_clean)
print(summary(fit1)$coefficients)
cat("R^2:", round(summary(fit1)$r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 2: + age + sex
# ----------------------------------------------------------------
cat("\n=== Test 2: log NfL ~ AXIS_METABOLIC + age + sex ===\n")
fit2 <- lm(logNfL ~ AXIS_METABOLIC + AGE_AT_SYMPTOM_ONSET + SEX,
           data = df_clean)
print(summary(fit2)$coefficients)
cat("Adj R^2:", round(summary(fit2)$adj.r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 3: + subtype — direct parallel to Test 3 in CCT2 analysis
# CCT2 attenuated 67% and went null. The question: does the
# pathway-level signal retain a direct effect after subtype?
# ----------------------------------------------------------------
cat("\n=== Test 3: log NfL ~ AXIS_METABOLIC + age + sex + subtype ===\n")
fit3 <- lm(logNfL ~ AXIS_METABOLIC + AGE_AT_SYMPTOM_ONSET + SEX + cluster_final,
           data = df_clean)
print(summary(fit3)$coefficients)
cat("Adj R^2:", round(summary(fit3)$adj.r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 4: Within-subtype regressions — the load-bearing test
# CCT2 was null in both subtypes (β ~ 0.09, p > 0.58).
# A significant within-subtype AXIS_METABOLIC effect = direct
# evidence pathway signal is graded where gene signal isn't.
# ----------------------------------------------------------------
cat("\n=== Test 4: Within-subtype AXIS_METABOLIC → log NfL ===\n")
for (sub in levels(df_clean$cluster_final)) {
  d_sub <- df_clean[df_clean$cluster_final == sub, ]
  cat("\n--", sub, "(n =", nrow(d_sub), ") --\n")
  f <- lm(logNfL ~ AXIS_METABOLIC + AGE_AT_SYMPTOM_ONSET + SEX, data = d_sub)
  print(summary(f)$coefficients["AXIS_METABOLIC", ])
}

# ----------------------------------------------------------------
# Test 5: Head-to-head — CCT2 vs AXIS_METABOLIC
# Both in one model with subtype controlled. Direct test of the
# gene-vs-pathway resolution claim.
# ----------------------------------------------------------------
cat("\n=== Test 5: Head-to-head CCT2 vs AXIS_METABOLIC ===\n")
fit_hth <- lm(logNfL ~ CCT2_expr + AXIS_METABOLIC + AGE_AT_SYMPTOM_ONSET +
                SEX + cluster_final, data = df_clean)
print(summary(fit_hth)$coefficients)
cat("Adj R^2:", round(summary(fit_hth)$adj.r.squared, 4), "\n")

# Collinearity check — CCT2, AXIS_METABOLIC, and cluster_final all
# carry overlapping signal. High VIFs would qualify the head-to-head.
if (requireNamespace("car", quietly = TRUE)) {
  cat("\nVIFs:\n")
  print(car::vif(fit_hth))
}

# ----------------------------------------------------------------
# Test 6: BIC model comparison — extends the CCT2 BIC table
# ----------------------------------------------------------------
cat("\n=== Test 6: Model comparison (BIC) ===\n")
m_sub      <- lm(logNfL ~ cluster_final + AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
m_cct2     <- lm(logNfL ~ CCT2_expr + AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
m_axis     <- lm(logNfL ~ AXIS_METABOLIC + AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
m_sub_axis <- lm(logNfL ~ AXIS_METABOLIC + cluster_final + AGE_AT_SYMPTOM_ONSET + SEX,
                 data = df_clean)
m_all      <- lm(logNfL ~ CCT2_expr + AXIS_METABOLIC + cluster_final +
                   AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
bic_tab <- data.frame(
  model  = c("subtype only", "CCT2 only", "AXIS only",
             "subtype + axis", "subtype + axis + CCT2"),
  BIC    = c(BIC(m_sub), BIC(m_cct2), BIC(m_axis), BIC(m_sub_axis), BIC(m_all)),
  adj_R2 = c(summary(m_sub)$adj.r.squared,
             summary(m_cct2)$adj.r.squared,
             summary(m_axis)$adj.r.squared,
             summary(m_sub_axis)$adj.r.squared,
             summary(m_all)$adj.r.squared))
bic_tab$dBIC <- bic_tab$BIC - min(bic_tab$BIC)
print(bic_tab[order(bic_tab$BIC), ])

# ----------------------------------------------------------------
# Test 7: Sensitivity — full cohort
# ----------------------------------------------------------------
cat("\n=== Test 7: Full cohort (all three subtypes) ===\n")
fit_full <- lm(logNfL ~ AXIS_METABOLIC + AGE_AT_SYMPTOM_ONSET + SEX,
               data = df_axis)
print(summary(fit_full)$coefficients)


# ============================================================
# AXIS_METABOLIC → ALSFRS-R progression slope on clean subset
# Tests whether the existing β = −0.72 (p = 0.035) finding
# survives transcriptomic subtype adjustment.
#
# Two scenarios:
#   1. Collapses under subtype → axis is a pure state marker;
#      single coherent story for the paper
#   2. Survives subtype → axis affects functional decline
#      through a pathway not mediated by axonal damage burden
# ============================================================

# --- Build the analysis frame ---
# ALSFRS-R slope lives on the clinical metadata frame — adapt the
# variable name if it's stored as ALSFRS_R_PROGRESSION_SLOPE elsewhere
df_alsfrs <- merge(
  ax_clean_dedup[, c("Participant_ID", "AXIS_METABOLIC",
                     "AXIS_CELL_CYCLE", "AXIS_INFLAMMATION",
                     "cluster_final")],
  sporadic_final_check[, c("Participant_ID", "ALSFRS_R_PROGRESSION_SLOPE",
                  "AGE_AT_SYMPTOM_ONSET", "SEX")],
  by = "Participant_ID")

# Bring CCT2 in for the parallel head-to-head
cct2_vec <- expr_clean_mat["CCT2", ]
pid_from_sample <- sub("-[0-9]+-T$", "", names(cct2_vec))
cct2_df <- data.frame(Participant_ID = pid_from_sample,
                      CCT2_expr = as.numeric(cct2_vec))
df_alsfrs <- merge(df_alsfrs, cct2_df, by = "Participant_ID")

# Drop rows missing the outcome (not all patients have longitudinal scores)
df_alsfrs <- df_alsfrs[!is.na(df_alsfrs$ALSFRS_R_PROGRESSION_SLOPE), ]

# Clean subset
df_clean <- subset(df_alsfrs, cluster_final %in% c("S-ALS-2","S-ALS-3"))
df_clean$cluster_final <- droplevels(factor(df_clean$cluster_final))
cat("Clean subset n with ALSFRS slope:", nrow(df_clean), "\n")
print(table(df_clean$cluster_final))

# ----------------------------------------------------------------
# Test 1: Univariate AXIS_METABOLIC → ALSFRS slope
# Replicates the published direction in the clean subset
# ----------------------------------------------------------------
cat("\n=== Test 1: ALSFRS slope ~ AXIS_METABOLIC (univariate) ===\n")
fit1 <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC, data = df_clean)
print(summary(fit1)$coefficients)
cat("R^2:", round(summary(fit1)$r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 2: + age + sex
# Standard adjustment — age in particular affects progression rate
# ----------------------------------------------------------------
cat("\n=== Test 2: ALSFRS slope ~ AXIS_METABOLIC + age + sex ===\n")
fit2 <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC +
             AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
print(summary(fit2)$coefficients)
cat("Adj R^2:", round(summary(fit2)$adj.r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 3: + subtype — the load-bearing test
# Direct parallel to the NfL Test 3. Does AXIS_METABOLIC retain
# a continuous effect on functional decline beyond subtype?
# ----------------------------------------------------------------
cat("\n=== Test 3: ALSFRS slope ~ AXIS_METABOLIC + age + sex + subtype ===\n")
fit3 <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC +
             AGE_AT_SYMPTOM_ONSET + SEX + cluster_final, data = df_clean)
print(summary(fit3)$coefficients)
cat("Adj R^2:", round(summary(fit3)$adj.r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 4: Within-subtype regressions
# Power-limited (S-ALS-3 likely n < 50 with ALSFRS), but
# directionally informative
# ----------------------------------------------------------------
cat("\n=== Test 4: Within-subtype AXIS_METABOLIC → ALSFRS slope ===\n")
for (sub in levels(df_clean$cluster_final)) {
  d_sub <- df_clean[df_clean$cluster_final == sub, ]
  cat("\n--", sub, "(n =", nrow(d_sub), ") --\n")
  f <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC +
            AGE_AT_SYMPTOM_ONSET + SEX, data = d_sub)
  print(summary(f)$coefficients["AXIS_METABOLIC", ])
}

# ----------------------------------------------------------------
# Test 5: BIC model comparison — same structure as NfL panel
# ----------------------------------------------------------------
cat("\n=== Test 5: Model comparison (BIC) ===\n")
m_sub      <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ cluster_final +
                   AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
m_axis     <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC +
                   AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
m_cct2     <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr +
                   AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
m_sub_axis <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC +
                   cluster_final + AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
bic_tab <- data.frame(
  model  = c("subtype only", "AXIS only", "CCT2 only", "subtype + axis"),
  BIC    = c(BIC(m_sub), BIC(m_axis), BIC(m_cct2), BIC(m_sub_axis)),
  adj_R2 = c(summary(m_sub)$adj.r.squared,
             summary(m_axis)$adj.r.squared,
             summary(m_cct2)$adj.r.squared,
             summary(m_sub_axis)$adj.r.squared))
bic_tab$dBIC <- bic_tab$BIC - min(bic_tab$BIC)
print(bic_tab[order(bic_tab$BIC), ])

# ----------------------------------------------------------------
# Test 6: Sensitivity — full cohort
# Replicates the original published β = −0.72 framing and
# checks whether restriction to clean subset matters
# ----------------------------------------------------------------
cat("\n=== Test 6: Full cohort (all three subtypes) ===\n")
fit_full <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC +
                 AGE_AT_SYMPTOM_ONSET + SEX + cluster_final,
               data = df_alsfrs)
print(summary(fit_full)$coefficients)

# CCT2 → ALSFRS slope (parallel structure to AXIS_METABOLIC)
fit1 <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr, data = df_clean)
fit2 <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr + AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
fit3 <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr + AGE_AT_SYMPTOM_ONSET + SEX + cluster_final, data = df_clean)

# Within-subtype
for (sub in levels(df_clean$cluster_final)) {
  d_sub <- df_clean[df_clean$cluster_final == sub, ]
  cat("\n--", sub, "(n =", nrow(d_sub), ") --\n")
  f <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr + AGE_AT_SYMPTOM_ONSET + SEX, data = d_sub)
  print(summary(f)$coefficients["CCT2_expr", ])
}

# Head-to-head
fit_hth <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr + AXIS_METABOLIC +
                AGE_AT_SYMPTOM_ONSET + SEX + cluster_final, data = df_clean)
summary(fit_hth)$coefficients





# ============================================================
# CCT2 → ALSFRS-R progression slope on clean subset
# Companion to AXIS_METABOLIC → ALSFRS panel.
#
# Test 5 in the AXIS panel showed CCT2-only beat subtype-only
# on BIC for this endpoint. This script runs the full ladder
# to determine whether that ranking reflects a genuine graded
# effect or a parametric artifact.
#
# Two scenarios:
#   1. CCT2 survives subtype adjustment → genuine endpoint
#      dissociation. CCT2 is a graded biomarker for functional
#      decline, categorical state marker for NfL.
#   2. CCT2 collapses like on NfL → simpler "subtype is the
#      unit" story holds across endpoints; BIC win is artifact.
# ============================================================

# Assumes df_clean already contains: CCT2_expr, AXIS_METABOLIC,
# cluster_final, ALSFRS_R_PROGRESSION_SLOPE, AGE_AT_SYMPTOM_ONSET, SEX
# (built in the previous AXIS_METABOLIC → ALSFRS script)

cat("Clean subset n with ALSFRS slope:", nrow(df_clean), "\n")
print(table(df_clean$cluster_final))

# ----------------------------------------------------------------
# Test 1: Univariate CCT2 → ALSFRS slope
# Direct parallel to AXIS Test 1 (β = −0.10, p = 0.049)
# ----------------------------------------------------------------
cat("\n=== Test 1: ALSFRS slope ~ CCT2 (univariate) ===\n")
fit1 <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr, data = df_clean)
print(summary(fit1)$coefficients)
cat("R^2:", round(summary(fit1)$r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 2: + age + sex
# AXIS attenuated to p = 0.118 here. Does CCT2 retain its
# univariate effect more robustly under standard adjustment?
# ----------------------------------------------------------------
cat("\n=== Test 2: ALSFRS slope ~ CCT2 + age + sex ===\n")
fit2 <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr +
             AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
print(summary(fit2)$coefficients)
cat("Adj R^2:", round(summary(fit2)$adj.r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 3: + subtype — the load-bearing test
# This is the test that decides between the two scenarios.
# CCT2 collapsed completely on NfL Test 3 (β 0.23 → 0.08, p 0.03 → 0.57).
# If it collapses here too, the BIC ranking from the AXIS panel
# is an artifact. If it survives, real endpoint dissociation.
# ----------------------------------------------------------------
cat("\n=== Test 3: ALSFRS slope ~ CCT2 + age + sex + subtype ===\n")
fit3 <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr +
             AGE_AT_SYMPTOM_ONSET + SEX + cluster_final, data = df_clean)
print(summary(fit3)$coefficients)
cat("Adj R^2:", round(summary(fit3)$adj.r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 4: Within-subtype regressions
# AXIS within-subtype was null in both groups but directionally
# negative. CCT2 within-subtype on NfL was null (β ~ 0.09, p > 0.58).
# If CCT2 within-subtype on ALSFRS shows graded signal anywhere,
# that's the strongest single piece of evidence for the
# endpoint-dissociation framing.
# ----------------------------------------------------------------
cat("\n=== Test 4: Within-subtype CCT2 → ALSFRS slope ===\n")
for (sub in levels(df_clean$cluster_final)) {
  d_sub <- df_clean[df_clean$cluster_final == sub, ]
  cat("\n--", sub, "(n =", nrow(d_sub), ") --\n")
  f <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr +
            AGE_AT_SYMPTOM_ONSET + SEX, data = d_sub)
  print(summary(f)$coefficients["CCT2_expr", ])
}

# ----------------------------------------------------------------
# Test 5: Head-to-head — CCT2 vs AXIS_METABOLIC vs subtype
# All three encodings of the same biology in one model.
# AXIS panel head-to-head (NfL outcome) had CCT2 p = 0.46,
# AXIS p = 0.60, subtype p = 0.056 — categorical winning.
# Run the same structure on ALSFRS to compare directly.
# ----------------------------------------------------------------
cat("\n=== Test 5: Head-to-head CCT2 vs AXIS_METABOLIC vs subtype ===\n")
fit_hth <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr + AXIS_METABOLIC +
                AGE_AT_SYMPTOM_ONSET + SEX + cluster_final, data = df_clean)
print(summary(fit_hth)$coefficients)
cat("Adj R^2:", round(summary(fit_hth)$adj.r.squared, 4), "\n")

if (requireNamespace("car", quietly = TRUE)) {
  cat("\nVIFs:\n")
  print(car::vif(fit_hth))
}

# ----------------------------------------------------------------
# Test 6: BIC model comparison — extends the AXIS panel BIC table
# with the missing rows (CCT2 + subtype, all three together).
# AXIS panel showed: CCT2-only (394.1) < AXIS-only (396.5)
# < subtype-only (397.4) < subtype + axis (401.7).
# Need: subtype + CCT2, and full model — to see if categorical
# information adds anything to CCT2 on this endpoint.
# ----------------------------------------------------------------
cat("\n=== Test 6: Model comparison (BIC) ===\n")
m_sub      <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ cluster_final +
                   AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
m_cct2     <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr +
                   AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
m_axis     <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC +
                   AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
m_sub_cct2 <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr + cluster_final +
                   AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
m_sub_axis <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC + cluster_final +
                   AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
m_all      <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr + AXIS_METABOLIC +
                   cluster_final + AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
bic_tab <- data.frame(
  model  = c("subtype only", "CCT2 only", "AXIS only",
             "subtype + CCT2", "subtype + axis", "all three"),
  BIC    = c(BIC(m_sub), BIC(m_cct2), BIC(m_axis),
             BIC(m_sub_cct2), BIC(m_sub_axis), BIC(m_all)),
  adj_R2 = c(summary(m_sub)$adj.r.squared,
             summary(m_cct2)$adj.r.squared,
             summary(m_axis)$adj.r.squared,
             summary(m_sub_cct2)$adj.r.squared,
             summary(m_sub_axis)$adj.r.squared,
             summary(m_all)$adj.r.squared))
bic_tab$dBIC <- bic_tab$BIC - min(bic_tab$BIC)
print(bic_tab[order(bic_tab$BIC), ])

# ----------------------------------------------------------------
# Test 7: Sensitivity — full cohort (include S-ALS-1)
# Confirms direction and rough magnitude aren't an artifact of
# restricting to the clean subset. AXIS Test 7 showed β = −0.12,
# p = 0.062 in full cohort — borderline. Compare directly.
# ----------------------------------------------------------------
cat("\n=== Test 7: Full cohort (all three subtypes) ===\n")
fit_full <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ CCT2_expr +
                 AGE_AT_SYMPTOM_ONSET + SEX + cluster_final, data = df_alsfrs)
print(summary(fit_full)$coefficients)




# ============================================================
# PRSet specificity matrix re-run within the clean subset
# Original full-cohort result: PRS_CELL_CYCLE → AXIS_CELL_CYCLE
# R²=0.018, p=0.012. Question: is the on-diagonal hit
# carried by the differentiation-confounded S-ALS-1 subset, or
# does it survive in the unconfounded half of the cohort?
# ============================================================

library(dplyr)

workdir <- file.path(ALSPROJ_BASE_DIR, "prset")

# Reload PRS scores (all three pathway PRS in one .best file)
prs_scores <- read.table(file.path(workdir, "prset_axis_metabolic.best"),
                         header = TRUE, stringsAsFactors = FALSE) %>%
  filter(In_Regression == "Yes") %>%
  select(IID,
         PRS_METABOLIC    = AXIS_METABOLIC,
         PRS_CELL_CYCLE   = AXIS_CELL_CYCLE,
         PRS_INFLAMMATION = AXIS_INFLAMMATION)

# Rebuild analysis_df with cluster_final + axis scores + CCT2
cct2_vec <- expr_clean_mat["CCT2", ]
pid_from_sample <- sub("-[0-9]+-T$", "", names(cct2_vec))
cct2_df <- data.frame(Participant_ID = pid_from_sample,
                      CCT2_expr = as.numeric(cct2_vec))

analysis_df <- prs_meta_clean %>%
  select(VCF_ID, Participant_ID, AXIS_METABOLIC, AXIS_CELL_CYCLE,
         AXIS_INFLAMMATION, cluster_final, SEX) %>%
  inner_join(prs_scores, by = c("VCF_ID" = "IID")) %>%
  inner_join(cct2_df,    by = "Participant_ID")

# Restrict to clean half — the load-bearing subset for the paper
df_clean <- subset(analysis_df,
                   cluster_final %in% c("S-ALS-2", "S-ALS-3"))
df_clean$cluster_final <- droplevels(factor(df_clean$cluster_final))
cat("Clean subset n with PRS + axis scores:", nrow(df_clean), "\n")
print(table(df_clean$cluster_final))

# ----------------------------------------------------------------
# Test 1: 3×3 specificity matrix — clean subset, sex-adjusted
# Direct comparison to the full-cohort matrix
# ----------------------------------------------------------------
prs_vars  <- c("PRS_METABOLIC", "PRS_CELL_CYCLE", "PRS_INFLAMMATION")
axis_vars <- c("AXIS_METABOLIC", "AXIS_CELL_CYCLE", "AXIS_INFLAMMATION")

r2_mat <- matrix(NA, 3, 3, dimnames = list(prs_vars, axis_vars))
b_mat  <- matrix(NA, 3, 3, dimnames = list(prs_vars, axis_vars))
p_mat  <- matrix(NA, 3, 3, dimnames = list(prs_vars, axis_vars))

for (prs in prs_vars) for (ax in axis_vars) {
  fit <- lm(as.formula(paste(ax, "~", prs, "+ SEX")), data = df_clean)
  s   <- summary(fit)
  r2_mat[prs, ax] <- s$r.squared
  b_mat[prs, ax]  <- coef(s)[prs, "Estimate"]
  p_mat[prs, ax]  <- coef(s)[prs, "Pr(>|t|)"]
}

cat("\n=== Test 1: R² matrix (sex-adjusted, clean subset) ===\n")
print(round(r2_mat, 4))
cat("\n=== p-value matrix ===\n")
print(round(p_mat, 4))

cat("\nDiagonal hits:\n")
for (i in 1:3) cat(sprintf("  %s -> %s: R^2 = %.4f, p = %.4f\n",
                           prs_vars[i], axis_vars[i],
                           r2_mat[i, i], p_mat[i, i]))
cat("\nMean off-diagonal R^2:",
    round(mean(r2_mat[row(r2_mat) != col(r2_mat)]), 4), "\n")

# ----------------------------------------------------------------
# Test 2: PRS_CELL_CYCLE → AXIS_CELL_CYCLE — direct, subtype-adjusted
# Within the clean subset — does the on-diagonal effect survive?
# ----------------------------------------------------------------
cat("\n=== Test 2: PRS_CELL_CYCLE → AXIS_CELL_CYCLE (clean subset) ===\n")
fit_uni <- lm(AXIS_CELL_CYCLE ~ PRS_CELL_CYCLE + SEX, data = df_clean)
fit_sub <- lm(AXIS_CELL_CYCLE ~ PRS_CELL_CYCLE + SEX + cluster_final,
              data = df_clean)
cat("\nUnivariate (sex-adjusted):\n")
print(summary(fit_uni)$coefficients)
cat("\n+ subtype:\n")
print(summary(fit_sub)$coefficients)

# ----------------------------------------------------------------
# Test 3: PRS_CELL_CYCLE → CCT2 — gene-level anchor, clean subset
# Original full-cohort hit: β = -5.05, p = 0.042
# Question: does it survive in the clean half?
# ----------------------------------------------------------------
cat("\n=== Test 3: PRS_CELL_CYCLE → CCT2 (clean subset) ===\n")
fit_cct2 <- lm(CCT2_expr ~ PRS_CELL_CYCLE + SEX, data = df_clean)
print(summary(fit_cct2)$coefficients)
cat("R^2:", round(summary(fit_cct2)$r.squared, 4), "\n")

cat("\n=== Test 3b: + subtype ===\n")
fit_cct2_sub <- lm(CCT2_expr ~ PRS_CELL_CYCLE + SEX + cluster_final,
                   data = df_clean)
print(summary(fit_cct2_sub)$coefficients)


# ============================================================
# PRS_CELL_CYCLE → AXIS_CELL_CYCLE under glial + DQ adjustment
# Test 2 of the sensitivity panel showed AXIS_CELL_CYCLE has
# 48.5% of variance explained by glial + DQ composites.
# Question: does the PRS effect (R²=0.018, p=0.012) survive
# adjustment for that composition signal?
# ============================================================

library(dplyr)

# ----------------------------------------------------------------
# Build composite scores
# These were computed in Test 2 of the sensitivity panel.
# If they're already on a workspace object, skip the build step
# and inner_join. If not, compute from expr_clean_mat.
# ----------------------------------------------------------------
astrocyte_markers   <- c("GFAP","AQP4","ALDH1L1","S100B","SLC1A3",
                         "SLC1A2","GJA1","AGT")
microglia_markers   <- c("AIF1","CX3CR1","P2RY12","TMEM119","CD68",
                         "ITGAM","TYROBP")
oligo_markers       <- c("OLIG2","MBP","MOG","PLP1","SOX10","GALC")
maturity_markers    <- c("MAP2","SYN1","NEFL","NEFM","SNAP25","RBFOX3",
                         "STMN2","CHAT","ISL1","MNX1","SLC18A3","TUBB3",
                         "ELAVL3","ELAVL4")
progenitor_markers  <- c("PAX6","SOX2","NES","MKI67","TOP2A","PCNA",
                         "VIM","SOX1","MSI1")

# Helper — mean z-scored expression per sample for a marker set
composite_score <- function(mat, markers) {
  m <- markers[markers %in% rownames(mat)]
  if (length(m) == 0) return(NULL)
  z <- t(scale(t(mat[m, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

ast_score <- composite_score(expr_clean_mat, astrocyte_markers)
mic_score <- composite_score(expr_clean_mat, microglia_markers)
oli_score <- composite_score(expr_clean_mat, oligo_markers)
mat_score <- composite_score(expr_clean_mat, maturity_markers)
prog_score <- composite_score(expr_clean_mat, progenitor_markers)
diff_quality <- mat_score - prog_score

composites <- data.frame(
  sample_id  = colnames(expr_clean_mat),
  astrocyte  = ast_score,
  microglia  = mic_score,
  oligo      = oli_score,
  maturity   = mat_score,
  progenitor = prog_score,
  diff_quality = diff_quality,
  stringsAsFactors = FALSE
)
composites$Participant_ID <- sub("-[0-9]+-T$", "", composites$sample_id)
composites <- composites %>%
  distinct(Participant_ID, .keep_all = TRUE)

# Merge with analysis_df (built in Script 1)
df_g <- analysis_df %>%
  inner_join(composites, by = "Participant_ID")

cat("n with PRS + axis + composites:", nrow(df_g), "\n")

# ----------------------------------------------------------------
# Test 1: PRS_CELL_CYCLE → AXIS_CELL_CYCLE — full cohort, progressive adjustment
# ----------------------------------------------------------------
cat("\n=== Test 1: AXIS_CELL_CYCLE ~ PRS_CELL_CYCLE, progressive adjustment ===\n")

m0 <- lm(AXIS_CELL_CYCLE ~ PRS_CELL_CYCLE + SEX, data = df_g)
m1 <- lm(AXIS_CELL_CYCLE ~ PRS_CELL_CYCLE + SEX +
           astrocyte + microglia + oligo, data = df_g)
m2 <- lm(AXIS_CELL_CYCLE ~ PRS_CELL_CYCLE + SEX +
           astrocyte + microglia + oligo + diff_quality, data = df_g)
m3 <- lm(AXIS_CELL_CYCLE ~ PRS_CELL_CYCLE + SEX +
           astrocyte + microglia + oligo + diff_quality +
           maturity + progenitor, data = df_g)

results <- data.frame(
  model = c("baseline (sex)", "+ glial panel", "+ glial + DQ",
            "+ glial + DQ + mat + prog"),
  beta_PRS = c(coef(m0)["PRS_CELL_CYCLE"],
               coef(m1)["PRS_CELL_CYCLE"],
               coef(m2)["PRS_CELL_CYCLE"],
               coef(m3)["PRS_CELL_CYCLE"]),
  p_PRS    = c(summary(m0)$coefficients["PRS_CELL_CYCLE", 4],
               summary(m1)$coefficients["PRS_CELL_CYCLE", 4],
               summary(m2)$coefficients["PRS_CELL_CYCLE", 4],
               summary(m3)$coefficients["PRS_CELL_CYCLE", 4]),
  adj_R2   = c(summary(m0)$adj.r.squared,
               summary(m1)$adj.r.squared,
               summary(m2)$adj.r.squared,
               summary(m3)$adj.r.squared)
)
print(results)

# ----------------------------------------------------------------
# Test 2: per-covariate effects in the fully-adjusted model
# Confirms each covariate is doing independent work
# ----------------------------------------------------------------
cat("\n=== Test 2: full-model coefficients ===\n")
print(summary(m3)$coefficients)

# ----------------------------------------------------------------
# Test 3: same progression on PRS_CELL_CYCLE → CCT2
# Closes the gene-level genomic anchor under composition adjustment
# ----------------------------------------------------------------
cat("\n=== Test 3: CCT2 ~ PRS_CELL_CYCLE, progressive adjustment ===\n")

c0 <- lm(CCT2_expr ~ PRS_CELL_CYCLE + SEX, data = df_g)
c1 <- lm(CCT2_expr ~ PRS_CELL_CYCLE + SEX +
           astrocyte + microglia + oligo, data = df_g)
c2 <- lm(CCT2_expr ~ PRS_CELL_CYCLE + SEX +
           astrocyte + microglia + oligo + diff_quality, data = df_g)
c3 <- lm(CCT2_expr ~ PRS_CELL_CYCLE + SEX +
           astrocyte + microglia + oligo + diff_quality +
           maturity + progenitor, data = df_g)

cct2_results <- data.frame(
  model = c("baseline (sex)", "+ glial panel", "+ glial + DQ",
            "+ glial + DQ + mat + prog"),
  beta_PRS = c(coef(c0)["PRS_CELL_CYCLE"],
               coef(c1)["PRS_CELL_CYCLE"],
               coef(c2)["PRS_CELL_CYCLE"],
               coef(c3)["PRS_CELL_CYCLE"]),
  p_PRS    = c(summary(c0)$coefficients["PRS_CELL_CYCLE", 4],
               summary(c1)$coefficients["PRS_CELL_CYCLE", 4],
               summary(c2)$coefficients["PRS_CELL_CYCLE", 4],
               summary(c3)$coefficients["PRS_CELL_CYCLE", 4]),
  adj_R2   = c(summary(c0)$adj.r.squared,
               summary(c1)$adj.r.squared,
               summary(c2)$adj.r.squared,
               summary(c3)$adj.r.squared)
)
print(cct2_results)


# ============================================================
# CCT2 → disease duration on deceased patients
# Parallel structure to the CCT2 → ALSFRS-R analysis.
#
# Survival bias caveat: deceased subset is enriched for fast
# progressors. NfL → duration was strong (β=−4.61, p=1e-9) in
# this same subset; previous subtype ANOVA was null (p=0.61).
# Question: does CCT2 predict disease duration as a graded
# trait beyond subtype, paralleling the ALSFRS-R result?
# ============================================================

library(dplyr)

# ----------------------------------------------------------------
# Build the analysis frame
# ----------------------------------------------------------------
cct2_vec <- expr_clean_mat["CCT2", ]
pid_from_sample <- sub("-[0-9]+-T$", "", names(cct2_vec))
cct2_df <- data.frame(Participant_ID = pid_from_sample,
                      CCT2_expr = as.numeric(cct2_vec))

# sporadic_final_check has AGE_AT_DEATH and AGE_AT_SYMPTOM_ONSET
df_dur <- sporadic_final_check %>%
  distinct(Participant_ID, .keep_all = TRUE) %>%
  filter(!is.na(AGE_AT_DEATH), !is.na(AGE_AT_SYMPTOM_ONSET)) %>%
  mutate(disease_duration = AGE_AT_DEATH - AGE_AT_SYMPTOM_ONSET)

# Join CCT2 + axis scores + cluster_final
df_dur <- df_dur %>%
  inner_join(cct2_df, by = "Participant_ID") %>%
  inner_join(ax_clean_dedup %>%
               select(Participant_ID, cluster_final),
             by = "Participant_ID")

cat("Deceased patients with CCT2 + cluster:", nrow(df_dur), "\n")
print(table(df_dur$cluster_final))

# Clean subset
df_clean <- subset(df_dur, cluster_final %in% c("S-ALS-2", "S-ALS-3"))
df_clean$cluster_final <- droplevels(factor(df_clean$cluster_final))
cat("\nClean subset (deceased S-ALS-2 + S-ALS-3):", nrow(df_clean), "\n")
print(table(df_clean$cluster_final))

# ----------------------------------------------------------------
# Test 1: Univariate CCT2 → duration
# Direction expected: high CCT2 (S-ALS-3 marker) → shorter duration
# ----------------------------------------------------------------
cat("\n=== Test 1: duration ~ CCT2 (univariate, clean subset) ===\n")
fit1 <- lm(disease_duration ~ CCT2_expr, data = df_clean)
print(summary(fit1)$coefficients)
cat("R^2:", round(summary(fit1)$r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 2: + age at onset + sex
# Age at onset is a known driver of survival in ALS
# ----------------------------------------------------------------
cat("\n=== Test 2: duration ~ CCT2 + age + sex ===\n")
fit2 <- lm(disease_duration ~ CCT2_expr + AGE_AT_SYMPTOM_ONSET + SEX,
           data = df_clean)
print(summary(fit2)$coefficients)
cat("Adj R^2:", round(summary(fit2)$adj.r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 3: + subtype — the load-bearing test
# Parallel to CCT2 → ALSFRS Test 3 (β=-0.60, p=0.075)
# If CCT2 retains direct effect, three-endpoint dissociation:
#   - NfL: subtype-mediated
#   - ALSFRS-R: CCT2 graded
#   - duration: ?
# ----------------------------------------------------------------
cat("\n=== Test 3: duration ~ CCT2 + age + sex + subtype ===\n")
fit3 <- lm(disease_duration ~ CCT2_expr + AGE_AT_SYMPTOM_ONSET + SEX +
             cluster_final, data = df_clean)
print(summary(fit3)$coefficients)
cat("Adj R^2:", round(summary(fit3)$adj.r.squared, 4), "\n")

# ----------------------------------------------------------------
# Test 4: Within-subtype regressions
# Power-limited (n likely ~30 in S-ALS-3 deceased)
# ----------------------------------------------------------------
cat("\n=== Test 4: Within-subtype CCT2 → duration ===\n")
for (sub in levels(df_clean$cluster_final)) {
  d_sub <- df_clean[df_clean$cluster_final == sub, ]
  cat("\n--", sub, "(n =", nrow(d_sub), ") --\n")
  if (nrow(d_sub) > 10) {
    f <- lm(disease_duration ~ CCT2_expr + AGE_AT_SYMPTOM_ONSET + SEX,
            data = d_sub)
    print(summary(f)$coefficients["CCT2_expr", ])
  } else {
    cat("Insufficient n for regression\n")
  }
}

# ----------------------------------------------------------------
# Test 5: Head-to-head — CCT2 vs AXIS_METABOLIC vs subtype
# Parallel to the CCT2 → ALSFRS head-to-head (CCT2 won)
# ----------------------------------------------------------------
cat("\n=== Test 5: Head-to-head ===\n")
fit_hth <- lm(disease_duration ~ CCT2_expr + AXIS_METABOLIC +
                AGE_AT_SYMPTOM_ONSET + SEX + cluster_final, data = df_clean)
print(summary(fit_hth)$coefficients)
cat("Adj R^2:", round(summary(fit_hth)$adj.r.squared, 4), "\n")

if (requireNamespace("car", quietly = TRUE)) {
  cat("\nVIFs:\n")
  print(car::vif(fit_hth))
}

# ----------------------------------------------------------------
# Test 6: BIC model comparison — same structure as ALSFRS panel
# ----------------------------------------------------------------
cat("\n=== Test 6: Model comparison (BIC) ===\n")
m_sub      <- lm(disease_duration ~ cluster_final + AGE_AT_SYMPTOM_ONSET + SEX,
                 data = df_clean)
m_cct2     <- lm(disease_duration ~ CCT2_expr + AGE_AT_SYMPTOM_ONSET + SEX,
                 data = df_clean)
m_axis     <- lm(disease_duration ~ AXIS_METABOLIC + AGE_AT_SYMPTOM_ONSET + SEX,
                 data = df_clean)
m_sub_cct2 <- lm(disease_duration ~ CCT2_expr + cluster_final +
                   AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)
m_all      <- lm(disease_duration ~ CCT2_expr + AXIS_METABOLIC +
                   cluster_final + AGE_AT_SYMPTOM_ONSET + SEX, data = df_clean)

bic_tab <- data.frame(
  model  = c("subtype only", "CCT2 only", "AXIS only",
             "subtype + CCT2", "all three"),
  BIC    = c(BIC(m_sub), BIC(m_cct2), BIC(m_axis),
             BIC(m_sub_cct2), BIC(m_all)),
  adj_R2 = c(summary(m_sub)$adj.r.squared,
             summary(m_cct2)$adj.r.squared,
             summary(m_axis)$adj.r.squared,
             summary(m_sub_cct2)$adj.r.squared,
             summary(m_all)$adj.r.squared))
bic_tab$dBIC <- bic_tab$BIC - min(bic_tab$BIC)
print(bic_tab[order(bic_tab$BIC), ])

# ----------------------------------------------------------------
# Test 7: Full cohort (include S-ALS-1)
# ----------------------------------------------------------------
cat("\n=== Test 7: Full cohort (all three subtypes) ===\n")
fit_full <- lm(disease_duration ~ CCT2_expr + AGE_AT_SYMPTOM_ONSET + SEX +
                 cluster_final, data = df_dur)
print(summary(fit_full)$coefficients)

# ----------------------------------------------------------------
# Test 8: Survival bias check — does CCT2 predict logNfL within
# the deceased subset? If NfL variance is compressed in this
# subset, that's why mediation analyses on duration get fuzzy.
# ----------------------------------------------------------------
df_dur_nfl <- merge(df_dur,
                    df_nfl_uniq %>% select(Participant_ID, logNfL),
                    by = "Participant_ID")
cat("\n=== Test 8: CCT2 → logNfL within deceased subset (n =",
    nrow(df_dur_nfl), ") ===\n")
fit_nfl <- lm(logNfL ~ CCT2_expr + AGE_AT_SYMPTOM_ONSET + SEX,
              data = df_dur_nfl)
print(summary(fit_nfl)$coefficients)


