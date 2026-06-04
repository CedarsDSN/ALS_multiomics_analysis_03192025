# ============================================================
# Refactored entrypoint for original script: LassoClassifier.R
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
if (isTRUE(ALS_ARGS$help)) { als_print_help("02_lasso_classifier"); quit(save = "no", status = 0) }
CFG <- als_load_config(ALS_ARGS$config)
ALSPROJ_BASE_DIR <- als_choose(ALS_ARGS$base_dir, CFG, "project.base_dir", getwd())
ALSPROJ_OUTDIR <- als_choose(ALS_ARGS$outdir, CFG, "outputs.root", file.path(getwd(), "results", "02_lasso_classifier"))
ALSPROJ_SCRIPT_ID <- "02_lasso_classifier"
als_init(ALSPROJ_SCRIPT_ID, ALSPROJ_BASE_DIR, ALSPROJ_OUTDIR, ALS_ARGS$redirect_outputs)
als_register_output_wrappers()
als_load_objects(ALS_ARGS$objects)
on.exit(als_write_session_info(ALSPROJ_SCRIPT_ID), add = TRUE)
# ---- End bootstrap ----

# ============================================================================
# LASSO ONE-VS-REST SUBTYPE CLASSIFIER — HALLMARK AXIS GENES
# ============================================================================
# Goal: minimal per-subtype gene panels using only Hallmark axis genes
# Design: 3 binary LASSOs (S-ALS-1 vs rest, S-ALS-2 vs rest, S-ALS-3 vs rest)
# Validation: nested CV — outer 10-fold for performance, inner CV for lambda
# ============================================================================

library(glmnet)
library(DESeq2)
library(msigdbr)
library(dplyr)
library(tibble)
library(pROC)
library(caret)

set.seed(42)  # reproducibility

# ---------------------------------------------------------------------------
# STEP 1 — Pull the Hallmark axis gene universe
# ---------------------------------------------------------------------------
# Use msigdbr to get the same gene sets used in your ssGSEA
hallmark <- msigdbr(species = "Homo sapiens", category = "H")

# These three names should match what you used for AXIS scoring
# Adjust if your axis-to-Hallmark mapping included multiple sets per axis
axis_sets <- list(
  AXIS_METABOLIC    = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                        "HALLMARK_FATTY_ACID_METABOLISM",
                        "HALLMARK_GLYCOLYSIS"),  # CHECK against your ssGSEA call
  AXIS_INFLAMMATION = c("HALLMARK_INFLAMMATORY_RESPONSE",
                        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                        "HALLMARK_INTERFERON_GAMMA_RESPONSE"),  # CHECK
  AXIS_CELL_CYCLE   = c("HALLMARK_E2F_TARGETS",
                        "HALLMARK_G2M_CHECKPOINT",
                        "HALLMARK_MYC_TARGETS_V1")  # CHECK
)
# IMPORTANT: replace these with the EXACT Hallmark sets that fed your ssGSEA.
# If you used a different mapping, this is the one place to fix it.

axis_genes <- hallmark |>
  filter(gs_name %in% unlist(axis_sets)) |>
  pull(gene_symbol) |>
  unique()

cat("Total unique Hallmark axis genes:", length(axis_genes), "\n")

# ---------------------------------------------------------------------------
# STEP 2 — Build expression matrix restricted to axis genes
# ---------------------------------------------------------------------------
# Assumes you have a DESeq2 dds object or VST matrix in memory
# If you have `dds_sporadic`, do:
# vst_mat <- assay(vst(dds_sporadic, blind = FALSE))
# Otherwise, swap in your existing VST matrix object name

# vst_mat: rows = genes (symbols), cols = sample IDs
# Make sure rownames are gene symbols, not Ensembl IDs (convert if needed)

axis_genes_present <- intersect(axis_genes, rownames(expr_clean_mat))
cat("Axis genes present in expression matrix:",
    length(axis_genes_present), "/", length(axis_genes), "\n")

X <- t(expr_clean_mat[axis_genes_present, ])  # samples x genes
cat("Expression matrix:", nrow(X), "samples x", ncol(X), "genes\n")

# ---------------------------------------------------------------------------
# STEP 3 — Align with subtype labels
# ---------------------------------------------------------------------------
# Use ax_clean_dedup (391 patients, one row per patient, earliest timepoint)
# Need to match sample IDs in X (rownames) to ax_clean_dedup

# Adjust the join key as needed — your sample_id has the -N-T suffix,
# while ax_clean_dedup is keyed on Participant_ID
sample_to_pid <- sub("-[0-9]+-T$", "", rownames(X))

label_df <- data.frame(
  sample_id = rownames(X),
  Participant_ID = sample_to_pid
) |>
  inner_join(ax_clean_dedup |>
               select(Participant_ID, cluster_final),
             by = "Participant_ID")

# Filter X to matched samples and order
X <- X[label_df$sample_id, ]
y <- factor(label_df$cluster_final)

cat("Final aligned dataset:", nrow(X), "samples\n")
print(table(y))


# ---------------------------------------------------------------------------
# STEP 3b — Restrict to one sample per patient (earliest timepoint)
# ---------------------------------------------------------------------------
# Match the same dedup logic used for ax_clean_dedup: lowest numeric suffix wins
sample_meta <- data.frame(
  sample_id = rownames(X),
  Participant_ID = sub("-[0-9]+-T$", "", rownames(X)),
  timepoint = as.integer(sub(".*-([0-9]+)-T$", "\\1", rownames(X)))
)

earliest_samples <- sample_meta |>
  group_by(Participant_ID) |>
  slice_min(timepoint, n = 1, with_ties = FALSE) |>
  ungroup() |>
  pull(sample_id)

cat("Samples before dedup:", nrow(X), "\n")
cat("Samples after earliest-timepoint dedup:", length(earliest_samples), "\n")

# Re-filter X and rebuild label_df
X <- X[earliest_samples, ]
label_df <- data.frame(
  sample_id = rownames(X),
  Participant_ID = sub("-[0-9]+-T$", "", rownames(X))
) |>
  inner_join(ax_clean_dedup |> select(Participant_ID, cluster_final),
             by = "Participant_ID")

X <- X[label_df$sample_id, ]
y <- factor(label_df$cluster_final)

cat("Final patient-level dataset:", nrow(X), "samples (= patients)\n")
print(table(y))
stopifnot(!any(duplicated(label_df$Participant_ID)))


# ============================================================================
# STEP 4 — Nested CV for honest performance + final model for gene panels
# ============================================================================

library(glmnet)
library(pROC)
library(caret)

set.seed(42)

subtypes <- levels(y)
n_outer <- 10
n_inner <- 10

# Storage for outer-fold predictions and the final-model panels
results <- list()

for (st in subtypes) {
  cat("\n========================================\n")
  cat("Fitting one-vs-rest LASSO for:", st, "\n")
  cat("========================================\n")
  
  # Binary outcome: 1 = this subtype, 0 = other two
  y_bin <- as.integer(y == st)
  
  # Class weights to handle imbalance — inversely proportional to class freq
  w_pos <- sum(y_bin == 0) / length(y_bin)
  w_neg <- sum(y_bin == 1) / length(y_bin)
  obs_weights <- ifelse(y_bin == 1, w_pos, w_neg)
  
  # ----- Outer CV: stratified 10-fold -----
  outer_folds <- createFolds(factor(y_bin), k = n_outer, list = TRUE,
                             returnTrain = FALSE)
  
  outer_preds <- numeric(length(y_bin))  # held-out predicted probabilities
  outer_lambdas <- numeric(n_outer)
  outer_n_genes <- integer(n_outer)
  
  for (i in seq_along(outer_folds)) {
    test_idx <- outer_folds[[i]]
    train_idx <- setdiff(seq_along(y_bin), test_idx)
    
    X_tr <- X[train_idx, ]; y_tr <- y_bin[train_idx]; w_tr <- obs_weights[train_idx]
    X_te <- X[test_idx, ]
    
    # Inner CV picks lambda — stratified so rare class is in every inner fold
    inner_foldid <- as.integer(
      createFolds(factor(y_tr), k = n_inner, list = FALSE)
    )
    
    cv_fit <- cv.glmnet(
      x = X_tr, y = y_tr,
      family = "binomial",
      alpha = 1,                 # pure LASSO
      weights = w_tr,
      foldid = inner_foldid,
      type.measure = "auc",      # tune lambda on AUC, robust to imbalance
      standardize = TRUE
    )
    
    # Predict held-out probabilities at lambda.1se (sparser, more conservative)
    # Switch to lambda.min if you want maximum AUC at the cost of more genes
    outer_preds[test_idx] <- as.numeric(
      predict(cv_fit, newx = X_te, s = "lambda.1se", type = "response")
    )
    outer_lambdas[i] <- cv_fit$lambda.1se
    outer_n_genes[i] <- sum(coef(cv_fit, s = "lambda.1se")[-1] != 0)
    
    cat(sprintf("  Fold %2d: lambda=%.4f, genes=%d\n",
                i, outer_lambdas[i], outer_n_genes[i]))
  }
  
  # ----- Outer-CV performance -----
  roc_obj <- roc(y_bin, outer_preds, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  # Find threshold maximizing balanced accuracy
  thresholds <- seq(0.01, 0.99, by = 0.01)
  bal_acc <- sapply(thresholds, function(th) {
    pred_bin <- as.integer(outer_preds >= th)
    sens <- sum(pred_bin == 1 & y_bin == 1) / max(sum(y_bin == 1), 1)
    spec <- sum(pred_bin == 0 & y_bin == 0) / max(sum(y_bin == 0), 1)
    (sens + spec) / 2
  })
  best_th <- thresholds[which.max(bal_acc)]
  pred_bin <- as.integer(outer_preds >= best_th)
  sens <- sum(pred_bin == 1 & y_bin == 1) / sum(y_bin == 1)
  spec <- sum(pred_bin == 0 & y_bin == 0) / sum(y_bin == 0)
  
  
  cat(sprintf("  >>> Median genes per fold: %.0f (range %d-%d)\n",
              median(outer_n_genes), min(outer_n_genes), max(outer_n_genes)))
  
  # ----- Final model on ALL data for the reportable gene panel -----
  final_foldid <- as.integer(createFolds(factor(y_bin), k = n_inner, list = FALSE))
  final_fit <- cv.glmnet(
    x = X, y = y_bin,
    family = "binomial",
    alpha = 1,
    weights = obs_weights,
    foldid = final_foldid,
    type.measure = "auc",
    standardize = TRUE
  )
  
  coefs <- coef(final_fit, s = "lambda.1se")
  panel <- data.frame(
    gene = rownames(coefs)[-1],
    coef = as.numeric(coefs)[-1]
  ) |>
    filter(coef != 0) |>
    arrange(desc(abs(coef)))
  
  # Just re-run the whole STEP 4 block with this one-line fix applied:
  cat(sprintf("  >>> Median genes per fold: %.0f (range %d-%d)\n",
              median(outer_n_genes), min(outer_n_genes), max(outer_n_genes)))
  
  results[[st]] <- list(
    outer_auc = auc_val,
    outer_sens = sens,
    outer_spec = spec,
    outer_bal_acc = (sens + spec) / 2,
    outer_threshold = best_th,
    outer_preds = outer_preds,
    outer_n_genes_per_fold = outer_n_genes,
    final_lambda = final_fit$lambda.1se,
    panel = panel
  )
}

# ----- Summary table -----
summary_df <- data.frame(
  subtype = subtypes,
  AUC = sapply(results, \(r) r$outer_auc),
  sensitivity = sapply(results, \(r) r$outer_sens),
  specificity = sapply(results, \(r) r$outer_spec),
  bal_accuracy = sapply(results, \(r) r$outer_bal_acc),
  panel_size = sapply(results, \(r) nrow(r$panel)),
  median_fold_genes = sapply(results, \(r) median(r$outer_n_genes_per_fold))
)
cat("\n\n========================================\n")
cat("SUMMARY — Per-subtype LASSO classifiers\n")
cat("========================================\n")
print(summary_df, row.names = FALSE)


# Refit S-ALS-3 with lambda.min instead of lambda.1se
y_bin <- as.integer(y == "S-ALS-3")
w_pos <- sum(y_bin == 0) / length(y_bin)
w_neg <- sum(y_bin == 1) / length(y_bin)
obs_weights <- ifelse(y_bin == 1, w_pos, w_neg)

set.seed(42)
foldid <- as.integer(createFolds(factor(y_bin), k = 10, list = FALSE))
fit_min <- cv.glmnet(X, y_bin, family = "binomial", alpha = 1,
                     weights = obs_weights, foldid = foldid,
                     type.measure = "auc", standardize = TRUE)

# Compare panel sizes and CV AUCs
n_1se <- sum(coef(fit_min, s = "lambda.1se")[-1] != 0)
n_min <- sum(coef(fit_min, s = "lambda.min")[-1] != 0)
auc_1se <- max(fit_min$cvm[fit_min$lambda >= fit_min$lambda.1se])
auc_min <- max(fit_min$cvm)
cat(sprintf("lambda.1se: %d genes, inner-CV AUC %.3f\n", n_1se, auc_1se))
cat(sprintf("lambda.min: %d genes, inner-CV AUC %.3f\n", n_min, auc_min))

# Look at the full coefficient path — at what panel size does AUC plateau?
plot(fit_min)  # lambda vs CV AUC


# Identify the dominant gene(s) in S-ALS-3
# At lambda.1se (1 gene)
coef_1se <- coef(fit_min, s = "lambda.1se")
panel_1se <- data.frame(
  gene = rownames(coef_1se)[-1],
  coef = as.numeric(coef_1se)[-1]
) |> filter(coef != 0) |> arrange(desc(abs(coef)))
cat("S-ALS-3 panel at lambda.1se:\n")
print(panel_1se)

# At lambda.min (46 genes), top 10
coef_min <- coef(fit_min, s = "lambda.min")
panel_min <- data.frame(
  gene = rownames(coef_min)[-1],
  coef = as.numeric(coef_min)[-1]
) |> filter(coef != 0) |> arrange(desc(abs(coef)))
cat("\nS-ALS-3 top 10 genes at lambda.min:\n")
print(head(panel_min, 10))

# Check which Hallmark axis the top gene belongs to
top_gene <- panel_1se$gene[1]
hallmark |> 
  filter(gene_symbol == top_gene, gs_name %in% unlist(axis_sets)) |>
  select(gs_name, gene_symbol)


# Repeat sensitivity check for the other two subtypes
for (st in c("S-ALS-1", "S-ALS-2")) {
  y_bin <- as.integer(y == st)
  w_pos <- sum(y_bin == 0) / length(y_bin)
  w_neg <- sum(y_bin == 1) / length(y_bin)
  obs_weights <- ifelse(y_bin == 1, w_pos, w_neg)
  
  set.seed(42)
  foldid <- as.integer(createFolds(factor(y_bin), k = 10, list = FALSE))
  fit_st <- cv.glmnet(X, y_bin, family = "binomial", alpha = 1,
                      weights = obs_weights, foldid = foldid,
                      type.measure = "auc", standardize = TRUE)
  
  n_1se <- sum(coef(fit_st, s = "lambda.1se")[-1] != 0)
  n_min <- sum(coef(fit_st, s = "lambda.min")[-1] != 0)
  auc_1se <- fit_st$cvm[fit_st$lambda == fit_st$lambda.1se]
  auc_min <- max(fit_st$cvm)
  
  cat(sprintf("\n%s: lambda.1se = %d genes (AUC %.3f), lambda.min = %d genes (AUC %.3f)\n",
              st, n_1se, auc_1se, n_min, auc_min))
}



# ============================================================================
# STEP 5 — CCT2 follow-up: stability, expression, NfL link, panel composition
# ============================================================================

library(glmnet)
library(caret)
library(dplyr)
library(ggplot2)
library(tibble)

set.seed(42)

# ---------------------------------------------------------------------------
# 5A — BOOTSTRAP STABILITY SELECTION
# ---------------------------------------------------------------------------
# Question: how often does CCT2 land in the S-ALS-3 panel under resampling?
# Method: 200 bootstrap resamples, refit lambda.1se LASSO each time,
#         count selection frequency for every gene
# ---------------------------------------------------------------------------

n_boot <- 200
y_bin <- as.integer(y == "S-ALS-3")
w_pos <- sum(y_bin == 0) / length(y_bin)
w_neg <- sum(y_bin == 1) / length(y_bin)
obs_weights <- ifelse(y_bin == 1, w_pos, w_neg)

# Track which genes get selected, how often, with what coef sign
gene_selection <- matrix(0, nrow = n_boot, ncol = ncol(X),
                         dimnames = list(NULL, colnames(X)))

cat("Running", n_boot, "bootstraps for S-ALS-3 stability selection...\n")
for (b in seq_len(n_boot)) {
  # Stratified bootstrap: resample within each class to preserve balance
  pos_idx <- which(y_bin == 1)
  neg_idx <- which(y_bin == 0)
  boot_idx <- c(sample(pos_idx, length(pos_idx), replace = TRUE),
                sample(neg_idx, length(neg_idx), replace = TRUE))
  
  X_b <- X[boot_idx, ]
  y_b <- y_bin[boot_idx]
  w_b <- obs_weights[boot_idx]
  
  # Inner CV with stratified folds
  foldid_b <- as.integer(createFolds(factor(y_b), k = 10, list = FALSE))
  
  fit_b <- tryCatch(
    cv.glmnet(X_b, y_b, family = "binomial", alpha = 1,
              weights = w_b, foldid = foldid_b,
              type.measure = "auc", standardize = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit_b)) next
  
  coefs_b <- as.numeric(coef(fit_b, s = "lambda.1se"))[-1]
  gene_selection[b, ] <- as.integer(coefs_b != 0)
  
  if (b %% 25 == 0) cat("  bootstrap", b, "/", n_boot, "\n")
}

selection_freq <- colMeans(gene_selection)
stability_df <- data.frame(
  gene = names(selection_freq),
  selection_freq = selection_freq
) |> arrange(desc(selection_freq))

cat("\n=== Top 15 most stably selected genes for S-ALS-3 ===\n")
print(head(stability_df, 15))

cct2_freq <- stability_df$selection_freq[stability_df$gene == "CCT2"]
cat(sprintf("\nCCT2 selection frequency: %.1f%% across %d bootstraps\n",
            100 * cct2_freq, n_boot))

# ---------------------------------------------------------------------------
# 5B — CCT2 EXPRESSION BY SUBTYPE
# ---------------------------------------------------------------------------

cct2_expr <- data.frame(
  sample_id = rownames(X),
  CCT2 = X[, "CCT2"],
  subtype = y
)

# Quick descriptive stats
cct2_summary <- cct2_expr |>
  group_by(subtype) |>
  summarise(n = n(),
            mean = mean(CCT2),
            median = median(CCT2),
            sd = sd(CCT2),
            .groups = "drop")
cat("\n=== CCT2 expression by subtype (VST-normalized) ===\n")
print(cct2_summary)

# ANOVA + pairwise t-tests
aov_cct2 <- aov(CCT2 ~ subtype, data = cct2_expr)
cat("\nANOVA for CCT2 ~ subtype:\n")
print(summary(aov_cct2))

cat("\nPairwise Tukey HSD:\n")
print(TukeyHSD(aov_cct2))

# Plot
p_cct2 <- ggplot(cct2_expr, aes(x = subtype, y = CCT2, fill = subtype)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.08, alpha = 0.4, size = 0.8) +
  scale_fill_manual(values = c("S-ALS-1" = "#E41A1C",
                               "S-ALS-2" = "#377EB8",
                               "S-ALS-3" = "#4DAF4A")) +
  labs(title = "CCT2 expression by sporadic ALS subtype",
       subtitle = "Single-gene classifier for S-ALS-3 (AUC 0.918)",
       x = NULL, y = "CCT2 expression (VST)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
print(p_cct2)
ggsave("CCT2_by_subtype.pdf", p_cct2, width = 5, height = 4)

# ---------------------------------------------------------------------------
# 5C — DOES CCT2 EXPRESSION PREDICT NfL?
# ---------------------------------------------------------------------------
# Mirror the chain you established for AXIS_METABOLIC -> NfL -> survival
# Adjust the NfL column name and metadata object as needed
# ---------------------------------------------------------------------------

# Assumes ax_clean_dedup contains an NfL column (or you can join from elsewhere)
# Common names: NfL, plasma_NfL, NEFL, log_NfL — check yours
nfl_col_candidates <- c("NfL", "plasma_NfL", "log_NfL", "NEFL")
nfl_col <- intersect(nfl_col_candidates, colnames(ax_clean_dedup))

if (length(nfl_col) == 0) {
  cat("\n!! NfL column not found in ax_clean_dedup. Skipping 5C.\n")
  cat("   Add the join manually if NfL lives in another object.\n")
} else {
  nfl_col <- nfl_col[1]
  cat(sprintf("\nUsing NfL column: %s\n", nfl_col))
  
  cct2_nfl <- cct2_expr |>
    mutate(Participant_ID = sub("-[0-9]+-T$", "", sample_id)) |>
    inner_join(ax_clean_dedup |> select(Participant_ID, all_of(nfl_col)),
               by = "Participant_ID") |>
    rename(NfL = !!nfl_col) |>
    filter(!is.na(NfL))
  
  cat(sprintf("Patients with both CCT2 expression and NfL: %d\n", nrow(cct2_nfl)))
  
  # Log-transform NfL if it's not already (typical for NfL distributions)
  if (max(cct2_nfl$NfL, na.rm = TRUE) > 100) {
    cct2_nfl$NfL_used <- log10(cct2_nfl$NfL + 1)
    nfl_label <- "log10(NfL + 1)"
  } else {
    cct2_nfl$NfL_used <- cct2_nfl$NfL
    nfl_label <- nfl_col
  }
  
  # Univariate regression
  fit_nfl <- lm(NfL_used ~ CCT2, data = cct2_nfl)
  cat("\n=== CCT2 -> NfL univariate regression ===\n")
  print(summary(fit_nfl))
  
  # Adjusted for subtype — does CCT2 add information beyond subtype membership?
  fit_nfl_adj <- lm(NfL_used ~ CCT2 + subtype, data = cct2_nfl)
  cat("\n=== CCT2 -> NfL, adjusted for subtype ===\n")
  print(summary(fit_nfl_adj))
  
  # Plot
  p_nfl <- ggplot(cct2_nfl, aes(x = CCT2, y = NfL_used, color = subtype)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                aes(group = 1), linewidth = 0.7) +
    scale_color_manual(values = c("S-ALS-1" = "#E41A1C",
                                  "S-ALS-2" = "#377EB8",
                                  "S-ALS-3" = "#4DAF4A")) +
    labs(title = "CCT2 expression vs plasma NfL",
         x = "CCT2 expression (VST)", y = nfl_label) +
    theme_minimal(base_size = 12)
  print(p_nfl)
  ggsave("CCT2_vs_NfL.pdf", p_nfl, width = 6, height = 4)
}

# ---------------------------------------------------------------------------
# 5D — PANEL COMPOSITION BY HALLMARK AXIS (paper-ready table)
# ---------------------------------------------------------------------------
# For each subtype's final panel, what fraction comes from each axis?
# ---------------------------------------------------------------------------

# Build gene -> axis mapping (one gene can be in multiple axes; handle that)
gene_to_axis <- hallmark |>
  filter(gs_name %in% unlist(axis_sets)) |>
  mutate(axis = case_when(
    gs_name %in% axis_sets$AXIS_METABOLIC ~ "METABOLIC",
    gs_name %in% axis_sets$AXIS_INFLAMMATION ~ "INFLAMMATION",
    gs_name %in% axis_sets$AXIS_CELL_CYCLE ~ "CELL_CYCLE"
  )) |>
  select(gene_symbol, axis) |>
  distinct()

# For each subtype, summarize panel composition
panel_composition <- lapply(subtypes, function(st) {
  panel_genes <- results[[st]]$panel$gene
  axis_assignments <- gene_to_axis |> filter(gene_symbol %in% panel_genes)
  
  # A gene may be in multiple axes — count membership across all
  axis_counts <- axis_assignments |>
    count(axis, name = "n_genes") |>
    mutate(subtype = st,
           pct_of_panel = round(100 * n_genes / length(panel_genes), 1))
  
  axis_counts
}) |> bind_rows()

cat("\n=== Panel composition by Hallmark axis ===\n")
print(panel_composition)

# Wide format for the paper table
comp_wide <- panel_composition |>
  select(subtype, axis, n_genes) |>
  tidyr::pivot_wider(names_from = axis, values_from = n_genes, values_fill = 0) |>
  left_join(data.frame(
    subtype = subtypes,
    total_panel_size = sapply(subtypes, \(st) nrow(results[[st]]$panel))
  ), by = "subtype")

cat("\n=== Panel composition (wide format for paper) ===\n")
print(comp_wide)

# Save outputs
write.csv(stability_df, "S-ALS-3_stability_selection.csv", row.names = FALSE)
write.csv(panel_composition, "panel_axis_composition.csv", row.names = FALSE)
saveRDS(results, "lasso_results.rds")

