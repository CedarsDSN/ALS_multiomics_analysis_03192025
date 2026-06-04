# ============================================================
# Refactored entrypoint for original script: C9SOD1Convergence.R
# Use --config, --base-dir, --outdir, and --objects. See docs/SCRIPT_SPECS_FULL.md.
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
if (isTRUE(ALS_ARGS$help)) { als_print_help("02_c9_sod1_convergence"); quit(save = "no", status = 0) }
CFG <- als_load_config(ALS_ARGS$config)
ALSPROJ_BASE_DIR <- als_choose(ALS_ARGS$base_dir, CFG, "project.base_dir", getwd())
ALSPROJ_OUTDIR <- als_choose(ALS_ARGS$outdir, CFG, "outputs.root", file.path(getwd(), "results", "02_c9_sod1_convergence"))
ALSPROJ_SCRIPT_ID <- "02_c9_sod1_convergence"
als_init(ALSPROJ_SCRIPT_ID, ALSPROJ_BASE_DIR, ALSPROJ_OUTDIR, ALS_ARGS$redirect_outputs)
als_register_output_wrappers()
als_load_objects(ALS_ARGS$objects)
on.exit(als_write_session_info(ALSPROJ_SCRIPT_ID), add = TRUE)
# ---- End bootstrap ----


# ---- Libraries ----
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(matrixStats)
  library(NMF)              # for NMF + rank metrics
})

mainDir <- file.path(ALSPROJ_BASE_DIR)
setwd(mainDir)

counts <- read.csv(
  "transcriptomics/4_matrix/AnswerALS-651-T-v1-release6_raw-counts.csv",
  check.names = FALSE
)
stopifnot("Geneid" %in% colnames(counts))
counts <- column_to_rownames(counts, var = "Geneid")
# keep only transcriptomics samples
counts <- counts[, grepl("-T$", colnames(counts)), drop=FALSE]
colnames(counts) <- trimws(colnames(counts))

sample_map_tx <- read.csv("transcriptomics/4_matrix/Sample_Mapping_Information/Sample_Mapping_File_Feb_2024.csv")
colnames(sample_map_tx) <- gsub("\\.", "_", colnames(sample_map_tx))
sample_map_tx <- sample_map_tx %>%
  filter(grepl("-T$", Sample_ID)) %>%
  mutate(Sample_ID = trimws(Sample_ID),
         Participant_ID = trimws(Participant_ID))

repeat_df <- read_csv("metadata/clinical/aals_dataportal_datatable.csv") %>%
  mutate(Participant_ID = trimws(Participant_ID))

annot <- left_join(sample_map_tx, repeat_df, by = "Participant_ID")

# retain only samples with counts present
annot <- annot %>% filter(Sample_ID %in% colnames(counts))

mutation_subtypes_with_ALS_control <- read.csv(file.path(ALSPROJ_BASE_DIR, "ALS_Mutation_Subtypes_with_ALS_Control.csv"))
annot_als_test <- annot %>%
  left_join(mutation_subtypes_with_ALS_control, by = "Participant_ID") 
stopifnot(nrow(annot) > 0)
annot_als_filtered <- annot_als_test %>%
  filter(SUBJECT_GROUP.y %in% c("ALS", "Healthy Control"))
annot_als_filtered <- annot_als_filtered %>%
  filter(Mutation_Subtype %in% c("C9orf72", "SOD1", "Sporadic", "Control")) %>%
  droplevels()
# order counts columns to match ALS metadata
counts_als <- counts[, annot_als_filtered$Sample_ID, drop=FALSE]
counts_als <- counts_als[, match(annot_als_filtered$Sample_ID, colnames(counts_als)), drop=FALSE]
stopifnot(identical(colnames(counts_als), annot_als_filtered$Sample_ID))

# Common approach: keep genes with at least 10 counts in at least 10% of samples
min_count <- 10
min_samples <- round(0.1 * ncol(counts_als))

keep <- rowSums(counts_als >= min_count) >= min_samples
counts_filtered <- counts_als[keep, ]

# Check how many genes retained
nrow(counts_als)  # Before
nrow(counts_filtered)  # After

# Create DESeq2 object (minimal design for now)
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = annot_als_filtered,
  design = ~ 1
)

# VST transform
vst_data <- vst(dds, blind = TRUE)
expr_norm <- assay(vst_data)

pca_result <- prcomp(t(expr_norm), scale. = TRUE)
pca_meta <- data.frame(
  pca_result$x[, 1:10],
  annot_als_filtered
)

# Get variance explained
var_explained <- round(100 * summary(pca_result)$importance[2, 1:10], 1)
var_explained


# Variables to test
tech_vars <- c("Batch.x", "Site", "Shipment")
demo_vars <- c("SEX", "Age_at_First_PBMC_Collection", "RACE")
bio_vars <- c("Mutation_Subtype", "Site_of_Onset")
ipsc_vars <- c("PCT_SMI32", "PCT_ISL1", "PCT_TUJ1")

all_vars <- c(tech_vars, demo_vars, bio_vars, ipsc_vars)

# Test each variable against PC1-PC5
results <- data.frame()

for (pc in paste0("PC", 1:5)) {
  for (var in all_vars) {
    if (!var %in% colnames(pca_meta)) next
    
    tryCatch({
      if (is.numeric(pca_meta[[var]])) {
        test <- cor.test(pca_meta[[pc]], pca_meta[[var]], use = "complete.obs")
        results <- rbind(results, data.frame(
          PC = pc, Variable = var, Type = "continuous",
          P_value = test$p.value, Estimate = test$estimate
        ))
      } else {
        test <- summary(aov(pca_meta[[pc]] ~ pca_meta[[var]]))
        results <- rbind(results, data.frame(
          PC = pc, Variable = var, Type = "categorical",
          P_value = test[[1]]$`Pr(>F)`[1], Estimate = NA
        ))
      }
    }, error = function(e) NULL)
  }
}

# View significant associations
results %>% 
  filter(P_value < 0.05) %>%
  arrange(PC, P_value)


library(sva)

# Model with your biological variables of interest
mod <- model.matrix(~ Mutation_Subtype + SEX, data = annot_als_filtered)
mod0 <- model.matrix(~ SEX, data = annot_als_filtered)

# Find surrogate variables
svobj <- sva(expr_norm, mod, mod0)
svobj$n.sv

svobj <- sva(expr_norm, mod, mod0, n.sv = 10)
# Check correlation of SVs with Batch
sv_batch_cor <- sapply(1:10, function(i) {
  summary(aov(svobj$sv[,i] ~ annot_als_filtered$Batch.x))[[1]]$`Pr(>F)`[1]
})
sum(sv_batch_cor < 0.05) 

sv_mutation_cor <- sapply(1:10, function(i) {
  summary(aov(svobj$sv[,i] ~ annot_als_filtered$Mutation_Subtype))[[1]]$`Pr(>F)`[1]
})
sum(sv_mutation_cor < 0.05) 
sv_mutation_cor

annot_als_filtered$SV1 <- svobj$sv[, 1]
annot_als_filtered$SV2 <- svobj$sv[, 2]
annot_als_filtered$SV3 <- svobj$sv[, 3]
annot_als_filtered$SV4 <- svobj$sv[, 4]
annot_als_filtered$SV5 <- svobj$sv[, 5]
annot_als_filtered$SV6 <- svobj$sv[, 6]
annot_als_filtered$SV7 <- svobj$sv[, 7]
annot_als_filtered$SV8 <- svobj$sv[, 8]
annot_als_filtered$SV9 <- svobj$sv[, 9]
annot_als_filtered$SV10 <- svobj$sv[, 10]

annot_als_filtered$group <- annot_als_filtered$Mutation_Subtype

table(annot_als_filtered$group)

rownames(annot_als_filtered) <- annot_als_filtered$Sample_ID

all(colnames(counts_als) == rownames(annot_als_filtered))  


dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = annot_als_filtered,
  design = ~ SEX + SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + group
)

dds <- DESeq(dds)

pdf("dispersion_fit_6sva.pdf", width = 6, height = 6)
plotDispEsts(dds)
dev.off()

dds_SV10 <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = annot_als_filtered,
  design = ~ SEX + SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + group
)
# More stringent filtering
keep <- rowSums(counts(dds_SV10) >= 10) >= 10  
# at least 10 counts in at least 10 samples
dds_SV10_filtered <- dds_SV10[keep, ]
dds_SV10_filtered <- DESeq(dds_SV10_filtered)

dds_SV10 <- DESeq(dds_SV10)
dds_SV10 <- estimateDispersions(dds_SV10)
dds_SV10 <- nbinomWaldTest(dds_SV10, maxit = 1000)
res_C9_sv10 <- results(dds_SV10, contrast = c("group", "C9orf72", "Control"))
summary(res_C9_simple)

# SOD1 vs Control  
res_SOD1_sv10 <- results(dds_SV10, contrast = c("group", "SOD1", "Control"))
summary(res_SOD1_simple)

c9_de_sv10 <- sum(
  res_C9_sv10$padj < 0.05 &
    abs(res_C9_sv10$log2FoldChange) > 0.58,
  na.rm = TRUE
)

sod1_de_sv10 <- sum(
  res_SOD1_sv10$padj < 0.05 &
    abs(res_SOD1_sv10$log2FoldChange) > 0.58,
  na.rm = TRUE
)

# Check convergence more carefully
converged <- mcols(dds_SV10)$betaConv

# Look at what you're working with
head(names(converged))
head(rownames(res_C9_sv10))

# How many total failed to converge
sum(!converged, na.rm = TRUE)  # should be your 49

# For DEGs specifically - match by position instead
deg_idx <- which(res_C9_sv10$padj < 0.05 & !is.na(res_C9_sv10$padj))
sum(!converged[deg_idx], na.rm = TRUE)
pdf("dispersion_fit_10sva.pdf", width = 6, height = 6)
plotDispEsts(dds_SV10)
dev.off()

vst_pca <- vst(dds, blind = TRUE)
vst_mat <- assay(vst_pca)

# Run PCA
pca_result <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)

# Extract PC scores
pca_data <- data.frame(
  Sample_ID = colnames(vst_mat),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  PC4 = pca_result$x[, 4]
) %>%
  left_join(annot_als_filtered, by = "Sample_ID")

# Calculate variance explained
var_explained <- round(100 * summary(pca_result)$importance[2, 1:4], 1)

# ===== PLOT 1: PC1 vs PC2 (uncorrected) =====
p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = 2) +  # 95% confidence ellipses
  scale_color_manual(values = c(
    "C9orf72" = "#1e3a8a",    # Deep blue
    "SOD1" = "#dc2626",        # Red
    "Control" = "#6b7280",     # Gray
    "Sporadic" = "#f59e0b"     # Orange
  )) +
  scale_shape_manual(values = c(
    "C9orf72" = 16,
    "SOD1" = 17,
    "Control" = 15,
    "Sporadic" = 18
  )) +
  labs(
    title = "PCA: C9orf72 vs SOD1 vs Control vs Sporadic",
    subtitle = "VST-normalized",
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )
print(p1)
ggsave("pca_with_6_svs.png", p1, width = 10, height = 7, dpi = 300)

# C9 vs Control
res_C9 <- results(dds, contrast = c("group", "C9orf72", "Control"))
res_C9 <- as.data.frame(res_C9)
res_C9$gene <- rownames(res_C9)

# SOD1 vs Control
res_SOD1 <- results(dds, contrast = c("group", "SOD1", "Control"))
res_SOD1 <- as.data.frame(res_SOD1)
res_SOD1$gene <- rownames(res_SOD1)

# Check results
summary(res_C9)
summary(res_SOD1)

# How many significant genes in each?
sum(res_C9$padj < 0.05, na.rm = TRUE)
sum(res_SOD1$padj < 0.05, na.rm = TRUE)

# Create summary
de_summary <- data.frame(
  Comparison = c("C9 vs Control", "SOD1 vs Control"),
  N_samples = c("46 vs 283", "13 vs 283"),
  DE_genes = c(
    sum(res_C9$padj < 0.05, na.rm = TRUE),
    sum(res_SOD1$padj < 0.05, na.rm = TRUE)
  ),
  Median_absLog2FC = c(
    round(median(abs(res_C9$log2FoldChange[res_C9$padj < 0.05]), na.rm = TRUE), 2),
    round(median(abs(res_SOD1$log2FoldChange[res_SOD1$padj < 0.05]), na.rm = TRUE), 2)
  )
)

de_summary

# Merge fold changes for divergence test
divergence <- data.frame(
  gene = rownames(res_C9),
  log2FC_C9 = res_C9$log2FoldChange,
  padj_C9 = res_C9$padj,
  log2FC_SOD1 = res_SOD1$log2FoldChange,
  padj_SOD1 = res_SOD1$padj
)

# Remove NAs
divergence_clean <- divergence %>%
  filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1))

# Test correlation (negative = divergence)
cor.test(divergence_clean$log2FC_C9, divergence_clean$log2FC_SOD1)

# Option A
nrow(divergence %>% filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1)) %>% 
       filter(padj_C9 < 0.05 | padj_SOD1 < 0.05))

# Option B  
nrow(divergence %>% filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1)) %>% 
       filter(padj_C9 < 0.05 & padj_SOD1 < 0.05))


lfc_threshold <- 0.5
# Genes significant AND meet log2FC threshold in C9 OR SOD1
divergence_sig_lfc <- divergence %>%
  filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1)) %>%
  filter(
    (padj_C9 < 0.05 & abs(log2FC_C9) > lfc_threshold) | 
      (padj_SOD1 < 0.05 & abs(log2FC_SOD1) > lfc_threshold)
  )
nrow(divergence_sig_lfc)

# Test correlation
cor.test(divergence_sig_lfc$log2FC_C9, divergence_sig_lfc$log2FC_SOD1)

library(ggplot2)

pdf("C9_and_SOD1_divergence_transcriptomic_plot.pdf", width = 6, height = 6)
ggplot(divergence_clean, aes(x = log2FC_C9, y = log2FC_SOD1)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_smooth(method = "lm", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = paste0("C9 vs SOD1 Transcriptomic Signatures (r = ", 
                   round(0.40, 2), ")"),
    x = "log2FC (C9 vs Control)",
    y = "log2FC (SOD1 vs Control)"
  ) +
  theme_minimal()
dev.off()



# Let's run simple model and see what we get
# Create DESeq2 object with simple design
dds_simple <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = annot_als_filtered,
  design = ~ SEX + group
)

# Run DESeq2
dds_simple <- DESeq(dds_simple)

# Check dispersion fit
pdf("dispersion_fit_simple_model.pdf", width = 6, height = 6)
plotDispEsts(dds_simple)
dev.off()

# C9 vs Control
res_C9_simple <- results(dds_simple, contrast = c("group", "C9orf72", "Control"))
summary(res_C9_simple)

# SOD1 vs Control  
res_SOD1_simple <- results(dds_simple, contrast = c("group", "SOD1", "Control"))
summary(res_SOD1_simple)

# Count DE genes
sum(res_C9_simple$padj < 0.05, na.rm = TRUE)
sum(res_SOD1_simple$padj < 0.05, na.rm = TRUE)

c9_de_simple <- sum(
  res_C9_simple$padj < 0.05 &
    abs(res_C9_simple$log2FoldChange) > 0.58,
  na.rm = TRUE
)

sod1_de_simple <- sum(
  res_SOD1_simple$padj < 0.05 &
    abs(res_SOD1_simple$log2FoldChange) > 0.58,
  na.rm = TRUE
)
divergence <- data.frame(
  gene = rownames(res_C9_simple),
  log2FC_C9 = res_C9_simple$log2FoldChange,
  padj_C9 = res_C9_simple$padj,
  log2FC_SOD1 = res_SOD1_simple$log2FoldChange,
  padj_SOD1 = res_SOD1_simple$padj
)

# Option A
nrow(divergence %>% filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1)) %>% 
       filter(padj_C9 < 0.05 | padj_SOD1 < 0.05))

# Option B  
nrow(divergence %>% filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1)) %>% 
       filter(padj_C9 < 0.05 & padj_SOD1 < 0.05))

# Remove NAs
divergence_clean <- divergence %>%
  filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1))

# Test correlation (negative = divergence)
cor.test(divergence_clean$log2FC_C9, divergence_clean$log2FC_SOD1)

divergence_sig <- divergence %>%
  filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1)) %>%
  filter(padj_C9 < 0.05 | padj_SOD1 < 0.05)

# Test correlation
cor.test(divergence_sig$log2FC_C9, divergence_sig$log2FC_SOD1)


lfc_threshold <- 0.5
# Genes significant AND meet log2FC threshold in C9 OR SOD1
divergence_sig_lfc <- divergence %>%
  filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1)) %>%
  filter(
    (padj_C9 < 0.05 & abs(log2FC_C9) > lfc_threshold) | 
      (padj_SOD1 < 0.05 & abs(log2FC_SOD1) > lfc_threshold)
  )
nrow(divergence_sig_lfc)

# Test correlation
cor.test(divergence_sig_lfc$log2FC_C9, divergence_sig_lfc$log2FC_SOD1)

# Using res_C9_simple and res_SOD1_simple
divergence_simple <- data.frame(
  gene = rownames(res_C9_simple),
  log2FC_C9 = res_C9_simple$log2FoldChange,
  padj_C9 = res_C9_simple$padj,
  log2FC_SOD1 = res_SOD1_simple$log2FoldChange,
  padj_SOD1 = res_SOD1_simple$padj
)

divergence_sig_lfc_simple <- divergence_simple %>%
  filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1)) %>%
  filter(
    (padj_C9 < 0.05 & abs(log2FC_C9) > 0.5) | 
      (padj_SOD1 < 0.05 & abs(log2FC_SOD1) > 0.5)
  )

nrow(divergence_sig_lfc_simple)
cor.test(divergence_sig_lfc_simple$log2FC_C9, divergence_sig_lfc_simple$log2FC_SOD1)

# Using res_C9 and res_SOD1 (from SVA model)
divergence_sva <- data.frame(
  gene = rownames(res_C9),
  log2FC_C9 = res_C9$log2FoldChange,
  padj_C9 = res_C9$padj,
  log2FC_SOD1 = res_SOD1$log2FoldChange,
  padj_SOD1 = res_SOD1$padj
)

divergence_sig_lfc_sva <- divergence_sva %>%
  filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1)) %>%
  filter(
    (padj_C9 < 0.05 & abs(log2FC_C9) > 0.5) | 
      (padj_SOD1 < 0.05 & abs(log2FC_SOD1) > 0.5)
  )

nrow(divergence_sig_lfc_sva)
cor.test(divergence_sig_lfc_sva$log2FC_C9, divergence_sig_lfc_sva$log2FC_SOD1)

# Leek method (default, aggressive)
n_leek <- num.sv(expr_norm, mod, method = "leek")
n_leek

# Buja-Eyuboglu method (more conservative)
n_be <- num.sv(expr_norm, mod, method = "be")
n_be

# Then use the estimate
svobj <- sva(expr_norm, mod, mod0, n.sv = n_be)


library(dplyr)
library(sva)
library(DESeq2)

# ----------------------------
# Helpers (remove redundancy)
# ----------------------------

pval_aov <- function(x, f) summary(aov(x ~ f))[[1]]$`Pr(>F)`[1]

sv_assoc_counts <- function(sv_mat, annot, alpha = 0.05) {
  n <- ncol(sv_mat)
  batch_p <- sapply(seq_len(n), function(i) pval_aov(sv_mat[, i], annot$Batch.x))
  group_p <- sapply(seq_len(n), function(i) pval_aov(sv_mat[, i], annot$group))
  list(
    batch_p = batch_p,
    group_p = group_p,
    n_batch_corr = sum(batch_p < alpha, na.rm = TRUE),
    n_group_corr = sum(group_p < alpha, na.rm = TRUE)
  )
}

add_svs <- function(annot, sv_mat, idx, prefix = "SV") {
  out <- annot
  if (length(idx) == 0) return(out)
  for (j in seq_along(idx)) out[[paste0(prefix, j)]] <- sv_mat[, idx[j]]
  out
}

build_design <- function(n_sv, mode = c("clean", "all")) {
  mode <- match.arg(mode)
  if (n_sv <= 0) return(~ SEX + group)
  sv_terms <- paste0("SV", seq_len(n_sv), collapse = " + ")
  as.formula(paste("~ SEX +", sv_terms, "+ group"))
}

run_deseq_metrics <- function(counts_filtered, annot_with_svs, design_formula,
                              lfc_thr = 0.5, alpha = 0.05, min_sig_n = 10) {
  dds <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData   = annot_with_svs,
    design    = design_formula
  )
  dds <- DESeq(dds, quiet = TRUE)
  
  res_c9   <- results(dds, contrast = c("group", "C9orf72", "Control"))
  res_sod1 <- results(dds, contrast = c("group", "SOD1", "Control"))
  
  c9_de   <- sum(res_c9$padj < alpha, na.rm = TRUE)
  sod1_de <- sum(res_sod1$padj < alpha, na.rm = TRUE)
  
  # correlation (all genes)
  ok_all <- !is.na(res_c9$log2FoldChange) & !is.na(res_sod1$log2FoldChange)
  r_all <- if (sum(ok_all) > 3) cor.test(res_c9$log2FoldChange[ok_all],
                                         res_sod1$log2FoldChange[ok_all])$estimate else NA_real_
  
  # correlation (sig genes + LFC threshold in either)
  ok_sig <- ok_all &
    ((res_c9$padj < alpha & abs(res_c9$log2FoldChange) > lfc_thr) |
       (res_sod1$padj < alpha & abs(res_sod1$log2FoldChange) > lfc_thr))
  
  n_sig <- sum(ok_sig, na.rm = TRUE)
  r_sig <- if (n_sig > min_sig_n) cor.test(res_c9$log2FoldChange[ok_sig],
                                           res_sod1$log2FoldChange[ok_sig])$estimate else NA_real_
  
  list(
    c9_de = c9_de,
    sod1_de = sod1_de,
    r_all = unname(r_all),
    r_sig = unname(r_sig),
    n_sig = n_sig
  )
}

# ----------------------------
# Main loop (one version)
# Choose mode = "clean" (filter group-associated SVs) OR "all" (include all SVs)
# ----------------------------

test_nsv <- c(2, 5, 10, 15, 20)
mode <- "clean"   # "clean" or "all"
alpha <- 0.05
lfc_thr <- 0.5

sensitivity_results <- bind_rows(lapply(test_nsv, function(n) {
  cat("\n=== Testing n.sv =", n, "===\n")
  
  svobj <- sva(expr_norm, mod, mod0, n.sv = n)
  assoc <- sv_assoc_counts(svobj$sv, annot_als_filtered, alpha = alpha)
  
  if (mode == "clean") {
    keep_idx <- which(assoc$group_p >= alpha)  # keep SVs not associated with group
  } else {
    keep_idx <- seq_len(n)                    # keep all SVs
  }
  
  temp_annot <- add_svs(annot_als_filtered, svobj$sv, keep_idx)
  design_formula <- build_design(length(keep_idx), mode = mode)
  
  cat("  Batch correlations:", assoc$n_batch_corr, "/", n, "\n")
  cat("  Group correlations:", assoc$n_group_corr, "/", n, "\n")
  cat("  SVs used in model:", length(keep_idx), "\n")
  cat("  Design:", deparse(design_formula), "\n")
  
  # If you filtered everything out in clean mode, fall back to ~ SEX + group
  metrics <- run_deseq_metrics(
    counts_filtered = counts_filtered,
    annot_with_svs  = temp_annot,
    design_formula  = design_formula,
    lfc_thr = lfc_thr,
    alpha   = alpha
  )
  
  data.frame(
    n_sv = n,
    mode = mode,
    n_batch_corr = assoc$n_batch_corr,
    n_group_corr = assoc$n_group_corr,
    n_used = length(keep_idx),
    c9_de_genes = metrics$c9_de,
    sod1_de_genes = metrics$sod1_de,
    cor_all_genes = round(metrics$r_all, 3),
    cor_sig_genes = round(metrics$r_sig, 3),
    n_sig_genes = metrics$n_sig
  )
}))

print(sensitivity_results)

# Save results
write.csv(sensitivity_results, "sva_sensitivity_all_svs.csv", row.names = FALSE)


# PCA check

library(ggplot2)
library(ggrepel)

# ===== PREPARE DATA =====

# Select samples: all C9, all SOD1, all controls, + random sporadic
set.seed(123)  # For reproducibility

# Get sample IDs for each group
c9_samples <- annot_als_filtered %>% filter(group == "C9orf72") %>% pull(Sample_ID)
sod1_samples <- annot_als_filtered %>% filter(group == "SOD1") %>% pull(Sample_ID)
control_samples <- annot_als_filtered %>% filter(group == "Control") %>% pull(Sample_ID)
sporadic_samples <- annot_als_filtered %>% filter(group == "Sporadic") %>% pull(Sample_ID)

# Sample ~50 random sporadic (or all if you want)
n_sporadic <- min(50, length(sporadic_samples))
sporadic_subset <- sample(sporadic_samples, n_sporadic)

# Combine all samples
pca_samples <- c(c9_samples, sod1_samples, control_samples, sporadic_subset)

# Subset the data
pca_annot <- annot_als_filtered %>% filter(Sample_ID %in% pca_samples)
pca_counts <- counts_filtered[, pca_samples]
rownames(pca_annot) <- pca_annot$Sample_ID
pca_annot <- pca_annot[colnames(pca_counts), ]
cat("Sample counts:\n")
cat("  C9orf72:", length(c9_samples), "\n")
cat("  SOD1:", length(sod1_samples), "\n")
cat("  Control:", length(control_samples), "\n")
cat("  Sporadic:", n_sporadic, "\n")
cat("  Total:", length(pca_samples), "\n")

# ===== OPTION 1: PCA on VST-normalized (no batch correction) =====
library(DESeq2)

dds_pca <- DESeqDataSetFromMatrix(
  countData = pca_counts,
  colData = pca_annot,
  design = ~ group
)

# VST transformation
vst_pca <- vst(dds_pca, blind = TRUE)
vst_mat <- assay(vst_pca)

# Run PCA
pca_result <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)

# Extract PC scores
pca_data <- data.frame(
  Sample_ID = colnames(vst_mat),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  PC4 = pca_result$x[, 4]
) %>%
  left_join(pca_annot, by = "Sample_ID")

# Calculate variance explained
var_explained <- round(100 * summary(pca_result)$importance[2, 1:4], 1)

# ===== PLOT 1: PC1 vs PC2 (uncorrected) =====
p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = 2) +  # 95% confidence ellipses
  scale_color_manual(values = c(
    "C9orf72" = "#1e3a8a",    # Deep blue
    "SOD1" = "#dc2626",        # Red
    "Control" = "#6b7280",     # Gray
    "Sporadic" = "#f59e0b"     # Orange
  )) +
  scale_shape_manual(values = c(
    "C9orf72" = 16,
    "SOD1" = 17,
    "Control" = 15,
    "Sporadic" = 18
  )) +
  labs(
    title = "PCA: C9orf72 vs SOD1 vs Control vs Sporadic",
    subtitle = "VST-normalized, no batch correction",
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )

print(p1)
ggsave("pca_no_batch_correction.png", p1, width = 10, height = 7, dpi = 300)

# ===== PLOT 2: PC1 vs PC2 with batch labels =====
p2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(aes(shape = group), size = 3, alpha = 0.7) +
  scale_color_manual(values = c(
    "C9orf72" = "#1e3a8a",
    "SOD1" = "#dc2626",
    "Control" = "#6b7280",
    "Sporadic" = "#f59e0b"
  )) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(
    title = "PCA colored by mutation, shaped by batch structure",
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)")
  ) +
  theme_bw(base_size = 12)

print(p2)

# ===== PLOT 3: PC2 vs PC3 (sometimes separation is clearer here) =====
p3 <- ggplot(pca_data, aes(x = PC2, y = PC3, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = 2) +
  scale_color_manual(values = c(
    "C9orf72" = "#1e3a8a",
    "SOD1" = "#dc2626",
    "Control" = "#6b7280",
    "Sporadic" = "#f59e0b"
  )) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(
    title = "PCA: PC2 vs PC3",
    subtitle = "Sometimes separation is clearer in later PCs",
    x = paste0("PC2 (", var_explained[2], "% variance)"),
    y = paste0("PC3 (", var_explained[3], "% variance)")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", legend.title = element_blank())

print(p3)
ggsave("pca_pc2_pc3.png", p3, width = 10, height = 7, dpi = 300)

# Plot with just C9 and SOD1 (remove controls and sporadic for clarity)
p_familial <- ggplot(pca_data %>% filter(group %in% c("C9orf72", "SOD1")), 
                     aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2, size = 1) +
  scale_color_manual(values = c("C9orf72" = "#1e3a8a", "SOD1" = "#dc2626")) +
  labs(
    title = "PCA: C9orf72 vs SOD1 Only",
    subtitle = paste0("Distance between centroids: ", round(dist_c9_sod1, 2)),
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)")
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right", legend.title = element_blank())

print(p_familial)
ggsave("pca_c9_vs_sod1_only.png", p_familial, width = 10, height = 7, dpi = 300)

p_familial <- ggplot(pca_data %>% filter(group %in% c("C9orf72", "Control")),
                     aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2, size = 1) +
  scale_color_manual(values = c("C9orf72" = "#1e3a8a", "Control" = "#dc2626")) +
  labs(
    title = "PCA: C9orf72 vs Control Only",
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)")
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right", legend.title = element_blank())

print(p_familial)
ggsave("pca_c9_vs_control_only.png", p_familial, width = 10, height = 7, dpi = 300)

# PC1 vs PC3 (sometimes clearer)
p_pc1_pc3 <- ggplot(pca_data, aes(x = PC1, y = PC3, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = 2) +
  scale_color_manual(values = c(
    "C9orf72" = "#1e3a8a",
    "SOD1" = "#dc2626",
    "Control" = "#6b7280",
    "Sporadic" = "#f59e0b"
  )) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(
    title = "PCA: PC1 vs PC3",
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC3 (", var_explained[3], "% variance)")
  ) +
  theme_bw(base_size = 12)

print(p_pc1_pc3)
ggsave("pca_pc1_pc3.png", p_pc1_pc3, width = 10, height = 7, dpi = 300)

# ===== CALCULATE DISTANCES BETWEEN GROUP CENTROIDS =====
cat("\n=== Group Centroid Distances (PC1-PC2 space) ===\n")

centroids <- pca_data %>%
  group_by(group) %>%
  summarise(
    PC1_mean = mean(PC1),
    PC2_mean = mean(PC2),
    n = n()
  )

# Calculate Euclidean distances between centroids
dist_c9_sod1 <- sqrt((centroids$PC1_mean[centroids$group == "C9orf72"] - 
                        centroids$PC1_mean[centroids$group == "SOD1"])^2 +
                       (centroids$PC2_mean[centroids$group == "C9orf72"] - 
                          centroids$PC2_mean[centroids$group == "SOD1"])^2)

dist_c9_control <- sqrt((centroids$PC1_mean[centroids$group == "C9orf72"] - 
                           centroids$PC1_mean[centroids$group == "Control"])^2 +
                          (centroids$PC2_mean[centroids$group == "C9orf72"] - 
                             centroids$PC2_mean[centroids$group == "Control"])^2)

dist_sod1_control <- sqrt((centroids$PC1_mean[centroids$group == "SOD1"] - 
                             centroids$PC1_mean[centroids$group == "Control"])^2 +
                            (centroids$PC2_mean[centroids$group == "SOD1"] - 
                               centroids$PC2_mean[centroids$group == "Control"])^2)

cat("C9orf72 <-> SOD1:", round(dist_c9_sod1, 2), "\n")
cat("C9orf72 <-> Control:", round(dist_c9_control, 2), "\n")
cat("SOD1 <-> Control:", round(dist_sod1_control, 2), "\n")

# Ratio: if >1, C9 and SOD1 are farther apart than either is from control
ratio <- dist_c9_sod1 / mean(c(dist_c9_control, dist_sod1_control))
cat("\nC9-SOD1 distance / avg(disease-control distance):", round(ratio, 2), "\n")
if (ratio < 0.5) {
  cat("→ C9 and SOD1 are VERY CLOSE (supports convergence)\n")
} else if (ratio < 1) {
  cat("→ C9 and SOD1 are closer to each other than to controls (moderate convergence)\n")
} else {
  cat("→ C9 and SOD1 are farther apart than from controls (supports divergence)\n")
}

# ===== OPTION 2: PCA on BATCH-CORRECTED DATA (using limma removeBatchEffect) =====
cat("\n=== Creating batch-corrected PCA ===\n")

library(limma)

# Use SVA surrogate variables for correction
# First run SVA on this subset
mod_pca <- model.matrix(~ group, data = pca_annot)
mod0_pca <- model.matrix(~ 1, data = pca_annot)

# VST normalize first
vst_for_sva <- vst(dds_pca, blind = TRUE)
expr_for_sva <- assay(vst_for_sva)

# Run SVA
library(sva)
svobj_pca <- sva(expr_for_sva, mod_pca, mod0_pca, n.sv = 10)

# Remove batch effects using SVs
expr_corrected <- removeBatchEffect(
  expr_for_sva,
  covariates = svobj_pca$sv
)

# Run PCA on corrected data
pca_corrected <- prcomp(t(expr_corrected), center = TRUE, scale. = FALSE)

# Extract PC scores
pca_data_corrected <- data.frame(
  Sample_ID = colnames(expr_corrected),
  PC1 = pca_corrected$x[, 1],
  PC2 = pca_corrected$x[, 2],
  PC3 = pca_corrected$x[, 3],
  PC4 = pca_corrected$x[, 4]
) %>%
  left_join(pca_annot, by = "Sample_ID")

var_explained_corr <- round(100 * summary(pca_corrected)$importance[2, 1:4], 1)

# ===== PLOT 4: Batch-corrected PCA =====
p4 <- ggplot(pca_data_corrected, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = 2) +
  scale_color_manual(values = c(
    "C9orf72" = "#1e3a8a",
    "SOD1" = "#dc2626",
    "Control" = "#6b7280",
    "Sporadic" = "#f59e0b"
  )) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(
    title = "PCA: C9orf72 vs SOD1 vs Control vs Sporadic",
    subtitle = "Batch-corrected using SVA (10 SVs)",
    x = paste0("PC1 (", var_explained_corr[1], "% variance)"),
    y = paste0("PC2 (", var_explained_corr[2], "% variance)")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", legend.title = element_blank())

print(p4)
ggsave("pca_batch_corrected.png", p4, width = 10, height = 7, dpi = 300)

# ===== PERMANOVA TEST: Statistical test for group separation =====
library(vegan)

# Test on corrected data
dist_matrix <- dist(pca_corrected$x[, 1:10])  # Use first 10 PCs
permanova_result <- adonis2(dist_matrix ~ group, data = pca_data_corrected, permutations = 999)

cat("\n=== PERMANOVA: Testing for group separation ===\n")
print(permanova_result)

# If p < 0.05, groups are significantly different overall
# R² tells you how much variance is explained by group membership


library(sva)
library(DESeq2)
library(dplyr)

# ===== STEP 1: RUN SVA WITH n=10 =====
cat("Running SVA with n.sv = 10...\n")
mod <- model.matrix(~ Mutation_Subtype + SEX, data = annot_als_filtered)
mod0 <- model.matrix(~ SEX, data = annot_als_filtered)
svobj <- sva(expr_norm, mod, mod0, n.sv = 10)

# ===== STEP 2: IDENTIFY CONFOUNDED SVs =====
cat("\nChecking which SVs correlate with mutation subtype...\n")

annot_als_filtered$group <- annot_als_filtered$Mutation_Subtype

sv_mutation <- sapply(1:10, function(i) {
  summary(aov(svobj$sv[,i] ~ annot_als_filtered$group))[[1]]$`Pr(>F)`[1]
})

sv_batch <- sapply(1:10, function(i) {
  summary(aov(svobj$sv[,i] ~ annot_als_filtered$Batch.x))[[1]]$`Pr(>F)`[1]
})

# Create summary table
sv_summary <- data.frame(
  SV = paste0("SV", 1:10),
  Batch_pval = round(sv_batch, 4),
  Mutation_pval = round(sv_mutation, 4),
  Batch_sig = ifelse(sv_batch < 0.05, "YES", "no"),
  Mutation_sig = ifelse(sv_mutation < 0.05, "YES", "no")
)

print(sv_summary)

# Identify clean SVs (NOT correlated with mutation)
confounded_svs <- which(sv_mutation < 0.05)
clean_svs <- which(sv_mutation >= 0.05)

cat("\nConfounded SVs (correlated with mutation):", confounded_svs, "\n")
cat("Clean SVs (not correlated with mutation):", clean_svs, "\n")
cat("Number of clean SVs to use:", length(clean_svs), "\n")

# ===== STEP 3: ADD ONLY CLEAN SVs TO ANNOTATION =====
cat("\nAdding clean SVs to annotation...\n")

annot_sva_filtered <- annot_als_filtered

# Add only the clean SVs, but rename them sequentially as SV1, SV2, etc.
for (i in seq_along(clean_svs)) {
  annot_sva_filtered[[paste0("SV", i)]] <- svobj$sv[, clean_svs[i]]
}

cat("Added", length(clean_svs), "clean SVs to the model\n")

# ===== STEP 4: BUILD MODEL WITH FILTERED SVs =====
# Build design formula
if (length(clean_svs) > 0) {
  sv_terms <- paste0("SV", 1:length(clean_svs), collapse = " + ")
  design_filtered <- as.formula(paste("~ SEX +", sv_terms, "+ group"))
} else {
  design_filtered <- ~ SEX + group
  cat("WARNING: No clean SVs found, using simple model\n")
}

cat("Design formula:", deparse(design_filtered), "\n")

# ===== STEP 5: RUN DESeq2 =====
cat("\nRunning DESeq2...\n")

dds_filtered <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = annot_sva_filtered,
  design = design_filtered
)

dds_filtered <- DESeq(dds_filtered)
dds_filtered <- nbinomWaldTest(dds_filtered, maxit = 500)
sum(!mcols(dds_filtered)$betaConv)

# Get results
res_C9_filtered <- results(dds_filtered, contrast = c("group", "C9orf72", "Control"))
res_SOD1_filtered <- results(dds_filtered, contrast = c("group", "SOD1", "Control"))
res_Sporadic_filtered <- results(dds_filtered, contrast = c("group", "Sporadic", "Control"))

cat("DESeq2 complete!\n")

# ===== STEP 6: CALCULATE DETAILED METRICS =====
cat("\n========================================\n")
cat("DETAILED METRICS: SVA FILTERED\n")
cat(paste0("(n.sv=10 total, ", length(clean_svs), " clean SVs used)\n"))
cat("========================================\n\n")

# Count DE genes
c9_de <- sum(
  res_C9_filtered$padj < 0.05 &
    abs(res_C9_filtered$log2FoldChange) > 0.58,
  na.rm = TRUE
)
sod1_de <- sum(
  res_SOD1_filtered$padj < 0.05 &
    abs(res_SOD1_filtered$log2FoldChange) > 0.58,
  na.rm = TRUE
)

cat("C9 DEGs (padj < 0.05):", c9_de, "\n")
cat("SOD1 DEGs (padj < 0.05):", sod1_de, "\n\n")

# Get significant genes
c9_sig_genes <- rownames(res_C9_filtered)[which(res_C9_filtered$padj < 0.05)]
sod1_sig_genes <- rownames(res_SOD1_filtered)[which(res_SOD1_filtered$padj < 0.05)]

# Overlap (both significant)
overlap_genes <- intersect(c9_sig_genes, sod1_sig_genes)
cat("Overlap (both sig):", length(overlap_genes), "genes\n")

# Same direction?
if (length(overlap_genes) > 0) {
  overlap_data <- data.frame(
    gene = overlap_genes,
    c9_lfc = res_C9_filtered[overlap_genes, "log2FoldChange"],
    sod1_lfc = res_SOD1_filtered[overlap_genes, "log2FoldChange"]
  )
  same_dir <- sum(sign(overlap_data$c9_lfc) == sign(overlap_data$sod1_lfc))
  cat("  Same direction:", same_dir, "(", round(100*same_dir/length(overlap_genes), 1), "%)\n\n")
}

# Create divergence dataframe
divergence <- data.frame(
  gene = rownames(res_C9_filtered),
  log2FC_C9 = res_C9_filtered$log2FoldChange,
  padj_C9 = res_C9_filtered$padj,
  log2FC_SOD1 = res_SOD1_filtered$log2FoldChange,
  padj_SOD1 = res_SOD1_filtered$padj
)

# Significant in EITHER
div_either <- divergence %>%
  filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1)) %>%
  filter(padj_C9 < 0.05 | padj_SOD1 < 0.05)

cat("Sig in either:", nrow(div_either), "genes\n\n")

# ===== CORRELATIONS =====
cat("--- Correlation Analysis ---\n")

# 1. All genes
div_all <- divergence %>%
  filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1))

cor_all <- cor.test(div_all$log2FC_C9, div_all$log2FC_SOD1)
cat("All genes: r =", round(cor_all$estimate, 3), 
    ", p =", format(cor_all$p.value, scientific = TRUE), 
    ", n =", nrow(div_all), "\n")

# 2. Sig in either (no LFC threshold)
if (nrow(div_either) > 10) {
  cor_either <- cor.test(div_either$log2FC_C9, div_either$log2FC_SOD1)
  cat("Sig in either: r =", round(cor_either$estimate, 3), 
      ", p =", format(cor_either$p.value, scientific = TRUE), 
      ", n =", nrow(div_either), "\n")
} else {
  cor_either <- list(estimate = NA)
}

# 3. Sig with |LFC| > 0.5
div_sig_lfc <- divergence %>%
  filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1)) %>%
  filter(
    (padj_C9 < 0.05 & abs(log2FC_C9) > 0.5) | 
      (padj_SOD1 < 0.05 & abs(log2FC_SOD1) > 0.5)
  )

if (nrow(div_sig_lfc) > 10) {
  cor_sig_lfc <- cor.test(div_sig_lfc$log2FC_C9, div_sig_lfc$log2FC_SOD1)
  cat("Sig + |LFC| > 0.5: r =", round(cor_sig_lfc$estimate, 3), 
      ", p =", format(cor_sig_lfc$p.value, scientific = TRUE), 
      ", n =", nrow(div_sig_lfc), "\n")
} else {
  cor_sig_lfc <- list(estimate = NA)
}

# ===== CREATE SUMMARY TABLE =====
sva_filtered_summary <- data.frame(
  Model = paste0("SVA filtered (", length(clean_svs), " clean SVs)"),
  n_total_sv = 10,
  n_clean_sv = length(clean_svs),
  C9_DEGs = c9_de,
  SOD1_DEGs = sod1_de,
  Overlap_both_sig = length(overlap_genes),
  Sig_in_either = nrow(div_either),
  Cor_all_genes = round(cor_all$estimate, 3),
  Cor_sig_either = round(cor_either$estimate, 3),
  Cor_sig_LFC05 = round(cor_sig_lfc$estimate, 3),
  n_genes_LFC05 = nrow(div_sig_lfc)
)

cat("\n========================================\n")
cat("SUMMARY TABLE\n")
cat("========================================\n")
print(sva_filtered_summary)

# ===== COMPARE TO SIMPLE MODEL =====
cat("\n========================================\n")
cat("COMPARISON: SIMPLE vs SVA FILTERED\n")
cat("========================================\n\n")

comparison <- data.frame(
  Model = c("Simple (SEX only)", paste0("SVA filtered (", length(clean_svs), " clean SVs)")),
  C9_DEGs = c(333, c9_de),
  SOD1_DEGs = c(28, sod1_de),
  Overlap = c(7, length(overlap_genes)),
  Sig_either = c(354, nrow(div_either)),
  Cor_all = c(0.258, round(cor_all$estimate, 3)),
  Cor_sig_either = c(0.362, round(cor_either$estimate, 3)),
  Cor_LFC05 = c(0.378, round(cor_sig_lfc$estimate, 3)),
  n_LFC05 = c(179, nrow(div_sig_lfc))
)

print(comparison)

# Calculate improvement
cat("\n--- Improvement over Simple Model ---\n")
cat("C9 DEGs:", c9_de, "vs 333 (", round(100*(c9_de-333)/333, 1), "% change)\n")
cat("SOD1 DEGs:", sod1_de, "vs 28 (", round(100*(sod1_de-28)/28, 1), "% change)\n")
cat("Correlation (LFC>0.5):", round(cor_sig_lfc$estimate, 3), "vs 0.378 (", 
    round(100*(cor_sig_lfc$estimate - 0.378)/0.378, 1), "% improvement)\n")

# Save results
write.csv(sva_filtered_summary, "sva_filtered_summary.csv", row.names = FALSE)
write.csv(comparison, "model_comparison_simple_vs_filtered.csv", row.names = FALSE)

cat("\n✓ Results saved!\n")

# ===== GET DETAILED METRICS FOR SVA ALL (n=10) =====

cat("\n========================================\n")
cat("DETAILED METRICS: SVA ALL (10 SVs)\n")
cat("========================================\n\n")

# Extract the saved results from sensitivity analysis
res_C9_all <- results_storage[["n10"]]$res_c9
res_SOD1_all <- results_storage[["n10"]]$res_sod1

# Count DE genes (we already have these)
c9_de_all <- sum(res_C9_all$padj < 0.05, na.rm = TRUE)
sod1_de_all <- sum(res_SOD1_all$padj < 0.05, na.rm = TRUE)

cat("C9 DEGs (padj < 0.05):", c9_de_all, "\n")
cat("SOD1 DEGs (padj < 0.05):", sod1_de_all, "\n\n")

# Get significant genes
c9_sig_genes_all <- rownames(res_C9_all)[which(res_C9_all$padj < 0.05)]
sod1_sig_genes_all <- rownames(res_SOD1_all)[which(res_SOD1_all$padj < 0.05)]

# Overlap (both significant)
overlap_genes_all <- intersect(c9_sig_genes_all, sod1_sig_genes_all)
cat("Overlap (both sig):", length(overlap_genes_all), "genes\n")

# Same direction?
if (length(overlap_genes_all) > 0) {
  overlap_data_all <- data.frame(
    gene = overlap_genes_all,
    c9_lfc = res_C9_all[overlap_genes_all, "log2FoldChange"],
    sod1_lfc = res_SOD1_all[overlap_genes_all, "log2FoldChange"]
  )
  same_dir_all <- sum(sign(overlap_data_all$c9_lfc) == sign(overlap_data_all$sod1_lfc))
  cat("  Same direction:", same_dir_all, "(", round(100*same_dir_all/length(overlap_genes_all), 1), "%)\n\n")
}

# Create divergence dataframe
divergence_all <- data.frame(
  gene = rownames(res_C9_all),
  log2FC_C9 = res_C9_all$log2FoldChange,
  padj_C9 = res_C9_all$padj,
  log2FC_SOD1 = res_SOD1_all$log2FoldChange,
  padj_SOD1 = res_SOD1_all$padj
)

# Significant in EITHER
div_either_all <- divergence_all %>%
  filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1)) %>%
  filter(padj_C9 < 0.05 | padj_SOD1 < 0.05)

cat("Sig in either:", nrow(div_either_all), "genes\n\n")

# ===== CORRELATIONS =====
cat("--- Correlation Analysis ---\n")

# 1. All genes (we already have this)
div_all_all <- divergence_all %>%
  filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1))

cor_all_all <- cor.test(div_all_all$log2FC_C9, div_all_all$log2FC_SOD1)
cat("All genes: r =", round(cor_all_all$estimate, 3), 
    ", p =", format(cor_all_all$p.value, scientific = TRUE), 
    ", n =", nrow(div_all_all), "\n")

# 2. Sig in either (THIS IS WHAT'S MISSING)
if (nrow(div_either_all) > 10) {
  cor_either_all <- cor.test(div_either_all$log2FC_C9, div_either_all$log2FC_SOD1)
  cat("Sig in either: r =", round(cor_either_all$estimate, 3), 
      ", p =", format(cor_either_all$p.value, scientific = TRUE), 
      ", n =", nrow(div_either_all), "\n")
} else {
  cor_either_all <- list(estimate = NA)
}

# 3. Sig with |LFC| > 0.5 (we already have this)
div_sig_lfc_all <- divergence_all %>%
  filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1)) %>%
  filter(
    (padj_C9 < 0.05 & abs(log2FC_C9) > 0.5) | 
      (padj_SOD1 < 0.05 & abs(log2FC_SOD1) > 0.5)
  )

if (nrow(div_sig_lfc_all) > 10) {
  cor_sig_lfc_all <- cor.test(div_sig_lfc_all$log2FC_C9, div_sig_lfc_all$log2FC_SOD1)
  cat("Sig + |LFC| > 0.5: r =", round(cor_sig_lfc_all$estimate, 3), 
      ", p =", format(cor_sig_lfc_all$p.value, scientific = TRUE), 
      ", n =", nrow(div_sig_lfc_all), "\n")
}

# ===== UPDATE COMPLETE COMPARISON TABLE =====
cat("\n========================================\n")
cat("COMPLETE THREE-MODEL COMPARISON\n")
cat("========================================\n\n")

complete_comparison <- data.frame(
  Model = c("Simple (SEX only)", 
            "SVA filtered (6 clean SVs)", 
            "SVA all (10 SVs)"),
  C9_DEGs = c(333, 812, 1893),
  SOD1_DEGs = c(28, 17, 312),
  Overlap_both = c(7, 7, length(overlap_genes_all)),
  Sig_either = c(354, 822, nrow(div_either_all)),
  Cor_all = c(0.258, 0.320, 0.377),
  Cor_sig_either = c(0.362, 0.634, round(cor_either_all$estimate, 3)),
  Cor_sig_LFC05 = c(0.378, 0.638, 0.661),
  n_genes_LFC05 = c(179, 226, 324)
)

print(complete_comparison)

# Pretty print
cat("\n--- Summary ---\n")
for (i in 1:nrow(complete_comparison)) {
  cat("\n", complete_comparison$Model[i], "\n")
  cat("  C9 DEGs:", complete_comparison$C9_DEGs[i], "\n")
  cat("  SOD1 DEGs:", complete_comparison$SOD1_DEGs[i], "\n")
  cat("  Overlap:", complete_comparison$Overlap_both[i], "\n")
  cat("  Sig either:", complete_comparison$Sig_either[i], "\n")
  cat("  Cor (all):", complete_comparison$Cor_all[i], "\n")
  cat("  Cor (sig either):", complete_comparison$Cor_sig_either[i], "\n")
  cat("  Cor (LFC>0.5):", complete_comparison$Cor_sig_LFC05[i], 
      " (", complete_comparison$n_genes_LFC05[i], " genes)\n")
}

# Save complete table
write.csv(complete_comparison, "complete_three_model_comparison.csv", row.names = FALSE)

cat("\n✓ Complete comparison table saved!\n")

# Create a clean presentation table
presentation_table <- data.frame(
  Model = c("Simple", "SVA Filtered", "SVA All"),
  `C9 DEGs` = c(333, 812, 1893),
  `SOD1 DEGs` = c(28, 17, 312),
  `Overlap` = c(7, 7, 63),
  `Total Sig` = c(354, 822, 2142),
  `Correlation` = c("0.378 (179)", "0.638 (226)", "0.661 (324)"),
  check.names = FALSE
)

cat("\n========================================\n")
cat("FINAL MODEL COMPARISON\n")
cat("========================================\n\n")
print(presentation_table)

cat("\n--- Key Findings ---\n")
cat("✓ Filtered model shows 2.4x more C9 genes than simple\n")
cat("✓ Filtered model is conservative for SOD1 (17 vs 312 in all-SV)\n")
cat("✓ Convergence robust: r = 0.638 (69% stronger than simple)\n")
cat("✓ All-SV model shows suspicious explosions (18x SOD1, 9x overlap)\n")
cat("\n→ Recommendation: Use SVA filtered (6 clean SVs) as final model\n")




# ===== SETUP =====
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(msigdbr)

# You should have these from your filtered SVA model:
# res_C9_filtered
# res_SOD1_filtered

# ===== STEP 1: PREPARE GENE LISTS =====

# Convert gene symbols to Entrez IDs
# (clusterProfiler needs Entrez IDs for most databases)

# Get all genes with both C9 and SOD1 results
all_genes <- data.frame(
  symbol = rownames(res_C9_filtered),
  log2FC_C9 = res_C9_filtered$log2FoldChange,
  padj_C9 = res_C9_filtered$padj,
  log2FC_SOD1 = res_SOD1_filtered$log2FoldChange,
  padj_SOD1 = res_SOD1_filtered$padj
) %>%
  filter(!is.na(log2FC_C9) & !is.na(log2FC_SOD1))

# Convert symbols to Entrez IDs
gene_mapping <- bitr(all_genes$symbol, 
                     fromType = "ENSEMBL",
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

# Add Entrez IDs
all_genes <- all_genes %>%
  left_join(gene_mapping, by = c("symbol" = "ENSEMBL"))

cat("Total genes with mapping:", sum(!is.na(all_genes$ENTREZID)), "\n")

# ===== STEP 2: DEFINE CONVERGENT GENE SETS =====

# Convergent genes: significant in EITHER, changing in SAME direction
convergent_genes <- all_genes %>%
  filter(padj_C9 < 0.05 | padj_SOD1 < 0.05) %>%
  filter(sign(log2FC_C9) == sign(log2FC_SOD1)) %>%
  filter(!is.na(ENTREZID))

cat("\nConvergent genes (same direction):", nrow(convergent_genes), "\n")
cat("  Upregulated in both:", sum(convergent_genes$log2FC_C9 > 0), "\n")
cat("  Downregulated in both:", sum(convergent_genes$log2FC_C9 < 0), "\n")

# High-confidence convergent (|LFC| > 0.5)
convergent_high_conf <- convergent_genes %>%
  filter(
    (padj_C9 < 0.05 & abs(log2FC_C9) > 0.5) | 
      (padj_SOD1 < 0.05 & abs(log2FC_SOD1) > 0.5)
  )

cat("High-confidence convergent:", nrow(convergent_high_conf), "\n")

# ===== STEP 3: GENE SET ENRICHMENT ANALYSIS (GSEA) =====
cat("\n========================================\n")
cat("RUNNING GSEA\n")
cat("========================================\n\n")

# Create ranked gene lists for GSEA
# Rank by log2FC (C9 and SOD1 separately)

# C9 ranked list
ranked_C9 <- all_genes %>%
  filter(!is.na(ENTREZID)) %>%
  arrange(desc(log2FC_C9))

gene_list_C9 <- ranked_C9$log2FC_C9
names(gene_list_C9) <- ranked_C9$ENTREZID

# SOD1 ranked list
ranked_SOD1 <- all_genes %>%
  filter(!is.na(ENTREZID)) %>%
  arrange(desc(log2FC_SOD1))

gene_list_SOD1 <- ranked_SOD1$log2FC_SOD1
names(gene_list_SOD1) <- ranked_SOD1$ENTREZID

# ===== FIX: HANDLE DUPLICATE ENTREZ IDs (CORRECTED) =====

# C9 ranked list
ranked_C9 <- all_genes %>%
  filter(!is.na(ENTREZID)) %>%
  group_by(ENTREZID) %>%
  # Keep gene with highest absolute log2FC if duplicates
  slice_max(abs(log2FC_C9), n = 1, with_ties = FALSE) %>%
  ungroup()

# Create named vector
gene_list_C9 <- ranked_C9$log2FC_C9
names(gene_list_C9) <- ranked_C9$ENTREZID

# Add tiny noise to break ties, THEN sort
set.seed(123)
gene_list_C9 <- gene_list_C9 + rnorm(length(gene_list_C9), 0, 1e-10)
gene_list_C9 <- sort(gene_list_C9, decreasing = TRUE)

cat("C9 list length:", length(gene_list_C9), "\n")
cat("C9 duplicates:", sum(duplicated(names(gene_list_C9))), "\n")
cat("C9 is sorted:", all(diff(gene_list_C9) <= 0), "\n")

# SOD1 ranked list
ranked_SOD1 <- all_genes %>%
  filter(!is.na(ENTREZID)) %>%
  group_by(ENTREZID) %>%
  slice_max(abs(log2FC_SOD1), n = 1, with_ties = FALSE) %>%
  ungroup()

gene_list_SOD1 <- ranked_SOD1$log2FC_SOD1
names(gene_list_SOD1) <- ranked_SOD1$ENTREZID

# Add tiny noise to break ties, THEN sort
set.seed(124)
gene_list_SOD1 <- gene_list_SOD1 + rnorm(length(gene_list_SOD1), 0, 1e-10)
gene_list_SOD1 <- sort(gene_list_SOD1, decreasing = TRUE)

cat("SOD1 list length:", length(gene_list_SOD1), "\n")
cat("SOD1 duplicates:", sum(duplicated(names(gene_list_SOD1))), "\n")
cat("SOD1 is sorted:", all(diff(gene_list_SOD1) <= 0), "\n\n")

# ===== NOW RUN GSEA =====
cat("Running GSEA for C9orf72...\n")
gsea_C9 <- gseGO(
  geneList = gene_list_C9,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

cat("Running GSEA for SOD1...\n")
gsea_SOD1 <- gseGO(
  geneList = gene_list_SOD1,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

cat("C9 enriched pathways:", nrow(gsea_C9), "\n")
cat("SOD1 enriched pathways:", nrow(gsea_SOD1), "\n")

# Run GSEA for both (using GO Biological Process)
cat("Running GSEA for C9orf72...\n")
gsea_C9 <- gseGO(
  geneList = gene_list_C9,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

cat("Running GSEA for SOD1...\n")
gsea_SOD1 <- gseGO(
  geneList = gene_list_SOD1,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

cat("C9 enriched pathways:", nrow(gsea_C9), "\n")
cat("SOD1 enriched pathways:", nrow(gsea_SOD1), "\n")

# Find SHARED pathways between C9 and SOD1
if (nrow(gsea_C9) > 0 & nrow(gsea_SOD1) > 0) {
  shared_pathways <- intersect(gsea_C9$ID, gsea_SOD1$ID)
  cat("Shared pathways (enriched in BOTH):", length(shared_pathways), "\n\n")
  
  # Show top shared pathways
  if (length(shared_pathways) > 0) {
    shared_results <- data.frame(
      Pathway = gsea_C9[gsea_C9$ID %in% shared_pathways, "Description"],
      C9_NES = gsea_C9[gsea_C9$ID %in% shared_pathways, "NES"],
      C9_padj = gsea_C9[gsea_C9$ID %in% shared_pathways, "p.adjust"],
      SOD1_NES = gsea_SOD1[gsea_SOD1$ID %in% shared_pathways, "NES"],
      SOD1_padj = gsea_SOD1[gsea_SOD1$ID %in% shared_pathways, "p.adjust"]
    ) %>%
      arrange(C9_padj)
    
    cat("Top 10 shared pathways:\n")
    print(head(shared_results, 10))
  }
}

# ===== STEP 4: OVER-REPRESENTATION ANALYSIS (ORA) =====
cat("\n========================================\n")
cat("OVER-REPRESENTATION ANALYSIS\n")
cat("========================================\n\n")

# ORA on convergent genes
cat("Analyzing convergent genes...\n")

# GO Biological Process
ora_bp <- enrichGO(
  gene = convergent_high_conf$ENTREZID,
  universe = all_genes$ENTREZID[!is.na(all_genes$ENTREZID)],
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

cat("GO BP terms enriched:", nrow(ora_bp), "\n")

# KEGG pathways
ora_kegg <- enrichKEGG(
  gene = convergent_high_conf$ENTREZID,
  universe = all_genes$ENTREZID[!is.na(all_genes$ENTREZID)],
  organism = 'hsa',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

cat("KEGG pathways enriched:", nrow(ora_kegg), "\n")

# Reactome pathways
library(ReactomePA)
ora_reactome <- enrichPathway(
  gene = convergent_high_conf$ENTREZID,
  universe = all_genes$ENTREZID[!is.na(all_genes$ENTREZID)],
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE
)

cat("Reactome pathways enriched:", nrow(ora_reactome), "\n\n")

# ===== STEP 5: CHECK FOR ALS-RELEVANT PATHWAYS =====
cat("========================================\n")
cat("ALS-RELEVANT PATHWAY ANALYSIS\n")
cat("========================================\n\n")

# Key ALS pathways to look for
als_keywords <- c(
  "RNA", "splicing", "transport", "processing",
  "protein", "ubiquitin", "proteasome", "autophagy",
  "stress", "endoplasmic reticulum", "ER",
  "mitochondri", "oxidative",
  "axon", "neuron", "synap",
  "inflammation", "immune",
  "apoptosis", "death"
)

# Search in GO results
if (nrow(ora_bp) > 0) {
  als_related <- ora_bp@result %>%
    filter(grepl(paste(als_keywords, collapse = "|"), Description, ignore.case = TRUE)) %>%
    select(Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
    arrange(p.adjust)
  
  cat("ALS-relevant GO terms found:", nrow(als_related), "\n")
  if (nrow(als_related) > 0) {
    cat("\nTop ALS-relevant pathways:\n")
    print(head(als_related, 15))
  }
}

# ===== STEP 6: VISUALIZATIONS =====
cat("\n========================================\n")
cat("CREATING VISUALIZATIONS\n")
cat("========================================\n\n")

# 1. Dotplot of top GO terms
if (nrow(ora_bp) > 0) {
  p1 <- dotplot(ora_bp, showCategory = 20, title = "GO BP: Convergent Genes")
  print(p1)
  ggsave("pathway_go_bp_dotplot.png", p1, width = 10, height = 8, dpi = 300)
}

# 2. Barplot of KEGG pathways
if (nrow(ora_kegg) > 0) {
  p2 <- barplot(ora_kegg, showCategory = 15, title = "KEGG: Convergent Genes")
  print(p2)
  ggsave("pathway_kegg_barplot.png", p2, width = 10, height = 6, dpi = 300)
}

# 3. GSEA enrichment plots for shared pathways
if (length(shared_pathways) > 0) {
  # Plot top 3 shared pathways
  top_shared <- head(shared_pathways, 3)
  
  for (i in 1:length(top_shared)) {
    pathway_id <- top_shared[i]
    pathway_name <- gsea_C9[gsea_C9$ID == pathway_id, "Description"]
    
    p_c9 <- gseaplot2(gsea_C9, geneSetID = pathway_id, 
                      title = paste0("C9: ", pathway_name))
    p_sod1 <- gseaplot2(gsea_SOD1, geneSetID = pathway_id,
                        title = paste0("SOD1: ", pathway_name))
    
    ggsave(paste0("gsea_shared_", i, "_C9.png"), p_c9, width = 8, height = 6, dpi = 300)
    ggsave(paste0("gsea_shared_", i, "_SOD1.png"), p_sod1, width = 8, height = 6, dpi = 300)
  }
}

# 4. Compare C9 vs SOD1 pathway NES scores
if (length(shared_pathways) >= 10) {
  comparison_plot_data <- data.frame(
    Pathway = gsea_C9[gsea_C9$ID %in% head(shared_pathways, 20), "Description"],
    C9_NES = gsea_C9[gsea_C9$ID %in% head(shared_pathways, 20), "NES"],
    SOD1_NES = gsea_SOD1[gsea_SOD1$ID %in% head(shared_pathways, 20), "NES"]
  )
  
  p4 <- ggplot(comparison_plot_data, aes(x = C9_NES, y = SOD1_NES)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(
      title = "Pathway Convergence: C9 vs SOD1",
      subtitle = "Normalized Enrichment Scores for Shared Pathways",
      x = "C9orf72 NES",
      y = "SOD1 NES"
    ) +
    theme_bw(base_size = 12)
  
  print(p4)
  ggsave("pathway_c9_vs_sod1_comparison.png", p4, width = 8, height = 6, dpi = 300)
  
  # Calculate correlation
  cor_nes <- cor.test(comparison_plot_data$C9_NES, comparison_plot_data$SOD1_NES)
  cat("\nCorrelation of pathway NES scores:\n")
  cat("r =", round(cor_nes$estimate, 3), ", p =", format(cor_nes$p.value, scientific = TRUE), "\n")
}

# ===== STEP 7: SAVE RESULTS =====
cat("\n========================================\n")
cat("SAVING RESULTS\n")
cat("========================================\n\n")

# Save convergent gene list
write.csv(convergent_high_conf %>% select(symbol, log2FC_C9, padj_C9, log2FC_SOD1, padj_SOD1),
          "convergent_genes_high_confidence.csv", row.names = FALSE)

# Save pathway results
if (nrow(ora_bp) > 0) {
  write.csv(ora_bp@result, "pathway_GO_BP_results.csv", row.names = FALSE)
}
if (nrow(ora_kegg) > 0) {
  write.csv(ora_kegg@result, "pathway_KEGG_results.csv", row.names = FALSE)
}
if (nrow(ora_reactome) > 0) {
  write.csv(ora_reactome@result, "pathway_Reactome_results.csv", row.names = FALSE)
}

# Save shared pathways
if (exists("shared_results")) {
  write.csv(shared_results, "pathway_shared_C9_SOD1.csv", row.names = FALSE)
}

cat("✓ All results saved!\n")


#Include batch in model
annot_als_filtered$group <- relevel(factor(annot_als_filtered$group), ref = "Control")
annot_als_filtered$Batch.x <- factor(annot_als_filtered$Batch.x)

dds_batch <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = annot_als_filtered,
  design    = ~ Batch.x + SEX + group
)

dds_batch <- DESeq(dds_batch)

res_c9_vs_ctrl_batch <- results(dds_batch, contrast = c("group", "C9orf72", "control"))




# Load the mutation data
sod1_mutations <- read.csv("sod1_sample_mutations.csv", stringsAsFactors=FALSE)

# Extract genomic ID from Sample_ID
extract_genomic_id <- function(sample_id) {
  parts <- strsplit(sample_id, "-")[[1]]
  paste(parts[1], parts[2], sep="-")
}

annot_als_filtered$genomic_id <- sapply(annot_als_filtered$Sample_ID, extract_genomic_id)

# Match mutations
annot_als_filtered$sod1_mutation <- NA
annot_als_filtered$sod1_genotype <- NA

for(i in 1:nrow(annot_als_filtered)) {
  match_idx <- which(sod1_mutations$sample_id == annot_als_filtered$genomic_id[i])
  if(length(match_idx) > 0) {
    annot_als_filtered$sod1_mutation[i] <- sod1_mutations$mutation_canonical[match_idx]
    annot_als_filtered$sod1_genotype[i] <- sod1_mutations$genotype[match_idx]
  }
}

# Check SOD1 patients using Mutation_Subtype
cat("\n=== YOUR SOD1 PATIENTS ===\n")
sod1_patients <- annot_als_filtered[annot_als_filtered$Mutation_Subtype == "SOD1", 
                                    c("Sample_ID", "genomic_id", "Mutation_Subtype", "sod1_mutation")]

print(sod1_patients)

cat("\n=== SUMMARY ===\n")
cat("Total SOD1 patients:", nrow(sod1_patients), "\n")
cat("With A4V/D90A/G93A/H46R:", sum(!is.na(sod1_patients$sod1_mutation)), "\n")
cat("Missing mutation info:", sum(is.na(sod1_patients$sod1_mutation)), "\n")

cat("\nMutation breakdown:\n")
print(table(sod1_patients$sod1_mutation, useNA="always"))

# Who's missing?
missing_patients <- sod1_patients[is.na(sod1_patients$sod1_mutation), ]
if(nrow(missing_patients) > 0) {
  cat("\n=== PATIENTS WITHOUT A4V/D90A/G93A/H46R ===\n")
  cat("These", nrow(missing_patients), "patients have other SOD1 mutations:\n")
  print(missing_patients)
  
  # Save their IDs
  writeLines(missing_patients$genomic_id, "missing_sod1_genomic_ids.txt")
  cat("\nSaved to: missing_sod1_genomic_ids.txt\n")
}

# Save updated metadata
write.csv(annot_als_filtered, "annot_als_filtered_with_sod1_mutations.csv", row.names=FALSE)
cat("\nSaved updated metadata to: annot_als_filtered_with_sod1_mutations.csv\n")

# Simple table of 13 identified SOD1 patients
sod1_identified <- annot_als_filtered[annot_als_filtered$Mutation_Subtype == "SOD1" & 
                                        !is.na(annot_als_filtered$sod1_mutation_detailed), 
                                      c("Sample_ID", "genomic_id", "sod1_mutation_detailed")]

cat("=== 13 IDENTIFIED SOD1 PATIENTS ===\n\n")

# Count by mutation
mutation_counts <- as.data.frame(table(sod1_identified$sod1_mutation_detailed))
colnames(mutation_counts) <- c("Mutation", "Count")
mutation_counts <- mutation_counts[order(-mutation_counts$Count), ]

# Add grouping
mutation_counts$Group <- NA
mutation_counts$Group[mutation_counts$Mutation == "p.A5V"] <- "A4V (Wong 2019)"
mutation_counts$Group[mutation_counts$Mutation %in% c("p.D91A", "p.G94A", "p.H49R")] <- "Wong 2019 common"
mutation_counts$Group[mutation_counts$Mutation == "p.I114T"] <- "I113T (largest group)"
mutation_counts$Group[is.na(mutation_counts$Group)] <- "Other pathogenic"

print(mutation_counts)

cat("\n=== SUMMARY ===\n")
cat("Wong 2019 mutations (A4V, D90A, G93A, H46R):", 
    sum(sod1_identified$sod1_mutation_detailed %in% c("p.A5V", "p.D91A", "p.G94A", "p.H49R")), "\n")
cat("I113T (your largest single group):", 
    sum(sod1_identified$sod1_mutation_detailed == "p.I114T"), "\n")
cat("Other mutations:", 
    sum(!sod1_identified$sod1_mutation_detailed %in% c("p.A5V", "p.D91A", "p.G94A", "p.H49R", "p.I114T")), "\n")
cat("TOTAL:", nrow(sod1_identified), "\n")



# Required packages
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
})

# Optional but strongly recommended
# install.packages(c("matrixStats"))
# BiocManager::install(c("apeglm", "sva"))
# install.packages(c("fgsea", "msigdbr"))
rownames(annot_als_filtered) <- annot_als_filtered$Sample_ID
stopifnot(identical(colnames(counts_filtered), rownames(annot_als_filtered)))
annot_als_filtered$group <- annot_als_filtered$Mutation_Subtype
# Make factors explicit + set reference
annot_als_filtered <- annot_als_filtered %>%
  mutate(
    SEX = factor(SEX),
    Batch.x = factor(Batch.x),
    group = factor(group, levels = c("Control", "C9orf72", "SOD1", "Sporadic"))
  )

# Subset to C9 vs control only for a clean fit
keep_c9_ctrl <- annot_als_filtered$group %in% c("Control", "C9orf72")
counts_cc <- counts_filtered[, keep_c9_ctrl]
col_cc <- annot_als_filtered[keep_c9_ctrl, ] %>%
  droplevels()

col_cc$group <- relevel(col_cc$group, ref = "Control")

table(col_cc$group)

dds_simple <- DESeqDataSetFromMatrix(
  countData = counts_cc,
  colData   = col_cc,
  design    = ~ SEX + group
)
dds_simple <- DESeq(dds_simple)


dds_batch <- DESeqDataSetFromMatrix(
  countData = counts_cc,
  colData   = col_cc,
  design    = ~ Batch.x + SEX + group
)
dds_batch <- DESeq(dds_batch)

library(sva)

# Model matrix for svaseq: full and null
mod  <- model.matrix(~ SEX + group, data = col_cc)
mod0 <- model.matrix(~ SEX, data = col_cc)

svobj <- svaseq(as.matrix(counts_cc), mod, mod0)  # estimates SVs
sv <- as.data.frame(svobj$sv)
colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
col_cc_sva <- bind_cols(col_cc, sv)

# Quick check: are SVs associated with group?
sv_p <- sapply(sv, function(x) summary(lm(x ~ col_cc_sva$group))$coef[2,4])
sv_p

# Choose a small number of SVs (e.g., 2–6). You can start with 3.
k <- min(3, ncol(sv))

design_sva <- as.formula(
  paste("~", paste(c(paste0("SV", 1:k), "SEX", "group"), collapse = " + "))
)

dds_sva <- DESeqDataSetFromMatrix(
  countData = counts_cc,
  colData   = col_cc_sva,
  design    = design_sva
)
dds_sva <- DESeq(dds_sva)


library(apeglm)

library(dplyr)
library(tibble)

get_shrunken <- function(dds, coef_name){
  res_sh <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  as.data.frame(res_sh) %>%
    rownames_to_column("gene") %>%
    select(any_of(c("gene","log2FoldChange","lfcSE","stat","pvalue","padj")))
}


coef_simple <- "group_C9orf72_vs_Control"
coef_batch  <- "group_C9orf72_vs_Control"
coef_sva    <- "group_C9orf72_vs_Control"

res_simple <- get_shrunken(dds_simple, coef_simple)
res_batch  <- get_shrunken(dds_batch,  coef_batch)
res_sva    <- get_shrunken(dds_sva,    coef_sva)


merge_lfc <- function(a, b, suffix_a, suffix_b){
  inner_join(a %>% select(gene, log2FoldChange),
             b %>% select(gene, log2FoldChange),
             by = "gene",
             suffix = c(paste0("_", suffix_a), paste0("_", suffix_b)))
}

ab <- merge_lfc(res_simple, res_batch, "simple", "batch")
ac <- merge_lfc(res_simple, res_sva,   "simple", "sva")
bc <- merge_lfc(res_batch,  res_sva,   "batch",  "sva")

cor_simple_batch <- cor(ab$log2FoldChange_simple, ab$log2FoldChange_batch, use="pairwise.complete.obs")
cor_simple_sva   <- cor(ac$log2FoldChange_simple, ac$log2FoldChange_sva,   use="pairwise.complete.obs")
cor_batch_sva    <- cor(bc$log2FoldChange_batch,  bc$log2FoldChange_sva,   use="pairwise.complete.obs")

c(cor_simple_batch = cor_simple_batch,
  cor_simple_sva   = cor_simple_sva,
  cor_batch_sva    = cor_batch_sva)


top_genes <- function(res, n=100){
  res %>%
    filter(!is.na(log2FoldChange)) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    slice_head(n = n) %>%
    pull(gene)
}

top100_simple <- top_genes(res_simple, 100)
top100_batch  <- top_genes(res_batch,  100)
top100_sva    <- top_genes(res_sva,    100)

overlap <- function(x, y) length(intersect(x, y))

c(overlap_simple_batch = overlap(top100_simple, top100_batch),
  overlap_simple_sva   = overlap(top100_simple, top100_sva),
  overlap_batch_sva    = overlap(top100_batch,  top100_sva))


library(msigdbr)
library(fgsea)

# Build Hallmark gene sets
msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol)

pathways <- split(msig_h$gene_symbol, msig_h$gs_name)

library(dplyr)
library(tibble)

get_stat_res <- function(dds, coef_name){
  results(dds, name = coef_name) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    select(gene, stat, pvalue, padj)
}

res_simple_stat <- get_stat_res(dds_simple, coef_simple)
res_batch_stat  <- get_stat_res(dds_batch,  coef_batch)
res_sva_stat    <- get_stat_res(dds_sva,    coef_sva)


run_fgsea <- function(res_stat){
  ranks <- res_stat %>%
    filter(!is.na(stat)) %>%
    select(gene, stat) %>%
    distinct() %>%
    deframe()
  
  fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500, nperm = 10000) %>%
    arrange(padj)
}

fg_simple <- run_fgsea(res_simple_stat)
fg_batch  <- run_fgsea(res_batch_stat)
fg_sva    <- run_fgsea(res_sva_stat)

library(fgsea)
# Install once if needed:
# BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(tibble)

res_batch_stat_sym <- res_batch_stat %>%
  mutate(symbol = mapIds(org.Hs.eg.db,
                         keys = gene,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first") %>% unname())

res_simple_stat_sym <- res_simple_stat %>%
  mutate(symbol = mapIds(org.Hs.eg.db,
                         keys = gene,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first") %>% unname())

res_sva_stat_sym <- res_sva_stat %>%
  mutate(symbol = mapIds(org.Hs.eg.db,
                         keys = gene,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first") %>% unname())

make_ranks <- function(res_sym){
  res_sym %>%
    filter(!is.na(stat), !is.na(symbol)) %>%
    distinct(symbol, .keep_all = TRUE) %>%
    dplyr::select(symbol, stat) %>%
    deframe()
}

ranks_simple <- make_ranks(res_simple_stat_sym)
ranks_batch  <- make_ranks(res_batch_stat_sym)
ranks_sva    <- make_ranks(res_sva_stat_sym)

# sanity check: should now be >0
length(intersect(names(ranks_batch), msig_syms))

library(fgsea)

run_fgsea <- function(ranks){
  fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500) %>%
    arrange(padj)
}

fg_simple <- run_fgsea(ranks_simple)
fg_batch  <- run_fgsea(ranks_batch)
fg_sva    <- run_fgsea(ranks_sva)

head(fg_batch, 15)


head(fg_batch, 15)  # look at top signals from batch-adjusted model

topN <- 20
stable_sets <- Reduce(intersect, list(head(fg_simple$pathway, topN),
                                      head(fg_batch$pathway,  topN),
                                      head(fg_sva$pathway,    topN)))
stable_sets

library(matrixStats)

# Variance-stabilized expression for scoring
vsd <- vst(dds_batch, blind = FALSE)
expr <- assay(vsd)

# Make sure rownames(expr) are gene symbols matching res_batch$gene
# If not, adjust mapping (depends on how you built counts_filtered)

sig_up <- res_batch %>% arrange(desc(log2FoldChange)) %>% filter(!is.na(log2FoldChange)) %>% slice_head(n=100) %>% pull(gene)
sig_dn <- res_batch %>% arrange(log2FoldChange)       %>% filter(!is.na(log2FoldChange)) %>% slice_head(n=100) %>% pull(gene)

sig_up <- intersect(sig_up, rownames(expr))
sig_dn <- intersect(sig_dn, rownames(expr))

c(length(sig_up), length(sig_dn))

score_c9 <- colMeans(expr[sig_up, , drop=FALSE]) - colMeans(expr[sig_dn, , drop=FALSE])
score_df <- data.frame(Sample_ID = colnames(expr), C9_signature_score = score_c9) %>%
  left_join(col_cc %>% dplyr::select(Sample_ID, group, SEX, Batch.x), by = "Sample_ID")


sod1 <- read.csv("sod1_sample_mutations.csv", stringsAsFactors = FALSE)

sod1_flags <- sod1 %>%
  group_by(sample_id) %>%
  summarise(
    SOD1_any = TRUE,
    SOD1_A4V = any(mutation_canonical == "A4V"),
    .groups = "drop"
  ) %>%
  rename(Sample_ID = sample_id)

score_df <- score_df %>%
  left_join(sod1_flags, by = "Sample_ID") %>%
  mutate(
    SOD1_any = ifelse(is.na(SOD1_any), FALSE, SOD1_any),
    SOD1_A4V = ifelse(is.na(SOD1_A4V), FALSE, SOD1_A4V)
  )

# Compare SOD1_any vs controls among the same C9/control subset
dat <- score_df %>% filter(group == "control" | SOD1_any)

with(dat, table(group, SOD1_any))

wilcox.test(C9_signature_score ~ SOD1_any, data = dat %>% filter(group == "control" | SOD1_any))

# A4V-only vs controls
dat_a4v <- score_df %>% filter(group == "control" | SOD1_A4V)
wilcox.test(C9_signature_score ~ SOD1_A4V, data = dat_a4v)


perm_test <- function(values, labels, nperm = 20000){
  obs <- abs(median(values[labels]) - median(values[!labels]))
  perm <- replicate(nperm, {
    labp <- sample(labels)
    abs(median(values[labp]) - median(values[!labp]))
  })
  mean(perm >= obs)
}

# SOD1_any vs control
vals <- dat$C9_signature_score
lab  <- dat$SOD1_any
perm_test(vals, lab)

# A4V vs control
vals2 <- dat_a4v$C9_signature_score
lab2  <- dat_a4v$SOD1_A4V
perm_test(vals2, lab2)


c(
  simple_NA = sum(is.na(res_simple$log2FoldChange)),
  batch_NA  = sum(is.na(res_batch$log2FoldChange)),
  sva6_NA   = sum(is.na(res_sva$log2FoldChange))
)

c(
  simple_nonfinite = sum(!is.finite(res_simple$log2FoldChange)),
  batch_nonfinite  = sum(!is.finite(res_batch$log2FoldChange)),
  sva6_nonfinite   = sum(!is.finite(res_sva$log2FoldChange))
)

c(
  simple_total = nrow(res_simple),
  batch_total  = nrow(res_batch),
  sva6_total   = nrow(res_sva)
)


## Run C9 signature score three ways (Simple / Batch / SVA6) and test SOD1 projection
## Assumes you already have: dds_simple, dds_batch, dds_sva (SVA6),
## counts_filtered, annot_als_filtered, and sod1_sample_mutations.csv

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(apeglm)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ----------------------------
# Helpers
# ----------------------------

map_ensg_to_symbol <- function(ensg_vec){
  sym <- mapIds(org.Hs.eg.db,
                keys = ensg_vec,
                column = "SYMBOL",
                keytype = "ENSEMBL",
                multiVals = "first")
  unname(sym)
}

get_signature_genes <- function(dds, coef_name, Nsig = 100){
  # shrunken LFC for stable gene selection
  res_sh <- lfcShrink(dds, coef = coef_name, type = "apeglm") %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    mutate(symbol = map_ensg_to_symbol(gene)) %>%
    filter(!is.na(symbol)) %>%
    distinct(symbol, .keep_all = TRUE)
  
  up <- res_sh %>% arrange(desc(log2FoldChange)) %>% slice_head(n = Nsig) %>% pull(symbol)
  dn <- res_sh %>% arrange(log2FoldChange)       %>% slice_head(n = Nsig) %>% pull(symbol)
  
  list(up = up, dn = dn, res_sh = res_sh)
}

collapse_expr_to_symbol <- function(expr_mat){
  # expr_mat rownames are ENSG
  sym_map <- map_ensg_to_symbol(rownames(expr_mat))
  expr_df <- as.data.frame(expr_mat)
  expr_df$symbol <- sym_map
  expr_df <- expr_df %>% filter(!is.na(symbol))
  
  expr_sym <- expr_df %>%
    group_by(symbol) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop") %>%
    column_to_rownames("symbol") %>%
    as.matrix()
  
  expr_sym
}

compute_signature_score <- function(expr_sym, sig_up, sig_dn){
  up2 <- intersect(sig_up, rownames(expr_sym))
  dn2 <- intersect(sig_dn, rownames(expr_sym))
  if (length(up2) < 10 || length(dn2) < 10) {
    warning(paste("Very small overlap with expression matrix:",
                  "up =", length(up2), "down =", length(dn2)))
  }
  colMeans(expr_sym[up2, , drop = FALSE]) - colMeans(expr_sym[dn2, , drop = FALSE])
}

perm_test_median <- function(values, labels, nperm = 20000, seed = 1){
  set.seed(seed)
  labels <- as.logical(labels)
  obs <- abs(median(values[labels]) - median(values[!labels]))
  perm <- replicate(nperm, {
    labp <- sample(labels)
    abs(median(values[labp]) - median(values[!labp]))
  })
  mean(perm >= obs)
}

test_sod1_projection <- function(score_df){
  # expects columns: C9_score, group, SOD1_any, SOD1_A4V
  # Compare against Controls only (as you requested earlier)
  out <- list()
  
  # SOD1_any vs Control
  dat_any <- score_df %>% filter(group == "Control" | SOD1_any)
  out$wilcox_any <- wilcox.test(C9_score ~ SOD1_any, data = dat_any)
  out$perm_any <- perm_test_median(dat_any$C9_score, dat_any$SOD1_any)
  
  # A4V vs Control
  dat_a4v <- score_df %>% filter(group == "Control" | SOD1_A4V)
  out$wilcox_a4v <- wilcox.test(C9_score ~ SOD1_A4V, data = dat_a4v)
  out$perm_a4v <- perm_test_median(dat_a4v$C9_score, dat_a4v$SOD1_A4V)
  
  out
}

# ----------------------------
# Prepare expression matrix (once)
# ----------------------------

# Use vst from a DESeq object that contains ALL samples. Batch is fine since it's your primary.
# If dds_batch only contains C9/control, then instead build a dds_all and vst(dds_all).
vsd_all <- vst(dds_batch, blind = FALSE)   # change to dds_all if needed
expr <- assay(vsd_all)

expr_sym <- collapse_expr_to_symbol(expr)

# ----------------------------
# Add SOD1 flags
# ----------------------------

sod1 <- read.csv("sod1_sample_mutations.csv", stringsAsFactors = FALSE)


sod1_flags_pid <- sod1 %>%
  group_by(sample_id) %>%            # sample_id == Participant_ID
  summarise(
    SOD1_any = TRUE,
    SOD1_A4V = any(mutation_canonical == "A4V"),
    .groups = "drop"
  ) %>%
  rename(Participant_ID = sample_id)


# Build a metadata df aligned to expression columns
meta <- annot_als_filtered %>%
  right_join(data.frame(Sample_ID = colnames(expr_sym)), by = "Sample_ID") %>%
  left_join(sod1_flags, by = "Sample_ID") %>%
  mutate(
    SOD1_any = ifelse(is.na(SOD1_any), FALSE, SOD1_any),
    SOD1_A4V = ifelse(is.na(SOD1_A4V), FALSE, SOD1_A4V)
  )

# ----------------------------
# Compute three signatures + scores
# ----------------------------

Nsig <- 100

coef_simple <- "group_C9orf72_vs_Control"
coef_batch  <- "group_C9orf72_vs_Control"
coef_sva6   <- "group_C9orf72_vs_Control"   # adjust if your resultsNames differs

sig_simple <- get_signature_genes(dds_simple, coef_simple, Nsig = Nsig)
sig_batch  <- get_signature_genes(dds_batch,  coef_batch,  Nsig = Nsig)
sig_sva6   <- get_signature_genes(dds_sva,    coef_sva6,   Nsig = Nsig)

score_simple <- compute_signature_score(expr_sym, sig_simple$up, sig_simple$dn)
score_batch  <- compute_signature_score(expr_sym, sig_batch$up,  sig_batch$dn)
score_sva6   <- compute_signature_score(expr_sym, sig_sva6$up,   sig_sva6$dn)

scores_long <- bind_rows(
  data.frame(Sample_ID = names(score_simple), model = "Simple", C9_score = as.numeric(score_simple)),
  data.frame(Sample_ID = names(score_batch),  model = "Batch",  C9_score = as.numeric(score_batch)),
  data.frame(Sample_ID = names(score_sva6),   model = "SVA6",   C9_score = as.numeric(score_sva6))
) %>%
  left_join(meta, by = "Sample_ID")

# Quick sanity check: C9 should generally separate from Controls
scores_long <- scores_long %>%
  dplyr::select(-SOD1_any, -SOD1_A4V) %>%             # remove the old (wrong-key) flags
  mutate(Participant_ID = as.character(Participant_ID)) %>%
  left_join(sod1_flags_pid, by = "Participant_ID") %>%
  mutate(
    SOD1_any = ifelse(is.na(SOD1_any), FALSE, SOD1_any),
    SOD1_A4V = ifelse(is.na(SOD1_A4V), FALSE, SOD1_A4V)
  )

table(scores_long$SOD1_any)
table(scores_long$SOD1_A4V)


table(scores_long$SOD1_any, useNA = "always")
table(scores_long$SOD1_A4V, useNA = "always")


# ----------------------------
# Run SOD1 projection tests per model
# ----------------------------

proj_results <- scores_long %>%
  group_split(model) %>%
  setNames(unique(scores_long$model)) %>%
  lapply(function(df) test_sod1_projection(df))

# Print results in a compact way
for (m in names(proj_results)) {
  cat("\n=== Model:", m, "===\n")
  cat("SOD1_any vs Control: Wilcox p =", proj_results[[m]]$wilcox_any$p.value,
      "; Perm p =", proj_results[[m]]$perm_any, "\n")
  cat("SOD1_A4V vs Control: Wilcox p =", proj_results[[m]]$wilcox_a4v$p.value,
      "; Perm p =", proj_results[[m]]$perm_a4v, "\n")
}
# The above doesn't work because SOD1 is a mix of all variants which is different from the paper

# Optional: save scores for plotting/slide tables
write.csv(scores_long, "C9_signature_scores_three_models.csv", row.names = FALSE)







library(dplyr)
library(tibble)
library(stringr)

sod1_rows <- annot_als_filtered %>%
  # If Sample_ID is already a column (it is), no need for rownames_to_column
  filter(group == "SOD1") %>%
  distinct(Participant_ID, Sample_ID, HAS_GENOMICS, Mutation_Subtype, ClinReport_Mutations_Details)

nrow(sod1_rows)                       # expect 13
n_distinct(sod1_rows$Participant_ID)   # expect 12
sod1_rows %>% count(Participant_ID, sort = TRUE)  # shows who has 2 samples
sod1_pids <- sod1_rows %>% distinct(Participant_ID) %>% pull(Participant_ID)
sod1_pids
length(sod1_pids)


library(dplyr)
library(stringr)

sod1_clin_pids <- annot_als_filtered %>%
  filter(str_detect(toupper(ClinReport_Mutations_Details), "SOD1")) %>%
  distinct(Participant_ID) %>%
  pull(Participant_ID)

length(sod1_clin_pids)
sod1_clin_pids

sod1_tx_pids <- annot_als_filtered %>%
  filter(group == "SOD1") %>%
  distinct(Participant_ID) %>%
  pull(Participant_ID)

sod1_wgs_pids <- c("CASE-NEUAD952KAZ","CASE-NEUBR543CRA","CASE-NEUCD502BFU",
                   "CASE-NEUCT842RJV","CASE-NEUMA836GGL","CASE-NEURG115XK8")

list(
  tx_only   = setdiff(sod1_tx_pids, union(sod1_clin_pids, sod1_wgs_pids)),
  clin_only = setdiff(sod1_clin_pids, union(sod1_tx_pids, sod1_wgs_pids)),
  wgs_only  = setdiff(sod1_wgs_pids, union(sod1_tx_pids, sod1_clin_pids)),
  tx_clin_overlap = intersect(sod1_tx_pids, sod1_clin_pids),
  tx_wgs_overlap  = intersect(sod1_tx_pids, sod1_wgs_pids),
  clin_wgs_overlap = intersect(sod1_clin_pids, sod1_wgs_pids)
)


# ----------------------------
# Core prep (shared)
# ----------------------------
prep_pair <- function(resA, resB,
                      gene_col = "gene",
                      lfc_col  = "log2FoldChange",
                      padj_col = "padj",
                      alpha = 0.05,
                      universe = c("intersection", "union")) {
  
  universe <- match.arg(universe)
  
  # Ensure gene column exists; if not, take rownames
  if (!gene_col %in% colnames(resA)) resA[[gene_col]] <- rownames(resA)
  if (!gene_col %in% colnames(resB)) resB[[gene_col]] <- rownames(resB)
  
  A <- resA[, c(gene_col, lfc_col, padj_col)]
  B <- resB[, c(gene_col, lfc_col, padj_col)]
  colnames(A) <- c("gene", "lfcA", "padjA")
  colnames(B) <- c("gene", "lfcB", "padjB")
  
  U <- if (universe == "intersection") intersect(A$gene, B$gene) else union(A$gene, B$gene)
  
  A <- A[match(U, A$gene), ]
  B <- B[match(U, B$gene), ]
  
  df <- data.frame(
    gene  = U,
    lfcA  = A$lfcA, padjA = A$padjA,
    lfcB  = B$lfcB, padjB = B$padjB
  )
  
  df$DEG_A <- !is.na(df$padjA) & df$padjA <= alpha
  df$DEG_B <- !is.na(df$padjB) & df$padjB <= alpha
  
  df
}

# ----------------------------
# Fisher overlap + direction (paper style)
# ----------------------------
paper_style_cooccur <- function(res_C9_filtered, res_SOD1_filtered,
                                gene_col="gene", lfc_col="log2FoldChange", padj_col="padj",
                                alpha=0.05, universe=c("intersection","union")) {
  
  df <- prep_pair(res_C9_filtered, res_SOD1_filtered,
                  gene_col=gene_col, lfc_col=lfc_col, padj_col=padj_col,
                  alpha=alpha, universe=universe)
  
  # (1) Overlap Fisher
  a <- sum(df$DEG_A & df$DEG_B)
  b <- sum(df$DEG_A & !df$DEG_B)
  c <- sum(!df$DEG_A & df$DEG_B)
  d <- sum(!df$DEG_A & !df$DEG_B)
  
  tab_overlap <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE,
                        dimnames=list(A=c("DEG","NotDEG"), B=c("DEG","NotDEG")))
  ft_overlap <- fisher.test(tab_overlap)
  
  # (2) Direction Fisher among overlapping DEGs
  df_int <- df[df$DEG_A & df$DEG_B & !is.na(df$lfcA) & !is.na(df$lfcB), ]
  if (nrow(df_int) == 0) {
    return(list(overlap_table=tab_overlap, overlap_fisher=ft_overlap,
                direction_table=NULL, direction_fisher=NULL,
                n_overlap_deg=0, df_overlap_deg=df_int))
  }
  
  df_int$dirA <- ifelse(df_int$lfcA >= 0, "Up", "Down")
  df_int$dirB <- ifelse(df_int$lfcB >= 0, "Up", "Down")
  
  tab_dir <- table(df_int$dirA, df_int$dirB)
  # ensure consistent 2x2 ordering
  tab_dir <- tab_dir[c("Up","Down"), c("Up","Down")]
  ft_dir <- fisher.test(tab_dir)
  
  list(overlap_table=tab_overlap, overlap_fisher=ft_overlap,
       direction_table=tab_dir, direction_fisher=ft_dir,
       n_overlap_deg=nrow(df_int), df_overlap_deg=df_int)
}

# ----------------------------
# LFC concordance (one function covers the "set" logic)
# ----------------------------
lfc_concordance <- function(res_c9, res_sod1,
                            gene_col="gene", lfc_col="log2FoldChange", padj_col="padj",
                            alpha=0.05,
                            set=c("all","sig_either","overlap_deg","from_C9","from_SOD1"),
                            method="spearman") {
  
  set <- match.arg(set)
  
  df <- prep_pair(res_c9, res_sod1,
                  gene_col=gene_col, lfc_col=lfc_col, padj_col=padj_col,
                  alpha=alpha, universe="intersection")
  
  # Subset choice
  if (set == "sig_either")  df <- df[df$DEG_A | df$DEG_B, ]
  if (set == "overlap_deg") df <- df[df$DEG_A & df$DEG_B, ]
  if (set == "from_C9")     df <- df[df$DEG_A, ]
  if (set == "from_SOD1")   df <- df[df$DEG_B, ]
  
  # Need both LFCs
  df <- df[!is.na(df$lfcA) & !is.na(df$lfcB), ]
  
  if (nrow(df) < 5) {
    return(list(n=nrow(df), cor=NA_real_, sign_agree=NA_real_, df=df))
  }
  
  cor_val <- suppressWarnings(cor(df$lfcA, df$lfcB, method=method))
  sign_agree <- mean(sign(df$lfcA) == sign(df$lfcB))
  
  list(n=nrow(df), cor=cor_val, sign_agree=sign_agree, df=df)
}

# ----------------------------
# Usage (same outputs, less code)
# ----------------------------
out <- paper_style_cooccur(
  res_C9_filtered   = as.data.frame(res_C9_filtered),
  res_SOD1_filtered = as.data.frame(res_SOD1_filtered),
  lfc_col  = "log2FoldChange",
  padj_col = "padj",
  alpha = 0.05,
  universe = "intersection"
)

out_all <- lfc_concordance(res_C9_filtered, res_SOD1_filtered, set="all")
out_sig <- lfc_concordance(res_C9_filtered, res_SOD1_filtered, set="sig_either")
out_ov  <- lfc_concordance(res_C9_filtered, res_SOD1_filtered, set="overlap_deg")
out_setC9   <- lfc_concordance(res_C9_filtered, res_SOD1_filtered, set="from_C9")
out_setSOD1 <- lfc_concordance(res_C9_filtered, res_SOD1_filtered, set="from_SOD1")
# ----------------------------


#Protein coding

library(biomaRt)

all_ids <- unique(c(rownames(res_C9_filtered), rownames(res_SOD1_filtered)))

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

ann_all <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "description"),
  filters    = "ensembl_gene_id",
  values     = all_ids,
  mart       = mart
)

# keep one row per ensembl id (if duplicates, keep first)
ann_all <- ann_all[!duplicated(ann_all$ensembl_gene_id), ]
head(ann_all)

pc_ids <- ann_all$ensembl_gene_id[ann_all$gene_biotype == "protein_coding"]

resC9_pc   <- res_C9_filtered[rownames(res_C9_filtered) %in% pc_ids, , drop=FALSE]
resSOD1_pc <- res_SOD1_filtered[rownames(res_SOD1_filtered) %in% pc_ids, , drop=FALSE]

# paper_style_cooccur_rownames (your earlier function)
out_pc <- paper_style_cooccur(
  res_C9_filtered   = resC9_pc,
  res_SOD1_filtered = resSOD1_pc,
  lfc_col  = "log2FoldChange",
  padj_col = "padj",
  alpha = 0.05
)

out_pc$overlap_table
out_pc$overlap_fisher$estimate
out_pc$overlap_fisher$p.value

out_pc$direction_table
out_pc$direction_fisher$estimate
out_pc$direction_fisher$p.value

out_pc$n_overlap_deg

resC9_pc2 <- resC9_pc;     resC9_pc2$gene <- rownames(resC9_pc2)
resSOD1_pc2 <- resSOD1_pc; resSOD1_pc2$gene <- rownames(resSOD1_pc2)

out_all_pc <- lfc_concordance(resC9_pc2, resSOD1_pc2,
                              gene_col="gene", lfc_col="log2FoldChange", padj_col="padj",
                              alpha=0.05, set="all", method="spearman")

out_setC9_pc <- lfc_concordance_set_from_one(resC9_pc2, resSOD1_pc2,
                                             gene_col="gene",
                                             lfc_col="log2FoldChange", padj_col="padj",
                                             alpha=0.05, define_from="C9")

c(n=out_all_pc$n,   cor=out_all_pc$cor,   sign=out_all_pc$sign_agree)
c(n=out_setC9_pc$n, cor=out_setC9_pc$cor, sign=out_setC9_pc$sign_agree)


als_symbols <- c(
  "C9orf72","SOD1","TARDBP","FUS","TBK1","OPTN","SQSTM1","NEK1","UNC13A",
  "ATXN2","KIF5A","VCP","SETX","TIA1","PFN1","MATR3","SPG11","ANG","ALS2",
  "DCTN1","HNRNPA1","HNRNPA2B1","FIG4","SPAST","UBQLN2","CHCHD10","CCNF"
)
writeLines(als_symbols, "ALS_gene_list.txt")

als_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "description"),
  filters    = "external_gene_name",
  values     = als_symbols,
  mart       = mart
)

# keep protein-coding hits first if duplicates
als_map <- als_map[order(als_map$external_gene_name, als_map$gene_biotype != "protein_coding"), ]
als_map <- als_map[!duplicated(als_map$external_gene_name), ]
als_map

als_direction_table <- function(resC9_df, resSOD1_df, als_map,
                                lfc_col="log2FoldChange", padj_col="padj") {
  
  out <- data.frame(
    symbol  = als_map$external_gene_name,
    ensembl = als_map$ensembl_gene_id,
    biotype = als_map$gene_biotype,
    stringsAsFactors = FALSE
  )
  
  out$lfc_C9    <- resC9_df[out$ensembl, lfc_col]
  out$padj_C9   <- resC9_df[out$ensembl, padj_col]
  out$lfc_SOD1  <- resSOD1_df[out$ensembl, lfc_col]
  out$padj_SOD1 <- resSOD1_df[out$ensembl, padj_col]
  
  out$dir_C9   <- ifelse(is.na(out$lfc_C9), NA, ifelse(out$lfc_C9 >= 0, "Up", "Down"))
  out$dir_SOD1 <- ifelse(is.na(out$lfc_SOD1), NA, ifelse(out$lfc_SOD1 >= 0, "Up", "Down"))
  
  out$concordant_dir <- ifelse(is.na(out$dir_C9) | is.na(out$dir_SOD1),
                               NA, out$dir_C9 == out$dir_SOD1)
  
  out$sig_C9   <- !is.na(out$padj_C9)   & out$padj_C9 <= 0.05
  out$sig_SOD1 <- !is.na(out$padj_SOD1) & out$padj_SOD1 <= 0.05
  out$sig_either <- out$sig_C9 | out$sig_SOD1
  
  # show significant and/or discordant first
  out <- out[order(out$sig_either, is.na(out$concordant_dir), !out$concordant_dir, decreasing=TRUE), ]
  rownames(out) <- NULL
  out
}

als_tab <- als_direction_table(res_C9_filtered, res_SOD1_filtered, als_map)
als_tab

subset_tab <- als_tab[!is.na(als_tab$concordant_dir), ]
table(subset_tab$concordant_dir)


library(readxl)
library(dplyr)

# --------
# 1) Read paper DEGs (Supplemental Table 1A)
# --------
paper_path <- "supplemental_table_1_ddz104.xlsx"

raw <- read_excel(paper_path, sheet = "1A", col_names = FALSE)

tolower_chr <- function(x) tolower(trimws(as.character(x)))

# Find row where first four meaningful headers appear (Gene ID etc.)
hdr_row <- which(apply(raw, 1, function(r) any(tolower_chr(r) == "gene id")))[1]
if (is.na(hdr_row)) stop("Could not find header row containing 'Gene ID'.")

# Get header values from that row
hdr <- as.character(unlist(raw[hdr_row, ]))

# Identify the columns for the DEG table (exact match)
needed <- c("Gene ID", "Gene name", "log2FoldChange", "padj")
keep_cols <- which(!is.na(hdr) & hdr %in% needed)

# Safety check
if (length(keep_cols) != 4) {
  print(hdr)
  stop("Did not find exactly 4 DEG columns (Gene ID/Gene name/log2FoldChange/padj).")
}

# Extract the table below the header row
paper_c9 <- raw[(hdr_row + 1):nrow(raw), keep_cols]
colnames(paper_c9) <- hdr[keep_cols]

# Clean + type cast
paper_c9 <- paper_c9 %>%
  rename(
    ensembl = `Gene ID`,
    symbol_paper = `Gene name`,
    lfc_paper = log2FoldChange,
    padj_paper = padj
  ) %>%
  filter(!is.na(ensembl)) %>%
  mutate(
    ensembl = sub("\\.\\d+$", "", as.character(ensembl)),
    symbol_paper = as.character(symbol_paper),
    lfc_paper = as.numeric(lfc_paper),
    padj_paper = as.numeric(padj_paper)
  )

head(paper_c9)

call_state <- function(lfc, padj, alpha = 0.10, lfc_min = 0.1) {
  if (is.na(lfc) || is.na(padj)) return("NoData")
  if (padj > alpha || abs(lfc) < lfc_min) return("NoChange")
  if (lfc > 0) "Up" else "Down"
}

alpha_paper <- 0.10   # paper often uses FDR<0.1 for NeuroLINCS
alpha_mine  <- 0.05
lfc_min     <- 0.10

als_compare <- als_tab %>%
  mutate(ensembl = sub("\\.\\d+$", "", ensembl)) %>%
  left_join(paper_c9, by = "ensembl") %>%
  mutate(
    state_paper = mapply(call_state, lfc_paper, padj_paper,
                         MoreArgs=list(alpha=alpha_paper, lfc_min=lfc_min)),
    state_mine  = mapply(call_state, lfc_C9, padj_C9,
                         MoreArgs=list(alpha=alpha_mine, lfc_min=lfc_min)),
    agree_dir = (state_paper == state_mine) & state_paper %in% c("Up","Down")
  )

# Summary
table(Paper = als_compare$state_paper, Mine = als_compare$state_mine)

n_comp <- sum(als_compare$state_paper %in% c("Up","Down") &
                als_compare$state_mine  %in% c("Up","Down"), na.rm=TRUE)
n_agree <- sum(als_compare$agree_dir, na.rm=TRUE)

c(n_dir_comparable = n_comp,
  n_dir_agree = n_agree,
  agreement_rate = ifelse(n_comp > 0, n_agree/n_comp, NA))

# Per-gene view
als_compare %>%
  dplyr::select(symbol, ensembl, state_paper, lfc_paper, padj_paper,
         state_mine, lfc_C9, padj_C9, agree_dir) %>%
  arrange(desc(agree_dir), symbol)





# familial score calculation

# ===== BROADER FAMILIAL SIGNATURE =====

# Genes significant in C9 at padj < 0.05
c9_up <- rownames(res_C9_filtered)[which(res_C9_filtered$padj < 0.05 & res_C9_filtered$log2FoldChange > 0)]
c9_down <- rownames(res_C9_filtered)[which(res_C9_filtered$padj < 0.05 & res_C9_filtered$log2FoldChange < 0)]

cat("C9 UP genes:", length(c9_up), "\n")
cat("C9 DOWN genes:", length(c9_down), "\n")

# Filter to those with concordant direction in SOD1 (even if not significant)
# And require at least nominal effect in SOD1 (|log2FC| > 0.2)
c9_up_concordant <- c9_up[
  !is.na(res_SOD1_filtered[c9_up, "log2FoldChange"]) &
    res_SOD1_filtered[c9_up, "log2FoldChange"] > 0.2
]

c9_down_concordant <- c9_down[
  !is.na(res_SOD1_filtered[c9_down, "log2FoldChange"]) &
    res_SOD1_filtered[c9_down, "log2FoldChange"] < -0.2
]

cat("\nBroader signature:\n")
cat("C9 UP + SOD1 concordant (LFC>0.2):", length(c9_up_concordant), "\n")
cat("C9 DOWN + SOD1 concordant (LFC<-0.2):", length(c9_down_concordant), "\n")
cat("Total signature genes:", length(c9_up_concordant) + length(c9_down_concordant), "\n")

# Save this signature
familial_signature_genes <- c(c9_up_concordant, c9_down_concordant)
familial_signature_direction <- c(
  rep("UP", length(c9_up_concordant)), 
  rep("DOWN", length(c9_down_concordant))
)

# Get details
sig_details <- data.frame(
  ensembl_id = familial_signature_genes,
  direction = familial_signature_direction,
  C9_log2FC = res_C9_filtered[familial_signature_genes, "log2FoldChange"],
  C9_padj = res_C9_filtered[familial_signature_genes, "padj"],
  SOD1_log2FC = res_SOD1_filtered[familial_signature_genes, "log2FoldChange"],
  SOD1_padj = res_SOD1_filtered[familial_signature_genes, "padj"]
)

sig_details <- sig_details[order(sig_details$direction, -abs(sig_details$C9_log2FC)), ]

cat("\nTop signature genes:\n")
print(head(sig_details, 20))

write.csv(sig_details, "familial_signature_broad.csv", row.names=FALSE)




library(biomaRt)

strip_ens <- function(x) sub("\\.\\d+$", "", x)

get_biotype <- function(ensembl_ids) {
  mart <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")
  ids <- unique(strip_ens(ensembl_ids))
  ann <- getBM(
    attributes = c("ensembl_gene_id","external_gene_name","gene_biotype"),
    filters = "ensembl_gene_id",
    values = ids,
    mart = mart
  )
  ann <- ann[!duplicated(ann$ensembl_gene_id), ]
  ann
}


make_concordant_signature <- function(resC9, resSOD1,
                                      alpha_c9 = 0.05,
                                      lfc_min = 0.2,
                                      sod1_lfc_min = 0.2,
                                      use_sod1_p = FALSE,
                                      sod1_p_max = 0.3,
                                      protein_coding_only = TRUE) {
  
  resC9 <- as.data.frame(resC9)
  resSOD1 <- as.data.frame(resSOD1)
  rownames(resC9) <- strip_ens(rownames(resC9))
  rownames(resSOD1) <- strip_ens(rownames(resSOD1))
  
  # Align to shared universe
  U <- intersect(rownames(resC9), rownames(resSOD1))
  resC9 <- resC9[U, , drop=FALSE]
  resSOD1 <- resSOD1[U, , drop=FALSE]
  
  # Optional protein-coding restriction
  ann <- NULL
  if (protein_coding_only) {
    ann <- get_biotype(U)
    pc <- ann$ensembl_gene_id[ann$gene_biotype == "protein_coding"]
    U <- intersect(U, pc)
    resC9 <- resC9[U, , drop=FALSE]
    resSOD1 <- resSOD1[U, , drop=FALSE]
  }
  
  # C9 DEG + effect threshold
  keep_c9 <- !is.na(resC9$padj) & resC9$padj <= alpha_c9 & abs(resC9$log2FoldChange) >= lfc_min
  
  # SOD1 directional guard: non-trivial effect + (optional) pvalue filter
  keep_sod1_eff <- !is.na(resSOD1$log2FoldChange) & abs(resSOD1$log2FoldChange) >= sod1_lfc_min
  
  if (use_sod1_p) {
    if (!("pvalue" %in% colnames(resSOD1))) stop("use_sod1_p=TRUE but resSOD1 has no 'pvalue' column")
    keep_sod1_eff <- keep_sod1_eff & !is.na(resSOD1$pvalue) & resSOD1$pvalue <= sod1_p_max
  }
  
  # Concordant direction
  c9_up  <- U[keep_c9 & resC9$log2FoldChange > 0]
  c9_dn  <- U[keep_c9 & resC9$log2FoldChange < 0]
  
  sod1_up <- U[keep_sod1_eff & resSOD1$log2FoldChange > 0]
  sod1_dn <- U[keep_sod1_eff & resSOD1$log2FoldChange < 0]
  
  up_conc <- intersect(c9_up, sod1_up)
  dn_conc <- intersect(c9_dn, sod1_dn)
  
  sig <- data.frame(
    ensembl_id = c(up_conc, dn_conc),
    direction = c(rep("UP", length(up_conc)), rep("DOWN", length(dn_conc))),
    C9_log2FC = resC9[c(up_conc, dn_conc), "log2FoldChange"],
    C9_padj   = resC9[c(up_conc, dn_conc), "padj"],
    SOD1_log2FC = resSOD1[c(up_conc, dn_conc), "log2FoldChange"],
    SOD1_padj   = resSOD1[c(up_conc, dn_conc), "padj"],
    stringsAsFactors = FALSE
  )
  
  if (!is.null(ann)) {
    sig <- merge(sig, ann, by.x="ensembl_id", by.y="ensembl_gene_id", all.x=TRUE)
  }
  
  sig <- sig[order(sig$direction, -abs(sig$C9_log2FC)), ]
  
  list(
    signature_table = sig,
    up = up_conc,
    down = dn_conc,
    params = list(alpha_c9=alpha_c9, lfc_min=lfc_min, sod1_lfc_min=sod1_lfc_min,
                  use_sod1_p=use_sod1_p, sod1_p_max=sod1_p_max,
                  protein_coding_only=protein_coding_only),
    universe_n = nrow(resC9)
  )
}

# Strict, defensible default
sig_strict <- make_concordant_signature(
  resC9 = res_C9_filtered,
  resSOD1 = res_SOD1_filtered,
  alpha_c9 = 0.05,
  lfc_min = 0.2,
  sod1_lfc_min = 0.2,
  use_sod1_p = FALSE,
  protein_coding_only = TRUE
)

cat("STRICT signature: UP", length(sig_strict$up), "DOWN", length(sig_strict$down),
    "TOTAL", nrow(sig_strict$signature_table), "Universe", sig_strict$universe_n, "\n")

write.csv(sig_strict$signature_table, "familial_signature_strict.csv", row.names=FALSE)

# Lenient (for sensitivity)
sig_lenient <- make_concordant_signature(
  resC9 = res_C9_filtered,
  resSOD1 = res_SOD1_filtered,
  alpha_c9 = 0.10,
  lfc_min = 0.2,
  sod1_lfc_min = 0.2,
  use_sod1_p = FALSE,
  protein_coding_only = FALSE
)

cat("LENIENT signature: UP", length(sig_lenient$up), "DOWN", length(sig_lenient$down),
    "TOTAL", nrow(sig_lenient$signature_table), "Universe", sig_lenient$universe_n, "\n")

write.csv(sig_lenient$signature_table, "familial_signature_lenient.csv", row.names=FALSE)

library(DESeq2)

meta <- as.data.frame(annot_als_filtered)
meta <- meta[colnames(counts_filtered), , drop=FALSE]
stopifnot(all(rownames(meta) == colnames(counts_filtered)))

dds0 <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_filtered)),
  colData   = meta,
  design    = ~ 1
)
vsd <- vst(dds0, blind=TRUE)
expr_vst <- assay(vsd)

score_from_sig <- function(expr_mat, up_genes, down_genes) {
  up <- intersect(up_genes, rownames(expr_mat))
  dn <- intersect(down_genes, rownames(expr_mat))
  if (length(up) < 5 || length(dn) < 5) stop("Too few genes present for scoring (check IDs/filters).")
  colMeans(expr_mat[up, , drop=FALSE], na.rm=TRUE) -
    colMeans(expr_mat[dn, , drop=FALSE], na.rm=TRUE)
}

meta$guardrail_score_strict  <- score_from_sig(expr_vst, sig_strict$up,  sig_strict$down)
meta$guardrail_score_lenient <- score_from_sig(expr_vst, sig_lenient$up, sig_lenient$down)

summary(meta$guardrail_score_strict)

df_prog <- meta[!is.na(meta$ALSFRS_R_PROGRESSION_SLOPE) & !is.na(meta$guardrail_score_strict), ]
df_prog$SEX <- factor(df_prog$SEX)

fit_prog <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ guardrail_score_strict + AGE_AT_SYMPTOM_ONSET + SEX,
               data = df_prog)
summary(fit_prog)

df_prog$Site_of_Onset <- factor(df_prog$Site_of_Onset)

fit_prog2 <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ guardrail_score_strict + AGE_AT_SYMPTOM_ONSET + SEX + Site_of_Onset,
                data = df_prog)
summary(fit_prog2)

library(survival)

df_surv <- meta[!is.na(meta$AGE_AT_SYMPTOM_ONSET) & !is.na(meta$guardrail_score_strict), ]
df_surv$time <- df_surv$AGE_AT_DEATH - df_surv$AGE_AT_SYMPTOM_ONSET
df_surv$status <- ifelse(is.na(df_surv$AGE_AT_DEATH), 0, 1)

df_surv <- df_surv[!is.na(df_surv$time) & df_surv$time > 0, ]
df_surv$SEX <- factor(df_surv$SEX)

cox_fit <- coxph(Surv(time, status) ~ guardrail_score_strict + AGE_AT_SYMPTOM_ONSET + SEX, data=df_surv)
summary(cox_fit)

# group should be in meta already (Control/C9orf72/SOD1/Sporadic)
table(meta$group)

tapply(meta$guardrail_score_strict, meta$group, summary)

# Nonparametric group differences
kruskal.test(guardrail_score_strict ~ group, data=meta)

# Pairwise comparisons vs Control
pairwise.wilcox.test(meta$guardrail_score_strict, meta$group, p.adjust.method="BH")

library(pROC)

sub <- meta[meta$group %in% c("C9orf72","Control") & !is.na(meta$guardrail_score_strict), ]
roc_obj <- roc(response = sub$group, predictor = sub$guardrail_score_strict, levels=c("Control","C9orf72"))
auc(roc_obj)

library(dplyr)

df_sp <- meta %>%
  filter(group == "Sporadic") %>%
  filter(!is.na(guardrail_score_strict))

q25 <- quantile(df_sp$guardrail_score_strict, 0.25)
q75 <- quantile(df_sp$guardrail_score_strict, 0.75)

df_sp$guardrail_bin <- ifelse(df_sp$guardrail_score_strict >= q75, "High",
                              ifelse(df_sp$guardrail_score_strict <= q25, "Low", NA))

table(df_sp$guardrail_bin, useNA="ifany")

# Keep only samples with complete covariates used in the design
keep_complete <- !is.na(meta_sp$AGE_AT_SYMPTOM_ONSET) &
  !is.na(meta_sp$SEX) &
  !is.na(meta_sp$guardrail_bin)

table(keep_complete)

meta_sp2 <- meta_sp[keep_complete, , drop=FALSE]
counts_sp2 <- counts_sp[, rownames(meta_sp2), drop=FALSE]

# Make sure factor levels are clean
# Ensure correct types
meta_sp2$SEX <- factor(meta_sp2$SEX)
meta_sp2$guardrail_bin <- factor(meta_sp2$guardrail_bin, levels=c("Low","High"))

# Center (and optionally scale) age to reduce collinearity
meta_sp2$AGE_ONSET_c <- scale(meta_sp2$AGE_AT_SYMPTOM_ONSET, center=TRUE, scale=FALSE)[,1]
# or if you prefer full z-score scaling:
# meta_sp2$AGE_ONSET_z <- scale(meta_sp2$AGE_AT_SYMPTOM_ONSET, center=TRUE, scale=TRUE)[,1]

dds_sp <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_sp2)),
  colData   = meta_sp2,
  design    = ~ SEX + AGE_ONSET_c + guardrail_bin
)

dds_sp0 <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_sp2)),
  colData   = meta_sp2,
  design    = ~ SEX + guardrail_bin
)
dds_sp0 <- DESeq(dds_sp0)
res_sp_guardrail0 <- results(dds_sp0, contrast=c("guardrail_bin","High","Low"))
summary(res_sp_guardrail0)


res_sp_df <- as.data.frame(res_sp_guardrail0)
res_sp_df$ensembl_id <- rownames(res_sp_df)
res_sp_df <- res_sp_df[order(res_sp_df$padj), ]

write.csv(res_sp_df, "sporadic_guardrail_High_vs_Low_DE.csv", row.names=FALSE)

library(clusterProfiler)
library(org.Hs.eg.db)

# rank by signed statistic (or log2FC)
ranks <- res_sp_df$stat
names(ranks) <- strip_ens(res_sp_df$ensembl_id)
ranks <- ranks[!is.na(ranks)]
ranks <- sort(ranks, decreasing=TRUE)

# convert Ensembl -> Entrez
conv <- bitr(names(ranks), fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
ranks2 <- ranks[conv$ENSEMBL]
names(ranks2) <- conv$ENTREZID
library(clusterProfiler)
library(org.Hs.eg.db)

# ranks: named vector (ENSEMBL -> stat), sorted decreasing
# conv: data.frame from bitr(ENSEMBL -> ENTREZID)

# 1) attach rank to mapping
conv$rank <- ranks[conv$ENSEMBL]

# 2) drop unmapped / missing
conv <- conv[!is.na(conv$ENTREZID) & !is.na(conv$rank), ]

# 3) if multiple ENSEMBL map to same ENTREZID, keep the strongest signal
# (largest absolute rank is a good default)
conv <- conv[order(abs(conv$rank), decreasing = TRUE), ]
conv <- conv[!duplicated(conv$ENTREZID), ]

# 4) build ENTREZ-ranked vector
ranks2 <- conv$rank
names(ranks2) <- conv$ENTREZID
ranks2 <- sort(ranks2, decreasing = TRUE)

stopifnot(!anyDuplicated(names(ranks2)))
length(ranks2)


gsea_go <- gseGO(
  geneList = ranks2,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05
)

# top pathways
head(as.data.frame(gsea_go), 20)
write.csv(as.data.frame(gsea_go), "GSEA_GO_BP_sporadic_guardrail.csv", row.names=FALSE)

meta0 <- as.data.frame(annot_als_filtered)

# Keep sporadic with genomics available (adjust if HAS_GENOMICS is "YES"/1/etc.)
sp <- meta0[meta0$group=="Sporadic" & toupper(meta0$HAS_GENOMICS)=="YES", , drop=FALSE]

# Build VCF sample IDs
sp$vcf_sample_id <- paste0(sp$NYGC_CGND_ID, "-b38")

# Drop missing IDs
sp <- sp[!is.na(sp$NYGC_CGND_ID) & sp$NYGC_CGND_ID != "", , drop=FALSE]

# Write sample list for bcftools
writeLines(unique(sp$vcf_sample_id), "sporadic_samples_vcf.txt")

length(unique(sp$vcf_sample_id))
head(unique(sp$vcf_sample_id), 20)



vcf_samples <- readLines("vcf_samples.txt")

strip_b38 <- function(x) sub("-b38$", "", x)

vcf_df <- data.frame(
  vcf_sample = vcf_samples,
  base_id = strip_b38(vcf_samples),
  stringsAsFactors = FALSE
)

# sanity: how many have -b38?
table(grepl("-b38$", vcf_df$vcf_sample))
head(vcf_df, 20)

meta0 <- as.data.frame(annot_als_filtered)

# pick sporadic with genomics
sp <- meta0[meta0$group=="Sporadic" & toupper(meta0$HAS_GENOMICS)=="YES", , drop=FALSE]
sp <- sp[!is.na(sp$NYGC_CGND_ID) & sp$NYGC_CGND_ID != "", , drop=FALSE]

sp_base <- unique(sp$NYGC_CGND_ID)   # base IDs (no -b38)
length(sp_base)
head(sp_base, 20)

# Keep only VCF samples whose base_id is in sporadic set
hit <- vcf_df[vcf_df$base_id %in% sp_base, ]

# If both "CGND-HDA-00001" and "CGND-HDA-00001-b38" exist, keep one.
# Prefer the -b38 version when available.
hit <- hit[order(hit$base_id, !grepl("-b38$", hit$vcf_sample)), ]  # TRUE sorts after FALSE, so invert
# Better: explicitly prioritize -b38
hit$priority <- ifelse(grepl("-b38$", hit$vcf_sample), 1, 0)
hit <- hit[order(hit$base_id, -hit$priority), ]
hit <- hit[!duplicated(hit$base_id), ]

# Final sample list to feed bcftools
sporadic_vcf_samples <- hit$vcf_sample
writeLines(sporadic_vcf_samples, "sporadic_samples_present.txt")

length(sporadic_vcf_samples)
head(sporadic_vcf_samples, 20)




tab_bg <- table(annot_als_filtered$Batch.x, annot_als_filtered$group)
addmargins(tab_bg)
batch_n_groups <- apply(tab_bg > 0, 1, sum)
table(batch_n_groups)
problem_batches <- names(batch_n_groups[batch_n_groups == 1])
length(problem_batches)
head(problem_batches, 20)

batch_sizes <- rowSums(tab_bg)
table(batch_sizes)
tiny_batches <- names(batch_sizes[batch_sizes < 3])
length(tiny_batches)
head(tiny_batches, 20)


marker_cols <- c("PCT_SMI32","PCT_ISL1","PCT_NKX61","PCT_TUJ1","PCT_S100b","PCT_Nestin")

# keep only samples with complete marker data
keep <- complete.cases(annot_als_filtered[, marker_cols])

annot2 <- annot_als_filtered[keep, ]
counts2 <- counts_filtered[, rownames(annot2)]

# PCA on scaled marker %s
M <- scale(annot2[, marker_cols])
pc <- prcomp(M, center=TRUE, scale.=TRUE)

annot2$MarkerPC1 <- pc$x[,1]
annot2$MarkerPC2 <- pc$x[,2]


drop_batches <- c("51","52","55")
keep_batch <- !(annot_als_filtered$Batch %in% drop_batches)

annotB <- annot_als_filtered[keep_batch, ]
countsB <- counts_filtered[, rownames(annotB)]

marker_cols <- c("PCT_SMI32","PCT_ISL1","PCT_NKX61","PCT_TUJ1","PCT_S100b","PCT_Nestin")

annotM <- annotB  # from earlier step (after optional batch filtering)

# median-impute missing marker %s
for (cc in marker_cols) {
  v <- annotM[[cc]]
  med <- median(v, na.rm=TRUE)
  annotM[[cc]][is.na(v)] <- med
}

# build marker matrix
M0 <- annotM[, marker_cols]
M0 <- as.data.frame(M0)

# coerce to numeric (important if any are character)
for (cc in colnames(M0)) M0[[cc]] <- as.numeric(M0[[cc]])

# identify and drop zero-variance columns
sds <- sapply(M0, sd, na.rm=TRUE)
zero_cols <- names(sds)[sds == 0 | is.na(sds)]

if (length(zero_cols) > 0) {
  cat("Dropping zero-variance marker columns:", paste(zero_cols, collapse=", "), "\n")
}

marker_cols_use <- setdiff(marker_cols, zero_cols)
M <- scale(M0[, marker_cols_use, drop=FALSE])

pc <- prcomp(M, center=TRUE, scale.=TRUE)

annotM$MarkerPC1 <- pc$x[,1]
annotM$MarkerPC2 <- pc$x[,2]

marker_cols <- c("PCT_SMI32","PCT_ISL1","PCT_NKX61","PCT_TUJ1","PCT_S100b","PCT_Nestin")

colSums(is.na(annot_als_filtered[, marker_cols]))

sapply(annot_als_filtered[, marker_cols], class)
head(annot_als_filtered[, marker_cols])

tmp <- annot_als_filtered[, marker_cols]
tmp2 <- as.data.frame(lapply(tmp, function(x) as.numeric(as.character(x))))
colSums(is.na(tmp2))
summary(tmp2)

marker_cols_all <- c("PCT_SMI32","PCT_ISL1","PCT_NKX61","PCT_TUJ1","PCT_S100b","PCT_Nestin")
marker_cols_use <- c("PCT_SMI32","PCT_ISL1","PCT_NKX61","PCT_TUJ1","PCT_S100b")  # drop Nestin

annotM <- annot_als_filtered  # start from the full metadata that you confirmed has variation

# optional: drop singleton batches 51/52/55 if you want
drop_batches <- c("51","52","55")
annotM <- annotM[!(annotM$Batch.x %in% drop_batches), ]
countsM <- counts_filtered[, rownames(annotM)]

# impute missing values for the markers we keep
for (cc in marker_cols_use) {
  v <- annotM[[cc]]
  annotM[[cc]][is.na(v)] <- median(v, na.rm=TRUE)
}

# build scaled matrix and drop any remaining zero-variance markers
M0 <- as.data.frame(annotM[, marker_cols_use])
sds <- sapply(M0, sd, na.rm=TRUE)
zero_cols <- names(sds)[sds == 0 | is.na(sds)]
if (length(zero_cols) > 0) {
  cat("Dropping zero-variance marker columns:", paste(zero_cols, collapse=", "), "\n")
}
marker_cols_final <- setdiff(marker_cols_use, zero_cols)

M <- scale(M0[, marker_cols_final, drop=FALSE])
pc <- prcomp(M, center=TRUE, scale.=TRUE)

annotM$MarkerPC1 <- pc$x[,1]
annotM$MarkerPC2 <- pc$x[,2]

# sanity check
summary(annotM$MarkerPC1)
summary(annotM$MarkerPC2)


# Compare sex, batch, marker PCs

library(sva)

# use the same annotM + countsM you just built (has MarkerPC1/2 and batch_col)
mod  <- model.matrix(as.formula(paste("~ SEX +", "Batch.x", "+ MarkerPC1 + MarkerPC2 + group")), data=annotM)
mod0 <- model.matrix(as.formula(paste("~ SEX +", "Batch.x", "+ MarkerPC1 + MarkerPC2")), data=annotM)

# expr_norm should be log-scale expression (genes x samples) for sva
# Use your VST expression:

library(DESeq2)

stopifnot(all(colnames(countsM) == rownames(annotM)))

dds0 <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(countsM)),
  colData   = as.data.frame(annotM),
  design    = ~ 1
)

vsd <- vst(dds0, blind=TRUE)
expr_vst <- assay(vsd)

svobj <- sva(expr_vst, mod, mod0, n.sv = 10)

sv_group_p <- sapply(1:ncol(svobj$sv), function(i) {
  summary(aov(svobj$sv[, i] ~ annotM$group))[[1]]$`Pr(>F)`[1]
})

clean_svs <- which(sv_group_p >= 0.05)
conf_svs  <- which(sv_group_p < 0.05)

cat("Confounded SVs:", conf_svs, "\n")
cat("Clean SVs:", clean_svs, "\n")


annot_sva2 <- annotM

for (j in seq_along(clean_svs)) {
  annot_sva2[[paste0("SV", j)]] <- svobj$sv[, clean_svs[j]]
}

# sanity
head(annot_sva2[, paste0("SV", 1:length(clean_svs))])

library(DESeq2)

annot_sva2$SEX <- factor(annot_sva2$SEX)
annot_sva2$group <- factor(annot_sva2$group)
annot_sva2[["Batch.x"]] <- factor(annot_sva2[["Batch.x"]])

sv_terms <- paste0("SV", 1:length(clean_svs), collapse=" + ")

design_sva <- as.formula(
  paste("~ SEX +", "Batch.x", "+ MarkerPC1 + MarkerPC2 +", sv_terms, "+ group")
)

# full-rank guard
X <- model.matrix(design_sva, data=annot_sva2)
stopifnot(qr(X)$rank == ncol(X))
annot_sva2$group <- relevel(annot_sva2$group, ref="Control")

dds_sva <- DESeqDataSetFromMatrix(
  countData = countsM,
  colData   = annot_sva2,
  design    = design_sva
)
dds_sva <- DESeq(dds_sva)

res_C9_sva   <- results(dds_sva, contrast=c("group","C9orf72","Control"))
res_SOD1_sva <- results(dds_sva, contrast=c("group","SOD1","Control"))

sig_sva_strict <- make_concordant_signature(
  resC9 = res_C9_sva,
  resSOD1 = res_SOD1_sva,
  alpha_c9 = 0.05,
  lfc_min = 0.2,
  sod1_lfc_min = 0.2,
  use_sod1_p = FALSE,
  protein_coding_only = TRUE
)

cat("SVA+markers STRICT signature:",
    "UP", length(sig_sva_strict$up),
    "DOWN", length(sig_sva_strict$down),
    "TOTAL", nrow(sig_sva_strict$signature_table), "\n")

write.csv(sig_sva_strict$signature_table, "familial_signature_strict_sva_markers.csv", row.names=FALSE)


annot_sva2$familial_score_sva <- score_from_sig(expr_vst, sig_sva_strict$up, sig_sva_strict$down)
summary(annot_sva2$familial_score_sva)

df_prog <- annot_sva2[!is.na(annot_sva2$ALSFRS_R_PROGRESSION_SLOPE) & !is.na(annot_sva2$familial_score_sva), ]
df_prog$Batch_use <- df_prog[[batch_col]]

fit_prog <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ familial_score_sva + AGE_AT_SYMPTOM_ONSET + SEX + MarkerPC1 + MarkerPC2 + Batch_use,
               data=df_prog)
summary(fit_prog)

# =========================
# 0) Libraries
# =========================
library(DESeq2)
library(sva)

# =========================
# 1) Helpers
# =========================
strip_ens <- function(x) sub("\\.\\d+$", "", x)

score_from_sig <- function(expr_mat, up_genes, down_genes) {
  up <- intersect(up_genes, rownames(expr_mat))
  dn <- intersect(down_genes, rownames(expr_mat))
  if (length(up) < 5 || length(dn) < 5) stop("Too few genes present for scoring (check IDs/filters).")
  colMeans(expr_mat[up, , drop=FALSE], na.rm=TRUE) -
    colMeans(expr_mat[dn, , drop=FALSE], na.rm=TRUE)
}

align_res <- function(resA, resB, lfc_col="log2FoldChange", padj_col="padj") {
  A <- as.data.frame(resA); B <- as.data.frame(resB)
  A$gene <- strip_ens(rownames(A))
  B$gene <- strip_ens(rownames(B))
  A <- A[, c("gene", lfc_col, padj_col)]
  B <- B[, c("gene", lfc_col, padj_col)]
  colnames(A) <- c("gene","lfcA","padjA")
  colnames(B) <- c("gene","lfcB","padjB")
  U <- intersect(A$gene, B$gene)
  A <- A[match(U, A$gene), ]
  B <- B[match(U, B$gene), ]
  data.frame(gene=U, lfcA=A$lfcA, padjA=A$padjA, lfcB=B$lfcB, padjB=B$padjB)
}

lfc_concordance_simple <- function(res_c9, res_sod1, alpha=0.05, set=c("all","sig_either","overlap_deg")) {
  set <- match.arg(set)
  df <- align_res(res_c9, res_sod1)
  df$DEG_A <- !is.na(df$padjA) & df$padjA <= alpha
  df$DEG_B <- !is.na(df$padjB) & df$padjB <= alpha
  if (set=="sig_either") df <- df[df$DEG_A | df$DEG_B, ]
  if (set=="overlap_deg") df <- df[df$DEG_A & df$DEG_B, ]
  df <- df[!is.na(df$lfcA) & !is.na(df$lfcB), ]
  if (nrow(df) < 5) return(list(n=nrow(df), cor=NA_real_, sign_agree=NA_real_))
  list(
    n = nrow(df),
    cor = cor(df$lfcA, df$lfcB, method="spearman"),
    sign_agree = mean(sign(df$lfcA) == sign(df$lfcB))
  )
}

# --- your existing signature builder (keep as-is if already defined in your env) ---
# expects: make_concordant_signature(resC9,resSOD1, alpha_c9, lfc_min, sod1_lfc_min, use_sod1_p, protein_coding_only)

# =========================
# 2) Inputs assumed to exist
# =========================
# counts_filtered: genes x samples
# annot_als_filtered: rownames = sample IDs, includes group, SEX, Batch.x, markers
# res_C9_filtered, res_SOD1_filtered: your OLD model results already computed

stopifnot(all(colnames(counts_filtered) %in% rownames(annot_als_filtered)))

# =========================
# 3) Build annotM + countsM (Batch filter + Marker PCs)
# =========================
marker_cols_use <- c("PCT_SMI32","PCT_ISL1","PCT_NKX61","PCT_TUJ1","PCT_S100b")  # drop Nestin

annotM <- annot_als_filtered

# pick batch column
batch_col <- if ("Batch" %in% colnames(annotM)) "Batch" else if ("Batch.x" %in% colnames(annotM)) "Batch.x" else stop("No Batch column found")
drop_batches <- c("51","52","55")
b <- as.character(annotM[[batch_col]])
keep_rows <- is.na(b) | !(b %in% drop_batches)

annotM <- annotM[keep_rows, , drop=FALSE]
countsM <- counts_filtered[, rownames(annotM), drop=FALSE]
stopifnot(all(colnames(countsM) == rownames(annotM)))

# impute marker NAs (median)
for (cc in marker_cols_use) annotM[[cc]][is.na(annotM[[cc]])] <- median(annotM[[cc]], na.rm=TRUE)

# PCA -> MarkerPC1/2
M <- scale(as.data.frame(annotM[, marker_cols_use, drop=FALSE]))
pc <- prcomp(M, center=TRUE, scale.=TRUE)
annotM$MarkerPC1 <- pc$x[,1]
annotM$MarkerPC2 <- pc$x[,2]

# factors + set Control as reference (do BEFORE DESeq objects)
annotM$SEX <- factor(annotM$SEX)
annotM$group <- relevel(factor(annotM$group), ref="Control")
annotM[[batch_col]] <- droplevels(factor(annotM[[batch_col]]))

# =========================
# 4) VST expression for SVA + scoring
# =========================
dds0 <- DESeqDataSetFromMatrix(countData = round(as.matrix(countsM)), colData = as.data.frame(annotM), design = ~ 1)
vsd <- vst(dds0, blind=TRUE)
expr_vst <- assay(vsd)

# =========================
# 5) SVA on top of known covariates
# =========================
mod  <- model.matrix(as.formula(paste("~ SEX +", batch_col, "+ MarkerPC1 + MarkerPC2 + group")), data=annotM)
mod0 <- model.matrix(as.formula(paste("~ SEX +", batch_col, "+ MarkerPC1 + MarkerPC2")), data=annotM)

stopifnot(qr(mod)$rank  == ncol(mod))
stopifnot(qr(mod0)$rank == ncol(mod0))

svobj <- sva(expr_vst, mod, mod0, n.sv=10)

# keep only SVs NOT associated with group
sv_group_p <- sapply(1:ncol(svobj$sv), function(i) summary(aov(svobj$sv[, i] ~ annotM$group))[[1]]$`Pr(>F)`[1])
clean_svs <- which(sv_group_p >= 0.05)
conf_svs  <- which(sv_group_p < 0.05)
cat("Confounded SVs:", conf_svs, "\n")
cat("Clean SVs:", clean_svs, "\n")

annot_sva2 <- annotM
for (j in seq_along(clean_svs)) annot_sva2[[paste0("SV", j)]] <- svobj$sv[, clean_svs[j]]

# =========================
# 6) DESeq2 NEW model: Batch + MarkerPCs + clean SVs + group
# =========================
sv_terms <- if (length(clean_svs) > 0) paste0("SV", seq_along(clean_svs), collapse=" + ") else NULL
design_sva <- if (!is.null(sv_terms)) {
  as.formula(paste("~ SEX +", batch_col, "+ MarkerPC1 + MarkerPC2 +", sv_terms, "+ group"))
} else {
  as.formula(paste("~ SEX +", batch_col, "+ MarkerPC1 + MarkerPC2 + group"))
}

X <- model.matrix(design_sva, data=annot_sva2)
stopifnot(qr(X)$rank == ncol(X))

dds_sva <- DESeqDataSetFromMatrix(countData = countsM, colData = as.data.frame(annot_sva2), design = design_sva)
dds_sva <- DESeq(dds_sva)

res_C9_sva   <- results(dds_sva, contrast=c("group","C9orf72","Control"))
res_SOD1_sva <- results(dds_sva, contrast=c("group","SOD1","Control"))

# =========================
# 7) Summaries: old vs new concordance + signature overlap
# =========================
alpha <- 0.05

cat("\n--- OLD model concordance (C9 vs SOD1) ---\n")
print(lfc_concordance_simple(res_C9_filtered, res_SOD1_filtered, alpha=alpha, set="all"))
print(lfc_concordance_simple(res_C9_filtered, res_SOD1_filtered, alpha=alpha, set="sig_either"))

cat("\n--- NEW model concordance (C9 vs SOD1) ---\n")
print(lfc_concordance_simple(res_C9_sva, res_SOD1_sva, alpha=alpha, set="all"))
print(lfc_concordance_simple(res_C9_sva, res_SOD1_sva, alpha=alpha, set="sig_either"))

cat("\n--- OLD vs NEW (C9) LFC stability ---\n")
print(lfc_concordance_simple(res_C9_filtered, res_C9_sva, alpha=alpha, set="all"))
print(lfc_concordance_simple(res_C9_filtered, res_C9_sva, alpha=alpha, set="sig_either"))

cat("\n--- OLD vs NEW (SOD1) LFC stability ---\n")
print(lfc_concordance_simple(res_SOD1_filtered, res_SOD1_sva, alpha=alpha, set="all"))
print(lfc_concordance_simple(res_SOD1_filtered, res_SOD1_sva, alpha=alpha, set="sig_either"))

# signatures (requires make_concordant_signature defined)
sig_old <- make_concordant_signature(res_C9_filtered, res_SOD1_filtered,
                                     alpha_c9=0.05, lfc_min=0.2, sod1_lfc_min=0.2,
                                     use_sod1_p=FALSE, protein_coding_only=TRUE)

sig_new <- make_concordant_signature(res_C9_sva, res_SOD1_sva,
                                     alpha_c9=0.05, lfc_min=0.2, sod1_lfc_min=0.2,
                                     use_sod1_p=FALSE, protein_coding_only=TRUE)

old_genes <- sig_old$signature_table$ensembl_id
new_genes <- sig_new$signature_table$ensembl_id

cat("\n--- Signature comparison ---\n")
cat("OLD signature:", length(old_genes), "NEW signature:", length(new_genes),
    "Overlap:", length(intersect(old_genes, new_genes)), "\n")





# sparse mutations
library(dplyr)

# -----------------------
# Inputs
# -----------------------
car_file <- "allchr.ALS_PTVsplice.gene_carriers.long.tsv"
wgs_list <- "sporadic_samples_present.txt"

# meta is already in memory and has:
# - Participant_ID
# - NYGC_CGND_ID
# - guardrail_score_strict
# - group labels column (Control/C9/SOD1/Sporadic) either in meta$group or similar

# -----------------------
# Helper: detect / standardize group column
# -----------------------
if (!("group" %in% names(meta))) {
  expected <- c("Control","Sporadic","C9","SOD1")
  cand <- names(meta)[sapply(meta, function(x) any(x %in% expected))]
  if (length(cand) == 0) {
    stop("Could not find a group column. Please set meta$group to Control/C9/SOD1/Sporadic.")
  }
  message("Using group column: ", cand[1])
  meta <- meta %>% rename(group = all_of(cand[1]))
}

# -----------------------
# Build CASE_ID and attach guardrail
# -----------------------
meta$CASE_ID <- paste0("CASE-", sub("^[^-]+-", "", meta$Participant_ID))
meta$guardrail <- meta$guardrail_score_strict

# -----------------------
# Restrict to WGS-present cohort
# -----------------------
wgs <- sub("-b38$", "", readLines(wgs_list))
meta$CGND_BASE <- sub("-b38$", "", meta$NYGC_CGND_ID)
meta2 <- meta %>% filter(CGND_BASE %in% wgs)

# -----------------------
# Load carriers and keep those present in meta2
# -----------------------
car_all <- read.delim(car_file, stringsAsFactors=FALSE) %>%
  filter(sample %in% meta2$CASE_ID)

# Make two versions: with and without ATXN2
car_with <- car_all
car_no   <- car_all %>% filter(gene != "ATXN2")

# -----------------------
# Core analysis function
# -----------------------
run_all <- function(car2, tag, min_n_gene=3) {
  
  # ---- Table A: sparse carriers vs controls (guardrail)
  controls <- meta2 %>% filter(group=="Control") %>% pull(CASE_ID)
  
  gene_counts <- car2 %>%
    group_by(gene) %>%
    summarize(n_carriers = n_distinct(sample), .groups="drop") %>%
    arrange(desc(n_carriers))
  
  genes_to_test <- gene_counts %>% filter(n_carriers >= min_n_gene) %>% pull(gene)
  
  sparse_vs_ctrl <- lapply(genes_to_test, function(g){
    carriers <- unique(car2$sample[car2$gene==g])
    wt <- wilcox.test(meta2$guardrail[meta2$CASE_ID %in% carriers],
                      meta2$guardrail[meta2$CASE_ID %in% controls])
    data.frame(
      gene=g,
      n_carriers=length(carriers),
      median_carrier=median(meta2$guardrail[meta2$CASE_ID %in% carriers], na.rm=TRUE),
      median_ctrl=median(meta2$guardrail[meta2$CASE_ID %in% controls], na.rm=TRUE),
      delta=median_carrier - median_ctrl,
      p=wt$p.value
    )
  }) |> bind_rows()
  
  if (nrow(sparse_vs_ctrl) > 0) {
    sparse_vs_ctrl <- sparse_vs_ctrl %>%
      mutate(padj = p.adjust(p, "BH")) %>%
      arrange(padj, p)
  }
  
  write.table(sparse_vs_ctrl, paste0("A_sparse_vs_control_guardrail_", tag, ".tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  
  # ---- Table B: sporadic burden vs guardrail
  sp <- meta2 %>% filter(group=="Sporadic", !is.na(guardrail))
  
  burden <- car2 %>%
    group_by(sample) %>%
    summarize(ptvsplice_burden = n_distinct(gene), .groups="drop")
  
  sp <- sp %>% left_join(burden, by=c("CASE_ID"="sample"))
  sp$ptvsplice_burden[is.na(sp$ptvsplice_burden)] <- 0
  
  rho <- suppressWarnings(cor.test(sp$ptvsplice_burden, sp$guardrail, method="spearman"))
  
  q <- quantile(sp$guardrail, probs=c(0.25,0.75), na.rm=TRUE)
  sp2 <- sp %>%
    mutate(bin = case_when(
      guardrail >= q[2] ~ "High",
      guardrail <= q[1] ~ "Low",
      TRUE ~ NA_character_
    )) %>% filter(!is.na(bin))
  
  wt_burden <- wilcox.test(ptvsplice_burden ~ bin, data=sp2)
  
  sporadic_burden_summary <- data.frame(
    tag=tag,
    n_sporadic=nrow(sp),
    spearman_rho=unname(rho$estimate),
    spearman_p=rho$p.value,
    high_n=sum(sp2$bin=="High"),
    low_n=sum(sp2$bin=="Low"),
    high_median_burden=median(sp2$ptvsplice_burden[sp2$bin=="High"]),
    low_median_burden=median(sp2$ptvsplice_burden[sp2$bin=="Low"]),
    wilcox_p=wt_burden$p.value
  )
  
  write.table(sporadic_burden_summary, paste0("B_sporadic_burden_guardrail_", tag, ".tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  
  # ---- Table C: gene enrichment High vs Low in sporadics (Fisher)
  genes <- sort(unique(car2$gene))
  
  enrich <- lapply(genes, function(g){
    carriers <- unique(car2$sample[car2$gene==g])
    high <- sp2$CASE_ID[sp2$bin=="High"]
    low  <- sp2$CASE_ID[sp2$bin=="Low"]
    
    a <- sum(high %in% carriers)
    b <- sum(!(high %in% carriers))
    c <- sum(low  %in% carriers)
    d <- sum(!(low  %in% carriers))
    
    ft <- fisher.test(matrix(c(a,b,c,d), nrow=2))
    data.frame(gene=g, high_carriers=a, low_carriers=c, OR=unname(ft$estimate), p=ft$p.value)
  }) |> bind_rows() %>%
    mutate(padj = p.adjust(p, "BH")) %>%
    arrange(padj, p)
  
  write.table(enrich, paste0("C_sporadic_highlow_gene_enrichment_", tag, ".tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  
  list(gene_counts=gene_counts,
       sparse_vs_ctrl=sparse_vs_ctrl,
       sporadic_burden=sporadic_burden_summary,
       enrich=enrich)
}

# -----------------------
# Run both analyses
# -----------------------
out_with <- run_all(car_with, "withATXN2")
out_no   <- run_all(car_no,   "noATXN2")

# Print quick summaries
print(out_with$gene_counts)
print(out_with$sporadic_burden)
print(head(out_with$enrich, 10))

print(out_no$gene_counts)
print(out_no$sporadic_burden)
print(head(out_no$enrich, 10))

run_sporadic_only <- function(car2, tag) {
  
  # Sporadic only (WGS-present)
  sp <- meta2 %>% filter(group=="Sporadic", !is.na(guardrail))
  
  # Burden = distinct ALS genes with >=1 PTV/splice
  burden <- car2 %>%
    group_by(sample) %>%
    summarize(ptvsplice_burden = n_distinct(gene), .groups="drop")
  
  sp <- sp %>% left_join(burden, by=c("CASE_ID"="sample"))
  sp$ptvsplice_burden[is.na(sp$ptvsplice_burden)] <- 0
  
  # Spearman burden vs guardrail
  rho <- suppressWarnings(cor.test(sp$ptvsplice_burden, sp$guardrail, method="spearman"))
  
  # High/Low bins (quartiles)
  q <- quantile(sp$guardrail, probs=c(0.25,0.75), na.rm=TRUE)
  sp2 <- sp %>%
    mutate(bin = case_when(
      guardrail >= q[2] ~ "High",
      guardrail <= q[1] ~ "Low",
      TRUE ~ NA_character_
    )) %>% filter(!is.na(bin))
  
  wt_burden <- wilcox.test(ptvsplice_burden ~ bin, data=sp2)
  
  sporadic_burden_summary <- data.frame(
    tag=tag,
    n_sporadic=nrow(sp),
    spearman_rho=unname(rho$estimate),
    spearman_p=rho$p.value,
    high_n=sum(sp2$bin=="High"),
    low_n=sum(sp2$bin=="Low"),
    high_median_burden=median(sp2$ptvsplice_burden[sp2$bin=="High"]),
    low_median_burden=median(sp2$ptvsplice_burden[sp2$bin=="Low"]),
    wilcox_p=wt_burden$p.value
  )
  
  write.table(sporadic_burden_summary, paste0("B_sporadic_burden_guardrail_", tag, ".tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  
  # Gene-level enrichment High vs Low (Fisher)
  genes <- sort(unique(car2$gene))
  enrich <- lapply(genes, function(g){
    carriers <- unique(car2$sample[car2$gene == g])
    high <- sp2$CASE_ID[sp2$bin=="High"]
    low  <- sp2$CASE_ID[sp2$bin=="Low"]
    
    a <- sum(high %in% carriers)
    b <- sum(!(high %in% carriers))
    c <- sum(low  %in% carriers)
    d <- sum(!(low  %in% carriers))
    
    ft <- fisher.test(matrix(c(a,b,c,d), nrow=2))
    data.frame(gene=g, high_carriers=a, low_carriers=c, OR=unname(ft$estimate), p=ft$p.value)
  }) |> bind_rows() %>%
    mutate(padj = p.adjust(p, "BH")) %>%
    arrange(padj, p)
  
  write.table(enrich, paste0("C_sporadic_highlow_gene_enrichment_", tag, ".tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  
  list(sporadic_burden=sporadic_burden_summary, enrich=enrich)
  # Using sp2 from your run_sporadic_only() function:
  # Add at end of run_sporadic_only()
  list(sporadic_burden=sporadic_burden_summary, enrich=enrich, sp2=sp2)
}

# car_all already filtered to meta2 samples earlier
car_with <- car_all
car_no   <- car_all %>% filter(gene != "ATXN2")

out_with <- run_sporadic_only(car_with, "withATXN2")
out_no   <- run_sporadic_only(car_no,   "noATXN2")

out_with$sporadic_burden
head(out_with$enrich, 10)

out_no$sporadic_burden
head(out_no$enrich, 10)

out_with <- run_sporadic_only(car_with, "withATXN2")
sp2 <- out_with$sp2

sp2$any_hit <- as.integer(sp2$ptvsplice_burden > 0)
tab <- table(sp2$bin, sp2$any_hit)
tab
prop.table(tab, 1)
fisher.test(tab)

sp2_no <- out_no$sp2
sp2_no$any_hit <- as.integer(sp2_no$ptvsplice_burden > 0)
tab_no <- table(sp2_no$bin, sp2_no$any_hit)
tab_no
prop.table(tab_no, 1)
fisher.test(tab_no)

# using sp2_no from your noATXN2 run
low_ids <- sp2_no$CASE_ID[sp2_no$bin=="Low"]

car_no_low <- car_no %>% filter(sample %in% low_ids)

car_no_low %>%
  distinct(sample, gene) %>%
  count(gene, sort=TRUE)

car_no_low %>%
  filter(gene=="CHCHD10") %>%
  distinct(sample)


#Pathway comparison

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(msigdbr)
})

# Hallmark gene sets keyed by Ensembl
TERM2GENE <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, ensembl_gene) %>%
  distinct() %>%
  filter(!is.na(ensembl_gene) & ensembl_gene != "")

strip_version <- function(x) sub("\\..*$", "", x)

make_rank_from_deseq <- function(res) {
  df <- as.data.frame(res)
  df$ensembl <- strip_version(rownames(df))
  
  # choose rank stat
  if ("stat" %in% colnames(df) && any(!is.na(df$stat))) {
    r <- df$stat
  } else {
    r <- sign(df$log2FoldChange) * (-log10(df$pvalue))
  }
  names(r) <- df$ensembl
  
  # keep finite values
  ok <- is.finite(r) & !is.na(names(r)) & names(r) != ""
  r <- r[ok]
  
  # dedupe Ensembl IDs and force numeric vector (not array)
  r2 <- tapply(r, names(r), function(x) x[which.max(abs(x))])
  r2 <- as.numeric(r2)
  names(r2) <- names(tapply(r, names(r), identity))  # restore names safely
  # simpler: keep names from r2 before numeric coercion:
  # nm <- names(r2); r2 <- as.numeric(r2); names(r2) <- nm
  
  r2 <- sort(r2, decreasing = TRUE)
  return(r2)
}

gsea_run <- function(rank_vec, tag) {
  # sanity: must be numeric vector
  stopifnot(is.numeric(rank_vec), !is.null(names(rank_vec)))
  
  g <- GSEA(
    geneList = rank_vec,
    TERM2GENE = TERM2GENE,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE,
    eps = 0  # addresses fgsea warning about tiny p-values
  )
  tab <- as.data.frame(g) %>%
    dplyr::select(ID, Description, NES, pvalue, p.adjust, core_enrichment)
  tab$contrast <- tag
  tab
}

rank_C9   <- make_rank_from_deseq(res_C9_filtered)
rank_SOD1 <- make_rank_from_deseq(res_SOD1_filtered)

# quick overlap check
cat("Ranked genes:", length(rank_C9), "\n")
cat("Mapped to Hallmark:", mean(names(rank_C9) %in% TERM2GENE$ensembl_gene), "\n")

C9_tab   <- gsea_run(rank_C9, "C9_vs_CTRL")
SOD1_tab <- gsea_run(rank_SOD1, "SOD1_vs_CTRL")
SPORADIC_tab <- gsea_run(make_rank_from_deseq(res_Sporadic_filtered), "SPORADIC_vs_CTRL")

library(dplyr)

C9_sig   <- C9_tab   %>% mutate(sig_C9 = p.adjust < 0.05) %>% dplyr::select(ID, NES_C9 = NES, FDR_C9 = p.adjust, sig_C9)
SOD1_sig <- SOD1_tab %>% mutate(sig_SOD1 = p.adjust < 0.05) %>% dplyr::select(ID, NES_SOD1 = NES, FDR_SOD1 = p.adjust, sig_SOD1)

comp <- full_join(C9_sig, SOD1_sig, by="ID") %>%
  mutate(
    shared_sig = (sig_C9 %in% TRUE) & (sig_SOD1 %in% TRUE),
    concordant = ifelse(shared_sig, sign(NES_C9) == sign(NES_SOD1), NA)
  )

# headline stats
comp %>% summarize(
  n_sig_C9 = sum(sig_C9, na.rm=TRUE),
  n_sig_SOD1 = sum(sig_SOD1, na.rm=TRUE),
  n_shared_sig = sum(shared_sig, na.rm=TRUE),
  frac_concordant_shared = mean(concordant, na.rm=TRUE),
  spearman_NES_all = cor(NES_C9, NES_SOD1, method="spearman", use="pairwise.complete.obs")
)

# top shared pathways for slides
top_shared <- comp %>%
  filter(shared_sig) %>%
  arrange(FDR_C9, FDR_SOD1) %>%
  dplyr::select(ID, NES_C9, FDR_C9, NES_SOD1, FDR_SOD1, concordant)

top_shared

plot_df <- comp %>% filter(!is.na(NES_C9) & !is.na(NES_SOD1))

plot(plot_df$NES_C9, plot_df$NES_SOD1,
     xlab="NES (C9 vs Control)", ylab="NES (SOD1 vs Control)",
     main="Hallmark pathway concordance: C9 vs SOD1")

abline(h=0, v=0, lty=2)
abline(lm(NES_SOD1 ~ NES_C9, data=plot_df))



suppressPackageStartupMessages({
  library(dplyr)
  library(DESeq2)
  library(clusterProfiler)
  library(msigdbr)
})

# ---- user inputs ----
# meta2 must have: Sample_ID, group, guardrail_score_strict, sex (and optionally SVs)
# dds_filtered: DESeqDataSet with colData(dds_filtered)$Sample_ID matching meta2$Sample_ID

# 1) Build High/Low bins within sporadics
sp <- meta2 %>%
  filter(group == "Sporadic", !is.na(guardrail_score_strict), !is.na(Sample_ID))

q <- quantile(sp$guardrail_score_strict, probs = c(0.25, 0.75), na.rm = TRUE)

sp <- sp %>%
  mutate(
    guard_bin = case_when(
      guardrail_score_strict <= q[1] ~ "Low",
      guardrail_score_strict >= q[2] ~ "High",
      TRUE ~ NA_character_
    )
  )

table(sp$guard_bin, useNA="ifany")

keep_ids <- sp %>% filter(guard_bin %in% c("High","Low")) %>% pull(Sample_ID)

# 2) Subset DESeq2 object to High/Low sporadics
dds_sp <- dds_filtered[, colData(dds_filtered)$Sample_ID %in% keep_ids]
dds_sp <- dds_sp[, !is.na(colData(dds_sp)$Sample_ID)]  # extra safety

# attach guard_bin into colData in correct order
bin_map <- sp %>%
  filter(Sample_ID %in% keep_ids) %>%
  dplyr::select(Sample_ID, guard_bin) %>%
  distinct()

colData(dds_sp)$guard_bin <- bin_map$guard_bin[match(colData(dds_sp)$Sample_ID, bin_map$Sample_ID)]
colData(dds_sp)$guard_bin <- relevel(factor(colData(dds_sp)$guard_bin), ref = "Low")

stopifnot(all(!is.na(colData(dds_sp)$guard_bin)))  # must be only High/Low

# 3) DESeq2 model
# Start simple; add SVs later once this works.
# Ensure 'sex' exists in colData(dds_sp) or remove it.
design(dds_sp) <- ~ SEX + guard_bin
dds_sp <- DESeq(dds_sp)

res_sp <- results(dds_sp, contrast = c("guard_bin", "High", "Low"))

# 4) Hallmark GSEA using ENSEMBL IDs
TERM2GENE <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, ensembl_gene) %>%
  distinct() %>%
  filter(!is.na(ensembl_gene) & ensembl_gene != "")

strip_version <- function(x) sub("\\..*$", "", x)

dedupe_maxabs <- function(r) {
  df <- data.frame(id = names(r), val = as.numeric(r), stringsAsFactors = FALSE)
  df <- df %>% group_by(id) %>% slice_max(order_by = abs(val), n = 1, with_ties = FALSE) %>% ungroup()
  out <- df$val
  names(out) <- df$id
  sort(out, decreasing = TRUE)
}

make_rank_from_deseq <- function(res) {
  df <- as.data.frame(res)
  df$ensembl <- strip_version(rownames(df))
  
  if ("stat" %in% colnames(df) && any(!is.na(df$stat))) {
    r <- df$stat
  } else {
    r <- sign(df$log2FoldChange) * (-log10(df$pvalue))
  }
  names(r) <- df$ensembl
  r <- r[is.finite(r) & !is.na(names(r)) & names(r) != ""]
  dedupe_maxabs(r)
}

gsea_run <- function(rank_vec, tag) {
  g <- GSEA(
    geneList = rank_vec,
    TERM2GENE = TERM2GENE,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE,
    eps = 0
  )
  tab <- as.data.frame(g) %>%
    dplyr::select(ID, Description, NES, pvalue, p.adjust, core_enrichment)
  tab$contrast <- tag
  tab
}

rank_sp <- make_rank_from_deseq(res_sp)

cat("Sporadic High/Low ranked genes:", length(rank_sp), "\n")
cat("Mapped to Hallmark:", mean(names(rank_sp) %in% TERM2GENE$ensembl_gene), "\n")

SP_tab <- gsea_run(rank_sp, "Sporadic_High_vs_Low")

dim(SP_tab)
head(SP_tab, 10)

# If you still have C9_tab and SOD1_tab from earlier:
comp_paths <- C9_tab %>% dplyr::select(ID, NES_C9 = NES, FDR_C9 = p.adjust) %>%
  full_join(SOD1_tab %>% dplyr::select(ID, NES_SOD1 = NES, FDR_SOD1 = p.adjust), by="ID") %>%
  full_join(SP_tab %>% dplyr::select(ID, NES_SP = NES, FDR_SP = p.adjust), by="ID")

# focus on pathways significant in both C9 and SOD1 (your shared set)
comp_paths %>%
  filter(FDR_C9 < 0.05, FDR_SOD1 < 0.05) %>%
  arrange(FDR_C9, FDR_SOD1) %>%
  dplyr::select(ID, NES_C9, FDR_C9, NES_SOD1, FDR_SOD1, NES_SP, FDR_SP) %>%
  head(20)

SP_tab %>% arrange(p.adjust) %>% head(10)


# data resources

table(annot_als_filtered$group)

pdf("mutation_subtype_distribution.pdf", width=6, height=4)
barplot(table(annot_als_filtered$group),
        col=c("darkgreen","red","blue"),
        main="Mutation Subtype Distribution",
        ylab="Number of Patients",
        xlab="Mutation Subtype")
dev.off()

pdf("sex_distribution_answerALs.pdf", width=6, height=4)
barplot(table(annot_als_filtered$SEX),
        col=c("steelblue","orange"),
        main="Sex Distribution",
        ylab="Number of Patients",
        xlab="Sex")
dev.off()

pdf("age_at_symptom_onset_histogram.pdf", width=6, height=4)
hist(annot_als_filtered$AGE_AT_SYMPTOM_ONSET,
     breaks=20,
     col="gray",
     main="Age at Symptom Onset",
     xlab="Age (years)")
dev.off()

pdf("alsfrs_progression_slope_histogram.pdf", width=6, height=4)
hist(annot_als_filtered$ALSFRS_R_PROGRESSION_SLOPE,
     col="lightgray",
     main="ALSFRS Progression Rate Distribution",
     xlab="ALSFRS-R slope")
dev.off()

pdf("sexvsmutationsubtype.pdf", width=6, height=4)
par(mar=c(5,4,4,8), xpd=TRUE)

barplot(table(annot_als_filtered$SEX, annot_als_filtered$Mutation_Subtype),
        beside=TRUE,
        col=c("steelblue","orange"),
        main="Sex Distribution Across Mutation Subtypes",
        ylab="Number of Participants",
        xlab="Mutation Subtype")

legend("topright",
       inset=c(0.15,0.15),
       legend=c("Female","Male"),
       fill=c("steelblue","orange"),
       bty="n")
dev.off()





#Standardized comparison

############################################
## ALS subtype convergence / divergence
## RNA + proteomics comparison framework
############################################

## Packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(fgsea)
  library(msigdbr)
})

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

############################################
## 0. USER SETTINGS
############################################

## Significance cutoffs
rna_padj_cutoff  <- 0.05
prot_padj_cutoff <- 0.1

## Optional effect-size cutoffs
rna_lfc_cutoff   <- 0
prot_lfc_cutoff  <- 0

p_col = "pval"
prot_pval_cutoff <- 0.05
p_cutoff = prot_pval_cutoff

## Column names in your DE result tables
## Change these if needed
feature_col_candidates <- c("gene", "gene_name", "symbol", "Gene", "SYMBOL", "external_gene_name")
lfc_col_candidates     <- c("log2FoldChange", "logFC", "log2FC", "LFC")
padj_col_candidates    <- c("padj", "adj.P.Val", "FDR", "qvalue", "p_adj")

############################################
## 1. INPUT OBJECTS
############################################

## RNA signatures
rna_list <- list(
  C9       = res_C9_filtered,
  SOD1     = res_SOD1_filtered,
  Sporadic = res_Sporadic_filtered
)

## Proteomics signatures
prot_list <- list(
  C9       = res_c9_ship,
  SOD1     = res_sod1_ship,   
  Sporadic = res_spor_ship
)

############################################
## 2. HELPER FUNCTIONS
############################################

## Find a column name safely
pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

## Standardize DE result object
standardize_rna_de <- function(df) {
  df <- as.data.frame(df)
  
  ## keep rownames if needed
  if (!"feature" %in% colnames(df)) {
    df <- df %>% rownames_to_column("feature")
  }
  
  ## remove version numbers if present
  df$feature_clean <- sub("\\..*$", "", df$feature)
  
  ## map Ensembl -> symbol
  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(df$feature_clean),
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  ) %>%
    distinct(ENSEMBL, .keep_all = TRUE)
  
  out <- df %>%
    left_join(gene_map, by = c("feature_clean" = "ENSEMBL")) %>%
    transmute(
      feature = feature_clean,
      symbol = SYMBOL,
      log2FC = as.numeric(log2FoldChange),
      padj = as.numeric(padj)
    ) %>%
    filter(!is.na(log2FC), !is.na(padj), !is.na(symbol), symbol != "") %>%
    mutate(
      feature_upper = toupper(symbol)
    ) %>%
    distinct(feature_upper, .keep_all = TRUE)
  
  out
}

standardize_prot_de <- function(df) {
  df <- as.data.frame(df)
  
  df %>%
    transmute(
      feature = gene,
      log2FC = as.numeric(logFC),
      padj = as.numeric(adj.P.Val),
      pval = as.numeric(P.Value),
      tstat = as.numeric(t)
    ) %>%
    filter(!is.na(feature), feature != "", !is.na(log2FC)) %>%
    mutate(
      feature_upper = toupper(feature)
    ) %>%
    distinct(feature_upper, .keep_all = TRUE)
}

standardize_prot <- function(df){
  
  df %>%
    rename(
      log2FC = logFC,
      pval = P.Value,
      padj = adj.P.Val,
      feature = gene
    ) %>%
    mutate(
      feature_upper = toupper(feature)
    )
}
## Call significant features
call_sig <- function(df, padj_cutoff = 0.05, lfc_cutoff = 0) {
  df %>%
    mutate(
      significant = !is.na(padj) & padj < padj_cutoff & abs(log2FC) > lfc_cutoff,
      direction   = case_when(
        log2FC > 0 ~ "up",
        log2FC < 0 ~ "down",
        TRUE ~ "zero"
      )
    )
}

call_sig2 <- function(df, p_col = c("padj","pval"), p_cutoff = 0.05, lfc_cutoff = 0) {
  
  p_col <- match.arg(p_col)
  
  df %>%
    mutate(
      p_use = if (p_col == "padj") padj else pval,
      significant = !is.na(p_use) & p_use < p_cutoff & abs(log2FC) > lfc_cutoff,
      direction = case_when(
        log2FC > 0 ~ "up",
        log2FC < 0 ~ "down",
        TRUE ~ "zero"
      )
    )
}

## Basic summary metrics for one signature
summarize_signature <- function(df, padj_cutoff = 0.05, lfc_cutoff = 0, label = NULL) {
  tmp <- call_sig(df, padj_cutoff, lfc_cutoff)
  
  tibble(
    signature = label,
    n_features_tested = nrow(tmp),
    n_sig = sum(tmp$significant, na.rm = TRUE),
    median_abs_log2FC_all = median(abs(tmp$log2FC), na.rm = TRUE),
    median_abs_log2FC_sig = ifelse(sum(tmp$significant, na.rm = TRUE) > 0,
                                   median(abs(tmp$log2FC[tmp$significant]), na.rm = TRUE),
                                   NA_real_)
  )
}

## Pairwise similarity between two subtype signatures within one omics layer
compare_two_signatures <- function(df1, df2, name1, name2,
                                   padj_cutoff = 0.05,
                                   lfc_cutoff = 0) {
  
  x1 <- call_sig(df1, padj_cutoff, lfc_cutoff)
  x2 <- call_sig(df2, padj_cutoff, lfc_cutoff)
  
  merged <- inner_join(
    x1 %>% select(feature_upper, feature1 = feature, log2FC_1 = log2FC, padj_1 = padj, sig_1 = significant),
    x2 %>% select(feature_upper, feature2 = feature, log2FC_2 = log2FC, padj_2 = padj, sig_2 = significant),
    by = "feature_upper"
  ) %>%
    mutate(
      sign_1 = sign(log2FC_1),
      sign_2 = sign(log2FC_2),
      sign_agree = sign_1 == sign_2
    )
  
  sig1 <- x1 %>% filter(significant) %>% pull(feature_upper)
  sig2 <- x2 %>% filter(significant) %>% pull(feature_upper)
  
  universe <- union(x1$feature_upper, x2$feature_upper)
  overlap_sig <- intersect(sig1, sig2)
  
  ## Fisher / hypergeometric-style overlap test
  a <- length(overlap_sig)
  b <- length(setdiff(sig1, sig2))
  c <- length(setdiff(sig2, sig1))
  d <- length(setdiff(universe, union(sig1, sig2)))
  
  fisher_mat <- matrix(c(a, b, c, d), nrow = 2,
                       dimnames = list(c(name1, "not_" %+% name1),
                                       c(name2, "not_" %+% name2)))
  
  fisher_res <- fisher.test(fisher_mat, alternative = "greater")
  
  tibble(
    comparison = paste(name1, "vs", name2),
    n_shared_features = nrow(merged),
    spearman_log2FC = suppressWarnings(cor(merged$log2FC_1, merged$log2FC_2,
                                           method = "spearman", use = "pairwise.complete.obs")),
    pearson_log2FC  = suppressWarnings(cor(merged$log2FC_1, merged$log2FC_2,
                                           method = "pearson", use = "pairwise.complete.obs")),
    sign_agreement_all = mean(merged$sign_agree, na.rm = TRUE),
    n_sig_1 = length(sig1),
    n_sig_2 = length(sig2),
    n_sig_overlap = a,
    jaccard_sig = ifelse(length(union(sig1, sig2)) > 0,
                         a / length(union(sig1, sig2)),
                         NA_real_),
    fisher_overlap_p = fisher_res$p.value,
    fisher_overlap_or = unname(fisher_res$estimate)
  )
}

`%+%` <- function(a, b) paste0(a, b)

############################################
## 3. PATHWAY SETS (Hallmark default)
############################################

hallmark_df <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

hallmark_list <- split(toupper(hallmark_df$gene_symbol), hallmark_df$gs_name)

############################################
## 4. FGSEA WRAPPER
############################################

run_fgsea_from_de <- function(df, pathways, stat_col = "log2FC", minSize = 10, maxSize = 500) {
  stopifnot(all(c("feature_upper", stat_col) %in% colnames(df)))
  
  stats_df <- df %>%
    select(feature_upper, stat = all_of(stat_col)) %>%
    filter(!is.na(feature_upper), !is.na(stat), is.finite(stat)) %>%
    group_by(feature_upper) %>%
    summarise(stat = mean(stat), .groups = "drop") %>%
    arrange(desc(stat))
  
  stats <- stats_df$stat
  names(stats) <- stats_df$feature_upper
  
  if (length(stats) == 0) {
    return(tibble(
      pathway = character(),
      pval = numeric(),
      padj = numeric(),
      ES = numeric(),
      NES = numeric(),
      size = integer(),
      leadingEdge = list(),
      direction = character(),
      signed_score = numeric()
    ))
  }
  
  fg <- tryCatch(
    fgsea::fgsea(
      pathways = pathways,
      stats = stats,
      minSize = minSize,
      maxSize = maxSize
    ),
    error = function(e) NULL
  )
  
  if (is.null(fg) || nrow(fg) == 0) {
    return(tibble(
      pathway = character(),
      pval = numeric(),
      padj = numeric(),
      ES = numeric(),
      NES = numeric(),
      size = integer(),
      leadingEdge = list(),
      direction = character(),
      signed_score = numeric()
    ))
  }
  
  fg %>%
    as_tibble() %>%
    mutate(
      pathway = as.character(pathway),
      direction = case_when(
        NES > 0 ~ "up",
        NES < 0 ~ "down",
        TRUE ~ "zero"
      ),
      signed_score = sign(NES) * -log10(pmax(pval, 1e-300))
    ) %>%
    arrange(padj)
}

## Summarize pathway enrichment
summarize_fgsea <- function(fg, label = NULL, padj_cutoff = 0.05) {
  tibble(
    signature = label,
    n_pathways_tested = nrow(fg),
    n_pathways_sig = sum(fg$padj < padj_cutoff, na.rm = TRUE),
    median_abs_NES_all = median(abs(fg$NES), na.rm = TRUE),
    median_abs_NES_sig = ifelse(sum(fg$padj < padj_cutoff, na.rm = TRUE) > 0,
                                median(abs(fg$NES[fg$padj < padj_cutoff]), na.rm = TRUE),
                                NA_real_)
  )
}

## Compare pathways between two subtype signatures within one omics layer
compare_two_fgsea <- function(fg1, fg2, name1, name2, padj_cutoff = 0.05) {
  x1 <- fg1 %>%
    select(pathway, NES_1 = NES, padj_1 = padj, signed_score_1 = signed_score) %>%
    mutate(sig_1 = padj_1 < padj_cutoff)
  
  x2 <- fg2 %>%
    select(pathway, NES_2 = NES, padj_2 = padj, signed_score_2 = signed_score) %>%
    mutate(sig_2 = padj_2 < padj_cutoff)
  
  merged <- inner_join(x1, x2, by = "pathway") %>%
    mutate(
      sign_1 = sign(NES_1),
      sign_2 = sign(NES_2),
      sign_agree = sign_1 == sign_2
    )
  
  sig1 <- x1 %>% filter(sig_1) %>% pull(pathway)
  sig2 <- x2 %>% filter(sig_2) %>% pull(pathway)
  
  universe <- union(x1$pathway, x2$pathway)
  overlap_sig <- intersect(sig1, sig2)
  
  a <- length(overlap_sig)
  b <- length(setdiff(sig1, sig2))
  c <- length(setdiff(sig2, sig1))
  d <- length(setdiff(universe, union(sig1, sig2)))
  
  fisher_mat <- matrix(c(a, b, c, d), nrow = 2)
  fisher_res <- fisher.test(fisher_mat, alternative = "greater")
  
  tibble(
    comparison = paste(name1, "vs", name2),
    n_shared_pathways = nrow(merged),
    spearman_NES = suppressWarnings(cor(merged$NES_1, merged$NES_2,
                                        method = "spearman", use = "pairwise.complete.obs")),
    pearson_NES = suppressWarnings(cor(merged$NES_1, merged$NES_2,
                                       method = "pearson", use = "pairwise.complete.obs")),
    spearman_signed_score = suppressWarnings(cor(merged$signed_score_1, merged$signed_score_2,
                                                 method = "spearman", use = "pairwise.complete.obs")),
    sign_agreement_all = mean(merged$sign_agree, na.rm = TRUE),
    n_sig_path_1 = length(sig1),
    n_sig_path_2 = length(sig2),
    n_sig_path_overlap = a,
    jaccard_sig_path = ifelse(length(union(sig1, sig2)) > 0,
                              a / length(union(sig1, sig2)),
                              NA_real_),
    fisher_path_overlap_p = fisher_res$p.value,
    fisher_path_overlap_or = unname(fisher_res$estimate)
  )
}

############################################
## 5. CROSS-OMICS COMPARISON FOR SAME SUBTYPE
############################################

compare_rna_protein_same_signature <- function(rna_df, prot_df, label,
                                               rna_padj_cutoff = 0.05,
                                               prot_padj_cutoff = 0.05,
                                               rna_lfc_cutoff = 0,
                                               prot_lfc_cutoff = 0) {
  
  rna_sig  <- call_sig(rna_df,  rna_padj_cutoff,  rna_lfc_cutoff)
  prot_sig <- call_sig(prot_df, prot_padj_cutoff, prot_lfc_cutoff)
  
  merged <- inner_join(
    rna_sig  %>% select(feature_upper, log2FC_rna = log2FC, padj_rna = padj, sig_rna = significant),
    prot_sig %>% select(feature_upper, log2FC_prot = log2FC, padj_prot = padj, sig_prot = significant),
    by = "feature_upper"
  ) %>%
    mutate(
      sign_rna  = sign(log2FC_rna),
      sign_prot = sign(log2FC_prot),
      sign_agree = sign_rna == sign_prot
    )
  
  sig_rna  <- rna_sig  %>% filter(significant) %>% pull(feature_upper)
  sig_prot <- prot_sig %>% filter(significant) %>% pull(feature_upper)
  
  universe <- union(rna_sig$feature_upper, prot_sig$feature_upper)
  overlap_sig <- intersect(sig_rna, sig_prot)
  
  a <- length(overlap_sig)
  b <- length(setdiff(sig_rna, sig_prot))
  c <- length(setdiff(sig_prot, sig_rna))
  d <- length(setdiff(universe, union(sig_rna, sig_prot)))
  
  fisher_mat <- matrix(c(a, b, c, d), nrow = 2)
  fisher_res <- fisher.test(fisher_mat, alternative = "greater")
  
  tibble(
    signature = label,
    n_shared_features = nrow(merged),
    spearman_log2FC = suppressWarnings(cor(merged$log2FC_rna, merged$log2FC_prot,
                                           method = "spearman", use = "pairwise.complete.obs")),
    pearson_log2FC = suppressWarnings(cor(merged$log2FC_rna, merged$log2FC_prot,
                                          method = "pearson", use = "pairwise.complete.obs")),
    sign_agreement_all = mean(merged$sign_agree, na.rm = TRUE),
    n_sig_rna = length(sig_rna),
    n_sig_prot = length(sig_prot),
    n_sig_overlap = a,
    jaccard_sig = ifelse(length(union(sig_rna, sig_prot)) > 0,
                         a / length(union(sig_rna, sig_prot)),
                         NA_real_),
    fisher_overlap_p = fisher_res$p.value,
    fisher_overlap_or = unname(fisher_res$estimate)
  )
}

compare_rna_protein_fgsea <- function(rna_fg, prot_fg, label, padj_cutoff = 0.05) {
  x1 <- rna_fg %>%
    select(pathway, NES_rna = NES, padj_rna = padj, signed_score_rna = signed_score) %>%
    mutate(sig_rna = padj_rna < padj_cutoff)
  
  x2 <- prot_fg %>%
    select(pathway, NES_prot = NES, padj_prot = padj, signed_score_prot = signed_score) %>%
    mutate(sig_prot = padj_prot < padj_cutoff)
  
  merged <- inner_join(x1, x2, by = "pathway") %>%
    mutate(
      sign_rna = sign(NES_rna),
      sign_prot = sign(NES_prot),
      sign_agree = sign_rna == sign_prot
    )
  
  sig1 <- x1 %>% filter(sig_rna) %>% pull(pathway)
  sig2 <- x2 %>% filter(sig_prot) %>% pull(pathway)
  
  universe <- union(x1$pathway, x2$pathway)
  overlap_sig <- intersect(sig1, sig2)
  
  a <- length(overlap_sig)
  b <- length(setdiff(sig1, sig2))
  c <- length(setdiff(sig2, sig1))
  d <- length(setdiff(universe, union(sig1, sig2)))
  
  fisher_mat <- matrix(c(a, b, c, d), nrow = 2)
  fisher_res <- fisher.test(fisher_mat, alternative = "greater")
  
  tibble(
    signature = label,
    n_shared_pathways = nrow(merged),
    spearman_NES = suppressWarnings(cor(merged$NES_rna, merged$NES_prot,
                                        method = "spearman", use = "pairwise.complete.obs")),
    spearman_signed_score = suppressWarnings(cor(merged$signed_score_rna, merged$signed_score_prot,
                                                 method = "spearman", use = "pairwise.complete.obs")),
    sign_agreement_all = mean(merged$sign_agree, na.rm = TRUE),
    n_sig_path_rna = length(sig1),
    n_sig_path_prot = length(sig2),
    n_sig_path_overlap = a,
    jaccard_sig_path = ifelse(length(union(sig1, sig2)) > 0,
                              a / length(union(sig1, sig2)),
                              NA_real_),
    fisher_path_overlap_p = fisher_res$p.value,
    fisher_path_overlap_or = unname(fisher_res$estimate)
  )
}

compare_two_fgsea <- function(fg1, fg2, name1, name2, padj_cutoff = 0.05) {
  x1 <- fg1 %>%
    select(pathway, NES_1 = NES, padj_1 = padj, signed_score_1 = signed_score) %>%
    mutate(sig_1 = padj_1 < padj_cutoff)
  
  x2 <- fg2 %>%
    select(pathway, NES_2 = NES, padj_2 = padj, signed_score_2 = signed_score) %>%
    mutate(sig_2 = padj_2 < padj_cutoff)
  
  merged <- inner_join(x1, x2, by = "pathway") %>%
    mutate(
      sign_1 = sign(NES_1),
      sign_2 = sign(NES_2),
      sign_agree = sign_1 == sign_2
    )
  
  sig1 <- x1 %>% filter(sig_1) %>% pull(pathway)
  sig2 <- x2 %>% filter(sig_2) %>% pull(pathway)
  
  universe <- union(x1$pathway, x2$pathway)
  overlap_sig <- intersect(sig1, sig2)
  
  a <- length(overlap_sig)
  b <- length(setdiff(sig1, sig2))
  c <- length(setdiff(sig2, sig1))
  d <- length(setdiff(universe, union(sig1, sig2)))
  
  fisher_p <- NA_real_
  fisher_or <- NA_real_
  
  if ((a + b + c + d) > 0) {
    fisher_out <- tryCatch(
      fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater"),
      error = function(e) NULL
    )
    if (!is.null(fisher_out)) {
      fisher_p <- fisher_out$p.value
      fisher_or <- unname(fisher_out$estimate)
    }
  }
  
  if (nrow(merged) == 0) {
    return(tibble(
      comparison = paste(name1, "vs", name2),
      n_shared_pathways = 0,
      spearman_NES = NA_real_,
      pearson_NES = NA_real_,
      spearman_signed_score = NA_real_,
      sign_agreement_all = NA_real_,
      n_sig_path_1 = length(sig1),
      n_sig_path_2 = length(sig2),
      n_sig_path_overlap = a,
      jaccard_sig_path = ifelse(length(union(sig1, sig2)) > 0,
                                a / length(union(sig1, sig2)),
                                NA_real_),
      fisher_path_overlap_p = fisher_p,
      fisher_path_overlap_or = fisher_or
    ))
  }
  
  tibble(
    comparison = paste(name1, "vs", name2),
    n_shared_pathways = nrow(merged),
    spearman_NES = suppressWarnings(cor(merged$NES_1, merged$NES_2, method = "spearman")),
    pearson_NES = suppressWarnings(cor(merged$NES_1, merged$NES_2, method = "pearson")),
    spearman_signed_score = suppressWarnings(cor(merged$signed_score_1, merged$signed_score_2, method = "spearman")),
    sign_agreement_all = mean(merged$sign_agree, na.rm = TRUE),
    n_sig_path_1 = length(sig1),
    n_sig_path_2 = length(sig2),
    n_sig_path_overlap = a,
    jaccard_sig_path = ifelse(length(union(sig1, sig2)) > 0,
                              a / length(union(sig1, sig2)),
                              NA_real_),
    fisher_path_overlap_p = fisher_p,
    fisher_path_overlap_or = fisher_or
  )
}
############################################
## 6. STANDARDIZE INPUTS
############################################

rna_std <- list(
  C9       = standardize_rna_de(res_C9_filtered),
  SOD1     = standardize_rna_de(res_SOD1_filtered),
  Sporadic = standardize_rna_de(res_Sporadic_filtered)
)
prot_std <- list(
  C9       = standardize_prot_de(res_c9_ship),
  SOD1     = standardize_prot_de(res_sod1_ship),
  Sporadic = standardize_prot_de(res_spor_ship)
)

prot_std <- list(
  C9 = standardize_prot(res_c9_ship) %>%
    call_sig2(p_col = "pval", p_cutoff = 0.05, lfc_cutoff = prot_lfc_cutoff),
  
  SOD1 = standardize_prot(res_sod1_ship) %>%
    call_sig2(p_col = "pval", p_cutoff = 0.05, lfc_cutoff = prot_lfc_cutoff),
  
  Sporadic = standardize_prot(res_spor_ship) %>%
    call_sig2(p_col = "pval", p_cutoff = 0.05, lfc_cutoff = prot_lfc_cutoff)
)

############################################
## 7. BASIC SIGNATURE SUMMARIES
############################################

rna_signature_summary <- imap_dfr(
  rna_std,
  ~ summarize_signature(.x, padj_cutoff = rna_padj_cutoff, lfc_cutoff = rna_lfc_cutoff, label = .y)
)

prot_signature_summary <- imap_dfr(
  prot_std,
  ~ summarize_signature(.x, padj_cutoff = prot_padj_cutoff, lfc_cutoff = prot_lfc_cutoff, label = .y)
)

cat("\nRNA signature summary:\n")
print(rna_signature_summary)

cat("\nProtein signature summary:\n")
print(prot_signature_summary)

############################################
## 8. PAIRWISE SUBTYPE COMPARISONS WITHIN RNA
############################################

pair_names <- list(
  c("C9", "Sporadic"),
  c("C9", "SOD1"),
  c("SOD1", "Sporadic")
)

rna_pairwise_feature <- map_dfr(pair_names, function(p) {
  compare_two_signatures(
    df1 = rna_std[[p[1]]],
    df2 = rna_std[[p[2]]],
    name1 = p[1],
    name2 = p[2],
    padj_cutoff = rna_padj_cutoff,
    lfc_cutoff = rna_lfc_cutoff
  )
})

cat("\nRNA pairwise feature-level comparisons:\n")
print(rna_pairwise_feature)

############################################
## 9. PAIRWISE SUBTYPE COMPARISONS WITHIN PROTEIN
############################################

prot_pairwise_feature <- map_dfr(pair_names, function(p) {
  compare_two_signatures(
    df1 = prot_std[[p[1]]],
    df2 = prot_std[[p[2]]],
    name1 = p[1],
    name2 = p[2],
    padj_cutoff = prot_padj_cutoff,
    lfc_cutoff = prot_lfc_cutoff
  )
})

cat("\nProtein pairwise feature-level comparisons:\n")
print(prot_pairwise_feature)

############################################
## 10. RUN PATHWAY ENRICHMENT
############################################

rna_fgsea <- purrr::imap(rna_std, ~ run_fgsea_from_de(.x, hallmark_list, stat_col = "log2FC"))
prot_fgsea <- purrr::imap(prot_std, ~ run_fgsea_from_de(.x, hallmark_list, stat_col = "tstat"))

############################################
## 11. PATHWAY SUMMARIES
############################################

rna_pathway_summary <- purrr::imap_dfr(
  rna_fgsea,
  ~ summarize_fgsea(.x, label = .y, padj_cutoff = 0.05)
)

prot_pathway_summary <- purrr::imap_dfr(
  prot_fgsea,
  ~ summarize_fgsea(.x, label = .y, padj_cutoff = 0.05)
)

print(rna_pathway_summary)
print(prot_pathway_summary)

cat("\nRNA pathway summary:\n")
print(rna_pathway_summary)

cat("\nProtein pathway summary:\n")
print(prot_pathway_summary)

############################################
## 12. PAIRWISE SUBTYPE COMPARISONS AT PATHWAY LEVEL
############################################

pair_names <- list(
  c("C9", "Sporadic"),
  c("C9", "SOD1"),
  c("SOD1", "Sporadic")
)

rna_pairwise_pathway <- purrr::map_dfr(pair_names, function(p) {
  compare_two_fgsea(
    fg1 = rna_fgsea[[p[1]]],
    fg2 = rna_fgsea[[p[2]]],
    name1 = p[1],
    name2 = p[2],
    padj_cutoff = 0.05
  )
})

prot_pairwise_pathway <- purrr::map_dfr(pair_names, function(p) {
  compare_two_fgsea(
    fg1 = prot_fgsea[[p[1]]],
    fg2 = prot_fgsea[[p[2]]],
    name1 = p[1],
    name2 = p[2],
    padj_cutoff = 0.05
  )
})

print(rna_pairwise_pathway)
print(prot_pairwise_pathway)

cat("\nRNA pairwise pathway-level comparisons:\n")
print(rna_pairwise_pathway)

cat("\nProtein pairwise pathway-level comparisons:\n")
print(prot_pairwise_pathway)

############################################
## 13. CROSS-OMICS COMPARISON FOR SAME SUBTYPE
############################################

cross_omics_feature <- map_dfr(names(rna_std), function(nm) {
  compare_rna_protein_same_signature(
    rna_df = rna_std[[nm]],
    prot_df = prot_std[[nm]],
    label = nm,
    rna_padj_cutoff = rna_padj_cutoff,
    prot_padj_cutoff = prot_padj_cutoff,
    rna_lfc_cutoff = rna_lfc_cutoff,
    prot_lfc_cutoff = prot_lfc_cutoff
  )
})

cross_omics_pathway <- map_dfr(names(rna_fgsea), function(nm) {
  compare_rna_protein_fgsea(
    rna_fg = rna_fgsea[[nm]],
    prot_fg = prot_fgsea[[nm]],
    label = nm,
    padj_cutoff = 0.05
  )
})

cat("\nCross-omics feature-level comparisons (same subtype):\n")
print(cross_omics_feature)

cat("\nCross-omics pathway-level comparisons (same subtype):\n")
print(cross_omics_pathway)

############################################
## 14. BUILD A CLEAN RESULTS MATRIX
############################################

results_matrix <- rna_signature_summary %>%
  rename(
    RNA_n_sig = n_sig,
    RNA_median_abs_log2FC_all = median_abs_log2FC_all,
    RNA_median_abs_log2FC_sig = median_abs_log2FC_sig
  ) %>%
  select(signature, RNA_n_sig, RNA_median_abs_log2FC_all, RNA_median_abs_log2FC_sig) %>%
  left_join(
    prot_signature_summary %>%
      rename(
        Protein_n_sig = n_sig,
        Protein_median_abs_log2FC_all = median_abs_log2FC_all,
        Protein_median_abs_log2FC_sig = median_abs_log2FC_sig
      ) %>%
      select(signature, Protein_n_sig, Protein_median_abs_log2FC_all, Protein_median_abs_log2FC_sig),
    by = "signature"
  ) %>%
  left_join(
    rna_pathway_summary %>%
      rename(
        RNA_n_path_sig = n_pathways_sig,
        RNA_median_abs_NES_all = median_abs_NES_all,
        RNA_median_abs_NES_sig = median_abs_NES_sig
      ) %>%
      select(signature, RNA_n_path_sig, RNA_median_abs_NES_all, RNA_median_abs_NES_sig),
    by = "signature"
  ) %>%
  left_join(
    prot_pathway_summary %>%
      rename(
        Protein_n_path_sig = n_pathways_sig,
        Protein_median_abs_NES_all = median_abs_NES_all,
        Protein_median_abs_NES_sig = median_abs_NES_sig
      ) %>%
      select(signature, Protein_n_path_sig, Protein_median_abs_NES_all, Protein_median_abs_NES_sig),
    by = "signature"
  ) %>%
  left_join(
    cross_omics_feature %>%
      select(signature, cross_omics_feature_spearman = spearman_log2FC,
             cross_omics_feature_sign_agreement = sign_agreement_all,
             cross_omics_feature_overlap = n_sig_overlap),
    by = "signature"
  ) %>%
  left_join(
    cross_omics_pathway %>%
      select(signature, cross_omics_pathway_spearman = spearman_NES,
             cross_omics_pathway_sign_agreement = sign_agreement_all,
             cross_omics_pathway_overlap = n_sig_path_overlap),
    by = "signature"
  )

cat("\nPer-signature integrated results matrix:\n")
print(results_matrix)

############################################
## 15. OPTIONAL: SAVE OUTPUTS
############################################

write.csv(rna_signature_summary, "rna_signature_summary.csv", row.names = FALSE)
write.csv(prot_signature_summary, "protein_signature_summary.csv", row.names = FALSE)
write.csv(rna_pairwise_feature, "rna_pairwise_feature.csv", row.names = FALSE)
write.csv(prot_pairwise_feature, "protein_pairwise_feature.csv", row.names = FALSE)
write.csv(rna_pathway_summary, "rna_pathway_summary.csv", row.names = FALSE)
write.csv(prot_pathway_summary, "protein_pathway_summary.csv", row.names = FALSE)
write.csv(rna_pairwise_pathway, "rna_pairwise_pathway.csv", row.names = FALSE)
write.csv(prot_pairwise_pathway, "protein_pairwise_pathway.csv", row.names = FALSE)
write.csv(cross_omics_feature, "cross_omics_feature.csv", row.names = FALSE)
write.csv(cross_omics_pathway, "cross_omics_pathway.csv", row.names = FALSE)
write.csv(results_matrix, "integrated_signature_results_matrix.csv", row.names = FALSE)


compare_two_signatures_pval <- function(df1, df2, name1, name2,
                                        p_cutoff = 0.05,
                                        lfc_cutoff = 0) {
  
  x1 <- call_sig2(df1, p_col = "pval", p_cutoff = p_cutoff, lfc_cutoff = lfc_cutoff)
  x2 <- call_sig2(df2, p_col = "pval", p_cutoff = p_cutoff, lfc_cutoff = lfc_cutoff)
  
  merged <- inner_join(
    x1 %>% select(feature_upper, feature1 = feature, log2FC_1 = log2FC, p_1 = p_use, sig_1 = significant),
    x2 %>% select(feature_upper, feature2 = feature, log2FC_2 = log2FC, p_2 = p_use, sig_2 = significant),
    by = "feature_upper"
  ) %>%
    mutate(
      sign_1 = sign(log2FC_1),
      sign_2 = sign(log2FC_2),
      sign_agree = sign_1 == sign_2
    )
  
  sig1 <- x1 %>% filter(significant) %>% pull(feature_upper)
  sig2 <- x2 %>% filter(significant) %>% pull(feature_upper)
  
  universe <- union(x1$feature_upper, x2$feature_upper)
  overlap_sig <- intersect(sig1, sig2)
  
  a <- length(overlap_sig)
  b <- length(setdiff(sig1, sig2))
  c <- length(setdiff(sig2, sig1))
  d <- length(setdiff(universe, union(sig1, sig2)))
  
  fisher_mat <- matrix(c(a, b, c, d), nrow = 2)
  fisher_res <- fisher.test(fisher_mat, alternative = "greater")
  
  tibble(
    comparison = paste(name1, "vs", name2),
    n_shared_features = nrow(merged),
    spearman_log2FC = suppressWarnings(cor(merged$log2FC_1, merged$log2FC_2,
                                           method = "spearman", use = "pairwise.complete.obs")),
    pearson_log2FC  = suppressWarnings(cor(merged$log2FC_1, merged$log2FC_2,
                                           method = "pearson", use = "pairwise.complete.obs")),
    sign_agreement_all = mean(merged$sign_agree, na.rm = TRUE),
    n_sig_1 = length(sig1),
    n_sig_2 = length(sig2),
    n_sig_overlap = a,
    jaccard_sig = ifelse(length(union(sig1, sig2)) > 0,
                         a / length(union(sig1, sig2)),
                         NA_real_),
    fisher_overlap_p = fisher_res$p.value,
    fisher_overlap_or = unname(fisher_res$estimate)
  )
}

prot_pairwise_feature_pval <- map_dfr(pair_names, function(p) {
  compare_two_signatures_pval(
    df1 = prot_std[[p[1]]],
    df2 = prot_std[[p[2]]],
    name1 = p[1],
    name2 = p[2],
    p_cutoff = prot_pval_cutoff,
    lfc_cutoff = prot_lfc_cutoff
  )
})

print(prot_pairwise_feature_pval)


compare_rna_protein_same_signature_mixed <- function(rna_df, prot_df, label,
                                                     rna_padj_cutoff = 0.05,
                                                     prot_pval_cutoff = 0.05,
                                                     rna_lfc_cutoff = 0,
                                                     prot_lfc_cutoff = 0) {
  
  rna_sig  <- call_sig2(rna_df,  p_col = "padj", p_cutoff = rna_padj_cutoff,  lfc_cutoff = rna_lfc_cutoff)
  prot_sig <- call_sig2(prot_df, p_col = "pval", p_cutoff = prot_pval_cutoff, lfc_cutoff = prot_lfc_cutoff)
  
  merged <- inner_join(
    rna_sig  %>% select(feature_upper, log2FC_rna = log2FC, sig_rna = significant),
    prot_sig %>% select(feature_upper, log2FC_prot = log2FC, sig_prot = significant),
    by = "feature_upper"
  ) %>%
    mutate(
      sign_rna = sign(log2FC_rna),
      sign_prot = sign(log2FC_prot),
      sign_agree = sign_rna == sign_prot
    )
  
  sig_rna  <- rna_sig  %>% filter(significant) %>% pull(feature_upper)
  sig_prot <- prot_sig %>% filter(significant) %>% pull(feature_upper)
  
  universe <- union(rna_sig$feature_upper, prot_sig$feature_upper)
  overlap_sig <- intersect(sig_rna, sig_prot)
  
  a <- length(overlap_sig)
  b <- length(setdiff(sig_rna, sig_prot))
  c <- length(setdiff(sig_prot, sig_rna))
  d <- length(setdiff(universe, union(sig_rna, sig_prot)))
  
  fisher_mat <- matrix(c(a, b, c, d), nrow = 2)
  fisher_res <- fisher.test(fisher_mat, alternative = "greater")
  
  tibble(
    signature = label,
    n_shared_features = nrow(merged),
    spearman_log2FC = suppressWarnings(cor(merged$log2FC_rna, merged$log2FC_prot,
                                           method = "spearman", use = "pairwise.complete.obs")),
    pearson_log2FC = suppressWarnings(cor(merged$log2FC_rna, merged$log2FC_prot,
                                          method = "pearson", use = "pairwise.complete.obs")),
    sign_agreement_all = mean(merged$sign_agree, na.rm = TRUE),
    n_sig_rna = length(sig_rna),
    n_sig_prot = length(sig_prot),
    n_sig_overlap = a,
    jaccard_sig = ifelse(length(union(sig_rna, sig_prot)) > 0,
                         a / length(union(sig_rna, sig_prot)),
                         NA_real_),
    fisher_overlap_p = fisher_res$p.value,
    fisher_overlap_or = unname(fisher_res$estimate)
  )
}

cross_omics_feature_pval <- map_dfr(names(rna_std), function(nm) {
  compare_rna_protein_same_signature_mixed(
    rna_df = rna_std[[nm]],
    prot_df = prot_std[[nm]],
    label = nm,
    rna_padj_cutoff = rna_padj_cutoff,
    prot_pval_cutoff = prot_pval_cutoff,
    rna_lfc_cutoff = rna_lfc_cutoff,
    prot_lfc_cutoff = prot_lfc_cutoff
  )
})

print(cross_omics_feature_pval)
