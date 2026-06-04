# ============================================================
# Script: CompareProteomicNormalization.R
# Project: ALS multi-omics analysis
# Source status: documented/sanitized for GitHub handoff
# ============================================================
# Documentation summary from docs/SCRIPT_SPECS_FULL.md
# 
# ## `CompareProteomicNormalization.R`
# 
# **Category:** Proteomics / normalization QC
# 
# **Purpose**
# 
# Compares proteomics normalization strategies: log2(x+1), quantile normalization, row-wise z-score normalization, and median-centering. Produces normalized matrices and a PCA comparison PDF.
# 
# **Inputs**
# 
# - `proteomics/4_matrix/AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt`.
# 
# **Outputs**
# 
# - `Proteomics_Log2Plus1.csv`.
# - `Proteomics_QuantileNormalized.csv`.
# - `Proteomics_ZscoreNormalized.csv`.
# - `Proteomics_MedianCentered.csv`.
# - `PCA_Normalization_Comparison.pdf`.
# 
# **Upstream dependencies**
# 
# - None. This is a QC/preprocessing comparison script.
# 
# **How to run**
# 
# Legacy documented version:
# 
# ```bash
# Rscript scripts/documented/proteomics/CompareProteomicNormalization.R
# ```
# 
# Refactored command-line version:
# 
# ```bash
# Rscript scripts/refactored/proteomics/01_compare_normalization.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/proteomics/01_compare_normalization
# ```
# 
# **Refactor notes**
# 
# - Parameterize the input matrix and output directory.
# - Use a consistent sample/protein orientation check and write orientation to log.
# - Move outputs to `results/proteomics/normalization_qc/`.
# ============================================================

# Load required libraries
library(preprocessCore)  # For quantile normalization
library(matrixStats)     # For rowMedians, if needed

mainDir <- file.path(ALSPROJ_BASE_DIR)
setwd(mainDir)
# ----------------------------
# 1. Load raw proteomics data
# ----------------------------
# Replace this path with your actual file
raw_data <- read.delim("proteomics/4_matrix/AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt", row.names = 1)

# If samples are in rows and proteins in columns, transpose
datExpr_raw <- as.matrix(raw_data)
if (nrow(datExpr_raw) < ncol(datExpr_raw)) {
  datExpr_raw <- t(datExpr_raw)
}

# ----------------------------
# 2. Log2 + 1 transformation
# ----------------------------
datExpr_log <- log2(datExpr_raw + 1)

# ----------------------------
# 3. Quantile Normalization
# ----------------------------
datExpr_quantile <- normalize.quantiles(datExpr_log)
rownames(datExpr_quantile) <- rownames(datExpr_log)
colnames(datExpr_quantile) <- colnames(datExpr_log)

# ----------------------------
# 4. Z-score Normalization (Row-wise)
# ----------------------------
datExpr_zscore <- t(scale(t(datExpr_log), center = TRUE, scale = TRUE))  # Row-wise

# ----------------------------
# 5. (Optional) Median Normalization
# ----------------------------
# Subtract median per sample
sample_medians <- apply(datExpr_log, 2, median, na.rm = TRUE)
datExpr_median_norm <- sweep(datExpr_log, 2, sample_medians, FUN = "-")

# ----------------------------
# 6. Save Normalized Outputs
# ----------------------------
write.csv(datExpr_log, "Proteomics_Log2Plus1.csv")
write.csv(datExpr_quantile, "Proteomics_QuantileNormalized.csv")
write.csv(datExpr_zscore, "Proteomics_ZscoreNormalized.csv")
write.csv(datExpr_median_norm, "Proteomics_MedianCentered.csv")

# ----------------------------
# 7. Quick PCA Visualization (Optional)
# ----------------------------
library(ggplot2)

pca_plot <- function(mat, title) {
  pca <- prcomp(t(mat), scale. = TRUE)
  df <- as.data.frame(pca$x)
  ggplot(df, aes(PC1, PC2)) +
    geom_point() +
    ggtitle(title) +
    theme_minimal()
}

pdf("PCA_Normalization_Comparison.pdf", width = 10, height = 10)
print(pca_plot(datExpr_log, "Log2+1"))
print(pca_plot(datExpr_quantile, "Quantile Normalization"))
print(pca_plot(datExpr_zscore, "Z-score Normalization"))
print(pca_plot(datExpr_median_norm, "Median Normalization"))
dev.off()
