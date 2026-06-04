# ============================================================
# Refactored entrypoint for original script: CompareProteomicNormalization.R
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
if (isTRUE(ALS_ARGS$help)) { als_print_help("01_compare_normalization"); quit(save = "no", status = 0) }
CFG <- als_load_config(ALS_ARGS$config)
ALSPROJ_BASE_DIR <- als_choose(ALS_ARGS$base_dir, CFG, "project.base_dir", getwd())
ALSPROJ_OUTDIR <- als_choose(ALS_ARGS$outdir, CFG, "outputs.root", file.path(getwd(), "results", "01_compare_normalization"))
ALSPROJ_SCRIPT_ID <- "01_compare_normalization"
als_init(ALSPROJ_SCRIPT_ID, ALSPROJ_BASE_DIR, ALSPROJ_OUTDIR, ALS_ARGS$redirect_outputs)
als_register_output_wrappers()
als_load_objects(ALS_ARGS$objects)
on.exit(als_write_session_info(ALSPROJ_SCRIPT_ID), add = TRUE)
# ---- End bootstrap ----

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
