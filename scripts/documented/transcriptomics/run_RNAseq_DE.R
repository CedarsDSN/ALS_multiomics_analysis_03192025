# ============================================================
# Script: run_RNAseq_DE.R
# Project: ALS multi-omics analysis
# Source status: documented/sanitized for GitHub handoff
# ============================================================
# Documentation summary from docs/SCRIPT_SPECS_FULL.md
# 
# ## `run_RNAseq_DE.R`
# 
# **Category:** Transcriptomics / QC / DE
# 
# **Purpose**
# 
# Initial RNA-seq differential-expression and QC workflow. Compares ALS progression classes using DESeq2, including a race-adjusted model and a simpler progression-only model, then produces PCA, sample-distance, lower-PC, and batch-visualization plots.
# 
# **Inputs**
# 
# - `transcriptomics/4_matrix/AnswerALS-651-T-v1-release6_raw-counts.csv` — raw gene counts with `Geneid` column.
# - `transcriptomics/4_matrix/Sample_Mapping_Information/Sample_Mapping_File_Feb_2024.csv` — sample-to-participant map.
# - `ALS_Progression_Rates_Full.csv` — progression class table.
# - `metadata/aals_dataportal_datatable.csv` — clinical metadata.
# 
# **Outputs**
# 
# - `DGE_Fast_vs_Slow_with_race.csv`
# - `DGE_Intermediate_vs_Slow_with_race.csv`
# - `DGE_Fast_vs_Intermediate_with_race.csv`
# - `DGE_Fast_vs_Slow_dispersion_model.csv`
# - `DGE_Intermediate_vs_Slow_dispersion_model.csv`
# - `DGE_Fast_vs_Intermediate_dispersion_model.csv`
# - `PCA_ProgressorClass.pdf`
# - `Sample_Distance_Heatmap.pdf`
# - `Top20_Loadings_<PC>.csv`
# - `Lower_PC_Plots_Top500.pdf`
# - `Batch_visualization.pdf`
# 
# **Upstream dependencies**
# 
# - None. This is an entry-point transcriptomics script, assuming all input files are present.
# 
# **How to run**
# 
# Legacy documented version:
# 
# ```bash
# Rscript scripts/documented/transcriptomics/run_RNAseq_DE.R
# ```
# 
# Refactored command-line version:
# 
# ```bash
# Rscript scripts/refactored/transcriptomics/01_progression_rnaseq_de.R --config config/config.yml --base-dir /path/to/private/ALS --outdir results/transcriptomics/01_progression_rnaseq_de
# ```
# 
# **Refactor notes**
# 
# - Replace `setwd(mainDir)` with config-driven paths.
# - Write outputs under `results/transcriptomics/progression_de/` and `figures/transcriptomics/qc/`.
# - Save the DESeq2 dataset and VST object if downstream QC is reused.
# ============================================================

# Load libraries
library(tidyverse)
library(DESeq2)

library(dplyr)

#mainDir <- file.path(ALSPROJ_BASE_DIR, "transcriptomics_de") # For running on cluster
mainDir <- file.path(ALSPROJ_BASE_DIR) # For running locally
setwd(mainDir)
# Set file paths
# /path/to/private/ALS/transcriptomics/4_matrix/
count_file <- "transcriptomics/4_matrix/AnswerALS-651-T-v1-release6_raw-counts.csv"
mapping_file <- "transcriptomics/4_matrix/Sample_Mapping_Information/Sample_Mapping_File_Feb_2024.csv"
metadata_file <- "ALS_Progression_Rates_Full.csv"
clinical_file <- "metadata/aals_dataportal_datatable.csv"

# Load files
counts <- read.csv(count_file, check.names = FALSE)
counts <- counts[, -1]  # drop the first column
counts <- column_to_rownames(counts, var = "Geneid")
mapping <- read.csv(mapping_file)
progression <- read.csv(metadata_file)
clinical_data <- read.csv(clinical_file)

# Clean up IDs
mapping$Sample_ID <- trimws(mapping$`Sample.ID`)
mapping$Participant_ID <- trimws(mapping$`Participant.ID`)
progression$Participant_ID <- trimws(progression$Participant_ID)
mapping$Participant.ID <- NULL
colnames(mapping)
mapping <- mapping %>%
  filter(grepl("-T$", Sample_ID))

# Step 2: Merge mapping with progression metadata
merged <- mapping %>%
  inner_join(progression, by = "Participant_ID") %>%
  filter(QC_Flag == "OK") %>%
  filter(toupper(Batch) != "BTC")  # Remove Batch Technical controls

# Step 3: Align to count matrix by Participant_ID
colnames(counts) <- gsub("-\\d+-T$", "", colnames(counts))  # If needed
valid_ids <- intersect(merged$Participant_ID, colnames(counts))
counts <- counts[, valid_ids]
merged <- merged[match(valid_ids, merged$Participant_ID), ]
merged$Progressor_Class <- factor(merged$Progressor_Class, levels = c("Slow", "Intermediate", "Fast"))
#merged$Batch <- factor(merged$Batch)

# Filter clinical_data to only relevant participants
clinical_subset <- clinical_data[clinical_data$Participant_ID %in% merged$Participant_ID, ]

# Remove duplicates (keep first if multiple rows per patient)
clinical_subset <- clinical_subset[!duplicated(clinical_subset$Participant_ID), ]

merged <- merge(
  merged,
  clinical_subset[, c("Participant_ID", "RACE")],
  by = "Participant_ID"
)
merged$RACE <- as.factor(merged$RACE)
# Standard model
# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = merged,
                              design = ~ RACE + Progressor_Class)

#
# Filter low count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds$Progressor_Class <- factor(dds$Progressor_Class, levels = c("Slow", "Intermediate", "Fast"))

# Run DESeq2
dds <- DESeq(dds)

# Pairwise comparisons
res_fast_vs_slow <- results(dds, contrast = c("Progressor_Class", "Fast", "Slow"))
res_inter_vs_slow <- results(dds, contrast = c("Progressor_Class", "Intermediate", "Slow"))
res_fast_vs_inter <- results(dds, contrast = c("Progressor_Class", "Fast", "Intermediate"))

# Save results
write.csv(as.data.frame(res_fast_vs_slow), "DGE_Fast_vs_Slow_with_race.csv")
write.csv(as.data.frame(res_inter_vs_slow), "DGE_Intermediate_vs_Slow_with_race.csv")
write.csv(as.data.frame(res_fast_vs_inter), "DGE_Fast_vs_Intermediate_with_race.csv")

# Manual dispersion model
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = merged,
                              design = ~ Progressor_Class)
# Step 2: Set factor levels
dds$Progressor_Class <- factor(dds$Progressor_Class, levels = c("Slow", "Intermediate", "Fast"))

# Step 3: Estimate size factors
dds <- estimateSizeFactors(dds)

# Step 4: Filter low count genes based on normalized counts
nc <- counts(dds, normalized = TRUE)
filter <- rowSums(nc >= 10) >= 4  # You can try >=3 or >=4 if needed
dds <- dds[filter, ]

# Step 5: Estimate dispersions
dds <- estimateDispersions(dds)

# Step 6: Run Wald test with increased max iterations
dds <- nbinomWaldTest(dds, maxit = 500)
dds <- dds[which(mcols(dds)$betaConv), ]

# Pairwise comparisons
res_fast_vs_slow <- results(dds, contrast = c("Progressor_Class", "Fast", "Slow"))
res_inter_vs_slow <- results(dds, contrast = c("Progressor_Class", "Intermediate", "Slow"))
res_fast_vs_inter <- results(dds, contrast = c("Progressor_Class", "Fast", "Intermediate"))

# Save results
write.csv(as.data.frame(res_fast_vs_slow), "DGE_Fast_vs_Slow_dispersion_model.csv")
write.csv(as.data.frame(res_inter_vs_slow), "DGE_Intermediate_vs_Slow_dispersion_model.csv")
write.csv(as.data.frame(res_fast_vs_inter), "DGE_Fast_vs_Intermediate_dispersion_model.csv")


# Conduct PCA and plot PCA and Sample distance heatmap
library(pheatmap)
library(RColorBrewer)

# 1. Variance-stabilizing transformation
vsd <- vst(dds, blind = TRUE)

# ------------------ 1. PCA Plot ------------------
pdf("PCA_ProgressorClass.pdf", width = 6, height = 5)
plotPCA(vsd, intgroup = "Progressor_Class") +
  ggtitle("PCA: ALS Progressor Classes") +
  theme_minimal()
dev.off()

plotPCA(vsd, intgroup = "Progressor_Class") +
  ggtitle("PCA of RNA-seq Samples by Progressor Class")


# ------------------ 2. Optional: Sample Distance Heatmap ------------------
if (ncol(vsd) <= 60) {
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(vsd)
  colnames(sampleDistMatrix) <- NULL
  
  pdf("Sample_Distance_Heatmap.pdf", width = 7, height = 6)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           main = "Sample-to-Sample Distance")
  dev.off()
}

#Sanity checks
print("dds progressor class")
print(table(colData(dds)$Progressor_Class))

print("Mean expression")
print(mean(rowMeans(counts(dds, normalized = TRUE))))

print("how many genes passed filtering")
print(dim(dds))

print("Levels of Progressor Class")
print(levels(dds$Progressor_Class))

print("How many samples per level")
print(table(dds$Progressor_Class))

# Convert results to data frames
res_fast_vs_slow_df <- as.data.frame(res_fast_vs_slow)
res_inter_vs_slow_df <- as.data.frame(res_inter_vs_slow)
res_fast_vs_inter_df <- as.data.frame(res_fast_vs_inter)

# Function to summarize LFCs
lfc_summary <- function(res_df, label) {
  cat("\n", label, "\n")
  print(summary(res_df$log2FoldChange))
  cat("Number of genes with padj < 0.05:", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
}

# Get summaries
lfc_summary(res_fast_vs_slow_df, "Fast vs Slow")
lfc_summary(res_inter_vs_slow_df, "Intermediate vs Slow")
lfc_summary(res_fast_vs_inter_df, "Fast vs Intermediate")


# Find which PC is best for separating progressor class
library(DESeq2)
library(ggplot2)
library(dplyr)

# Step 1: VST transform
vsd <- vst(dds, blind = TRUE)

# Step 2: Run PCA manually on transposed assay
mat <- t(assay(vsd))
pca <- prcomp(mat, scale. = TRUE)

# Step 3: Create PCA data frame
pca_data <- as.data.frame(pca$x)
pca_data$Progressor_Class <- colData(vsd)$Progressor_Class

# Step 4: Run ANOVA for PCs 1 to 5
anova_results <- sapply(1:5, function(i) {
  pc_col <- paste0("PC", i)
  formula <- as.formula(paste0(pc_col, " ~ Progressor_Class"))
  fit <- aov(formula, data = pca_data)
  summary(fit)[[1]][["Pr(>F)"]][1]
})
names(anova_results) <- paste0("PC", 1:5)
print("ANOVA p-values by PC:")
print(anova_results)

# Step 5: Plot the most significant PC
best_pc <- names(anova_results)[which.min(anova_results)]
cat("\nMost significant PC for separating classes:", best_pc, "\n")

ggplot(pca_data, aes_string(x = "Progressor_Class", y = best_pc, fill = "Progressor_Class")) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  theme_minimal() +
  labs(title = paste(best_pc, "vs ALS Progressor Class"),
       y = paste0(best_pc, " Score"),
       x = "Progressor Class")

# Step 6: Save top contributing genes for that PC
loadings <- as.data.frame(pca$rotation[, best_pc, drop = FALSE])
colnames(loadings) <- best_pc
loadings$Gene <- rownames(loadings)

top_genes <- loadings %>%
  arrange(desc(abs(.data[[best_pc]]))) %>%
  head(20)

write.csv(top_genes, paste0("Top20_Loadings_", best_pc, ".csv"), row.names = FALSE)


# Get lower PCs
# Step 1: Extract top 500 most variable genes
top500_var_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)

# Step 2: Run PCA on transposed assay matrix for those genes
pca_top500 <- prcomp(t(assay(vsd)[top500_var_genes, ]), scale. = TRUE)

# Step 3: Create PCA data frame and add metadata
pca_data_500 <- as.data.frame(pca_top500$x)
pca_data_500$Progressor_Class <- colData(vsd)$Progressor_Class

# Define lower PC pairs
pc_pairs <- list(
  c("PC1", "PC3"),
  c("PC1", "PC4"),
  c("PC2", "PC3"),
  c("PC3", "PC4")
)

# Save all to PDF
pdf("Lower_PC_Plots_Top500.pdf", width = 7, height = 6)
for (pair in pc_pairs) {
  x_pc <- pair[1]
  y_pc <- pair[2]
  
  p <- ggplot(pca_data, aes_string(x = x_pc, y = y_pc, color = "Progressor_Class")) +
    geom_point(alpha = 0.7, size = 2) +
    labs(title = paste(x_pc, "vs", y_pc, "(Top 500 genes)")) +
    theme_minimal()
  
  print(p)
}
dev.off()


# Find confounding variables - either with top 500 or full pca object above
pca_data_500$Participant_ID <- rownames(pca_data_500)

merged_pca_data <- merge(pca_data_500, clinical_data, by = "Participant_ID")
rownames(merged_pca_data) <- merged_pca_data$Participant_ID  # Optional: restore as rownames

# Step 4: Select variables you want to screen
# Define variable groups
categorical_vars <- c("Progressor_Class", "SEX", "RACE", "Site_of_Onset")
numeric_vars <- c("AGE_AT_SYMPTOM_ONSET")

# Coerce types
merged_pca_data <- merged_pca_data %>%
  mutate(across(all_of(categorical_vars), as.factor)) %>%
  mutate(across(all_of(numeric_vars), as.numeric))

# Step 5: Function to test association of each variable with a principal component
test_pc_association <- function(pc, data, vars) {
  
  results <- lapply(vars, function(var) {
    x <- data[[var]]
    
    # Skip variables with no variance
    if (is.factor(x) && length(unique(x)) < 2) {
      message("Skipping ", var, ": only one level.")
      return(data.frame(variable = var, p_value = NA))
    }
    if (is.numeric(x) && length(unique(na.omit(x))) < 2) {
      message("Skipping ", var, ": numeric column has no variation.")
      return(data.frame(variable = var, p_value = NA))
    }
    
    # Build formula and fit model
    formula <- as.formula(paste(pc, "~", var))
    model <- tryCatch(lm(formula, data = data), error = function(e) NULL)
    
    if (is.null(model)) {
      message("Model failed for ", var)
      return(data.frame(variable = var, p_value = NA))
    }
    
    # Run ANOVA and extract p-value
    pval <- tryCatch({
      aov_result <- anova(model)
      round(aov_result[["Pr(>F)"]][1], 5)
    }, error = function(e) NA)
    
    return(data.frame(variable = var, p_value = pval))
  })
  
  results_df <- do.call(rbind, results)
  rownames(results_df) <- NULL
  return(results_df)
}

variables_to_test <- c(categorical_vars, numeric_vars)
# Step 6: Run tests for PC1 and PC2
results_pc1 <- test_pc_association("PC1", merged_pca_data, variables_to_test)
results_pc2 <- test_pc_association("PC2", merged_pca_data, variables_to_test)

# Step 7: Display results
cat("\n📊 Confounder Screening Results:\n")
cat("\n--- PC1 ---\n"); print(results_pc1)
cat("\n--- PC2 ---\n"); print(results_pc2)


#Sanity check to check if getting same value as above
model_check <- lm(PC2 ~ Progressor_Class, data = pca_data)
summary(model_check)
anova(model_check)

# Visualize batch to see if there are batch effects
colnames(mapping)[colnames(mapping) == "Participant.ID"] <- "Participant_ID"
pca_data_500 <- merge(pca_data_500, mapping[, c("Participant_ID", "Batch")], by = "Participant_ID")
pca_data_500$Batch <- as.factor(pca_data_500$Batch)

pdf("Batch_visualization.pdf", width = 7, height = 6)
ggplot(pca_data_500, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA: PC1 vs PC2 Colored by Batch", x = "PC1", y = "PC2") +
  theme_minimal()
dev.off()
