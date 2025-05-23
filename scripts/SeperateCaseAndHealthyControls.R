# Analysis to find if healthy and control can be seperated

# Load required libraries
library(DESeq2)
library(ggplot2)
library(readr)
library(tidyr)
library(tidyverse)

mainDir <- "~/Documents/DSN/ALS/" # For running locally
setwd(mainDir)

# ==== STEP 1: Load Data ====
# RNA-seq count matrix
count_file <- "transcriptomics/4_matrix/AnswerALS-651-T-v1-release6_raw-counts.csv"
counts <- read.csv(count_file, check.names = FALSE)
counts <- counts[, -1]  # drop the first column
counts <- column_to_rownames(counts, var = "Geneid")

# Participant metadata
participants <- read.csv("metadata/aals_participants.csv")

# ==== STEP 2: Filter to samples with RNA-seq data ====
sample_ids <- colnames(counts)
participant_ids <- sub("^((CTRL|CASE)-[A-Z0-9]+)-.*$", "\\1", sample_ids)
# Add Sample_ID and Participant_ID to metadata
metadata <- participants
metadata$Participant_ID <- as.character(metadata$Participant_ID)

# Create a metadata frame matching the counts
count_metadata <- data.frame(
  Sample_ID = sample_ids,
  Participant_ID = participant_ids,
  stringsAsFactors = FALSE
)

# Merge metadata from participants file onto count metadata
full_metadata <- merge(count_metadata, metadata, by = "Participant_ID")

# Ensure sample order matches counts
counts <- counts[, full_metadata$Sample_ID]
stopifnot(all(colnames(counts) == full_metadata$Sample_ID))


# ==== STEP 3: Create DESeq2 object and normalize ====
dds <- DESeqDataSetFromMatrix(countData = counts, colData = full_metadata, design = ~ Case_Control)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)

# ==== STEP 4: PCA ====
pca <- prcomp(t(assay(vsd)), scale. = TRUE)
pca_data <- as.data.frame(pca$x)
# 1. Extract Participant_ID from rownames of pca_data
pca_data$Participant_ID <- sub("^((CTRL|CASE)-[A-Z0-9]+)-.*$", "\\1", rownames(pca_data))

# 2. Make sure Participant_IDs in metadata are character
metadata$Participant_ID <- as.character(metadata$Participant_ID)

# 3. Match Case_Control using Participant_ID
pca_data$Case_Control <- metadata$Case_Control[match(pca_data$Participant_ID, metadata$Participant_ID)]

# 4. Check it worked
table(pca_data$Case_Control, useNA = "ifany")


# ==== STEP 5: Plot PCA ====
pdf("PCA_CasevsControl.pdf", width = 10, height = 12)
plotPCA(vsd, intgroup = "Case_Control") +
  ggtitle("PCA: ALS vs control") +
  theme_minimal()
dev.off()

# ==== Optional: Save variance explained ====
explained <- summary(pca)$importance[2, 1:5]
write.csv(explained, "PCA_Variance_Explained.csv")


pvals <- sapply(1:10, function(i) {
  pc_col <- paste0("PC", i)
  model <- lm(pca_data[[pc_col]] ~ pca_data$Case_Control)
  summary(model)$coefficients[2, 4]  # p-value
})
names(pvals) <- paste0("PC", 1:10)
print(round(pvals, 5))

# Sanity check with ggplot
pca_data <- plotPCA(vsd, intgroup = "Case_Control", returnData = TRUE)

percentVar <- round(attr(pca_data, "percentVar") * 100, 1)

pdf("PC1vsPC2_casevscontrol.pdf", width = 10, height = 12)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Case_Control)) +
  geom_point(size = 2, alpha = 0.8) +
  xlab(paste0("PC1 (", percentVar[1], "% variance)")) +
  ylab(paste0("PC2 (", percentVar[2], "% variance)")) +
  ggtitle("PCA: Case vs Control (PC1 vs PC2)") +  # ✅ Add your title here
  theme_minimal()
dev.off()

# Plot lower dimensional PCs
library(matrixStats)

# Step 1: Top 500 most variable genes
ntop <- 500
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[1:ntop]

# Step 2: Create the matrix used by plotPCA()
mat <- t(assay(vsd)[select, ])

# Step 3: Run PCA manually on this matrix
pca_result <- prcomp(mat, center = TRUE, scale. = TRUE)

# Step 4: Extract % variance
percentVar <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)

# Step 5: Build PCA data frame
pca_data <- as.data.frame(pca_result$x)
pca_data$Case_Control <- colData(vsd)$Case_Control

# Step 6: Plot PC1 vs PC3
pdf("PC1vsPC3_casevscontrol.pdf", width = 10, height = 12)
ggplot(pca_data, aes(x = PC1, y = PC3, color = Case_Control)) +
  geom_point(size = 2, alpha = 0.8) +
  xlab(paste0("PC1 (", percentVar[1], "% variance)")) +
  ylab(paste0("PC3 (", percentVar[3], "% variance)")) +
  ggtitle("PCA: Case vs Control (PC1 vs PC3)") +  # ✅ Add your title here
  theme_minimal()
dev.off()


pdf("PC1vsPC4_casevscontrol.pdf", width = 10, height = 12)
ggplot(pca_data, aes(x = PC1, y = PC4, color = Case_Control)) +
  geom_point(size = 2, alpha = 0.8) +
  xlab(paste0("PC1 (", percentVar[1], "% variance)")) +
  ylab(paste0("PC4 (", percentVar[4], "% variance)")) +
  ggtitle("PCA: Case vs Control (PC1 vs PC4)") +  # ✅ Add your title here
  theme_minimal()
dev.off()

pdf("PC2vsPC3_casevscontrol.pdf", width = 10, height = 12)
ggplot(pca_data, aes(x = PC2, y = PC3, color = Case_Control)) +
  geom_point(size = 2, alpha = 0.8) +
  xlab(paste0("PC2 (", percentVar[2], "% variance)")) +
  ylab(paste0("PC3 (", percentVar[3], "% variance)")) +
  ggtitle("PCA: Case vs Control (PC2 vs PC3)") +  # ✅ Add your title here
  theme_minimal()
dev.off()

pdf("PC3vsPC4_casevscontrol.pdf", width = 10, height = 12)
ggplot(pca_data, aes(x = PC3, y = PC4, color = Case_Control)) +
  geom_point(size = 2, alpha = 0.8) +
  xlab(paste0("PC3 (", percentVar[3], "% variance)")) +
  ylab(paste0("PC4 (", percentVar[4], "% variance)")) +
  ggtitle("PCA: Case vs Control (PC3 vs PC4)") +  # ✅ Add your title here
  theme_minimal()
dev.off()


pdf("PC1vsPC2_sex.pdf", width = 10, height = 12)
plotPCA(vsd, intgroup = "Sex") +
  ggtitle("PCA: male vs female") +
  theme_minimal()
dev.off()


# PCA plot for smoker vs non-smoker
# Load packages
library(DESeq2)
library(matrixStats)
library(ggplot2)

# Step 1: Read clinical metadata
clinical <- read.csv("metadata/clinical/Environmental_Questionnaire.csv")  # Change to your file path
clinical$Participant_ID <- as.character(clinical$Participant_ID)

# Step 2: Extract expression matrix
expr <- assay(vsd)  # or vst
ntop <- 500
rv <- rowVars(expr)
top500 <- order(rv, decreasing = TRUE)[1:ntop]
mat <- t(expr[top500, ])

# Step 3: Extract sample IDs
# Assuming your colnames are like "CASE-NEUABC123-XXXX-T"
sample_ids <- sub("^((CTRL|CASE)-[A-Z0-9]+)-.*$", "\\1", colnames(vsd))

# Step 4: Match metadata
clinical_matched <- clinical[match(sample_ids, clinical$Participant_ID), ]
rownames(clinical_matched) <- colnames(vsd)  # Now same order as expression matrix

# Step 5: Run PCA
pca <- prcomp(mat, scale. = TRUE)
pca <- prcomp(t(assay(vsd)), scale. = TRUE)
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
pca_df <- as.data.frame(pca$x)
pca_df$Participant_ID <- sub("^((CTRL|CASE)-[A-Z0-9]+)-.*$", "\\1", rownames(pca_df))

# 2. Make sure Participant_IDs in metadata are character
metadata$Participant_ID <- as.character(metadata$Participant_ID)

# 3. Match Case_Control using Participant_ID
pca_df$Case_Control <- metadata$Case_Control[match(pca_df$Participant_ID, metadata$Participant_ID)]

#pca_df$Case_Control <- clinical_matched$Case_Control
pca_df$Case_Control <- metadata$Case_Control[match(pca_df$Participant_ID, metadata$Participant_ID)]
pca_df$smoker <- clinical$smokerb[match(pca_df$Participant_ID, clinical$Participant_ID)]
full_metadata$smoker <- clinical$smokerb[match(rownames(full_metadata), clinical$Participant_ID)]

# Convert smoker column to factor with clear labels
pca_df$smoker <- factor(pca_df$smoker, levels = c(1, 2), labels = c("Smoker", "Non-smoker"))

# Make the PCA plot
pdf("PCA_Smoker_vs_NonSmoker.pdf", width = 10, height = 12)
ggplot(pca_df, aes(x = PC1, y = PC2, color = smoker)) +
  geom_point(size = 2) +
  labs(
    title = "PCA: Smoker vs Non-smoker",
    x = paste0("PC1 (", percentVar[1], "% variance)"),
    y = paste0("PC2 (", percentVar[2], "% variance)"),
    color = "Smoking Status"
  ) +
  theme_minimal()
dev.off()
