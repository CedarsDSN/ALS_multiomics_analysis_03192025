# Analysis to find if healthy and control can be seperated

# Load required libraries
library(DESeq2)
library(ggplot2)
library(readr)
library(tidyr)
library(tidyverse)
library(dplyr)
library(edgeR)
library(matrixStats)
library(biomaRt)

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
mapping_file <- "transcriptomics/4_matrix/Sample_Mapping_Information/Sample_Mapping_File_Feb_2024.csv"
mapping <- read.csv(mapping_file, stringsAsFactors = FALSE)

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
full_metadata$Case_Control <- factor(full_metadata$Case_Control, levels = c("Healthy Control", "Case"))
full_metadata <- left_join(full_metadata, mapping, by = c("Sample_ID" = "Sample.ID"))

# Ensure sample order matches counts
counts <- counts[, full_metadata$Sample_ID]
stopifnot(all(colnames(counts) == full_metadata$Sample_ID))


# ==== STEP 3: Create DESeq2 object and normalize ====
dds <- DESeqDataSetFromMatrix(countData = counts, colData = full_metadata, design = ~ Case_Control)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)

# ==== STEP 5: Plot PCA ====
pdf("PCA_PC1vsPC2_CaseandControl.pdf", width = 10, height = 12)
plotPCA(vsd, intgroup = "Case_Control") +
  ggtitle("PCA: ALS vs control") +
  theme_minimal()
dev.off()

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

# Color by sex - Top 500 variable genes
pdf("PCA_plotPCA_top_500_genes.pdf", width = 10, height = 12)
plotPCA(vsd, intgroup = "Sex") +
  ggtitle("PCA: male vs female (top 500 variable genes)") +
  theme_minimal()
dev.off()
# Color by sex - Top 1000 variable genes  
pdf(file = "PCA_plotPCA_top_1000_genes.pdf", width = 10, height = 12)
plotPCA(vsd, intgroup = "Sex", ntop = 1000)+
  ggtitle("PCA: male vs female (top 1000 genes)") +
  theme_minimal()
dev.off()
# Color by sex - All genes
pdf("PC1vsPC2_sex_all_genes.pdf", width = 10, height = 12)
plotPCA(vsd, intgroup = "Sex", ntop = 50772) +
  ggtitle("PCA: male vs female (all 50772 genes)") +
  theme_minimal()
dev.off()

#Batch in PCA
# Color by Batch - all genes
pdf("PCA_batch_all_genes.pdf", width = 10, height = 12)
plotPCA(vsd, intgroup = "Batch", ntop = "50772") +
  ggtitle("PCA: Colored by Batch") +
  theme_minimal()
dev.off()

# Try a different function to confirm pattern in PC1 vs PC2

# Step 1: VST transformation
vsd <- vst(dds, blind = TRUE)

# Step 2: Select top 500 most variable genes
ntop <- 500
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(ntop)]
mat <- t(assay(vsd)[select, ])

# Step 3: PCA
pca_500 <- prcomp(mat, scale. = FALSE)
percentVar <- round(100 * (pca_500$sdev^2 / sum(pca_500$sdev^2)), 1)

# Step 4: Prepare data
pca_data <- as.data.frame(pca_500$x)
pca_data$Sex <- colData(vsd)$Sex  # Replace with your grouping variable
pca_data$Case_Control <- colData(vsd)$Case_Control  # Add your case/control column
pca_data$Sample <- rownames(pca_data)

pdf("PC1vsPC2_CaseControl_prcomp_500.pdf", width = 10, height = 12)
# Step 5: Plot (styled like plotPCA)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: PC1 vs PC2 (Top 500 Genes)") +
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
dev.off()

cat("Total genes in dataset:", nrow(vsd), "\n")
cat("Top 500 genes: 500\n")
cat("Difference:", nrow(vsd) - 500, "additional genes\n")

# PCA with all genes (50772 genes)
# VST transformation
vsd <- vst(dds, blind = TRUE)

# Use ALL genes (no selection step)
mat_all <- t(assay(vsd))  # Transpose: samples as rows, genes as columns

# PCA (matching DESeq2's method - center but don't scale)
pca_all <- prcomp(mat_all, center = TRUE, scale. = FALSE)

# Calculate variance explained
percentVar_all <- round(100 * (pca_all$sdev^2 / sum(pca_all$sdev^2)), 1)

# Prepare data for plotting
pca_data_all <- as.data.frame(pca_all$x)
pca_data_all$Sex <- colData(vsd)$Sex
pca_data_all$Case_Control <- colData(vsd)$Case_Control  # Add your case/control column
pca_data_all$Sample <- rownames(pca_data_all)

# Plot
pdf(file = "PCA_All_Genes_prcomp.pdf", width = 10, height = 12)
ggplot(pca_data_all, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_all[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_all[2], "% variance")) +
  ggtitle("PCA: All Genes") +
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
dev.off()



# Stringent filtering - Does this change the PCA pattern in PC1 vs PC2?
dds <- DESeqDataSetFromMatrix(countData = counts, colData = full_metadata, design = ~ Case_Control)
smallestGroupSize <- 4 
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize 
dds <- dds[keep,] 
nrow(dds)
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)

pdf("PCA_Sex_stringent.pdf", width = 10, height = 12)
plotPCA(vsd, intgroup = "Sex") +
  ggtitle("PCA: ALS vs control") +
  theme_minimal()
dev.off()

# PCA with edgeR log-CPM transformation - Does a different normalization method change the PCA pattern in PC1 vs PC2?

# Step 1: Normalize and log-transform
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

# Log-CPM transformation
logCPM <- cpm(dge, log = TRUE, prior.count = 1)

# Step 2: Get top 500 most variable genes
vars_cpm <- apply(logCPM, 1, var)
top500_cpm <- head(order(vars_cpm, decreasing = TRUE), 500)
logCPM_top <- logCPM[top500_cpm, ]

# Step 3: Run PCA
pca_cpm <- prcomp(t(logCPM_top), scale. = FALSE)
pca_df_cpm <- as.data.frame(pca_cpm$x)

# Step 4: Add metadata
pca_df_cpm$Sample <- rownames(pca_df_cpm)
pca_df_cpm <- merge(pca_df_cpm, full_metadata, by.x = "Sample", by.y = "Sample_ID")

# Step 5: Calculate percent variance
percentVar_cpm <- round(100 * (pca_cpm$sdev^2 / sum(pca_cpm$sdev^2))[1:2], 1)

# Step 6: Plot with percent variance
pdf("PCA_Top500_Variable_Genes_edgeR_logCPM.pdf", width = 10, height = 12)
ggplot(pca_df_cpm, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: PC1 vs PC2 (Top 500 Genes)") +
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
dev.off()



# Empirical evidence to show PCA differs with top 500 variable genes and all genes
summary(pca_500)$importance["Proportion of Variance", 1:5]
summary(pca_all)$importance["Proportion of Variance", 1:5]

plot(100 * pca_500$sdev^2 / sum(pca_500$sdev^2), type = "b", col = "blue", main = "Top 500 Genes", ylab = "% Variance")
plot(100 * pca_all$sdev^2 / sum(pca_all$sdev^2), type = "b", col = "red", main = "All Genes", ylab = "% Variance")

pdf("PCA_Top500_vs_All_Genes_Variance.pdf", width = 15, height = 12)
cum_var_500 <- cumsum(pca_500$sdev^2) / sum(pca_500$sdev^2)
cum_var_all <- cumsum(pca_all$sdev^2) / sum(pca_all$sdev^2)

plot(cum_var_500[1:50], type = "l", col = "blue", lwd = 2,
     ylab = "Cumulative Variance Explained", xlab = "PC", main = "Top 500 vs All Genes")
lines(cum_var_all[1:50], col = "red", lwd = 2)
legend("bottomright", legend = c("Top 500", "All Genes"), col = c("blue", "red"), lwd = 2)
dev.off()

cor(pca_500$x[, 1], pca_all$x[, 1])  # How similar are PC1 coordinates?
cor(pca_500$x[, 2], pca_all$x[, 2])


#Conduct DE analysis
dds <- DESeqDataSetFromMatrix(countData = counts, colData = full_metadata, design = ~ Sex + Case_Control)
smallestGroupSize <- 4 
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize 
dds <- dds[keep,] 
nrow(dds)
dds <- DESeq(dds)
res_case_vs_control <- results(dds, contrast = c("Case_Control", "Case", "Healthy Control"))

#res_case_vs_control <- results(dds, name = "Case_Control_Healthy.Control_vs_Case")
res_sig <- res_case_vs_control[which(res_case_vs_control$padj < 0.05), ]
summary(res_sig)

# Remove NAs first
res_filtered <- res_case_vs_control[!is.na(res_case_vs_control$padj), ]

# Apply both padj and log2 fold change thresholds
res_sig_strict <- res_filtered[
  res_filtered$padj < 0.05 & abs(res_filtered$log2FoldChange) > 1, ]

# View number of genes passing both filters
nrow(res_sig_strict)
summary(res_sig_strict)


#Get PCA loadings
loadings <- pca_500$rotation   # genes × PCs

# Top contributors to PC1 and PC2
top_PC1_genes <- names(sort(abs(loadings[, 1]), decreasing = TRUE)[1:20])
top_PC2_genes <- names(sort(abs(loadings[, 2]), decreasing = TRUE)[1:20])

# Check how many genes in the PCA loadings are from sex chromosomes
# Connect to Ensembl (human dataset)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Your Ensembl gene IDs
gene_ids <- unique(c(top_PC1_genes, top_PC2_genes))

# Query Ensembl for chromosome and gene name
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = ensembl
)

# Optional: highlight sex chromosome genes
gene_info$is_sex_chromosome <- gene_info$chromosome_name %in% c("X", "Y")

# View
print(gene_info)

# Remove sex genes and check PCA

# Get chromosome info
gene_chr <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name"),
  filters = "ensembl_gene_id",
  values = rownames(dds),
  mart = ensembl
)

# Filter genes on chrX or chrY
sex_chr_genes <- gene_chr$ensembl_gene_id[gene_chr$chromosome_name %in% c("X", "Y")]

# Remove sex chromosome genes
non_sex_genes <- setdiff(rownames(dds), sex_chr_genes)
dds_no_sex <- dds[non_sex_genes, ]
vsd_no_sex <- vst(dds_no_sex, blind = TRUE)

# Same top 500 variance selection
ntop = 500
rv2 <- rowVars(assay(vsd_no_sex))
select2 <- order(rv2, decreasing = TRUE)[seq_len(ntop)]
mat_no_sex <- t(assay(vsd_no_sex)[select2, ])

pca_no_sex <- prcomp(mat_no_sex, scale. = FALSE)

percentVar <- round(100 * (pca_no_sex$sdev^2 / sum(pca_no_sex$sdev^2)), 1)

# Step 4: Prepare data
pca_no_sex <- as.data.frame(pca_no_sex$x)
pca_no_sex$Sex <- colData(vsd)$Sex  # Replace with your grouping variable
pca_no_sex$Case_Control <- colData(vsd)$Case_Control  # Add your case/control column
pca_no_sex$Sample <- rownames(pca_no_sex)
#pca_no_sex <- merge(pca_no_sex, full_metadata, by.x = "Sample", by.y = "Sample_ID")

pdf("PC1vsPC2_no_sex_chromosomes.pdf", width = 10, height = 12)
# Step 5: Plot (styled like plotPCA)
ggplot(pca_no_sex, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: PC1 vs PC2 (Top 500 Genes)") +
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
dev.off()


# PCA association study with binary variables
pca_df <- pca_no_sex
# Merge with Participant_ID from your full_metadata
pca_df$Participant_ID <- full_metadata$Participant_ID[match(pca_df$Sample, full_metadata$Sample_ID)]

# --- 2. Load Variable Map ---
binary_vars <- read_csv("binary_variables.csv")
names(binary_vars) <- make.names(names(binary_vars))
binary_vars <- binary_vars %>%
  mutate(
    Form.Name = gsub(" ", "_", trimws(Form.Name)),
    Field.Name = trimws(Field.Name),
    Field.Label = trimws(Field.Label)
  )
data_dir <- path.expand("~/Documents/DSN/ALS/metadata/clinical")


test_var_association <- function(form_name, field_name, field_label, pca_df, data_dir) {
  file_path <- file.path(data_dir, paste0(form_name, ".csv"))
  
  if (!file.exists(file_path)) {
    message("[MISSING FILE] ", file_path)
    return(NULL)
  }
  
  table_df <- tryCatch(
    read.csv(file_path, stringsAsFactors = FALSE),
    error = function(e) {
      message("[READ ERROR] ", file_path, " -- ", e$message)
      return(NULL)
    }
  )
  if (is.null(table_df)) return(NULL)
  
  if (!("Participant_ID" %in% colnames(table_df)) || !(field_name %in% colnames(table_df))) {
    message("[MISSING COLUMN] ", field_name, " in ", form_name)
    return(NULL)
  }
  
  # Merge with PCA
  merged <- left_join(pca_df, table_df[, c("Participant_ID", field_name)], by = "Participant_ID")
  x <- merged[[field_name]]
  
  if (all(is.na(x))) {
    message("[ALL NA] ", field_name, " in ", form_name)
    return(NULL)
  }
  
  if (length(unique(na.omit(x))) < 2) {
    message("[NO VARIABILITY] ", field_name, " in ", form_name)
    return(NULL)
  }
  
  # Initialize result
  res <- list(
    Variable = field_label,
    Form = form_name,
    Field = field_name,
    Type = class(x),
    p_PC1 = NA,
    p_PC2 = NA
  )
  
  # Run association tests
  if (is.numeric(x)) {
    res$p_PC1 <- tryCatch(summary(lm(PC1 ~ x, data = merged))$coefficients[2, 4], error = function(e) NA)
    res$p_PC2 <- tryCatch(summary(lm(PC2 ~ x, data = merged))$coefficients[2, 4], error = function(e) NA)
  } else {
    x <- factor(x)
    if (length(unique(x)) > 1) {
      res$p_PC1 <- tryCatch(summary(aov(PC1 ~ x, data = merged))[[1]][["Pr(>F)"]][1], error = function(e) NA)
      res$p_PC2 <- tryCatch(summary(aov(PC2 ~ x, data = merged))[[1]][["Pr(>F)"]][1], error = function(e) NA)
    }
  }
  
  return(as.data.frame(res))
}

# --- 5. Loop Over All Variables and Collect Results ---
results_list <- lapply(1:nrow(binary_vars), function(i) {
  row <- binary_vars[i, ]
  test_var_association(
    form_name = row$Form.Name,
    field_name = row$Field.Name,
    field_label = row$Field.Label,
    pca_df = pca_df,
    data_dir = data_dir
  )
})

# --- 6. Combine, Rank, and Save ---
results_df <- bind_rows(results_list) %>%
  mutate(min_p = pmin(p_PC1, p_PC2, na.rm = TRUE)) %>%
  arrange(min_p)

# Optional: Adjust for FDR
results_df$FDR_PC1 <- p.adjust(results_df$p_PC1, method = "fdr")
results_df$FDR_PC2 <- p.adjust(results_df$p_PC2, method = "fdr")

# Save
write_csv(results_df, "PCA_PC1_PC2_variable_associations.csv")

# Top 4 variables with strongest association
top_vars <- results_df %>%
  filter(!is.na(min_p)) %>%
  arrange(min_p) %>%
  slice(1:4)

# Bar plot of top 15 associations
top_bar <- results_df %>%
  filter(!is.na(min_p)) %>%
  arrange(min_p) %>%
  slice(1:15) %>%
  mutate(Variable_Label = paste0(Variable, "\n(", Form, "::", Field, ")"))

pdf("Top_15_Variable_Associations_PCA.pdf", width = 10, height = 8)
ggplot(top_bar, aes(x = reorder(Variable_Label, min_p), y = -log10(min_p))) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 15 Variables Associated with PC1 or PC2",
    x = "Variable (Form::Field)",
    y = expression(-log[10](p-value))
  ) +
  theme_minimal(base_size = 14)
dev.off()
