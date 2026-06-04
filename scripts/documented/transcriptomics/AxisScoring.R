# ============================================================
# Script: AxisScoring.R
# Project: ALS multi-omics analysis
# Source status: documented/sanitized for GitHub handoff
# ============================================================
# Documentation summary from scripts/README.md
# 
# ## `AxisScoring.R`
# 
# **Category:** Transcriptomics / axis scoring
# 
# **Purpose**
# 
# Computes Hallmark-derived metabolic, inflammation, and cell-cycle axis scores using transformed RNA-seq expression, removes mutation-confounded surrogate variables, and tests axis associations with mutation subtype, sporadic ALS, NfL, disease duration, ALSFRS-R progression, sex, and other clinical/marker covariates.
# 
# **Inputs**
# 
# - `transcriptomics/4_matrix/AnswerALS-651-T-v1-release6_raw-counts.csv`.
# - `transcriptomics/4_matrix/Sample_Mapping_Information/Sample_Mapping_File_Feb_2024.csv`.
# - `metadata/clinical/aals_dataportal_datatable.csv`.
# - `ALS_Mutation_Subtypes_with_ALS_Control.csv`.
# - `answerals/NfL/NfL_Seimens_Atellica_Immunoassay_Proteomics_2022_07.csv`.
# - MSigDB Hallmark sets via `msigdbr`.
# 
# **Outputs**
# 
# - Figures: `axis_boxplots_adjusted.pdf`, `metabolic_axis_vs_progression.pdf`, `metabolic_quartiles_vs_progression.pdf`, `metabolic_partial_effect_vs_progression_annotated.pdf`, `axis_metabolic_sex.pdf`, plus one covariate-specific plot created through `outfile`.
# - Tables: `axis_covariate_screen_continuous_ALS_only.csv`, `axis_covariate_screen_clinical_ALS_only.csv`, `axis_covariate_screen_marker_ALS_only.csv`, `axis_covariate_screen_top_hits_table.csv`, `axis_covariate_screen_interactions_ALS_only.csv`.
# - Recommended new intermediates: `axis_scores.rds`, `df_adj.rds`, `df_nfl_uniq.rds`, `expr_symbol.rds`, `expr_clean_mat.rds`, `counts_filtered.rds`, `annot_als_filtered.rds`, `svobj.rds`, `clean_svs.rds`, `axis_sets.rds`.
# 
# **Upstream dependencies**
# 
# - None beyond data files. This should become the main upstream producer for sporadic subtyping, PRSet, eQTL, LASSO, IPA, and sensitivity scripts.
# 
# **How to run**
# 
# Legacy documented version:
# 
# ```bash
# Rscript scripts/documented/transcriptomics/AxisScoring.R
# ```
# 
# Refactored command-line version:
# 
# ```bash
# Rscript scripts/refactored/transcriptomics/03_axis_scoring.R --config config/config.yml --base-dir /path/to/private/ALS --outdir analysis/batch_03_axis_scoring_sporadic_subtyping/results/axis_scoring
# ```
# 
# **Refactor notes**
# 
# - Add explicit `saveRDS()` calls for all shared objects.
# - Move ssGSEA axis gene definitions into `R/axis_definitions.R` so LASSO/genomics/PRS use exactly the same gene sets.
# - Create one table `results/axis_scores/axis_scores_with_metadata.csv` for collaborators and one RDS for downstream code.
# ============================================================

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

library(sva)
library(DESeq2)
library(dplyr)

# Model with your biological variables of interest
mod <- model.matrix(~ Mutation_Subtype + SEX, data = annot_als_filtered)
mod0 <- model.matrix(~ SEX, data = annot_als_filtered)
# ===== STEP 1: RUN SVA WITH n=10 =====
cat("Running SVA with n.sv = 10...\n")

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

stopifnot(exists("dds_filtered"))
stopifnot("group" %in% colnames(colData(dds_filtered)))  # change if your column is named differently

# Ensure group is a factor
colData(dds_filtered)$group <- as.factor(colData(dds_filtered)$group)

# -----------------------------
# 1) VST expression matrix
# -----------------------------
vsd <- vst(dds_filtered, blind = FALSE)
expr <- assay(vsd)  # genes x samples

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# expr: genes x samples (rownames are ENSEMBL IDs)
ens <- rownames(expr)

# If your ENSEMBL IDs sometimes have version suffix like ENSG... .12
ens_clean <- sub("\\..*$", "", ens)

# Map
map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ens_clean,
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
)

# Keep 1:1 mapping (drop duplicates + NAs)
map <- map[!is.na(map$SYMBOL), ]
map <- map[!duplicated(map$ENSEMBL), ]

# Subset expr to mapped genes and replace rownames with SYMBOL
keep_ens <- map$ENSEMBL
expr_mapped <- expr[match(keep_ens, ens_clean), , drop = FALSE]
rownames(expr_mapped) <- map$SYMBOL

# If multiple ENSEMBL map to same SYMBOL, collapse by mean (recommended)
expr_df <- as.data.frame(expr_mapped)
expr_df$SYMBOL <- rownames(expr_mapped)

expr_symbol <- aggregate(. ~ SYMBOL, data = expr_df, FUN = mean)
rownames(expr_symbol) <- expr_symbol$SYMBOL
expr_symbol$SYMBOL <- NULL

expr_symbol <- as.matrix(expr_symbol)  # genes(symbol) x samples

suppressPackageStartupMessages({
  library(DESeq2)
  library(msigdbr)
  library(dplyr)
  library(GSVA)
  library(ggplot2)
  library(tidyr)
})
# -----------------------------
# 2) Build Hallmark gene sets
# -----------------------------
hallmark_tbl <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  distinct()

hallmark_list <- split(hallmark_tbl$gene_symbol, hallmark_tbl$gs_name)

# -----------------------------
# 3) Define axis_sets (THIS answers your question)
# -----------------------------
# Each axis is a union of multiple Hallmark gene sets.
axis_sets <- list(
  AXIS_METABOLIC = unique(c(
    hallmark_list[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],
    hallmark_list[["HALLMARK_FATTY_ACID_METABOLISM"]],
    hallmark_list[["HALLMARK_GLYCOLYSIS"]]  # optional but useful
  )),
  
  AXIS_CELL_CYCLE = unique(c(
    hallmark_list[["HALLMARK_E2F_TARGETS"]],
    hallmark_list[["HALLMARK_G2M_CHECKPOINT"]],
    hallmark_list[["HALLMARK_MYC_TARGETS_V1"]],
    hallmark_list[["HALLMARK_DNA_REPAIR"]]  # optional; aligns with your slides
  )),
  
  AXIS_INFLAMMATION = unique(c(
    hallmark_list[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
    hallmark_list[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
    hallmark_list[["HALLMARK_COMPLEMENT"]]
  ))
)

# -----------------------------
# 4) Coverage check (recommended)
# -----------------------------
coverage <- sapply(axis_sets, function(gs) mean(gs %in% rownames(expr_symbol)))
print(round(coverage, 3))


library(GSVA)

# Make sure it's a numeric matrix
expr_mat <- as.matrix(expr_symbol)
storage.mode(expr_mat) <- "double"

# Build ssGSEA parameter object
param <- ssgseaParam(expr_mat, axis_sets, normalize = TRUE)

# Run ssGSEA
axis_scores <- gsva(param)
axis_df <- as.data.frame(t(axis_scores))  # samples x axes
axis_df <- as.data.frame(scale(axis_df))  # z-score each axis across samples
axis_df$sample_id <- rownames(axis_df)

meta <- as.data.frame(colData(dds_filtered))
meta$sample_id <- rownames(meta)

axis_merged <- merge(meta, axis_df, by = "sample_id")

table(axis_merged$group)
print(head(axis_merged[, c("sample_id", "group", "SEX", "AXIS_METABOLIC", "AXIS_CELL_CYCLE", "AXIS_INFLAMMATION")]))

par(mfrow = c(1,3))
boxplot(AXIS_METABOLIC ~ group, data = axis_merged, main="Metabolic axis (z)", las=2)
boxplot(AXIS_CELL_CYCLE ~ group, data = axis_merged, main="Cell-cycle axis (z)", las=2)
boxplot(AXIS_INFLAMMATION ~ group, data = axis_merged, main="Inflammation axis (z)", las=2)
par(mfrow = c(1,1))


plot(axis_merged$AXIS_METABOLIC, axis_merged$AXIS_CELL_CYCLE,
     col = as.integer(axis_merged$group), pch = 16,
     xlab = "Metabolic axis (z)", ylab = "Cell-cycle axis (z)")
legend("topright", legend = levels(axis_merged$group),
       col = seq_along(levels(axis_merged$group)), pch = 16, bty = "n")

summary(aov(AXIS_METABOLIC ~ group, data = axis_merged))
summary(aov(AXIS_CELL_CYCLE ~ group, data = axis_merged))
summary(aov(AXIS_INFLAMMATION ~ group, data = axis_merged))

als_only <- subset(axis_merged, group != "Control")

summary(aov(AXIS_METABOLIC ~ group, data = als_only))
summary(aov(AXIS_CELL_CYCLE ~ group, data = als_only))
summary(aov(AXIS_INFLAMMATION ~ group, data = als_only))


# axis_scores is axes x samples from gsva(param)
axis_raw <- as.data.frame(t(axis_scores))   # samples x axes
axis_raw$sample_id <- rownames(axis_raw)

meta <- as.data.frame(colData(dds_filtered))
meta$sample_id <- rownames(meta)

df <- merge(meta, axis_raw, by="sample_id")

# Compute delta vs control mean (per axis)
axis_names <- c("AXIS_METABOLIC","AXIS_CELL_CYCLE","AXIS_INFLAMMATION")

ctrl_means <- colMeans(df[df$group=="Control", axis_names], na.rm=TRUE)

for (ax in axis_names) {
  df[[paste0(ax, "_delta")]] <- df[[ax]] - ctrl_means[[ax]]
}

summary(aov(AXIS_METABOLIC_delta ~ group, data=df))
summary(aov(AXIS_CELL_CYCLE_delta ~ group, data=df))
summary(aov(AXIS_INFLAMMATION_delta ~ group, data=df))


# Compare variance of delta scores across groups (simple)
tapply(df$AXIS_METABOLIC_delta, df$group, var, na.rm=TRUE)
tapply(df$AXIS_CELL_CYCLE_delta, df$group, var, na.rm=TRUE)
tapply(df$AXIS_INFLAMMATION_delta, df$group, var, na.rm=TRUE)

var.test(df$AXIS_METABOLIC_delta[df$group=="Sporadic"],
         df$AXIS_METABOLIC_delta[df$group=="Control"])

# Control-based thresholds
axis_names <- c("AXIS_METABOLIC","AXIS_CELL_CYCLE","AXIS_INFLAMMATION")

thr <- setNames(
  sapply(axis_names, function(ax) quantile(df[df$group=="Control", ax], 0.90, na.rm = TRUE)),
  axis_names
)

thr

for (ax in axis_names) {
  df[[paste0(ax, "_high")]] <- df[[ax]] > thr[[ax]]
}

# How many "high" in each group?
with(df, table(group, AXIS_METABOLIC_high))
with(df, table(group, AXIS_CELL_CYCLE_high))
with(df, table(group, AXIS_INFLAMMATION_high))

# Use delta or raw; delta is easier to interpret
lm_met <- lm(AXIS_METABOLIC_delta ~ group + SEX + SV1 + SV2 + SV3 + SV4 + SV5 + SV6, data=df)
summary(lm_met)

lm_cell <- lm(AXIS_CELL_CYCLE_delta ~ group + SEX + SV1 + SV2 + SV3 + SV4 + SV5 + SV6, data=df)
summary(lm_cell)

lm_inf <- lm(AXIS_INFLAMMATION_delta ~ group + SEX + SV1 + SV2 + SV3 + SV4 + SV5 + SV6, data=df)
summary(lm_inf)

als <- df[df$group != "Control", ]

plot(als$AXIS_METABOLIC_delta, als$AXIS_CELL_CYCLE_delta,
     pch=16, col=ifelse(als$AXIS_METABOLIC_high, "red", "black"),
     xlab="Metabolic axis (delta vs Control)", ylab="Cell-cycle axis (delta vs Control)")
legend("topright", legend=c("Metabolic-high","Not high"), pch=16, col=c("red","black"), bty="n")


# Restrict to Sporadic vs Control
sub_sc <- droplevels(subset(df, group %in% c("Sporadic","Control")))

# Metabolic-high enrichment
tab_met <- table(sub_sc$group, sub_sc$AXIS_METABOLIC_high)
tab_met
fisher.test(tab_met)

# Cell-cycle-high enrichment
tab_cell <- table(sub_sc$group, sub_sc$AXIS_CELL_CYCLE_high)
tab_cell
fisher.test(tab_cell)

# Inflammation-high enrichment
tab_inf <- table(sub_sc$group, sub_sc$AXIS_INFLAMMATION_high)
tab_inf
fisher.test(tab_inf)

# Build covariate design (do NOT include group)
covar <- model.matrix(~ SEX + SV1 + SV2 + SV3 + SV4 + SV5 + SV6, data = meta)

# expr_symbol: genes x samples
# Residualize per gene: remove covariate effects
expr_resid <- t(apply(expr_symbol, 1, function(y) {
  fit <- lm(y ~ covar - 1)   # -1 because covar already includes intercept
  resid(fit)
}))

# Keep same dimnames
rownames(expr_resid) <- rownames(expr_symbol)
colnames(expr_resid) <- colnames(expr_symbol)

library(GSVA)

expr_resid_mat <- as.matrix(expr_resid)
storage.mode(expr_resid_mat) <- "double"

param2 <- ssgseaParam(expr_resid_mat, axis_sets, normalize = TRUE)
axis_scores_adj <- gsva(param2)

axis_adj_df <- as.data.frame(t(axis_scores_adj))
axis_adj_df$sample_id <- rownames(axis_adj_df)

df_adj <- merge(meta, axis_adj_df, by="sample_id")

pdf("axis_boxplots_adjusted.pdf", width=12, height=4)
par(mfrow = c(1,3))
boxplot(AXIS_METABOLIC ~ group, data = df_adj, main="Metabolic axis (z)", las=2)
boxplot(AXIS_CELL_CYCLE ~ group, data = df_adj, main="Cell-cycle axis (z)", las=2)
boxplot(AXIS_INFLAMMATION ~ group, data = df_adj, main="Inflammation axis (z)", las=2)
par(mfrow = c(1,1))
dev.off()

plot(df_adj$AXIS_METABOLIC, df_adj$AXIS_CELL_CYCLE,
     col = as.integer(df_adj$group), pch = 16,
     xlab = "Metabolic axis (z)", ylab = "Cell-cycle axis (z)")
legend("topright", legend = levels(axis_merged$group),
       col = seq_along(levels(df_adj$group)), pch = 16, bty = "n")

summary(aov(AXIS_METABOLIC ~ group, data = df_adj))
summary(aov(AXIS_CELL_CYCLE ~ group, data = df_adj))
summary(aov(AXIS_INFLAMMATION ~ group, data = df_adj))

nfl <- read.csv("answerals/NfL/NfL_Seimens_Atellica_Immunoassay_Proteomics_2022_07.csv", stringsAsFactors = FALSE)

# Keep only baseline
nfl_baseline <- subset(nfl, `Sample.Timepoint` == 1)

# Keep only necessary columns
nfl_baseline <- nfl_baseline[, c("SubjID", "Result.Value")]

# Rename cleanly
colnames(nfl_baseline) <- c("SubjectUID", "NfL_pgml")
# Remove commas if present
nfl_baseline$NfL_pgml <- gsub(",", "", nfl_baseline$NfL_pgml)

# Convert to numeric
nfl_baseline$NfL_pgml <- as.numeric(nfl_baseline$NfL_pgml)

# Check
summary(nfl_baseline$NfL_pgml)
# Log-transform (important)
nfl_baseline$logNfL <- log10(nfl_baseline$NfL_pgml)

hist(nfl_baseline$NfL_pgml, breaks=50, main="Raw NfL")
hist(nfl_baseline$logNfL, breaks=50, main="log10 NfL")


df_adj$SubjectUID <- sub("^CASE-", "", df_adj$Participant_ID)
length(intersect(df_adj$SubjectUID, nfl_baseline$SubjectUID))
df_nfl <- merge(df_adj, nfl_baseline, by = "SubjectUID")
dim(df_nfl)
table(df_nfl$group)
summary(df_nfl$logNfL)

dup_ids <- df_nfl$SubjectUID[duplicated(df_nfl$SubjectUID)]
df_nfl[df_nfl$SubjectUID %in% dup_ids,
       c("SubjectUID","Participant_ID","sample_id","group","NfL_pgml","logNfL")]

df_nfl_uniq <- df_nfl[!duplicated(df_nfl$SubjectUID), ]

dim(df_nfl_uniq)
table(df_nfl_uniq$group)

als_nfl <- subset(df_nfl_uniq, group != "Control")

m_all <- lm(
  logNfL ~ AXIS_CELL_CYCLE + AXIS_INFLAMMATION + AXIS_METABOLIC + SEX + group,
  data = als_nfl
)

summary(m_all)

m_all_age <- lm(
  logNfL ~ AXIS_CELL_CYCLE + AXIS_INFLAMMATION + AXIS_METABOLIC +
    AGE_AT_SYMPTOM_ONSET + SEX + group,
  data = als_nfl
)

summary(m_all_age)

plot(als_nfl$AXIS_INFLAMMATION,
     als_nfl$logNfL,
     pch=16,
     col=as.integer(als_nfl$group),
     xlab="Inflammation axis (SV-adjusted)",
     ylab="log10 NfL")

abline(lm(logNfL ~ AXIS_INFLAMMATION, data=als_nfl), col="red", lwd=2)

m_slope <- lm(
  ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_CELL_CYCLE +
    AXIS_INFLAMMATION +
    AXIS_METABOLIC +
    AGE_AT_SYMPTOM_ONSET +
    SEX + group,
  data = df_nfl_uniq
)

summary(m_slope)

m_int <- lm(
  logNfL ~ AXIS_INFLAMMATION * group +
    AGE_AT_SYMPTOM_ONSET + SEX,
  data = als_nfl
)

summary(m_int)

summary(df_nfl_uniq$ALSFRS_R_PROGRESSION_SLOPE)

df_plot <- subset(df_nfl_uniq,
                  complete.cases(AXIS_METABOLIC,
                                 ALSFRS_R_PROGRESSION_SLOPE,
                                 group))

grp <- as.character(df_plot$group)

pdf("metabolic_axis_vs_progression.pdf",
    width = 10,
    height = 6)

# add extra right margin
par(mar = c(5,5,2,9))
cols <- c(
  "group1" = "steelblue",
  "group2" = "firebrick",
  "group3" = "darkgreen"
)
plot(df_plot$AXIS_METABOLIC,
     df_plot$ALSFRS_R_PROGRESSION_SLOPE,
     pch = 19,
     cex = 1.2,
     col = cols[grp],
     xlab = "Metabolic axis (SV-adjusted)",
     ylab = "ALSFRS progression slope")

abline(lm(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC,
          data = df_plot),
       col = "black",
       lwd = 2)

# allow drawing in margin
par(xpd = NA)

legend("topright",
       inset = c(-0.25, 0),
       legend = c("C9orf72","SOD1","Sporadic"),
       col = cols[c("C9orf72","SOD1","Sporadic")],
       pch = 19,
       pt.cex = 1.2,
       cex = 1,
       title = "Mutation subtype",
       bty = "n")

dev.off()


df_plot <- subset(df_nfl_uniq,
                  complete.cases(AXIS_METABOLIC,
                                 ALSFRS_R_PROGRESSION_SLOPE,
                                 group))

# Quartiles of metabolic axis
df_plot$AXIS_METABOLIC_Q <- cut(df_plot$AXIS_METABOLIC,
                                breaks = quantile(df_plot$AXIS_METABOLIC,
                                                  probs = seq(0, 1, 0.25),
                                                  na.rm = TRUE),
                                include.lowest = TRUE,
                                labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"))

pdf("metabolic_quartiles_vs_progression.pdf", width = 8, height = 5)

par(mar = c(5,5,2,2))

qtab <- table(df_plot$AXIS_METABOLIC_Q)
xlabs <- paste0(names(qtab), "\n(n=", as.integer(qtab), ")")

boxplot(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC_Q,
        data = df_plot,
        outline = FALSE,
        xlab = "Metabolic axis quartile (SV-adjusted)",
        ylab = "ALSFRS progression slope",
        names = xlabs)

# Overlay jittered points
stripchart(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC_Q,
           data = df_plot,
           vertical = TRUE,
           method = "jitter",
           pch = 19,
           cex = 0.6,
           add = TRUE)

# Add trend of quartile means
q_means <- tapply(df_plot$ALSFRS_R_PROGRESSION_SLOPE, df_plot$AXIS_METABOLIC_Q, mean, na.rm = TRUE)
lines(x = 1:4, y = q_means, lwd = 2)

dev.off()


df_plot <- subset(df_nfl_uniq,
                  complete.cases(AXIS_METABOLIC,
                                 ALSFRS_R_PROGRESSION_SLOPE,
                                 AGE_AT_SYMPTOM_ONSET,
                                 SEX,
                                 group))

fit <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC + AGE_AT_SYMPTOM_ONSET + SEX + group,
          data = df_plot)

pr <- residuals(fit) + coef(fit)["AXIS_METABOLIC"] * df_plot$AXIS_METABOLIC

beta <- unname(coef(fit)["AXIS_METABOLIC"])
pval <- summary(fit)$coefficients["AXIS_METABOLIC", "Pr(>|t|)"]
n <- nobs(fit)

pdf("metabolic_partial_effect_vs_progression_annotated.pdf", width = 8, height = 5)

par(mar = c(5, 5, 2.5, 2))

plot(df_plot$AXIS_METABOLIC,
     pr,
     pch = 19,
     cex = 0.9,
     xlab = "Metabolic axis (SV-adjusted)",
     ylab = "ALSFRS slope (partial residual; adjusted)")

abline(a = 0, b = beta, lwd = 2)

mtext(sprintf("Adjusted model: beta = %.3f, p = %.3g, n = %d", beta, pval, n),
      side = 3, line = 0.5, adj = 0, cex = 0.9)

dev.off()


fit <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC + AGE_AT_SYMPTOM_ONSET + SEX + group,
          data = df_plot)

summary(fit)$coefficients["AXIS_METABOLIC", ]
confint(fit, "AXIS_METABOLIC")
nobs(fit)


fit_full <- lm(ALSFRS_R_PROGRESSION_SLOPE ~
                 AXIS_METABOLIC + AXIS_INFLAMMATION + AXIS_CELL_CYCLE +
                 AGE_AT_SYMPTOM_ONSET + SEX + group,
               data = df_plot)
summary(fit_full)$coefficients[c("AXIS_METABOLIC","AXIS_INFLAMMATION","AXIS_CELL_CYCLE"), ]


df_std <- df_plot
df_std$AXIS_METABOLIC_z <- scale(df_std$AXIS_METABOLIC)

fit_z <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC_z + AGE_AT_SYMPTOM_ONSET + SEX + group,
            data = df_std)
summary(fit_z)$coefficients["AXIS_METABOLIC_z", ]


axes <- df_plot[, c("AXIS_METABOLIC","AXIS_INFLAMMATION","AXIS_CELL_CYCLE")]
cor(axes, use="complete.obs")



make_partial_plot <- function(axis_var, outfile){
  dfp <- subset(df_plot, complete.cases(df_plot[, c(axis_var,
                                                    "ALSFRS_R_PROGRESSION_SLOPE",
                                                    "AXIS_METABOLIC",
                                                    "AXIS_INFLAMMATION",
                                                    "AXIS_CELL_CYCLE",
                                                    "AGE_AT_SYMPTOM_ONSET",
                                                    "SEX",
                                                    "group")]))
  
  fml <- as.formula(paste(
    "ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC + AXIS_INFLAMMATION + AXIS_CELL_CYCLE + AGE_AT_SYMPTOM_ONSET + SEX + group"
  ))
  fit <- lm(fml, data = dfp)
  
  beta <- coef(fit)[axis_var]
  pr <- residuals(fit) + beta * dfp[[axis_var]]
  pval <- summary(fit)$coefficients[axis_var, "Pr(>|t|)"]
  n <- nobs(fit)
  
  pdf(outfile, width=8, height=5)
  par(mar=c(5,5,2.5,2))
  plot(dfp[[axis_var]], pr,
       pch=19, cex=0.9,
       xlab=paste0(axis_var, " (SV-adjusted)"),
       ylab="ALSFRS slope (partial residual; adjusted)")
  abline(a=0, b=beta, lwd=2)
  mtext(sprintf("Adjusted model: beta = %.3f, p = %.3g, n = %d", beta, pval, n),
        side=3, line=0.5, adj=0, cex=0.9)
  dev.off()
}

make_partial_plot("AXIS_INFLAMMATION", "inflammation_partial_vs_progression.pdf")
make_partial_plot("AXIS_CELL_CYCLE",   "cellcycle_partial_vs_progression.pdf")


met_resid <- resid(lm(AXIS_METABOLIC ~ AXIS_INFLAMMATION + AXIS_CELL_CYCLE, data=df_plot))

fit_orth <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ met_resid + AGE_AT_SYMPTOM_ONSET + SEX + group,
               data=df_plot)

summary(fit_orth)$coefficients["met_resid", ]
confint(fit_orth, "met_resid")

pdf("axis_metabolic_sex.pdf", width=8, height=5)
boxplot(AXIS_METABOLIC ~ SEX,
        data=df_plot,
        col=c("steelblue","orange"),
        main="Metabolic Axis by Sex")
dev.off()

met_resid <- resid(lm(AXIS_METABOLIC ~ AXIS_INFLAMMATION + AXIS_CELL_CYCLE, data=df_plot))

# Let's test other variables

candidate_vars <- c(
  "ALSFRS_R_PROGRESSION_SLOPE",
  "COGNITIVE_PROGRESSION_SLOPE",
  "ALSFRS_R_Baseline_Value",
  "ALSFRS_R_Latest_Value",
  "CBS_Baseline_Value",
  "CBS_Latest_Value",
  "AGE_AT_SYMPTOM_ONSET",
  "C9orf72_repeat_length",
  "ATXN2_repeat_length",
  "PCT_SMI32",
  "PCT_ISL1",
  "PCT_NKX61",
  "PCT_TUJ1",
  "PCT_S100b",
  "PCT_Nestin",
  "Bulbar_Onset",
  "Axial_Onset",
  "Limb_Onset",
  "Other_Onset",
  "HasIPSC",
  "HasDIMN",
  "HasVariant",
  "EH_C9orf72",
  "EH_ATXN2",
  "Site_of_Onset",
  "REVISED_EL_ESCORIAL_CRITERIA"
)

qc_table <- data.frame(
  variable = candidate_vars,
  present = candidate_vars %in% colnames(df_nfl_uniq),
  non_missing = sapply(candidate_vars, function(v) {
    if (!v %in% colnames(df_nfl_uniq)) return(NA)
    sum(!is.na(df_nfl_uniq[[v]]) & df_nfl_uniq[[v]] != "")
  }),
  n_unique = sapply(candidate_vars, function(v) {
    if (!v %in% colnames(df_nfl_uniq)) return(NA)
    length(unique(na.omit(df_nfl_uniq[[v]])))
  })
)

qc_table[order(-qc_table$non_missing), ]



suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})

# =========================================================
# 0) SETUP
# =========================================================

axis_vars <- c("AXIS_METABOLIC", "AXIS_INFLAMMATION", "AXIS_CELL_CYCLE")
sv_vars   <- grep("^SV[0-9]+$", colnames(df_nfl_uniq), value = TRUE)

# Restrict to ALS only for clinical analyses
df_als_only <- subset(df_nfl_uniq, group != "Control")

# Make sure key variables are the right type
df_als_only$group <- as.factor(df_als_only$group)
df_als_only$SEX   <- as.factor(df_als_only$SEX)

# Continuous variables that look usable from your QC table
continuous_vars_priority <- c(
  "ALSFRS_R_PROGRESSION_SLOPE",
  "ALSFRS_R_Baseline_Value",
  "ALSFRS_R_Latest_Value",
  "COGNITIVE_PROGRESSION_SLOPE",
  "CBS_Baseline_Value",
  "CBS_Latest_Value",
  "AGE_AT_SYMPTOM_ONSET",
  "PCT_SMI32",
  "PCT_ISL1",
  "PCT_NKX61",
  "PCT_TUJ1",
  "PCT_S100b",
  "PCT_Nestin"
)

continuous_vars_priority <- continuous_vars_priority[
  continuous_vars_priority %in% colnames(df_als_only)
]

# =========================================================
# 1) HELPER: QC TABLE FOR MODELABLE CONTINUOUS VARIABLES
# =========================================================

qc_cont <- data.frame(
  variable = continuous_vars_priority,
  non_missing = sapply(continuous_vars_priority, function(v) {
    sum(!is.na(df_als_only[[v]]))
  }),
  n_unique = sapply(continuous_vars_priority, function(v) {
    length(unique(na.omit(df_als_only[[v]])))
  }),
  stringsAsFactors = FALSE
)

qc_cont <- qc_cont %>% arrange(desc(non_missing))
print(qc_cont)

# Keep only variables with enough data and variation
usable_continuous <- qc_cont %>%
  filter(non_missing >= 25, n_unique >= 5) %>%
  pull(variable)

cat("\nUsable continuous outcomes:\n")
print(usable_continuous)

# =========================================================
# 2) SCREENING FUNCTION FOR CONTINUOUS OUTCOMES
# =========================================================

run_continuous_screen <- function(data, outcomes, axis_vars,
                                  adjust_vars = c("SEX", "AGE_AT_SYMPTOM_ONSET", "group")) {
  out <- list()
  
  for (y in outcomes) {
    if (!y %in% colnames(data)) next
    
    # do not adjust for the outcome itself
    covars_y <- setdiff(adjust_vars, y)
    needed <- unique(c(y, axis_vars, covars_y))
    needed <- needed[needed %in% colnames(data)]
    
    df_sub <- data[, needed, drop = FALSE]
    df_sub <- df_sub[complete.cases(df_sub), , drop = FALSE]
    
    if (nrow(df_sub) < 20) next
    if (!is.numeric(df_sub[[y]])) next
    if (length(unique(df_sub[[y]])) < 5) next
    
    # -------------------------
    # Joint model: all axes together
    # -------------------------
    joint_terms <- c(axis_vars, covars_y)
    joint_terms <- joint_terms[joint_terms %in% colnames(df_sub)]
    
    fml_joint <- as.formula(
      paste(y, "~", paste(joint_terms, collapse = " + "))
    )
    
    fit_joint <- lm(fml_joint, data = df_sub)
    sm_joint  <- summary(fit_joint)
    coef_joint <- sm_joint$coefficients
    
    for (ax in axis_vars) {
      if (ax %in% rownames(coef_joint)) {
        out[[paste(y, ax, "joint", sep = "__")]] <- data.frame(
          outcome   = y,
          predictor = ax,
          model     = "joint",
          beta      = coef_joint[ax, "Estimate"],
          se        = coef_joint[ax, "Std. Error"],
          stat      = coef_joint[ax, "t value"],
          p         = coef_joint[ax, "Pr(>|t|)"],
          conf_low  = confint(fit_joint, ax)[1],
          conf_high = confint(fit_joint, ax)[2],
          adj_r2    = sm_joint$adj.r.squared,
          n         = nobs(fit_joint),
          stringsAsFactors = FALSE
        )
      }
    }
    
    # -------------------------
    # Single-axis models
    # -------------------------
    for (ax in axis_vars) {
      single_terms <- c(ax, covars_y)
      single_terms <- single_terms[single_terms %in% colnames(df_sub)]
      
      fml_single <- as.formula(
        paste(y, "~", paste(single_terms, collapse = " + "))
      )
      
      fit_single <- lm(fml_single, data = df_sub)
      sm_single  <- summary(fit_single)
      coef_single <- sm_single$coefficients
      
      if (ax %in% rownames(coef_single)) {
        out[[paste(y, ax, "single", sep = "__")]] <- data.frame(
          outcome   = y,
          predictor = ax,
          model     = "single",
          beta      = coef_single[ax, "Estimate"],
          se        = coef_single[ax, "Std. Error"],
          stat      = coef_single[ax, "t value"],
          p         = coef_single[ax, "Pr(>|t|)"],
          conf_low  = confint(fit_single, ax)[1],
          conf_high = confint(fit_single, ax)[2],
          adj_r2    = sm_single$adj.r.squared,
          n         = nobs(fit_single),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  res <- bind_rows(out)
  
  if (nrow(res) > 0) {
    res$FDR_global <- p.adjust(res$p, method = "BH")
    res <- res %>%
      group_by(model) %>%
      mutate(FDR_within_model = p.adjust(p, method = "BH")) %>%
      ungroup()
  }
  
  res
}

# =========================================================
# 3) RUN MAIN SCREEN
# =========================================================

adjust_vars_main <- unique(c("SEX", "AGE_AT_SYMPTOM_ONSET", "group", sv_vars))
adjust_vars_main <- adjust_vars_main[adjust_vars_main %in% colnames(df_als_only)]

cont_res <- run_continuous_screen(
  data = df_als_only,
  outcomes = usable_continuous,
  axis_vars = axis_vars,
  adjust_vars = adjust_vars_main
)

cont_res <- cont_res %>%
  arrange(FDR_global, p)

print(cont_res)

# Save full results
write_csv(cont_res, "axis_covariate_screen_continuous_ALS_only.csv")

# =========================================================
# 4) CLINICAL-ONLY SUBSET
# =========================================================

clinical_vars <- c(
  "ALSFRS_R_PROGRESSION_SLOPE",
  "ALSFRS_R_Baseline_Value",
  "ALSFRS_R_Latest_Value",
  "COGNITIVE_PROGRESSION_SLOPE",
  "CBS_Baseline_Value",
  "CBS_Latest_Value",
  "AGE_AT_SYMPTOM_ONSET"
)

clinical_vars <- intersect(clinical_vars, usable_continuous)

clinical_res <- cont_res %>%
  filter(outcome %in% clinical_vars) %>%
  arrange(FDR_global, p)

print(clinical_res)

write_csv(clinical_res, "axis_covariate_screen_clinical_ALS_only.csv")

# =========================================================
# 5) MARKER-ONLY SUBSET
# =========================================================

marker_vars <- c(
  "PCT_SMI32",
  "PCT_ISL1",
  "PCT_NKX61",
  "PCT_TUJ1",
  "PCT_S100b",
  "PCT_Nestin"
)

marker_vars <- intersect(marker_vars, usable_continuous)

marker_res <- cont_res %>%
  filter(outcome %in% marker_vars) %>%
  arrange(FDR_global, p)

print(marker_res)

write_csv(marker_res, "axis_covariate_screen_marker_ALS_only.csv")

# =========================================================
# 6) EASY SUMMARY TABLE FOR SLIDES / REPORT
# =========================================================

top_hits <- cont_res %>%
  mutate(
    beta_round = round(beta, 3),
    p_fmt = signif(p, 3),
    FDR_fmt = signif(FDR_global, 3),
    ci = paste0("(", round(conf_low, 3), ", ", round(conf_high, 3), ")")
  ) %>%
  dplyr::select(outcome, predictor, model, beta_round, ci, p_fmt, FDR_fmt, adj_r2, n) %>%
  arrange(FDR_fmt, p_fmt)

print(top_hits)

write_csv(top_hits, "axis_covariate_screen_top_hits_table.csv")

# =========================================================
# 7) OPTIONAL: INTERACTION TESTS WITH SUBTYPE
#    Useful if you want to know whether axis effects differ
#    by mutation subtype
# =========================================================

run_interaction_tests <- function(data, outcomes, axis_vars,
                                  adjust_vars = c("SEX", "AGE_AT_SYMPTOM_ONSET")) {
  out <- list()
  
  for (y in outcomes) {
    if (!y %in% colnames(data)) next
    
    covars_y <- setdiff(adjust_vars, y)
    needed <- unique(c(y, axis_vars, "group", covars_y))
    needed <- needed[needed %in% colnames(data)]
    
    df_sub <- data[, needed, drop = FALSE]
    df_sub <- df_sub[complete.cases(df_sub), , drop = FALSE]
    
    if (nrow(df_sub) < 30) next
    if (!is.numeric(df_sub[[y]])) next
    if (length(unique(df_sub[[y]])) < 5) next
    
    for (ax in axis_vars) {
      fml <- as.formula(
        paste(y, "~", paste(c(ax, "group", paste0(ax, ":group"), covars_y), collapse = " + "))
      )
      
      fit <- lm(fml, data = df_sub)
      an <- anova(fit)
      
      int_rows <- grep(paste0("^", ax, ":group$|^group:", ax), rownames(an))
      if (length(int_rows) == 0) {
        int_rows <- grep(":", rownames(an))
      }
      
      if (length(int_rows) > 0) {
        out[[paste(y, ax, sep = "__")]] <- data.frame(
          outcome = y,
          predictor = ax,
          interaction_p = an[int_rows[1], "Pr(>F)"],
          n = nobs(fit),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  res <- bind_rows(out)
  if (nrow(res) > 0) {
    res$FDR <- p.adjust(res$interaction_p, method = "BH")
  }
  res
}

interaction_res <- run_interaction_tests(
  data = df_als_only,
  outcomes = clinical_vars,
  axis_vars = axis_vars,
  adjust_vars = c("SEX", "AGE_AT_SYMPTOM_ONSET", sv_vars)
)

print(interaction_res)
write_csv(interaction_res, "axis_covariate_screen_interactions_ALS_only.csv")

# Axes are heavily influenced by cell differentiation markers

# Correlate each gene with axis score
get_axis_genes <- function(expr, axis_scores, axis_name, top_n = 50) {
  ax <- axis_scores[axis_name, ]
  
  cors <- apply(expr, 1, function(g) {
    cor(g, ax, use = "complete.obs")
  })
  
  df <- data.frame(
    gene = rownames(expr),
    cor = cors
  )
  
  df <- df[order(-abs(df$cor)), ]
  head(df, top_n)
}


met_genes <- get_axis_genes(expr_symbol, axis_scores, "AXIS_METABOLIC", 50)
inf_genes <- get_axis_genes(expr_symbol, axis_scores, "AXIS_INFLAMMATION", 50)
cc_genes  <- get_axis_genes(expr_symbol, axis_scores, "AXIS_CELL_CYCLE", 50)

marker_vars <- c("PCT_S100b","PCT_TUJ1","PCT_Nestin","PCT_NKX61")

remove_marker_driven <- function(gene_list, expr, df, threshold = 0.4) {
  keep <- sapply(gene_list$gene, function(g) {
    cors <- sapply(marker_vars, function(m) {
      if (m %in% colnames(df)) {
        cor(expr[g, ], df[[m]], use="complete.obs")
      } else NA
    })
    max(abs(cors), na.rm = TRUE) < threshold
  })
  
  gene_list[keep, ]
}

common_ids <- intersect(colnames(expr_symbol), df_als_only$sample_id)

expr_als <- expr_symbol[, common_ids, drop = FALSE]

df_als_aligned <- df_als_only[match(common_ids, df_als_only$sample_id), ]
inf_genes_clean <- remove_marker_driven(inf_genes, expr_als, df_als_aligned)
met_genes_clean <- remove_marker_driven(met_genes, expr_als, df_als_aligned)
cc_genes_clean  <- remove_marker_driven(cc_genes,  expr_als, df_als_aligned)

annot_als_filtered$DISEASE_DURATION_DEATH <- 
  annot_als_filtered$AGE_AT_DEATH - annot_als_filtered$AGE_AT_SYMPTOM_ONSET

# Merge axis scores with survival data
surv_axis <- merge(
  annot_als_filtered,
  ax_clean_df[, c("sample_id",
                  "AXIS_METABOLIC",
                  "AXIS_INFLAMMATION",
                  "AXIS_CELL_CYCLE")],
  by.x = "Sample_ID",
  by.y = "sample_id"
)
surv_axis$DISEASE_DURATION_DEATH <- as.numeric(as.character(surv_axis$DISEASE_DURATION_DEATH))
# Continuous correlation — what Eshima actually did
cor.test(surv_axis$AXIS_METABOLIC, surv_axis$DISEASE_DURATION_DEATH, 
         method = "pearson")
cor.test(surv_axis$AXIS_INFLAMMATION, surv_axis$DISEASE_DURATION_DEATH,
         method = "pearson")
cor.test(surv_axis$AXIS_CELL_CYCLE, surv_axis$DISEASE_DURATION_DEATH,
         method = "pearson")

# Also try regression adjusting for age at onset
lm(DISEASE_DURATION_DEATH ~ AXIS_METABOLIC + AXIS_INFLAMMATION + 
     AXIS_CELL_CYCLE + AGE_AT_SYMPTOM_ONSET + SEX, 
   data = surv_axis) |> summary()


library(mediation)

# Build the mediation dataset
# Need: cluster_final, logNfL, disease_duration, AGE_AT_SYMPTOM_ONSET, SEX
surv_med <- surv_axis |>
  filter(!is.na(logNfL), !is.na(disease_duration)) |>
  mutate(is_sals3 = as.integer(cluster_final == "S-ALS-3"))

cat("n for mediation:", nrow(surv_med), "\n")
table(surv_med$cluster_final)

# Step 1: mediator model (S-ALS-3 → NfL)
med_model <- lm(logNfL ~ is_sals3 + AGE_AT_SYMPTOM_ONSET + SEX,
                data = surv_med)
summary(med_model)

# Step 2: outcome model (NfL + S-ALS-3 → duration)
out_model <- lm(disease_duration ~ logNfL + is_sals3 + 
                  AGE_AT_SYMPTOM_ONSET + SEX,
                data = surv_med)
summary(out_model)

# Step 3: mediation decomposition
set.seed(42)
med_result <- mediate(med_model, out_model,
                      treat = "is_sals3", mediator = "logNfL",
                      boot = TRUE, sims = 1000)
summary(med_result)

# What to look for:
# ACME (average causal mediation effect) = indirect effect through NfL
# ADE  (average direct effect)           = direct S-ALS-3 → duration
# Total effect                           = ACME + ADE
# If ACME significant + ADE non-significant = full mediation

