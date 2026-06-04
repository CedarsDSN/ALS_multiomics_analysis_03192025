# ============================================================
# Refactored entrypoint for original script: sporadic_subtyping.R
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
if (isTRUE(ALS_ARGS$help)) { als_print_help("01_sporadic_subtyping"); quit(save = "no", status = 0) }
CFG <- als_load_config(ALS_ARGS$config)
ALSPROJ_BASE_DIR <- als_choose(ALS_ARGS$base_dir, CFG, "project.base_dir", getwd())
ALSPROJ_OUTDIR <- als_choose(ALS_ARGS$outdir, CFG, "outputs.root", file.path(getwd(), "results", "01_sporadic_subtyping"))
ALSPROJ_SCRIPT_ID <- "01_sporadic_subtyping"
als_init(ALSPROJ_SCRIPT_ID, ALSPROJ_BASE_DIR, ALSPROJ_OUTDIR, ALS_ARGS$redirect_outputs)
als_register_output_wrappers()
als_load_objects(ALS_ARGS$objects)
on.exit(als_write_session_info(ALSPROJ_SCRIPT_ID), add = TRUE)
# ---- End bootstrap ----

# ============================================================
# SPORADIC ALS SUBTYPING
# Picks up from: df_adj, axis_scores, df_nfl_uniq, expr_symbol
# ============================================================

suppressPackageStartupMessages({
  library(ConsensusClusterPlus)   # install.packages("ConsensusClusterPlus") or BiocManager
  library(umap)
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(pheatmap)
  library(survival)
  library(survminer)
})

mainDir <- file.path(ALSPROJ_BASE_DIR)
setwd(mainDir)
# ============================================================
# STEP 0: BUILD SPORADIC-ONLY AXIS SCORE MATRIX
# ============================================================
# df_adj comes directly from your pipeline (meta + axis scores)
# axis_scores is the raw GSVA output (axes x samples)

sporadic_ids <- df_adj$sample_id[df_adj$group == "Sporadic"]

# Axis score matrix: samples x axes, sporadic only
ax_sporadic <- df_adj[df_adj$group == "Sporadic",
                      c("sample_id",
                        "AXIS_METABOLIC",
                        "AXIS_INFLAMMATION",
                        "AXIS_CELL_CYCLE")]

rownames(ax_sporadic) <- ax_sporadic$sample_id
ax_sporadic$sample_id <- NULL

# Z-score each axis across sporadic samples only
# (re-center within group for unbiased clustering)
ax_scaled <- as.data.frame(scale(ax_sporadic))

# ============================================================
# STEP 1: CONSENSUS CLUSTERING
# Why: determines stable k (2–5) without overfitting;
# used in TCGA, ALS (Swindell 2021), and MS subtype papers.
# ============================================================

# ConsensusClusterPlus needs samples x features matrix
m_for_cc <- t(as.matrix(ax_scaled))   # features x samples

set.seed(42)
cc_res <- ConsensusClusterPlus(
  m_for_cc,
  maxK    = 5,      # test k = 2 to 5
  reps    = 1000,   # bootstrap iterations
  pItem   = 0.8,    # 80% of samples per run
  pFeature= 1.0,    # use all axes (only 3, so keep all)
  clusterAlg = "hc",
  distance   = "euclidean",
  seed       = 42,
  plot       = "pdf",
  title      = "consensus_clustering_sporadic"
)

# ---- Assess stability ----
# Delta area plot: look for where the curve flattens
icl <- calcICL(cc_res, plot = "pdf", title = "consensus_icl")

# Print PAC (proportion of ambiguous clustering) per k
# Lower PAC = more stable clusters
pac <- sapply(2:5, function(k) {
  m <- cc_res[[k]]$consensusMatrix
  # PAC = fraction of entries between 0.1 and 0.9
  sum(m > 0.1 & m < 0.9) / (nrow(m) * ncol(m))
})

pac_df <- data.frame(k = 2:5, PAC = round(pac, 4))
print(pac_df)
cat("Recommended k (lowest PAC):", pac_df$k[which.min(pac_df$PAC)], "\n")

# ============================================================
# STEP 2: ASSIGN CLUSTER LABELS (use your chosen k)
# ============================================================

# Replace k_chosen with whichever k PAC analysis recommends
# (most likely 2 or 3 based on 3 axes)
k_chosen <- 3

cluster_labels <- cc_res[[k_chosen]]$consensusClass
names(cluster_labels) <- colnames(m_for_cc)

# Add to sporadic metadata
ax_scaled$cluster <- factor(
  cluster_labels[rownames(ax_scaled)],
  levels = 1:k_chosen,
  labels = paste0("S-ALS-", 1:k_chosen)   # S-ALS-1, S-ALS-2, S-ALS-3
)

cat("\nCluster sizes:\n")
print(table(ax_scaled$cluster))

# ============================================================
# STEP 3: CHARACTERIZE EACH CLUSTER BIOLOGICALLY
# ============================================================

# 3a) Mean axis score per cluster
cluster_means <- ax_scaled %>%
  group_by(cluster) %>%
  summarise(
    Metabolic   = mean(AXIS_METABOLIC),
    Inflammation = mean(AXIS_INFLAMMATION),
    Cell_Cycle  = mean(AXIS_CELL_CYCLE),
    n = n()
  )

print(cluster_means)

# 3b) Heatmap of axis scores (publication-quality)
hm_mat <- t(as.matrix(ax_scaled[, c("AXIS_METABOLIC",
                                    "AXIS_INFLAMMATION",
                                    "AXIS_CELL_CYCLE")]))

ann_col <- data.frame(
  Subtype = ax_scaled$cluster,
  row.names = rownames(ax_scaled)
)

subtype_colors <- setNames(
  c("#E64B35", "#4DBBD5", "#00A087"),
  paste0("S-ALS-", 1:k_chosen)
)

pdf("sporadic_subtype_heatmap.pdf", width = 10, height = 3)
pheatmap(
  hm_mat,
  annotation_col   = ann_col,
  annotation_colors = list(Subtype = subtype_colors),
  cluster_rows     = FALSE,
  cluster_cols     = TRUE,
  show_colnames    = FALSE,
  color = colorRampPalette(c("#2166AC","white","#D6604D"))(100),
  main = "Axis scores by sporadic ALS subtype"
)
dev.off()

# 3c) ANOVA: are axes significantly different across clusters?
ax_scaled_long <- ax_scaled %>%
  tibble::rownames_to_column("sample_id") %>%
  tidyr::pivot_longer(
    cols = c(AXIS_METABOLIC, AXIS_INFLAMMATION, AXIS_CELL_CYCLE),
    names_to = "axis", values_to = "score"
  )

anova_results <- ax_scaled_long %>%
  group_by(axis) %>%
  summarise(
    F_stat = summary(aov(score ~ cluster))[[1]]$`F value`[1],
    p_val  = summary(aov(score ~ cluster))[[1]]$`Pr(>F)`[1],
    eta_sq = {
      a <- aov(score ~ cluster)
      ss <- summary(a)[[1]]$`Sum Sq`
      ss[1] / sum(ss)    # eta-squared effect size
    }
  ) %>%
  mutate(FDR = p.adjust(p_val, method = "BH"))

print(anova_results)

# ============================================================
# STEP 4: UMAP VISUALIZATION
# ============================================================
# Why UMAP over PCA: 3 axes → PCA is trivial; UMAP finds
# non-linear structure and gives a cleaner figure for a paper.
# For 3 axes specifically, PCA and UMAP will be similar, but
# UMAP is more visually convincing.

set.seed(42)
umap_res <- umap(
  as.matrix(ax_scaled[, c("AXIS_METABOLIC",
                          "AXIS_INFLAMMATION",
                          "AXIS_CELL_CYCLE")]),
  config = umap.defaults
)

umap_df <- data.frame(
  UMAP1  = umap_res$layout[, 1],
  UMAP2  = umap_res$layout[, 2],
  cluster = ax_scaled$cluster,
  sample_id = rownames(ax_scaled)
)

p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = subtype_colors) +
  theme_classic(base_size = 14) +
  labs(
    title  = "Sporadic ALS molecular subtypes",
    color  = "Subtype",
    caption = paste0("n=", nrow(umap_df), " sporadic ALS patients; k=", k_chosen)
  )

ggsave("sporadic_umap.pdf", p_umap, width = 6, height = 5)

# ============================================================
# STEP 5: DIFFERENTIAL GENE EXPRESSION BETWEEN SUBTYPES
# Why: gives biological meaning to each cluster beyond axes.
# If S-ALS-1 is "metabolic-high" we can ask: which genes drive it?
# ============================================================

# Add cluster labels back to full metadata
meta_cluster <- as.data.frame(colData(dds_filtered)) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(
    tibble::rownames_to_column(ax_scaled, "sample_id")[, c("sample_id","cluster")],
    by = "sample_id"
  )

# Sporadic-only DESeq2 object for inter-subtype DE
sporadic_samples <- meta_cluster$sample_id[
  !is.na(meta_cluster$cluster)
]

counts_sporadic <- counts_filtered[, sporadic_samples]
meta_sporadic   <- meta_cluster[meta_cluster$sample_id %in% sporadic_samples, ]
meta_sporadic   <- meta_sporadic[match(sporadic_samples, meta_sporadic$sample_id), ]
rownames(meta_sporadic) <- meta_sporadic$sample_id

# Build design: cluster as the variable of interest
# Still include SEX and clean SVs
sv_terms_sub <- paste0("SV", 1:length(clean_svs), collapse = " + ")
design_sub <- as.formula(paste("~ SEX +", sv_terms_sub, "+ cluster"))

dds_sub <- DESeqDataSetFromMatrix(
  countData = counts_sporadic,
  colData   = meta_sporadic,
  design    = design_sub
)

dds_sub <- DESeq(dds_sub)

# Pairwise contrasts between subtypes
res_sub12 <- results(dds_sub, contrast = c("cluster", "S-ALS-1", "S-ALS-2"))
res_sub13 <- results(dds_sub, contrast = c("cluster", "S-ALS-1", "S-ALS-3"))
res_sub23 <- results(dds_sub, contrast = c("cluster", "S-ALS-2", "S-ALS-3"))

# Significant DEGs per contrast (padj < 0.05, |LFC| > 0.5)
sig_filter <- function(res, label) {
  df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5) %>%
    arrange(padj) %>%
    mutate(contrast = label)
  df
}

degs_12 <- sig_filter(res_sub12, "S-ALS-1 vs S-ALS-2")
degs_13 <- sig_filter(res_sub13, "S-ALS-1 vs S-ALS-3")
degs_23 <- sig_filter(res_sub23, "S-ALS-2 vs S-ALS-3")

cat("DEGs S-ALS-1 vs S-ALS-2:", nrow(degs_12), "\n")
cat("DEGs S-ALS-1 vs S-ALS-3:", nrow(degs_13), "\n")
cat("DEGs S-ALS-2 vs S-ALS-3:", nrow(degs_23), "\n")

write.csv(degs_12, "degs_sporadic_subtype1_vs_2.csv", row.names=FALSE)
write.csv(degs_13, "degs_sporadic_subtype1_vs_3.csv", row.names=FALSE)
write.csv(degs_23, "degs_sporadic_subtype2_vs_3.csv", row.names=FALSE)

# ============================================================
# STEP 6: CLINICAL ASSOCIATION — NfL and ALSFRS-R
# ============================================================

# Merge cluster labels into clinical dataframe
df_clinical_cluster <- df_nfl_uniq %>%
  left_join(
    tibble::rownames_to_column(ax_scaled, "sample_id")[, c("sample_id","cluster")],
    by = "sample_id"
  ) %>%
  filter(!is.na(cluster))

cat("\nClinical data per subtype:\n")
print(table(df_clinical_cluster$cluster))

# --- 6a) NfL by subtype ---
m_nfl_sub <- lm(
  logNfL ~ cluster + AGE_AT_SYMPTOM_ONSET + SEX,
  data = df_clinical_cluster
)
summary(m_nfl_sub)

# Pairwise Tukey post-hoc
aov_nfl <- aov(logNfL ~ cluster + AGE_AT_SYMPTOM_ONSET + SEX,
               data = df_clinical_cluster)
tukey_nfl <- TukeyHSD(aov_nfl, which = "cluster")
print(tukey_nfl)

# Plot NfL by subtype
pdf("nfl_by_sporadic_subtype.pdf", width = 6, height = 5)
par(mar = c(5, 5, 2, 2))
boxplot(logNfL ~ cluster,
        data    = df_clinical_cluster,
        col     = subtype_colors[levels(df_clinical_cluster$cluster)],
        las     = 2,
        ylab    = "log10 NfL (pg/mL)",
        xlab    = "",
        main    = "NfL by sporadic ALS subtype",
        outline = FALSE)
stripchart(logNfL ~ cluster,
           data     = df_clinical_cluster,
           vertical = TRUE, method = "jitter",
           pch = 19, cex = 0.6,
           col = "black", add = TRUE)
dev.off()

# --- 6b) ALSFRS-R progression by subtype ---
df_slope <- df_clinical_cluster %>%
  filter(!is.na(ALSFRS_R_PROGRESSION_SLOPE))

m_slope_sub <- lm(
  ALSFRS_R_PROGRESSION_SLOPE ~ cluster + AGE_AT_SYMPTOM_ONSET + SEX,
  data = df_slope
)
summary(m_slope_sub)

aov_slope <- aov(ALSFRS_R_PROGRESSION_SLOPE ~ cluster + AGE_AT_SYMPTOM_ONSET + SEX,
                 data = df_slope)
tukey_slope <- TukeyHSD(aov_slope, which = "cluster")
print(tukey_slope)

pdf("alsfrs_progression_by_subtype.pdf", width = 6, height = 5)
par(mar = c(5, 5, 2, 2))
boxplot(ALSFRS_R_PROGRESSION_SLOPE ~ cluster,
        data    = df_slope,
        col     = subtype_colors[levels(df_slope$cluster)],
        las     = 2,
        ylab    = "ALSFRS-R progression slope",
        xlab    = "",
        main    = "Disease progression by sporadic ALS subtype",
        outline = FALSE)
stripchart(ALSFRS_R_PROGRESSION_SLOPE ~ cluster,
           data     = df_slope,
           vertical = TRUE, method = "jitter",
           pch = 19, cex = 0.6,
           col = "black", add = TRUE)
dev.off()

# --- 6c) Site of onset by subtype (chi-square) ---
if ("Site_of_Onset" %in% colnames(df_clinical_cluster)) {
  onset_tab <- table(df_clinical_cluster$cluster,
                     df_clinical_cluster$Site_of_Onset)
  print(onset_tab)
  print(chisq.test(onset_tab))
}

# --- 6d) Age at symptom onset by subtype ---
m_age <- lm(AGE_AT_SYMPTOM_ONSET ~ cluster + SEX,
            data = df_clinical_cluster)
summary(m_age)

# ============================================================
# STEP 7: ROBUSTNESS CHECKS
# ============================================================

# 7a) Silhouette width (how well-separated are clusters?)
library(cluster)
dist_mat <- dist(as.matrix(ax_scaled[, c("AXIS_METABOLIC",
                                         "AXIS_INFLAMMATION",
                                         "AXIS_CELL_CYCLE")]),
                 method = "euclidean")

sil_width <- silhouette(
  as.integer(ax_scaled$cluster),
  dist_mat
)

cat("\nMean silhouette width:", mean(sil_width[, "sil_width"]), "\n")
# > 0.5 is good; > 0.25 is acceptable for biological data

# 7b) Stability: re-run with k-1 and k+1, check if core of results holds
# This is qualitative — check that dominant cluster biology doesn't flip

# 7c) Sex-stratified sensitivity
for (s in unique(df_adj$SEX)) {
  cat("\n--- Sex:", s, "---\n")
  sub <- ax_scaled[df_adj$SEX[match(rownames(ax_scaled), df_adj$sample_id)] == s, ]
  if (nrow(sub) >= 10) print(table(sub$cluster))
}

# ============================================================
# STEP 8: LABEL CLUSTERS BIOLOGICALLY
# ============================================================
# After running the code, you will interpret cluster_means to name them.
# Example if k=3:
# S-ALS-1 high metabolic, low inflammation → "Metabolic-predominant"
# S-ALS-2 high inflammation, low cell-cycle → "Inflammatory-predominant"
# S-ALS-3 high cell-cycle, moderate metabolic → "Cell-cycle-predominant"
#
# Rename for figures:
bio_labels <- c(
  "S-ALS-1" = "Metabolic",       # REPLACE after inspecting cluster_means
  "S-ALS-2" = "Inflammatory",
  "S-ALS-3" = "Cell-cycle"
)

ax_scaled$bio_subtype <- bio_labels[as.character(ax_scaled$cluster)]
cat("\nBiological subtype assignments:\n")
print(table(ax_scaled$bio_subtype))



print(cluster_means)
# Also run this to confirm axis dominance per cluster
apply(ax_scaled[, c("AXIS_METABOLIC","AXIS_INFLAMMATION","AXIS_CELL_CYCLE")], 
      2, function(x) tapply(x, ax_scaled$cluster, mean))



suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

# Function to run ORA on a DEG list
run_ora <- function(deg_df, label, universe_genes) {
  
  # Strip ENSEMBL version suffixes (e.g. ENSG00000123456.7 -> ENSG00000123456)
  deg_df$gene_clean    <- sub("\\..*$", "", deg_df$gene)
  universe_clean       <- sub("\\..*$", "", universe_genes)
  
  up_genes   <- deg_df$gene_clean[deg_df$log2FoldChange > 0]
  down_genes <- deg_df$gene_clean[deg_df$log2FoldChange < 0]
  
  # Hallmark gene sets (uses Entrez IDs internally)
  msig <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene) %>%
    mutate(entrez_gene = as.character(entrez_gene))
  
  # Map universe ENSEMBL -> Entrez once (reused for both directions)
  universe_entrez <- bitr(
    universe_clean,
    fromType = "ENSEMBL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )$ENTREZID
  
  results <- list()
  
  for (direction in c("up", "down")) {
    genes <- if (direction == "up") up_genes else down_genes
    
    entrez <- bitr(
      genes,
      fromType = "ENSEMBL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )$ENTREZID
    
    if (length(entrez) < 5) {
      cat("Skipping", direction, "for", label, "— too few mapped genes\n")
      next
    }
    
    ora <- enricher(
      gene          = as.character(entrez),
      TERM2GENE     = msig,
      universe      = as.character(universe_entrez),
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    )
    
    if (!is.null(ora) && nrow(ora@result) > 0) {
      results[[direction]] <- ora@result %>%
        filter(p.adjust < 0.05) %>%
        mutate(contrast = label, direction = direction)
    }
  }
  
  bind_rows(results)
}

# Universe = all ENSEMBL IDs in your count matrix
universe <- rownames(counts_filtered)   # ENSEMBL IDs with version suffixes — cleaned inside function

ora_12 <- run_ora(degs_12, "S-ALS-1 vs S-ALS-2", universe)
ora_13 <- run_ora(degs_13, "S-ALS-1 vs S-ALS-3", universe)
ora_23 <- run_ora(degs_23, "S-ALS-2 vs S-ALS-3", universe)

all_ora <- bind_rows(ora_12, ora_13, ora_23)
write.csv(all_ora, "sporadic_subtype_pathway_enrichment.csv", row.names = FALSE)


# Dotplot per contrast — one figure per pairwise comparison
for (contrast_label in unique(all_ora$contrast)) {
  sub <- all_ora %>%
    filter(contrast == contrast_label) %>%
    mutate(
      neg_log10_padj = -log10(p.adjust),
      Description = gsub("HALLMARK_", "", Description)
    ) %>%
    arrange(desc(neg_log10_padj)) %>%
    head(20)
  
  p <- ggplot(sub, aes(x = reorder(Description, neg_log10_padj),
                       y = neg_log10_padj,
                       fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("up"="#D6604D","down"="#4393C3")) +
    theme_classic(base_size = 12) +
    labs(title = contrast_label,
         x = "", y = "-log10 adjusted p-value",
         fill = "Direction")
  
  fname <- paste0("ora_", gsub(" vs ", "_vs_", gsub("-","",contrast_label)), ".pdf")
  ggsave(fname, p, width=9, height=6)
}


# ============================================================
# STEP 1: CHECK HOW MUCH THE SUBTYPES CORRELATE WITH
#         DIFFERENTIATION MARKERS
# ============================================================

marker_vars <- c("PCT_SMI32", "PCT_ISL1", "PCT_NKX61",
                 "PCT_TUJ1", "PCT_S100b", "PCT_Nestin")

# Merge markers into axis_scaled df
marker_data <- df_nfl_uniq[, c("sample_id", marker_vars[marker_vars %in% colnames(df_nfl_uniq)])]
ax_with_markers <- merge(
  tibble::rownames_to_column(ax_scaled, "sample_id"),
  marker_data, by = "sample_id", all.x = TRUE
)

# Correlate each axis with each marker
avail_markers <- marker_vars[marker_vars %in% colnames(ax_with_markers)]
cat("Axis-marker correlations:\n")
for (ax in c("AXIS_METABOLIC","AXIS_INFLAMMATION","AXIS_CELL_CYCLE")) {
  for (m in avail_markers) {
    r <- cor(ax_with_markers[[ax]], ax_with_markers[[m]],
             use = "complete.obs")
    if (abs(r) > 0.2) cat(sprintf("  %s ~ %s: r = %.3f\n", ax, m, r))
  }
}

# ============================================================
# STEP 2: RESIDUALIZE EXPRESSION AGAINST DIFFERENTIATION MARKERS
# ============================================================

# Build marker covariate matrix (only available markers)
# Use sporadic samples only
sporadic_meta <- df_adj[df_adj$group == "Sporadic", ]
sporadic_meta <- merge(sporadic_meta, marker_data, by="sample_id", all.x=TRUE)
rownames(sporadic_meta) <- sporadic_meta$sample_id

# Sporadic expression matrix
expr_sporadic <- expr_symbol[, sporadic_meta$sample_id]

# Which markers have enough non-missing data?
usable_markers <- avail_markers[
  sapply(avail_markers, function(m) {
    sum(!is.na(sporadic_meta[[m]])) > 0.8 * nrow(sporadic_meta)
  })
]
cat("Using markers for residualization:", paste(usable_markers, collapse=", "), "\n")

# Build covariate matrix: existing SVs + differentiation markers
sv_cols <- grep("^SV[0-9]+$", colnames(sporadic_meta), value=TRUE)
all_covars <- c("SEX", sv_cols, usable_markers)
all_covars <- all_covars[all_covars %in% colnames(sporadic_meta)]

# Impute missing marker values with column mean before building model matrix
sporadic_meta_imp <- sporadic_meta
for (m in usable_markers) {
  sporadic_meta_imp[[m]][is.na(sporadic_meta_imp[[m]])] <-
    mean(sporadic_meta_imp[[m]], na.rm=TRUE)
}

covar_mat <- model.matrix(
  as.formula(paste("~", paste(all_covars, collapse=" + "))),
  data = sporadic_meta_imp
)

# Residualize each gene
cat("Residualizing expression matrix against differentiation markers...\n")
expr_clean <- t(apply(expr_sporadic, 1, function(y) {
  fit <- lm(y ~ covar_mat - 1)
  resid(fit)
}))
rownames(expr_clean) <- rownames(expr_sporadic)
colnames(expr_clean) <- colnames(expr_sporadic)

# ============================================================
# STEP 3: RE-COMPUTE AXIS SCORES ON CLEAN EXPRESSION
# ============================================================
library(msigdbr)
library(dplyr)

hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol)

axis_pathways <- axis_sets

axis_gene_sets <- lapply(axis_pathways, function(pathways) {
  hallmark %>%
    filter(gs_name %in% pathways) %>%
    pull(gene_symbol) %>%
    unique()
})

expr_clean_mat <- as.matrix(expr_clean)
storage.mode(expr_clean_mat) <- "double"

param_clean <- ssgseaParam(expr_clean_mat, axis_gene_sets, normalize = TRUE)
axis_scores_clean <- gsva(param_clean)

ax_clean_df <- as.data.frame(t(axis_scores_clean))
ax_clean_df <- as.data.frame(scale(ax_clean_df))
ax_clean_df$sample_id <- rownames(ax_clean_df)

# ============================================================
# STEP 4: RE-CLUSTER ON CLEAN AXIS SCORES
# ============================================================

ax_clean_mat <- t(as.matrix(ax_clean_df[,
                                        c("AXIS_METABOLIC","AXIS_INFLAMMATION","AXIS_CELL_CYCLE")]))

set.seed(42)
cc_clean <- ConsensusClusterPlus(
  ax_clean_mat,
  maxK = 5, reps = 1000,
  pItem = 0.8, pFeature = 1.0,
  clusterAlg = "hc", distance = "euclidean",
  seed = 42, plot = "pdf",
  title = "consensus_clustering_sporadic_clean"
)

# PAC for clean clustering
pac_clean <- sapply(2:5, function(k) {
  m <- cc_clean[[k]]$consensusMatrix
  sum(m > 0.1 & m < 0.9) / (nrow(m) * ncol(m))
})
cat("PAC (clean):\n")
print(data.frame(k=2:5, PAC=round(pac_clean,4)))

# Assign labels with same k
k_chosen <- 3   # or whatever PAC recommends
cluster_clean <- cc_clean[[k_chosen]]$consensusClass
names(cluster_clean) <- colnames(ax_clean_mat)

ax_clean_df$cluster_clean <- factor(
  cluster_clean[ax_clean_df$sample_id],
  levels = 1:k_chosen,
  labels = paste0("S-ALS-", 1:k_chosen)
)

# ============================================================
# STEP 5: MEASURE LABEL STABILITY vs ORIGINAL CLUSTERING
# ============================================================

# Merge old and new cluster labels
stability_df <- data.frame(
  sample_id     = ax_clean_df$sample_id,
  cluster_orig  = ax_scaled$cluster[match(ax_clean_df$sample_id, rownames(ax_scaled))],
  cluster_clean = ax_clean_df$cluster_clean
)

# Agreement table
print(table(stability_df$cluster_orig, stability_df$cluster_clean))

# Adjusted Rand Index: >0.6 means substantially stable
library(aricode)   # install.packages("aricode")
ari <- ARI(
  as.integer(stability_df$cluster_orig),
  as.integer(stability_df$cluster_clean)
)
cat(sprintf("\nAdjusted Rand Index (original vs clean): %.3f\n", ari))
cat("Interpretation: >0.6 = stable; 0.4-0.6 = moderate; <0.4 = unstable\n")

# ============================================================
# STEP 6: RE-RUN DESeq2 AND ORA ON CLEAN SUBTYPES
# ============================================================

# Rebuild DESeq2 with clean cluster labels
sporadic_meta_clean <- sporadic_meta_imp
sporadic_meta_clean$cluster_clean <- ax_clean_df$cluster_clean[
  match(sporadic_meta_clean$sample_id, ax_clean_df$sample_id)
]
sporadic_meta_clean <- sporadic_meta_clean[!is.na(sporadic_meta_clean$cluster_clean), ]
rownames(sporadic_meta_clean) <- sporadic_meta_clean$sample_id

counts_sp_clean <- counts_filtered[, rownames(sporadic_meta_clean)]

sv_terms_clean <- paste0("SV", seq_along(clean_svs), collapse=" + ")
design_clean <- as.formula(paste("~ SEX +", sv_terms_clean,
                                 "+ cluster_clean"))

dds_clean <- DESeqDataSetFromMatrix(
  countData = counts_sp_clean,
  colData   = sporadic_meta_clean,
  design    = design_clean
)
dds_clean <- DESeq(dds_clean)

# Pairwise contrasts
res_clean_12 <- results(dds_clean, contrast=c("cluster_clean","S-ALS-1","S-ALS-2"))
res_clean_13 <- results(dds_clean, contrast=c("cluster_clean","S-ALS-1","S-ALS-3"))
res_clean_23 <- results(dds_clean, contrast=c("cluster_clean","S-ALS-2","S-ALS-3"))

sig_filter <- function(res, label) {
  as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5) %>%
    mutate(contrast = label)
}

degs_clean_12 <- sig_filter(res_clean_12, "S-ALS-1 vs S-ALS-2")
degs_clean_13 <- sig_filter(res_clean_13, "S-ALS-1 vs S-ALS-3")
degs_clean_23 <- sig_filter(res_clean_23, "S-ALS-2 vs S-ALS-3")

cat("Clean DEGs 1v2:", nrow(degs_clean_12), "\n")
cat("Clean DEGs 1v3:", nrow(degs_clean_13), "\n")
cat("Clean DEGs 2v3:", nrow(degs_clean_23), "\n")

# ORA on clean DEGs (reuse run_ora function from before)
universe <- rownames(counts_filtered)
ora_clean_12 <- run_ora(degs_clean_12, "S-ALS-1 vs S-ALS-2 (clean)", universe)
ora_clean_13 <- run_ora(degs_clean_13, "S-ALS-1 vs S-ALS-3 (clean)", universe)
ora_clean_23 <- run_ora(degs_clean_23, "S-ALS-2 vs S-ALS-3 (clean)", universe)

all_ora_clean <- bind_rows(ora_clean_12, ora_clean_13, ora_clean_23)
write.csv(all_ora_clean, "sporadic_subtype_pathway_enrichment_clean.csv",
          row.names=FALSE)


# Check residual S100b correlation with clean inflammation axis
ax_clean_df_merged <- merge(
  ax_clean_df,
  marker_data,
  by = "sample_id", all.x = TRUE
)

cat("Residual correlations after cleaning:\n")
for (ax in c("AXIS_METABOLIC","AXIS_INFLAMMATION","AXIS_CELL_CYCLE")) {
  for (m in avail_markers) {
    r <- cor(ax_clean_df_merged[[ax]],
             ax_clean_df_merged[[m]],
             use = "complete.obs")
    cat(sprintf("  %s ~ %s: r = %.3f\n", ax, m, r))
  }
}
# If AXIS_INFLAMMATION ~ PCT_S100b drops from 0.696 to < 0.2, cleaning worked.
# If it remains > 0.4, you need to add S100b directly into the DESeq2 model.

# ============================================================
# STEP 1: IDENTIFY S100b-CORRELATED GENES IN THE INFLAMMATION
#         GENE SETS
# ============================================================

# Get S100b values matched to sporadic samples
s100b_vec <- sporadic_meta[[  "PCT_S100b.x"  ]]
names(s100b_vec) <- sporadic_meta$sample_id

# Restrict to samples with S100b data
s100b_clean <- s100b_vec[!is.na(s100b_vec)]
shared_samples <- intersect(names(s100b_clean), colnames(expr_sporadic))

expr_for_cor <- expr_sporadic[, shared_samples]
s100b_aligned <- s100b_clean[shared_samples]

# Correlate every gene in inflammation sets with S100b
inflam_genes_all <- unique(c(
  axis_sets[["AXIS_INFLAMMATION"]]
))

# Keep only genes present in expression matrix
inflam_genes_present <- inflam_genes_all[inflam_genes_all %in% rownames(expr_for_cor)]

cat("Inflammation genes in expression matrix:", length(inflam_genes_present), "\n")

# Compute correlation with S100b for each gene
inflam_s100b_cors <- sapply(inflam_genes_present, function(g) {
  cor(expr_for_cor[g, ], s100b_aligned, use = "complete.obs")
})

# Distribution
cat("\nS100b correlation distribution (inflammation genes):\n")
print(quantile(inflam_s100b_cors, probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1),
               na.rm=TRUE))

# How many genes are highly correlated with S100b?
thresholds <- c(0.2, 0.3, 0.4, 0.5)
for (t in thresholds) {
  cat(sprintf("  Genes |r| > %.1f with S100b: %d / %d (%.1f%%)\n",
              t,
              sum(abs(inflam_s100b_cors) > t, na.rm=TRUE),
              length(inflam_s100b_cors),
              100 * mean(abs(inflam_s100b_cors) > t, na.rm=TRUE)))
}

# ============================================================
# STEP 2: REBUILD INFLAMMATION GENE SET EXCLUDING S100b-DRIVEN GENES
# ============================================================

# Use |r| < 0.3 as cutoff — adjust based on distribution above
s100b_threshold <- 0.3

inflam_genes_clean <- names(inflam_s100b_cors)[
  !is.na(inflam_s100b_cors) & abs(inflam_s100b_cors) < s100b_threshold
]

cat(sprintf("\nInflammation genes retained after S100b filtering: %d / %d\n",
            length(inflam_genes_clean),
            length(inflam_genes_present)))

# If fewer than 15 genes remain, raise the threshold slightly
if (length(inflam_genes_clean) < 15) {
  warning("Too few genes remain — consider threshold of 0.4")
}

# Rebuild axis_sets with cleaned inflammation gene set
axis_sets_clean <- axis_sets
axis_sets_clean[["AXIS_INFLAMMATION"]] <- inflam_genes_clean

# ============================================================
# STEP 3: RECOMPUTE AXIS SCORES ON FULL COHORT
#         USING CLEAN GENE SETS
# ============================================================

# Use the full covariate-residualized expression (expr_resid from your pipeline)
# which already corrects for SVs and sex — run on full cohort not just sporadic
expr_resid_mat_full <- as.matrix(expr_resid)
storage.mode(expr_resid_mat_full) <- "double"

param_final <- ssgseaParam(expr_resid_mat_full, axis_sets_clean, normalize=TRUE)
axis_scores_final <- gsva(param_final)

axis_final_df <- as.data.frame(t(axis_scores_final))
axis_final_df <- as.data.frame(scale(axis_final_df))
axis_final_df$sample_id <- rownames(axis_final_df)

# Merge with metadata
meta_final <- as.data.frame(colData(dds_filtered))
meta_final$sample_id <- rownames(meta_final)
df_final <- merge(meta_final, axis_final_df, by="sample_id")

# ============================================================
# STEP 4: VERIFY THAT S100b CORRELATION IS NOW RESOLVED
# ============================================================

sporadic_final <- df_final[df_final$group == "Sporadic", ]
sporadic_final <- merge(sporadic_final, marker_data, by="sample_id", all.x=TRUE)

# Merge markers into sporadic_final first
sporadic_final_check <- merge(
  sporadic_final,
  marker_data,
  by = "sample_id",
  all.x = TRUE
)

# Now check residual correlations
avail_markers <- c("PCT_SMI32","PCT_ISL1","PCT_NKX61",
                   "PCT_TUJ1","PCT_S100b","PCT_Nestin")
avail_markers <- avail_markers[avail_markers %in% colnames(sporadic_final_check)]

cat("Residual correlations after cleaning:\n")
for (ax in c("AXIS_METABOLIC","AXIS_INFLAMMATION","AXIS_CELL_CYCLE")) {
  for (m in avail_markers) {
    r <- cor(sporadic_final_check[[ax]],
             sporadic_final_check[[m]],
             use = "complete.obs")
    if (abs(r) > 0.15) {
      cat(sprintf("  %s ~ %s: r = %.3f  %s\n",
                  ax, m, r,
                  ifelse(abs(r) > 0.3, "<-- still high", "")))
    }
  }
}

# ============================================================
# STEP 5: RECLUSTER SPORADIC SAMPLES ON FINAL CLEAN AXES
# ============================================================

ax_final_sporadic <- sporadic_final[,
                                    c("sample_id","AXIS_METABOLIC","AXIS_INFLAMMATION","AXIS_CELL_CYCLE")]
rownames(ax_final_sporadic) <- ax_final_sporadic$sample_id
ax_final_sporadic$sample_id <- NULL
ax_final_sporadic <- as.data.frame(scale(ax_final_sporadic))

m_final_cc <- t(as.matrix(ax_final_sporadic))

set.seed(42)
cc_final <- ConsensusClusterPlus(
  m_final_cc,
  maxK=5, reps=1000, pItem=0.8, pFeature=1.0,
  clusterAlg="hc", distance="euclidean",
  seed=42, plot="pdf",
  title="consensus_clustering_sporadic_final"
)

pac_final <- sapply(2:5, function(k) {
  m <- cc_final[[k]]$consensusMatrix
  sum(m > 0.1 & m < 0.9) / (nrow(m)*ncol(m))
})
cat("\nFinal PAC values:\n")
print(data.frame(k=2:5, PAC=round(pac_final,4)))

k_final <- 3  # update after checking PAC
cluster_final <- cc_final[[k_final]]$consensusClass
names(cluster_final) <- colnames(m_final_cc)

ax_final_sporadic$cluster <- factor(
  cluster_final[rownames(ax_final_sporadic)],
  levels=1:k_final, labels=paste0("S-ALS-", 1:k_final)
)

# Cluster means — this time the biological names will be real
cluster_means_final <- ax_final_sporadic %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    Metabolic    = round(mean(AXIS_METABOLIC), 3),
    Inflammation = round(mean(AXIS_INFLAMMATION), 3),
    Cell_Cycle   = round(mean(AXIS_CELL_CYCLE), 3)
  )
print(cluster_means_final)

# ARI vs original clustering
ari_final <- ARI(
  as.integer(ax_scaled$cluster[match(rownames(ax_final_sporadic),
                                     rownames(ax_scaled))]),
  as.integer(ax_final_sporadic$cluster)
)
cat(sprintf("\nARI vs original clustering: %.3f\n", ari_final))




# ============================================================
# DECISION: First-round marker-corrected clusters are FINAL
# ARI(original, clean_round1) = 0.791 → stable
# Second-round cleaning over-corrected → do not use
#
# Objects assumed to exist:
#   ax_clean_df       : first-round clean axis scores + cluster_clean
#   degs_clean_12/13/23 : DEGs on clean cluster labels
#   df_nfl_uniq       : metadata + NfL + ALSFRS-R slope
#   counts_filtered   : raw count matrix
#   marker_data       : PCT markers per sample
# ============================================================

# ============================================================
# STEP 0: CONFIRM FIRST-ROUND CLEANING IS SUFFICIENT
# Verify residual S100b correlation is acceptable
# ============================================================

ax_check <- merge(
  ax_clean_df,
  marker_data,
  by = "sample_id", all.x = TRUE
)

cat("=== Residual marker correlations after first-round cleaning ===\n")
for (ax in c("AXIS_METABOLIC", "AXIS_INFLAMMATION", "AXIS_CELL_CYCLE")) {
  for (m in c("PCT_S100b.x", "PCT_TUJ1", "PCT_NKX61", "PCT_Nestin")) {
    if (m %in% colnames(ax_check)) {
      r <- cor(ax_check[[ax]], ax_check[[m]], use = "complete.obs")
      flag <- ifelse(abs(r) > 0.35, " <-- DISCLOSE", "")
      cat(sprintf("  %s ~ %s: r = %.3f%s\n", ax, m, r, flag))
    }
  }
}
# Expected: AXIS_INFLAMMATION ~ PCT_S100b.x drops from 0.696 to ~-0.267 ✓
# If any |r| > 0.35 remains, disclose as limitation — do NOT re-clean

# ============================================================
# STEP 1: LOCK IN FINAL CLUSTER LABELS
# Merge clean cluster labels into clinical dataframe
# ============================================================

cluster_key <- ax_clean_df[, c("sample_id", "cluster_clean")]

df_clinical_clean <- df_nfl_uniq %>%
  left_join(cluster_key, by = "sample_id") %>%
  filter(!is.na(cluster_clean))

cat("\n=== Final clean cluster sizes ===\n")
print(table(df_clinical_clean$cluster_clean))

cat("\nClean cluster means (z-scored axes):\n")
ax_clean_df %>%
  group_by(cluster_clean) %>%
  summarise(
    n            = n(),
    Metabolic    = round(mean(AXIS_METABOLIC), 3),
    Inflammation = round(mean(AXIS_INFLAMMATION), 3),
    Cell_Cycle   = round(mean(AXIS_CELL_CYCLE), 3)
  ) %>%
  print()

# ============================================================
# STEP 2: ORA ON MARKER-CORRECTED DEGs
# Biologically names each subtype
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
})

run_ora <- function(deg_df, label, universe_genes) {
  deg_df$gene_clean  <- sub("\\..*$", "", deg_df$gene)
  universe_clean     <- sub("\\..*$", "", universe_genes)
  
  msig <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene) %>%
    mutate(entrez_gene = as.character(entrez_gene))
  
  universe_entrez <- suppressMessages(
    bitr(universe_clean, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  )$ENTREZID
  
  results <- list()
  for (direction in c("up", "down")) {
    genes <- if (direction == "up") {
      deg_df$gene_clean[deg_df$log2FoldChange > 0]
    } else {
      deg_df$gene_clean[deg_df$log2FoldChange < 0]
    }
    
    entrez <- suppressMessages(
      bitr(genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
    )$ENTREZID
    
    if (length(entrez) < 5) next
    
    ora <- enricher(
      gene          = as.character(entrez),
      TERM2GENE     = msig,
      universe      = as.character(universe_entrez),
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    )
    
    if (!is.null(ora) && nrow(ora@result) > 0) {
      results[[direction]] <- ora@result %>%
        filter(p.adjust < 0.05) %>%
        mutate(contrast = label, direction = direction)
    }
  }
  bind_rows(results)
}

universe <- rownames(counts_filtered)

ora_c12 <- run_ora(degs_clean_12, "S-ALS-1 vs S-ALS-2", universe)
ora_c13 <- run_ora(degs_clean_13, "S-ALS-1 vs S-ALS-3", universe)
ora_c23 <- run_ora(degs_clean_23, "S-ALS-2 vs S-ALS-3", universe)

all_ora_clean <- bind_rows(ora_c12, ora_c13, ora_c23) %>%
  mutate(Description = gsub("HALLMARK_", "", Description))

write.csv(all_ora_clean, "sporadic_subtype_pathway_enrichment_FINAL.csv",
          row.names = FALSE)

cat("\n=== Top pathways per contrast (clean subtypes) ===\n")
all_ora_clean %>%
  group_by(contrast, direction) %>%
  slice_min(p.adjust, n = 5) %>%
  select(contrast, direction, Description, p.adjust, Count) %>%
  print(n = 60)

# Dotplot for each pairwise contrast
subtype_colors <- c("S-ALS-1" = "#E64B35",
                    "S-ALS-2" = "#4DBBD5",
                    "S-ALS-3" = "#00A087")

for (cl in unique(all_ora_clean$contrast)) {
  sub <- all_ora_clean %>%
    filter(contrast == cl) %>%
    mutate(neg_log10_padj = -log10(p.adjust)) %>%
    arrange(desc(neg_log10_padj)) %>%
    head(20)
  
  if (nrow(sub) == 0) next
  
  p <- ggplot(sub, aes(x = reorder(Description, neg_log10_padj),
                       y = neg_log10_padj,
                       fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("up" = "#D6604D", "down" = "#4393C3")) +
    theme_classic(base_size = 12) +
    labs(title  = paste("Clean subtypes:", cl),
         x = "", y = "-log10 adjusted p-value",
         fill = "Direction")
  
  fname <- paste0("ora_clean_",
                  gsub(" vs ", "_v_", gsub("S-ALS-", "S", cl)), ".pdf")
  ggsave(fname, p, width = 9, height = 6)
  cat("Saved:", fname, "\n")
}

# ============================================================
# STEP 3: CLINICAL ASSOCIATIONS — NfL (CLEAN LABELS)
# ============================================================

df_nfl_clean <- df_clinical_clean %>%
  filter(!is.na(logNfL))

cat("\n=== NfL sample counts per clean subtype ===\n")
print(table(df_nfl_clean$cluster_clean))

# Linear model: logNfL ~ subtype + age + sex
m_nfl_clean <- lm(
  logNfL ~ cluster_clean + AGE_AT_SYMPTOM_ONSET + SEX,
  data = df_nfl_clean
)
cat("\n=== NfL model (clean labels) ===\n")
print(summary(m_nfl_clean))

# Tukey post-hoc
aov_nfl_clean <- aov(logNfL ~ cluster_clean + AGE_AT_SYMPTOM_ONSET + SEX,
                     data = df_nfl_clean)
tukey_nfl_clean <- TukeyHSD(aov_nfl_clean, which = "cluster_clean")
cat("\n=== Tukey post-hoc: NfL by clean subtype ===\n")
print(tukey_nfl_clean)

# Box + jitter plot
pdf("nfl_by_subtype_FINAL.pdf", width = 6, height = 5)
par(mar = c(5, 5, 2, 2))
boxplot(logNfL ~ cluster_clean,
        data    = df_nfl_clean,
        col     = subtype_colors[levels(df_nfl_clean$cluster_clean)],
        las     = 2,
        ylab    = "log10 NfL (pg/mL)",
        xlab    = "",
        main    = "NfL by sporadic ALS subtype (marker-corrected)",
        outline = FALSE)
stripchart(logNfL ~ cluster_clean,
           data     = df_nfl_clean,
           vertical = TRUE, method = "jitter",
           pch = 19, cex = 0.5, col = "black", add = TRUE)
dev.off()

# ============================================================
# STEP 4: CLINICAL ASSOCIATIONS — ALSFRS-R (CLEAN LABELS)
# ============================================================

df_slope_clean <- df_clinical_clean %>%
  filter(!is.na(ALSFRS_R_PROGRESSION_SLOPE))

cat("\n=== ALSFRS-R: samples with slope data per subtype ===\n")
print(table(df_slope_clean$cluster_clean))
cat("Total with slope:", nrow(df_slope_clean), "of", nrow(df_clinical_clean), "\n")

m_slope_clean <- lm(
  ALSFRS_R_PROGRESSION_SLOPE ~ cluster_clean + AGE_AT_SYMPTOM_ONSET + SEX,
  data = df_slope_clean
)
cat("\n=== ALSFRS-R progression model (clean labels) ===\n")
print(summary(m_slope_clean))

aov_slope_clean <- aov(ALSFRS_R_PROGRESSION_SLOPE ~ cluster_clean + AGE_AT_SYMPTOM_ONSET + SEX,
                       data = df_slope_clean)
tukey_slope_clean <- TukeyHSD(aov_slope_clean, which = "cluster_clean")
cat("\n=== Tukey post-hoc: ALSFRS-R by clean subtype ===\n")
print(tukey_slope_clean)

# ============================================================
# STEP 5: ALSFRS-R MISSINGNESS SENSITIVITY ANALYSIS
# Is missing slope data random, or subtype-biased?
# ============================================================

df_clinical_clean$has_slope <- !is.na(df_clinical_clean$ALSFRS_R_PROGRESSION_SLOPE)

cat("\n=== Slope data availability by subtype ===\n")
miss_tab <- df_clinical_clean %>%
  group_by(cluster_clean) %>%
  summarise(
    n_total    = n(),
    n_with_slope = sum(has_slope),
    pct_complete = round(100 * mean(has_slope), 1)
  )
print(miss_tab)

# Chi-square: is missingness independent of subtype?
miss_chi <- chisq.test(
  table(df_clinical_clean$cluster_clean, df_clinical_clean$has_slope)
)
cat("\n=== Chi-square: missingness vs subtype ===\n")
print(miss_chi)
# p < 0.05 → missingness is subtype-biased → results must be interpreted carefully
# p ≥ 0.05 → missingness is effectively random → result is reportable with n caveat

# Are samples with missing slope data different on key covariates?
cat("\n=== Covariate comparison: slope-present vs slope-missing ===\n")
miss_covar <- df_clinical_clean %>%
  group_by(has_slope) %>%
  summarise(
    n                    = n(),
    mean_age_onset       = round(mean(AGE_AT_SYMPTOM_ONSET, na.rm=TRUE), 1),
    pct_male             = round(100 * mean(SEX == "Male", na.rm=TRUE), 1),
    mean_logNfL          = round(mean(logNfL, na.rm=TRUE), 3)
  )
print(miss_covar)

# t-test: age at onset
t_age <- t.test(AGE_AT_SYMPTOM_ONSET ~ has_slope, data = df_clinical_clean)
cat(sprintf("\nAge at onset: slope-present=%.1f, slope-missing=%.1f, p=%.3f\n",
            mean(df_clinical_clean$AGE_AT_SYMPTOM_ONSET[df_clinical_clean$has_slope],  na.rm=TRUE),
            mean(df_clinical_clean$AGE_AT_SYMPTOM_ONSET[!df_clinical_clean$has_slope], na.rm=TRUE),
            t_age$p.value))

# NfL difference: is missingness confounded with disease burden?
t_nfl <- t.test(logNfL ~ has_slope, data = df_clinical_clean)
cat(sprintf("log NfL:       slope-present=%.3f, slope-missing=%.3f, p=%.3f\n",
            mean(df_clinical_clean$logNfL[df_clinical_clean$has_slope],  na.rm=TRUE),
            mean(df_clinical_clean$logNfL[!df_clinical_clean$has_slope], na.rm=TRUE),
            t_nfl$p.value))

# ============================================================
# STEP 6: SUMMARISE WHAT IS NOW REPORTABLE
# ============================================================

cat("\n\n========================================\n")
cat("REPORTABLE RESULTS SUMMARY\n")
cat("========================================\n")

cat("\n--- Subtype stability ---\n")
cat("ARI (original vs marker-corrected): 0.791\n")
cat("Label stability: 94% of samples retained assignment\n")
cat("S-ALS-3 most marker-influenced (23 reassigned)\n")

cat("\n--- NfL (use Tukey p-values from Step 3 output) ---\n")
cat("S-ALS-3 carries highest NfL burden\n")
cat("Check: does direction hold with clean labels?\n")

cat("\n--- ALSFRS-R (underpowered — see Step 5 missingness result) ---\n")
cat("Report as trend only if p > 0.05 remains after clean labels\n")
cat("Missingness chi-square p determines whether result is interpretable\n")

cat("\n--- Limitations to disclose ---\n")
cat("AXIS_CELL_CYCLE ~ PCT_S100b: residual correlation after cleaning (disclose)\n")
cat("ALSFRS-R: n=~186 of ~400 sporadic — sensitivity analysis above\n")
cat("SOD1 batch confound suppresses DEG counts\n")



# ============================================================
# PART A: DIAGNOSE EMT PERSISTENCE
# Is EMT driven by a differentiation state not captured
# by the PCT markers? If yes → add EMT score as covariate.
# ============================================================

suppressPackageStartupMessages({
  library(msigdbr)
  library(GSVA)
  library(dplyr)
  library(ggplot2)
})

# Get EMT gene set
emt_set <- msigdbr(species = "Homo sapiens", category = "H") %>%
  filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>%
  pull(gene_symbol)

# Map to your expression matrix rownames (gene symbols assumed)
emt_genes_present <- emt_set[emt_set %in% rownames(expr_symbol)]
cat("EMT genes in expression matrix:", length(emt_genes_present), "\n")

# Compute per-sample EMT score (mean expression of EMT genes)
# Using sporadic samples only
sporadic_ids <- ax_clean_df$sample_id
expr_sp <- expr_symbol[emt_genes_present, sporadic_ids]
emt_score <- colMeans(expr_sp, na.rm = TRUE)

# Merge with marker data
emt_df <- data.frame(
  sample_id = names(emt_score),
  EMT_score = emt_score
)
emt_df <- merge(emt_df, marker_data, by = "sample_id", all.x = TRUE)

# Correlate EMT score with all available PCT markers
cat("\n=== EMT score ~ differentiation marker correlations ===\n")
pct_cols <- grep("^PCT_", colnames(emt_df), value = TRUE)
for (m in pct_cols) {
  r <- cor(emt_df$EMT_score, emt_df[[m]], use = "complete.obs")
  flag <- ifelse(abs(r) > 0.3, " <-- ARTIFACT SIGNAL", "")
  cat(sprintf("  EMT_score ~ %s: r = %.3f%s\n", m, r, flag))
}

# Also: does EMT score differ by cluster?
emt_df2 <- merge(emt_df, ax_clean_df[, c("sample_id", "cluster_clean")],
                 by = "sample_id")
cat("\n=== Mean EMT score by clean cluster ===\n")
print(emt_df2 %>%
        group_by(cluster_clean) %>%
        summarise(mean_EMT = round(mean(EMT_score, na.rm=TRUE), 3), n=n()))

# Boxplot
pdf("emt_score_by_cluster.pdf", width=5, height=4)
boxplot(EMT_score ~ cluster_clean, data=emt_df2,
        col=c("#E64B35","#4DBBD5","#00A087"),
        ylab="Mean EMT gene expression", main="EMT score by clean subtype",
        outline=FALSE)
stripchart(EMT_score ~ cluster_clean, data=emt_df2,
           vertical=TRUE, method="jitter", pch=19, cex=0.4,
           col="black", add=TRUE)
dev.off()

# ============================================================
# DECISION POINT AFTER PART A:
#
# If EMT_score correlates strongly (|r|>0.3) with ANY PCT marker:
#   → The signal is artifact. Add EMT_score as covariate in DESeq2
#     (see Part B below), then rerun ORA.
#
# If EMT_score does NOT correlate with PCT markers:
#   → EMT signal may be real mesenchymal-like biology in S-ALS-1.
#     Report but acknowledge uncertainty (iPSC artifact cannot
#     be fully excluded without additional marker data).
#     Proceed to Part C (naming + clinical writeup).
# ============================================================

# ============================================================
# PART B: IF EMT IS ARTIFACT — ADD EMT SCORE AS COVARIATE
# Rerun DESeq2 with EMT_score as additional covariate
# ============================================================

# Add EMT score to sporadic metadata
sporadic_meta_emt <- merge(
  as.data.frame(colData(dds_filtered))[sporadic_ids, ],
  emt_df[, c("sample_id", "EMT_score")],
  by.x = "row.names", by.y = "sample_id",
  all.x = TRUE
)
rownames(sporadic_meta_emt) <- sporadic_meta_emt$Row.names

# Add cluster labels
sporadic_meta_emt <- merge(
  sporadic_meta_emt,
  ax_clean_df[, c("sample_id", "cluster_clean")],
  by.x = "Row.names", by.y = "sample_id",
  all.x = TRUE
)
rownames(sporadic_meta_emt) <- sporadic_meta_emt$Row.names
sporadic_meta_emt <- sporadic_meta_emt[!is.na(sporadic_meta_emt$cluster_clean), ]

# Impute missing EMT scores (should be minimal)
sporadic_meta_emt$EMT_score[is.na(sporadic_meta_emt$EMT_score)] <-
  mean(sporadic_meta_emt$EMT_score, na.rm = TRUE)

counts_sp_emt <- counts_filtered[, rownames(sporadic_meta_emt)]

sv_terms <- paste0("SV", seq_along(clean_svs), collapse = " + ")
design_emt <- as.formula(
  paste("~ SEX +", sv_terms, "+ EMT_score + cluster_clean")
)

dds_emt <- DESeqDataSetFromMatrix(
  countData = counts_sp_emt,
  colData   = sporadic_meta_emt,
  design    = design_emt
)
dds_emt <- DESeq(dds_emt)

# Pairwise contrasts with EMT covariate
res_emt_12 <- results(dds_emt, contrast = c("cluster_clean","S-ALS-1","S-ALS-2"))
res_emt_13 <- results(dds_emt, contrast = c("cluster_clean","S-ALS-1","S-ALS-3"))
res_emt_23 <- results(dds_emt, contrast = c("cluster_clean","S-ALS-2","S-ALS-3"))

sig_filter <- function(res, label) {
  as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5) %>%
    mutate(contrast = label)
}

degs_emt_12 <- sig_filter(res_emt_12, "S-ALS-1 vs S-ALS-2 (EMT-corrected)")
degs_emt_13 <- sig_filter(res_emt_13, "S-ALS-1 vs S-ALS-3 (EMT-corrected)")
degs_emt_23 <- sig_filter(res_emt_23, "S-ALS-2 vs S-ALS-3 (EMT-corrected)")

cat("\n=== DEG counts after EMT correction ===\n")
cat("1v2:", nrow(degs_emt_12), "\n")
cat("1v3:", nrow(degs_emt_13), "\n")
cat("2v3:", nrow(degs_emt_23), "\n")

# Rerun ORA on EMT-corrected DEGs
ora_emt_12 <- run_ora(degs_emt_12, "S-ALS-1 vs S-ALS-2", universe)
ora_emt_13 <- run_ora(degs_emt_13, "S-ALS-1 vs S-ALS-3", universe)
ora_emt_23 <- run_ora(degs_emt_23, "S-ALS-2 vs S-ALS-3", universe)

all_ora_emt <- bind_rows(ora_emt_12, ora_emt_13, ora_emt_23) %>%
  mutate(Description = gsub("HALLMARK_", "", Description))

cat("\n=== Top pathways after EMT correction ===\n")
all_ora_emt %>%
  group_by(contrast, direction) %>%
  slice_min(p.adjust, n = 5) %>%
  select(contrast, direction, Description, p.adjust, Count) %>%
  print(n = 60)

write.csv(all_ora_emt, "sporadic_subtype_pathway_enrichment_EMTcorrected.csv",
          row.names = FALSE)

# ============================================================
# PART C: BIOLOGICAL NAMING OF CLEAN CLUSTERS
# Based on axis means — confirm after ORA
# ============================================================

cat("\n=== Clean cluster axis profile (for naming) ===\n")
ax_clean_df %>%
  group_by(cluster_clean) %>%
  summarise(
    n            = n(),
    Metabolic    = round(mean(AXIS_METABOLIC), 3),
    Inflammation = round(mean(AXIS_INFLAMMATION), 3),
    Cell_Cycle   = round(mean(AXIS_CELL_CYCLE), 3)
  ) %>%
  mutate(
    dominant_axes = case_when(
      Metabolic > 0.5 & Inflammation > 0.5 & Cell_Cycle > 0.5 ~ "Pan-activated",
      Metabolic < -0.3 & Inflammation < -0.3 & Cell_Cycle < -0.3 ~ "Quiescent",
      Metabolic > 0.5 & Inflammation < 0 ~ "Metabolic-predominant",
      Inflammation > 0.5 & Metabolic < 0 ~ "Inflammatory-predominant",
      TRUE ~ "Mixed"
    )
  ) %>%
  print()

# Current clean cluster interpretation:
# S-ALS-1 (n=96): Met=+0.74, Inf=+1.27, CC=+1.20 → Pan-activated / High-stress
# S-ALS-2 (n=227): Met=-0.63, Inf=-0.41, CC=-0.67 → Quiescent / Low-activity
# S-ALS-3 (n=77): Met=+0.93, Inf=-0.37, CC=+0.48 → Metabolic-predominant

# ============================================================
# PART D: VERIFY ARI IS NOT MASKING LABEL PERMUTATION
# Compare cluster profile similarity, not just label matching
# ============================================================

# Original cluster means (from summary doc):
# S-ALS-1 orig: Met=+0.35, Inf=-1.12, CC=+0.43  → Inflammation-suppressed
# S-ALS-2 orig: Met=+0.93, Inf=+0.51, CC=+0.96  → Pan-activated
# S-ALS-3 orig: Met=-0.69, Inf=+0.01, CC=-0.73  → Quiescent

# Clean cluster means:
# S-ALS-1 clean: Met=+0.74, Inf=+1.27, CC=+1.20 → Pan-activated
# S-ALS-2 clean: Met=-0.63, Inf=-0.41, CC=-0.67 → Quiescent
# S-ALS-3 clean: Met=+0.93, Inf=-0.37, CC=+0.48 → Metabolic

# The labels have PERMUTED. Original S-ALS-2 ≈ Clean S-ALS-1 (Pan-activated)
# Original S-ALS-3 ≈ Clean S-ALS-2 (Quiescent)
# This is expected — ConsensusClusterPlus assigns labels arbitrarily.
# ARI handles this correctly (it's permutation-invariant).
# But you must remap labels for presentation consistency.

# Compute confusion table to find the best label mapping
orig_labels <- ax_scaled$cluster[match(ax_clean_df$sample_id, rownames(ax_scaled))]
clean_labels <- ax_clean_df$cluster_clean
confusion <- table(Original = orig_labels, Clean = clean_labels)
cat("\n=== Confusion table: original vs clean cluster labels ===\n")
print(confusion)
cat("\nUse this to remap clean labels to match original naming convention.\n")

# Find optimal label mapping (max overlap per row)
label_map <- apply(confusion, 1, which.max)
cat("\nSuggested label remapping (Original → Clean):\n")
print(label_map)
# e.g., if Original S-ALS-2 maps to Clean S-ALS-1,
# rename Clean S-ALS-1 → S-ALS-2 for consistency

# Apply remapping
remap <- setNames(names(label_map), paste0("S-ALS-", label_map))
ax_clean_df$cluster_final <- factor(
  remap[as.character(ax_clean_df$cluster_clean)],
  levels = paste0("S-ALS-", 1:3)
)
cat("\nFinal cluster sizes (label-remapped):\n")
print(table(ax_clean_df$cluster_final))

# ============================================================
# PART E: CONSOLIDATED REPORTABLE RESULTS SUMMARY
# ============================================================

cat("\n\n========================================\n")
cat("RESULTS SUMMARY — WHAT IS REPORTABLE\n")
cat("========================================\n")

cat("\n1. SUBTYPE STABILITY\n")
cat("   ARI(original, marker-corrected) = 0.791 → strong\n")
cat("   94% label stability; S-ALS-3 most affected (23 reassigned)\n")
cat("   NOTE: labels permuted in clean solution — remap before figures\n")

cat("\n2. NfL RESULT (reportable)\n")
cat("   S-ALS-3 > S-ALS-2: Tukey p = 0.008\n")
cat("   S-ALS-3 vs S-ALS-1: p = 0.125 (directional, ns)\n")
cat("   Overall model p = 4.55e-08, adj R² = 0.115\n")
cat("   CAVEAT: weaker than original — S-ALS-3 > S-ALS-2 only\n")

cat("\n3. ALSFRS-R (null result — report honestly)\n")
cat("   No subtype differences after clean labels\n")
cat("   Missingness is subtype-random (chi-sq p=0.71) — GOOD\n")
cat("   BUT slope-missing patients have higher NfL (p<0.001)\n")
cat("   → Severity-related dropout: more aggressive cases lost to follow-up\n")
cat("   → Null result may reflect survivor bias, not true absence of effect\n")

cat("\n4. ORA / PATHWAY BIOLOGY\n")
cat("   EMT dominates — must check Part A before interpreting\n")
cat("   Non-EMT signals: Inflammatory response, Coagulation,\n")
cat("   Interferon alpha/gamma (S-ALS-1 vs S-ALS-2/3)\n")
cat("   Cholesterol homeostasis DOWN in S-ALS-2 vs S-ALS-3\n")

cat("\n5. LIMITATIONS TO DISCLOSE\n")
cat("   EMT persistence: possible residual differentiation artifact\n")
cat("   AXIS_CELL_CYCLE ~ PCT_S100b: residual r after cleaning\n")
cat("   ALSFRS-R: 189/302 (63%) complete; severity-related dropout\n")
cat("   SOD1 batch confound: DEG counts suppressed\n")






# ============================================================
# STEP 1: FIX DESeq2 CONVERGENCE — SCALE EMT_score
# 106 non-converged genes were due to unscaled numeric covariate
# ============================================================

sporadic_meta_emt$EMT_score_scaled <- scale(sporadic_meta_emt$EMT_score)[,1]

design_emt_scaled <- as.formula(
  paste("~ SEX +", sv_terms, "+ EMT_score_scaled + cluster_clean")
)

dds_emt2 <- DESeqDataSetFromMatrix(
  countData = counts_sp_emt,
  colData   = sporadic_meta_emt,
  design    = design_emt_scaled
)
dds_emt2 <- DESeq(dds_emt2)
dds_emt2 <- nbinomWaldTest(dds_emt2, maxit = 500)

# Check convergence
n_nonconv <- sum(!mcols(dds_emt2)$betaConv, na.rm = TRUE)
cat("Non-converged genes after scaling:", n_nonconv, "\n")
# Should drop to 0 or near-0

# Rerun contrasts on converged model
res_final_12 <- results(dds_emt2, contrast = c("cluster_clean","S-ALS-1","S-ALS-2"))
res_final_13 <- results(dds_emt2, contrast = c("cluster_clean","S-ALS-1","S-ALS-3"))
res_final_23 <- results(dds_emt2, contrast = c("cluster_clean","S-ALS-2","S-ALS-3"))

sig_filter <- function(res, label) {
  as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5) %>%
    arrange(padj) %>%
    mutate(contrast = label)
}

degs_final_12 <- sig_filter(res_final_12, "S-ALS-1 vs S-ALS-2")
degs_final_13 <- sig_filter(res_final_13, "S-ALS-1 vs S-ALS-3")
degs_final_23 <- sig_filter(res_final_23, "S-ALS-2 vs S-ALS-3")

cat("Final DEGs 1v2:", nrow(degs_final_12), "\n")
cat("Final DEGs 1v3:", nrow(degs_final_13), "\n")
cat("Final DEGs 2v3:", nrow(degs_final_23), "\n")

# Final ORA
universe <- rownames(counts_filtered)
ora_final_12 <- run_ora(degs_final_12, "S-ALS-1 vs S-ALS-2", universe)
ora_final_13 <- run_ora(degs_final_13, "S-ALS-1 vs S-ALS-3", universe)
ora_final_23 <- run_ora(degs_final_23, "S-ALS-2 vs S-ALS-3", universe)

all_ora_final <- bind_rows(ora_final_12, ora_final_13, ora_final_23) %>%
  mutate(Description = gsub("HALLMARK_", "", Description))

cat("\n=== FINAL top pathways (scaled EMT covariate) ===\n")
all_ora_final %>%
  group_by(contrast, direction) %>%
  slice_min(p.adjust, n = 5) %>%
  select(contrast, direction, Description, p.adjust, Count) %>%
  print(n = 60)

write.csv(all_ora_final, "sporadic_subtype_pathway_enrichment_FINAL_scaled.csv",
          row.names = FALSE)

# ============================================================
# STEP 2: CLINICAL ASSOCIATIONS WITH cluster_final LABELS
# ============================================================

# Add cluster_final to clinical dataframe
cluster_final_key <- ax_clean_df[, c("sample_id", "cluster_final")]

df_clinical_final <- df_nfl_uniq %>%
  left_join(cluster_final_key, by = "sample_id") %>%
  filter(!is.na(cluster_final))

subtype_colors <- c("S-ALS-1" = "#E64B35",
                    "S-ALS-2" = "#4DBBD5",
                    "S-ALS-3" = "#00A087")

cat("\n=== Final cluster sizes in clinical data ===\n")
print(table(df_clinical_final$cluster_final))

# --- NfL ---
df_nfl_final <- df_clinical_final %>% filter(!is.na(logNfL))
cat("\n=== NfL counts per final subtype ===\n")
print(table(df_nfl_final$cluster_final))

m_nfl_final <- lm(
  logNfL ~ cluster_final + AGE_AT_SYMPTOM_ONSET + SEX,
  data = df_nfl_final
)
cat("\n=== NfL model (final labels) ===\n")
print(summary(m_nfl_final))

aov_nfl_final <- aov(logNfL ~ cluster_final + AGE_AT_SYMPTOM_ONSET + SEX,
                     data = df_nfl_final)
tukey_nfl_final <- TukeyHSD(aov_nfl_final, which = "cluster_final")
cat("\n=== Tukey NfL (final labels) ===\n")
print(tukey_nfl_final)

# NfL boxplot
pdf("nfl_by_subtype_FINAL_v2.pdf", width = 6, height = 5)
par(mar = c(5, 5, 2, 2))
boxplot(logNfL ~ cluster_final,
        data    = df_nfl_final,
        col     = subtype_colors[levels(df_nfl_final$cluster_final)],
        las     = 2,
        ylab    = "log10 NfL (pg/mL)",
        xlab    = "",
        main    = "NfL by sporadic ALS subtype",
        outline = FALSE)
stripchart(logNfL ~ cluster_final,
           data     = df_nfl_final,
           vertical = TRUE, method = "jitter",
           pch = 19, cex = 0.5, col = "black", add = TRUE)

# Add Tukey significance bracket for the reportable comparison
tukey_p <- tukey_nfl_final$cluster_final
sig_pairs <- which(tukey_p[,"p adj"] < 0.05)
if (length(sig_pairs) > 0) {
  cat("\nSignificant NfL pairwise comparisons:\n")
  print(tukey_p[sig_pairs, , drop=FALSE])
}
dev.off()

# --- ALSFRS-R ---
df_slope_final <- df_clinical_final %>%
  filter(!is.na(ALSFRS_R_PROGRESSION_SLOPE))

m_slope_final <- lm(
  ALSFRS_R_PROGRESSION_SLOPE ~ cluster_final + AGE_AT_SYMPTOM_ONSET + SEX,
  data = df_slope_final
)
cat("\n=== ALSFRS-R model (final labels) ===\n")
print(summary(m_slope_final))

# ============================================================
# STEP 3: BIOLOGICAL NAMING — LOCK IN SUBTYPE LABELS
# ============================================================

cat("\n=== FINAL SUBTYPE BIOLOGICAL NAMES ===\n")
cat("Based on axis profile + EMT-corrected ORA:\n\n")
cat("S-ALS-1 (n=96):  Met=+0.74, Inf=+1.27, CC=+1.20\n")
cat("  → TOP ORA: G2M_CHECKPOINT (p~1e-10), E2F_TARGETS, MITOTIC_SPINDLE\n")
cat("  → BIOLOGICAL NAME: 'Proliferative stress' or 'Cell-cycle activated'\n\n")
cat("S-ALS-2 (n=227): Met=-0.63, Inf=-0.41, CC=-0.67\n")
cat("  → TOP ORA: No strong UP enrichment; low across all axes\n")
cat("  → BIOLOGICAL NAME: 'Quiescent' or 'Low-activity'\n\n")
cat("S-ALS-3 (n=77):  Met=+0.93, Inf=-0.37, CC=+0.48\n")
cat("  → TOP ORA (from 2v3): CHOLESTEROL_HOMEOSTASIS down in S-ALS-2 vs S-ALS-3\n")
cat("     meaning CHOLESTEROL_HOMEOSTASIS enriched in S-ALS-3\n")
cat("  → BIOLOGICAL NAME: 'Metabolic-predominant' or 'Lipid-metabolic'\n\n")

# ============================================================
# STEP 4: SUMMARY FIGURE — AXIS HEATMAP WITH FINAL LABELS
# ============================================================

library(pheatmap)

hm_mat <- t(as.matrix(ax_clean_df[, c("AXIS_METABOLIC",
                                      "AXIS_INFLAMMATION",
                                      "AXIS_CELL_CYCLE")]))
colnames(hm_mat) <- ax_clean_df$sample_id

ann_col <- data.frame(
  Subtype = ax_clean_df$cluster_final,
  row.names = ax_clean_df$sample_id
)

# Order samples by subtype
sample_order <- order(ax_clean_df$cluster_final)
hm_mat_ordered <- hm_mat[, sample_order]
ann_col_ordered <- ann_col[sample_order, , drop=FALSE]

pdf("sporadic_subtype_heatmap_FINAL.pdf", width = 12, height = 3)
pheatmap(
  hm_mat_ordered,
  annotation_col    = ann_col_ordered,
  annotation_colors = list(Subtype = subtype_colors),
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,       # ordered by subtype
  show_colnames     = FALSE,
  color = colorRampPalette(c("#2166AC","white","#D6604D"))(100),
  breaks = seq(-3, 3, length.out = 101),
  main = "Sporadic ALS molecular subtypes — axis scores"
)
dev.off()

# ============================================================
# RESPONSE TO TEAM COMMENTS
# ============================================================

library(dplyr)
library(ggplot2)
library(ConsensusClusterPlus)
library(pheatmap)

subtype_colors <- c("S-ALS-1"="#E64B35",
                    "S-ALS-2"="#4DBBD5",
                    "S-ALS-3"="#00A087")

# ============================================================
# COMMENT 1: Is S-ALS-2 actually two clusters?
# ============================================================

# --- 1a: PCA of ssGSEA axis scores ---
# Use the clean axis scores (ax_clean_df)

ax_pca_input <- ax_clean_df[, c("AXIS_METABOLIC",
                                "AXIS_INFLAMMATION",
                                "AXIS_CELL_CYCLE")]
rownames(ax_pca_input) <- ax_clean_df$sample_id

pca_res <- prcomp(ax_pca_input, scale. = TRUE)

pca_df <- data.frame(
  sample_id = ax_clean_df$sample_id,
  PC1       = pca_res$x[, 1],
  PC2       = pca_res$x[, 2],
  PC3       = pca_res$x[, 3],
  cluster   = ax_clean_df$cluster_final
)

# Variance explained
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
cat("Variance explained: PC1=", var_exp[1],
    "% PC2=", var_exp[2],
    "% PC3=", var_exp[3], "%\n")

# PCA plot coloured by subtype
p_pca <- ggplot(pca_df, aes(PC1, PC2, color=cluster)) +
  geom_point(size=2, alpha=0.8) +
  scale_color_manual(values=subtype_colors) +
  theme_classic(base_size=13) +
  labs(
    x     = paste0("PC1 (", var_exp[1], "%)"),
    y     = paste0("PC2 (", var_exp[2], "%)"),
    color = "Subtype",
    title = "PCA of ssGSEA axis scores — sporadic ALS"
  )

ggsave("pca_ssgsea_subtypes.pdf", p_pca, width=6, height=5)

# --- 1b: Focus on S-ALS-2 only — does it look bimodal? ---
als2_df <- ax_pca_input[ax_clean_df$cluster_final == "S-ALS-2", ]
cat("\nS-ALS-2 n =", nrow(als2_df), "\n")

# PCA within S-ALS-2
pca_als2 <- prcomp(als2_df, scale.=TRUE)
pca_als2_df <- data.frame(
  sample_id = rownames(als2_df),
  PC1       = pca_als2$x[,1],
  PC2       = pca_als2$x[,2]
)

# Visual check for bimodality
p_als2 <- ggplot(pca_als2_df, aes(PC1, PC2)) +
  geom_point(size=2, alpha=0.7, color="#4DBBD5") +
  geom_density_2d(color="black", alpha=0.5) +
  theme_classic(base_size=13) +
  labs(title="S-ALS-2 internal structure (PCA)",
       subtitle="Look for two density peaks")

ggsave("pca_sals2_internal.pdf", p_als2, width=5, height=4)

# Statistical test for bimodality: Hartigan's dip test on PC1
library(diptest)
dip_result <- dip.test(pca_als2_df$PC1)
cat("\nHartigan's dip test (S-ALS-2 PC1):\n")
print(dip_result)
# p < 0.05 → evidence of bimodality

# --- 1c: Force k=4 and k=5 on S-ALS-2 only ---
# If bimodality exists, consensus clustering will find it
als2_mat <- t(as.matrix(als2_df))

set.seed(42)
cc_als2 <- ConsensusClusterPlus(
  als2_mat,
  maxK=4, reps=500,
  pItem=0.8, pFeature=1.0,
  clusterAlg="hc", distance="euclidean",
  seed=42, plot="pdf",
  title="consensus_sals2_substructure"
)

pac_als2 <- sapply(2:4, function(k) {
  m <- cc_als2[[k]]$consensusMatrix
  sum(m > 0.1 & m < 0.9) / (nrow(m)*ncol(m))
})
cat("\nPAC within S-ALS-2:\n")
print(data.frame(k=2:4, PAC=round(pac_als2, 4)))
# If PAC for k=2 is low → there really are 2 subclusters
# If PAC stays high → S-ALS-2 is genuinely one homogeneous group

# --- 1d: Check axis score distributions within S-ALS-2 ---
cat("\nS-ALS-2 axis score distributions:\n")
als2_full <- ax_clean_df[ax_clean_df$cluster_final == "S-ALS-2", ]
for (ax in c("AXIS_METABOLIC","AXIS_INFLAMMATION","AXIS_CELL_CYCLE")) {
  vals <- als2_full[[ax]]
  cat(sprintf("  %s: mean=%.3f, sd=%.3f, skew=%.3f\n",
              ax, mean(vals), sd(vals),
              mean((vals - mean(vals))^3) / sd(vals)^3))
}

# ============================================================
# COMMENT 2: Remove S100b-correlated genes from cell cycle axis
# ============================================================

# --- 2a: Identify S100b-correlated genes within cell cycle sets ---

sporadic_ids <- ax_clean_df$sample_id
expr_sp <- expr_symbol[, sporadic_ids]

# S100b values for sporadic samples
sporadic_meta_sub <- df_nfl_uniq[df_nfl_uniq$sample_id %in% sporadic_ids, ]
s100b_vec <- setNames(sporadic_meta_sub$PCT_S100b,
                      sporadic_meta_sub$sample_id)
s100b_vec <- s100b_vec[sporadic_ids]  # align order
s100b_vec <- s100b_vec[!is.na(s100b_vec)]

shared <- intersect(names(s100b_vec), colnames(expr_sp))
expr_for_cor <- expr_sp[, shared]
s100b_aligned <- s100b_vec[shared]

# Cell cycle genes
cc_genes_all <- unique(c(
  hallmark_list[["HALLMARK_E2F_TARGETS"]],
  hallmark_list[["HALLMARK_G2M_CHECKPOINT"]],
  hallmark_list[["HALLMARK_MYC_TARGETS_V1"]],
  hallmark_list[["HALLMARK_DNA_REPAIR"]]
))
cc_genes_present <- cc_genes_all[cc_genes_all %in% rownames(expr_for_cor)]

cc_s100b_cors <- sapply(cc_genes_present, function(g) {
  cor(expr_for_cor[g,], s100b_aligned, use="complete.obs")
})

cat("\n=== Cell cycle gene ~ S100b correlation distribution ===\n")
print(quantile(cc_s100b_cors, probs=c(0,.1,.25,.5,.75,.9,1), na.rm=TRUE))

for (t in c(0.2, 0.3, 0.4)) {
  cat(sprintf("  |r| > %.1f with S100b: %d / %d (%.1f%%)\n",
              t,
              sum(abs(cc_s100b_cors) > t, na.rm=TRUE),
              length(cc_s100b_cors),
              100*mean(abs(cc_s100b_cors) > t, na.rm=TRUE)))
}

# --- 2b: Build S100b-filtered cell cycle gene set ---
# Use same threshold as inflammation axis (0.3)
s100b_threshold <- 0.3

cc_genes_filtered <- names(cc_s100b_cors)[
  !is.na(cc_s100b_cors) & abs(cc_s100b_cors) < s100b_threshold
]

cat(sprintf("\nCell cycle genes retained after S100b filtering: %d / %d\n",
            length(cc_genes_filtered),
            length(cc_genes_present)))

# --- 2c: Recompute axis scores with filtered gene set ---
axis_sets_s100b_clean <- axis_sets
axis_sets_s100b_clean[["AXIS_CELL_CYCLE"]] <- cc_genes_filtered

# Run on sporadic expression (residualized)
expr_sp_resid <- expr_resid[, sporadic_ids]
expr_sp_resid_mat <- as.matrix(expr_sp_resid)
storage.mode(expr_sp_resid_mat) <- "double"

param_s100b <- ssgseaParam(expr_sp_resid_mat,
                           axis_sets_s100b_clean,
                           normalize=TRUE)
axis_scores_s100b <- gsva(param_s100b)

ax_s100b_df <- as.data.frame(t(axis_scores_s100b))
ax_s100b_df <- as.data.frame(scale(ax_s100b_df))
ax_s100b_df$sample_id <- rownames(ax_s100b_df)

# --- 2d: Recluster and check ARI vs original ---
ax_s100b_mat <- t(as.matrix(ax_s100b_df[,
                                        c("AXIS_METABOLIC","AXIS_INFLAMMATION","AXIS_CELL_CYCLE")]))

set.seed(42)
cc_s100b <- ConsensusClusterPlus(
  ax_s100b_mat,
  maxK=5, reps=1000,
  pItem=0.8, pFeature=1.0,
  clusterAlg="hc", distance="euclidean",
  seed=42, plot="pdf",
  title="consensus_s100b_filtered_cellcycle"
)

pac_s100b <- sapply(2:5, function(k) {
  m <- cc_s100b[[k]]$consensusMatrix
  sum(m > 0.1 & m < 0.9) / (nrow(m)*ncol(m))
})
cat("\nPAC (S100b-filtered cell cycle gene set):\n")
print(data.frame(k=2:5, PAC=round(pac_s100b,4)))

k_s100b <- 3
cluster_s100b <- cc_s100b[[k_s100b]]$consensusClass
names(cluster_s100b) <- colnames(ax_s100b_mat)

ax_s100b_df$cluster_s100b <- factor(
  cluster_s100b[ax_s100b_df$sample_id],
  levels=1:k_s100b, labels=paste0("S-ALS-",1:k_s100b)
)

# ARI vs primary clusters
library(aricode)
orig_labels <- ax_clean_df$cluster_final[
  match(ax_s100b_df$sample_id, ax_clean_df$sample_id)
]

ari_s100b <- ARI(as.integer(orig_labels),
                 as.integer(ax_s100b_df$cluster_s100b))
cat(sprintf("\nARI (primary clusters vs S100b-filtered): %.3f\n", ari_s100b))
# If ARI > 0.7 → clusters are robust to S100b gene removal
# Report this as a sensitivity analysis regardless of result

# Check residual correlation
ax_s100b_check <- merge(ax_s100b_df, marker_data, by="sample_id", all.x=TRUE)
cat("\nResidual correlations after S100b gene filtering:\n")
for (ax in c("AXIS_METABOLIC","AXIS_INFLAMMATION","AXIS_CELL_CYCLE")) {
  r <- cor(ax_s100b_check[[ax]], ax_s100b_check[["PCT_S100b"]],
           use="complete.obs")
  cat(sprintf("  %s ~ PCT_S100b: r = %.3f\n", ax, r))
}

# ============================================================
# COMMENT 3: Downsample S-ALS-2 to balance DEG comparisons
# ============================================================

# The concern: S-ALS-1 (n=96) vs S-ALS-3 (n=77) gives 338 DEGs
#              S-ALS-1 (n=96) vs S-ALS-2 (n=227) gives 3240 DEGs
# Is this biology or sample size?

# Approach: downsample S-ALS-2 to n=96 (size of S-ALS-1)
# Repeat 10 times, check whether DEG counts are stable

n_target    <- min(
  sum(ax_clean_df$cluster_final == "S-ALS-1"),
  sum(ax_clean_df$cluster_final == "S-ALS-3")
)
cat("Downsampling S-ALS-2 to n =", n_target, "\n")

# Samples per subtype
s1_ids <- ax_clean_df$sample_id[ax_clean_df$cluster_final == "S-ALS-1"]
s2_ids <- ax_clean_df$sample_id[ax_clean_df$cluster_final == "S-ALS-2"]
s3_ids <- ax_clean_df$sample_id[ax_clean_df$cluster_final == "S-ALS-3"]

n_iter <- 10
downsample_results <- list()

for (i in 1:n_iter) {
  set.seed(i * 100)
  
  # Downsample S-ALS-2 to n_target
  s2_sampled <- sample(s2_ids, n_target)
  
  # Build metadata for this iteration
  keep_ids  <- c(s1_ids, s2_sampled, s3_ids)
  meta_iter <- sporadic_meta_emt[sporadic_meta_emt$Row.names %in% keep_ids, ]
  
  counts_iter <- counts_filtered[, rownames(meta_iter)]
  
  sv_terms <- paste0("SV", seq_along(clean_svs), collapse=" + ")
  design_iter <- as.formula(
    paste("~ SEX +", sv_terms, "+ EMT_score_scaled + cluster_clean")
  )
  
  dds_iter <- DESeqDataSetFromMatrix(
    countData = counts_iter,
    colData   = meta_iter,
    design    = design_iter
  )
  
  # Use existing size factors and dispersions to save time
  suppressMessages(dds_iter <- estimateSizeFactors(dds_iter))
  suppressMessages(dds_iter <- estimateDispersions(dds_iter))
  suppressMessages(dds_iter <- nbinomWaldTest(dds_iter))
  
  res_12 <- results(dds_iter, contrast=c("cluster_clean","S-ALS-1","S-ALS-2"))
  res_13 <- results(dds_iter, contrast=c("cluster_clean","S-ALS-1","S-ALS-3"))
  res_23 <- results(dds_iter, contrast=c("cluster_clean","S-ALS-2","S-ALS-3"))
  
  count_degs <- function(res) {
    sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 0.5)
  }
  
  downsample_results[[i]] <- data.frame(
    iter    = i,
    n_s2    = n_target,
    degs_12 = count_degs(res_12),
    degs_13 = count_degs(res_13),
    degs_23 = count_degs(res_23)
  )
  
  cat(sprintf("Iter %d: 1v2=%d, 1v3=%d, 2v3=%d\n",
              i,
              downsample_results[[i]]$degs_12,
              downsample_results[[i]]$degs_13,
              downsample_results[[i]]$degs_23))
}

ds_df <- bind_rows(downsample_results)

cat("\n=== Downsampled DEG counts (S-ALS-2 matched to n =", n_target, ") ===\n")
cat("Full analysis (original n):  1v2=3240, 1v3=338, 2v3=1883\n\n")
print(ds_df %>%
        summarise(
          mean_12 = round(mean(degs_12)),
          sd_12   = round(sd(degs_12)),
          mean_13 = round(mean(degs_13)),
          sd_13   = round(sd(degs_13)),
          mean_23 = round(mean(degs_23)),
          sd_23   = round(sd(degs_23))
        ))

# Key interpretation:
# If downsampled 1v2 drops dramatically toward 1v3 levels
#   → DEG counts were partly sample-size driven
# If downsampled 1v2 stays substantially higher than 1v3
#   → the biology is genuinely different between these pairs

write.csv(ds_df, "deseq2_downsample_sensitivity.csv", row.names=FALSE)

# ============================================================
# COMMENT 4: Genomic subtyping — independent of transcriptomics
# ============================================================

# This requires the burden_all object from the genomic pipeline
# Run that first, then come back here

# Assuming burden_all exists with columns:
# cgnd_id, total_damaging, total_lof, total_nonsyn
# and cluster_cgnd bridges cgnd_id to cluster_final

# --- 4a: Does genomic burden profile cluster independently? ---

# Build per-sample feature matrix for genomic clustering
# Use: total_damaging, total_lof, total_nonsyn (per chromosome if available)
# or just genome-wide counts

burden_sporadic <- burden_all %>%
  inner_join(cluster_cgnd[, c("cgnd_id","cluster_final")],
             by="cgnd_id")

# Scale features
burden_mat <- as.matrix(
  scale(burden_sporadic[, c("total_damaging","total_lof","total_nonsyn")])
)
rownames(burden_mat) <- burden_sporadic$cgnd_id

cat("\nGenomic burden clustering input: n =", nrow(burden_mat), "\n")

# Consensus clustering on genomic burden
set.seed(42)
cc_genomic <- ConsensusClusterPlus(
  t(burden_mat),
  maxK=5, reps=500,
  pItem=0.8, pFeature=1.0,
  clusterAlg="hc", distance="euclidean",
  seed=42, plot="pdf",
  title="consensus_genomic_burden"
)

pac_genomic <- sapply(2:5, function(k) {
  m <- cc_genomic[[k]]$consensusMatrix
  sum(m > 0.1 & m < 0.9) / (nrow(m)*ncol(m))
})
cat("\nPAC (genomic burden clustering):\n")
print(data.frame(k=2:5, PAC=round(pac_genomic,4)))

# If PAC is uniformly high (>0.3 for all k)
#   → genomic burden alone does not form stable clusters
#   → report as null finding: transcriptomic subtypes are
#      not explained by aggregate mutation burden

# Assign labels at best k
k_genomic <- which.min(pac_genomic) + 1
cat("Best k by PAC:", k_genomic, "\n")

cluster_genomic <- cc_genomic[[k_genomic]]$consensusClass
names(cluster_genomic) <- rownames(burden_mat)

burden_sporadic$cluster_genomic <- factor(
  cluster_genomic[burden_sporadic$cgnd_id],
  levels=1:k_genomic,
  labels=paste0("G-",1:k_genomic)
)

# --- 4b: Does genomic subtype agree with transcriptomic subtype? ---
agreement_tab <- table(
  Transcriptomic = burden_sporadic$cluster_final,
  Genomic        = burden_sporadic$cluster_genomic
)
cat("\n=== Transcriptomic vs Genomic subtype agreement ===\n")
print(agreement_tab)

ari_genomic <- ARI(
  as.integer(burden_sporadic$cluster_final),
  as.integer(burden_sporadic$cluster_genomic)
)
cat(sprintf("\nARI (transcriptomic vs genomic): %.3f\n", ari_genomic))

# Low ARI: genomic burden and transcriptomics capture different axes of variation
# This is actually an interesting result either way:
# High ARI → genomic architecture shapes transcriptomic state
# Low ARI  → transcriptomic subtypes reflect post-genomic biology

# --- 4c: Clinical associations of genomic subtypes ---
burden_clinical_g <- burden_sporadic %>%
  left_join(
    df_nfl_uniq %>%
      left_join(crosswalk, by=c("sample_id"="expr_sample_id")) %>%
      dplyr::select(cgnd_id, logNfL, ALSFRS_R_PROGRESSION_SLOPE,
                    AGE_AT_SYMPTOM_ONSET, SEX),
    by="cgnd_id"
  ) %>%
  filter(!is.na(cluster_genomic))

# NfL by genomic subtype
m_nfl_g <- lm(
  logNfL ~ cluster_genomic + AGE_AT_SYMPTOM_ONSET + SEX,
  data = burden_clinical_g %>% filter(!is.na(logNfL))
)
cat("\n=== NfL ~ genomic subtype ===\n")
print(summary(m_nfl_g)$coefficients)

# Compare: if transcriptomic subtypes predict NfL but genomic don't
# → transcriptomics adds information beyond genomic architecture


