# 03b_marker_rescue_rerun.R
# Rescue the 7 marker genes missing from the backbone symbol matrix
# (GFAP; AIF1, CX3CR1, P2RY12, CD68, TYROBP; MOG) from the full VST
# matrix, then rerun ONLY the marker-dependent sensitivity analyses
# (S1 S100B ARI, S3 broad-panel R2 + ARI, S4 CCT2 ladder) with the
# complete panels. S2 uses maturity/progenitor panels that were already
# 14/14 and 9/9, so its numbers do not change, but it is recomputed and
# printed for completeness.
#
# Why the genes are missing: the backbone symbol matrix keeps 18,606
# genes after filtering and symbol mapping; these 7 are expressed low
# in iPSC-MN cultures (glial genes in a neuronal prep) and fell below
# the filter or lost their mapping. Rescuing them from the source VST
# reproduces the old pipeline's fuller panels.
#
# Inputs : cohort_inputs_FINAL.rds (for the source VST), backbone_FINAL.rds,
#          patient_level_FINAL.rds, diff_quality_sensitivity_FINAL.rds (03 run)
# Outputs: overwrites diff_quality_sensitivity_FINAL.rds (03 version renamed
#          _panels_as_backbone), sensitivity_summary.csv updated

suppressPackageStartupMessages({
  library(dplyr)
  library(mclust)
  library(ConsensusClusterPlus)
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::mutate,
                               dplyr::summarise, dplyr::slice,
                               base::intersect, base::setdiff, base::union,
                               base::unname, .quiet = TRUE)
}
set.seed(42)
setwd("~/Documents/DSN/ALS")

## ---- Load everything ----
b   <- readRDS("backbone_FINAL.rds")
pl  <- readRDS("patient_level_FINAL.rds")
s03 <- readRDS("diff_quality_sensitivity_FINAL.rds")
ci  <- readRDS("~/Documents/cohort_inputs_FINAL.rds")

expr_sym <- if (is.list(b) && "expr_symbol" %in% names(b)) b$expr_symbol else b
meta     <- s03$meta        # carries all scores + subtype labels from 03
AX <- c("AXIS_METABOLIC", "AXIS_CELL_CYCLE", "AXIS_INFLAMMATION")
meta$cluster_final <- factor(meta$cluster_final, levels = paste0("S-ALS-", 1:3))
stopifnot(nrow(meta) == 391)

## ---- Rescue the 7 genes from the source counts (VST them fresh) ----
missing_syms <- c("GFAP", "AIF1", "CX3CR1", "P2RY12", "CD68", "TYROBP", "MOG")

## Map symbols -> ENSEMBL, find rows in the corrected counts matrix
ens_map <- AnnotationDbi::select(org.Hs.eg.db, keys = missing_syms,
                                 keytype = "SYMBOL", columns = "ENSEMBL") %>%
  filter(!is.na(ENSEMBL))
counts_rows <- sub("\\..*$", "", rownames(ci$counts_filtered_corrected))

found <- ens_map %>% filter(ENSEMBL %in% counts_rows) %>%
  distinct(SYMBOL, .keep_all = TRUE)
cat("Rescuable from corrected counts:", nrow(found), "of 7:",
    paste(found$SYMBOL, collapse = ", "), "\n")
still_missing <- setdiff(missing_syms, found$SYMBOL)
if (length(still_missing) > 0)
  cat("Not in counts either (below the count filter; will stay absent):",
      paste(still_missing, collapse = ", "), "\n")

## VST the corrected counts once, pull the rescued rows, subset to the
## 391 patient samples in meta
dds_r <- DESeqDataSetFromMatrix(ci$counts_filtered_corrected,
                                ci$annot_corrected, design = ~ 1)
vst_r <- assay(vst(dds_r, blind = TRUE))
rownames(vst_r) <- sub("\\..*$", "", rownames(vst_r))

sample_key <- intersect(c("sample_id", "Sample_ID"), names(meta))[1]
stopifnot(all(meta[[sample_key]] %in% colnames(vst_r)))

rescue_mat <- vst_r[found$ENSEMBL, meta[[sample_key]], drop = FALSE]
rownames(rescue_mat) <- found$SYMBOL
rm(dds_r, vst_r); gc()

## ---- Rebuild the full-panel scores ----
expr_sym <- expr_sym[, meta[[sample_key]], drop = FALSE]
stopifnot(identical(colnames(expr_sym), meta[[sample_key]]))

## Combined matrix: backbone genes + rescued rows (backbone wins on overlap)
add_rows <- setdiff(rownames(rescue_mat), rownames(expr_sym))
expr_full <- rbind(expr_sym, rescue_mat[add_rows, , drop = FALSE])
cat("Expression matrix for panels:", nrow(expr_full), "genes",
    "(", length(add_rows), "rescued )\n")

astrocyte_markers  <- c("GFAP","AQP4","ALDH1L1","S100B","SLC1A3",
                        "SLC1A2","GJA1","AGT")
microglia_markers  <- c("AIF1","CX3CR1","P2RY12","TMEM119","CD68",
                        "ITGAM","TYROBP")
oligodendrocyte_mk <- c("OLIG2","MBP","MOG","PLP1","SOX10","GALC")

present <- function(g) intersect(g, rownames(expr_full))
for (nm in c("astrocyte_markers","microglia_markers","oligodendrocyte_mk")) {
  g <- get(nm)
  cat(nm, ":", length(present(g)), "/", length(g),
      " missing:", paste(setdiff(g, present(g)), collapse = ","), "\n")
}

zrow  <- function(x) (x - mean(x)) / sd(x)
expr_z <- t(apply(expr_full, 1, zrow))

meta$astro_score <- colMeans(expr_z[present(astrocyte_markers), , drop = FALSE])
meta$micro_score <- colMeans(expr_z[present(microglia_markers), , drop = FALSE])
meta$oligo_score <- colMeans(expr_z[present(oligodendrocyte_mk), , drop = FALSE])
meta$S100B_expr  <- as.numeric(expr_z["S100B", ])
## diff_quality, maturity, progenitor, prolif scores carry over from 03
## unchanged (their panels were complete)

recluster_ari <- function(mat3xN, title) {
  while (!is.null(grDevices::dev.list())) try(grDevices::dev.off(), silent = TRUE)
  ccp <- ConsensusClusterPlus(mat3xN, maxK = 3, reps = 1000, pItem = 0.8,
                              clusterAlg = "km", distance = "euclidean",
                              seed = 42, title = title, plot = "png")
  cl <- ccp[[3]]$consensusClass
  adjustedRandIndex(cl[meta[[sample_key]]], as.integer(meta$cluster_final))
}

## ---- S1 rerun: S100B residualisation (unchanged inputs, rerun for a
##      single consistent output file) ----
res_s100 <- sapply(AX, function(ax) residuals(lm(meta[[ax]] ~ meta$S100B_expr)))
m1 <- t(scale(res_s100)); colnames(m1) <- meta[[sample_key]]
ari_s100b <- recluster_ari(m1, "ccp_S1_s100b_fullpanel")
cat("\nS1 ARI (S100B):", round(ari_s100b, 3),
    " [03 run: 0.565 | old pipeline: 0.791]\n")

## ---- S3 rerun: broad panel with full glial coverage ----
r2 <- c()
res_broad <- sapply(AX, function(ax) {
  fit <- lm(meta[[ax]] ~ astro_score + micro_score + oligo_score +
              diff_quality, data = meta)
  r2[[ax]] <<- summary(fit)$adj.r.squared
  cat(" ", ax, "adjusted R^2:", round(summary(fit)$adj.r.squared, 3), "\n")
  residuals(fit)
})
cat("  [03 run: 0.197 / 0.398 / 0.605 | old: 0.257 / 0.485 / 0.659]\n")
m3 <- t(scale(res_broad)); colnames(m3) <- meta[[sample_key]]
ari_broad <- recluster_ari(m3, "ccp_S3_broad_fullpanel")
cat("S3 ARI (broad):", round(ari_broad, 3),
    " [03 run: 0.273 | old: 0.244]\n")

## ---- S4 rerun: CCT2 ladder with full-panel glial scores ----
m23 <- meta %>% filter(cluster_final %in% c("S-ALS-2", "S-ALS-3")) %>%
  mutate(cluster_final = droplevels(cluster_final))
f1 <- lm(CCT2 ~ cluster_final, data = m23)
f2 <- lm(CCT2 ~ cluster_final + prolif_score + maturity_score +
           diff_quality, data = m23)
f3 <- lm(CCT2 ~ cluster_final + prolif_score + maturity_score +
           diff_quality + astro_score + micro_score + oligo_score, data = m23)
cat("\nS4 CCT2 S-ALS-3 coefficient (unadj / composition / full glial):\n")
for (f in list(f1, f2, f3)) print(summary(f)$coefficients["cluster_finalS-ALS-3", ])

## ---- S2 (unchanged; reprint) ----
aov1 <- aov(diff_quality ~ cluster_final, data = meta)
ss <- summary(aov1)[[1]][["Sum Sq"]]
eta_sq <- ss[1] / sum(ss)
cat("\nS2 eta-squared (unchanged):", round(eta_sq, 3), "\n")

## ---- Save: full-panel version becomes FINAL ----
file.rename("diff_quality_sensitivity_FINAL.rds",
            "diff_quality_sensitivity_panels_as_backbone.rds")
sensitivity_summary <- data.frame(
  metric = c("ari_s100b", "eta_sq_diff_quality", "ari_broad", paste0("r2_", AX)),
  value  = c(ari_s100b, eta_sq, ari_broad, unlist(r2)),
  run03_backbone_panels = c(0.565, 0.382, 0.273, 0.197, 0.398, 0.605),
  old_pipeline = c(0.791, 0.346, 0.244, 0.257, 0.485, 0.659)
)
print(sensitivity_summary)
write.csv(sensitivity_summary, "sensitivity_summary.csv", row.names = FALSE)
write.csv(meta, "diff_quality_sensitivity_results.csv", row.names = FALSE)
saveRDS(list(summary = sensitivity_summary, meta = meta,
             cct2_models = list(unadjusted = summary(f1)$coefficients,
                                composition = summary(f2)$coefficients,
                                full = summary(f3)$coefficients),
             rescued_genes = found$SYMBOL,
             unrescuable = still_missing,
             panels = "full panels after VST rescue of low-expression glial genes",
             ccp_params = "maxK3 reps1000 pItem0.8 km euclidean seed42"),
        "diff_quality_sensitivity_FINAL.rds")
cat("\nSaved diff_quality_sensitivity_FINAL.rds (full-panel version)\n")
sessionInfo()
