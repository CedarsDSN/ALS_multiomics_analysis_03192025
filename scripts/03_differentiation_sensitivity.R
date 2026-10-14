# 03_differentiation_sensitivity.R
# Composition / differentiation-quality sensitivity, rerun on the FINAL
# rebuilt labels. File-based port of SensitivityAnalysis.R (March era):
# panels and ConsensusClusterPlus parameters verified IDENTICAL to the
# original (reps 1000, pItem 0.8, km, euclidean), inputs now read from
# pinned rds instead of session objects.
#
# Replaces these manuscript numbers (old pipeline values in brackets):
#   S2 eta-squared for diff quality by subtype        [0.346]
#   S1 ARI after S100B residualisation                [0.791]
#   S3 ARI after broad glial + diff-quality resid     [0.244]
#   S3 axis variance explained                        [25.7 / 48.5 / 65.9%]
# Plus S4: CCT2 S-ALS-3 vs S-ALS-2 under increasing adjustment.
#
# Inputs : backbone_FINAL.rds, patient_level_FINAL.rds
# Outputs: diff_quality_sensitivity_FINAL.rds,
#          diff_quality_sensitivity_results.csv, sensitivity_summary.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(mclust)                # adjustedRandIndex
  library(ConsensusClusterPlus)
})
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::mutate,
                               dplyr::summarise, dplyr::slice,
                               base::intersect, base::setdiff, base::union,
                               base::unname, .quiet = TRUE)
}
set.seed(42)
setwd("~/Documents/DSN/ALS")

## Preserve the stale Aug-13 output before writing anything
if (file.exists("diff_quality_sensitivity_FINAL.rds") &&
    !file.exists("diff_quality_sensitivity_OLD_aug13.rds")) {
  file.rename("diff_quality_sensitivity_FINAL.rds",
              "diff_quality_sensitivity_OLD_aug13.rds")
  cat("Renamed stale Aug-13 rds to _OLD_aug13.\n")
}

## ---- Inputs ----
b  <- readRDS("backbone_FINAL.rds")
pl <- readRDS("patient_level_FINAL.rds")
expr_sym <- if (is.list(b) && "expr_symbol" %in% names(b)) b$expr_symbol else b
meta     <- if (is.data.frame(pl)) pl else pl$pat

AX <- c("AXIS_METABOLIC", "AXIS_CELL_CYCLE", "AXIS_INFLAMMATION")
sub_col <- intersect(c("subtype", "cluster_final", "cluster"), names(meta))[1]
stopifnot(!is.na(sub_col), all(AX %in% names(meta)),
          "sample_id" %in% names(meta))
meta$cluster_final <- factor(meta[[sub_col]],
                             levels = paste0("S-ALS-", 1:3))
stopifnot(all(meta$sample_id %in% colnames(expr_sym)))
expr_sym <- expr_sym[, meta$sample_id, drop = FALSE]
stopifnot(identical(colnames(expr_sym), meta$sample_id))
cat("matrix:", paste(dim(expr_sym), collapse = " x "), "\n")
print(table(meta$cluster_final))          # expect 85 / 204 / 102

## ---- Marker panels (verified identical to SensitivityAnalysis.R) ----
maturity_markers   <- c("MAP2","SYN1","NEFL","NEFM","SNAP25","RBFOX3",
                        "STMN2","CHAT","ISL1","MNX1","SLC18A3",
                        "TUBB3","ELAVL3","ELAVL4")
progenitor_markers <- c("PAX6","SOX2","NES","MKI67","TOP2A","PCNA",
                        "VIM","SOX1","MSI1")
astrocyte_markers  <- c("GFAP","AQP4","ALDH1L1","S100B","SLC1A3",
                        "SLC1A2","GJA1","AGT")
microglia_markers  <- c("AIF1","CX3CR1","P2RY12","TMEM119","CD68",
                        "ITGAM","TYROBP")
oligodendrocyte_mk <- c("OLIG2","MBP","MOG","PLP1","SOX10","GALC")
prolif_markers     <- c("MKI67","PCNA","TOP2A","MCM2","CCNB1","CDK1","BIRC5")

present <- function(g) intersect(g, rownames(expr_sym))
for (nm in c("maturity_markers","progenitor_markers","astrocyte_markers",
             "microglia_markers","oligodendrocyte_mk","prolif_markers")) {
  g <- get(nm)
  cat(nm, ":", length(present(g)), "/", length(g),
      " missing:", paste(setdiff(g, present(g)), collapse = ","), "\n")
}

zrow  <- function(x) (x - mean(x)) / sd(x)
expr_z <- t(apply(expr_sym, 1, zrow))

meta$maturity_score   <- colMeans(expr_z[present(maturity_markers), , drop = FALSE])
meta$progenitor_score <- colMeans(expr_z[present(progenitor_markers), , drop = FALSE])
meta$diff_quality     <- meta$maturity_score - meta$progenitor_score
meta$astro_score      <- colMeans(expr_z[present(astrocyte_markers), , drop = FALSE])
meta$micro_score      <- colMeans(expr_z[present(microglia_markers), , drop = FALSE])
meta$oligo_score      <- colMeans(expr_z[present(oligodendrocyte_mk), , drop = FALSE])
meta$prolif_score     <- colMeans(expr_z[present(prolif_markers), , drop = FALSE])
stopifnot("S100B" %in% rownames(expr_z))
meta$S100B_expr       <- as.numeric(expr_z["S100B", ])

## CCT2 for S4 (from meta if present, else the expression matrix)
if (!"CCT2" %in% names(meta)) {
  stopifnot("CCT2" %in% rownames(expr_sym))
  meta$CCT2 <- as.numeric(expr_sym["CCT2", ])
}

recluster_ari <- function(mat3xN, title) {
  ccp <- ConsensusClusterPlus(mat3xN, maxK = 3, reps = 1000, pItem = 0.8,
                              clusterAlg = "km", distance = "euclidean",
                              seed = 42, title = title, plot = "png")
  cl <- ccp[[3]]$consensusClass
  adjustedRandIndex(cl[meta$sample_id], as.integer(meta$cluster_final))
}

## ---- S1: S100B residualisation -> recluster -> ARI ----
cat("\n===== S1: S100B residualisation =====\n")
res_s100 <- sapply(AX, function(ax)
  residuals(lm(meta[[ax]] ~ meta$S100B_expr)))
m1 <- t(scale(res_s100)); colnames(m1) <- meta$sample_id
ari_s100b <- recluster_ari(m1, "ccp_S1_s100b")
cat("ARI vs primary labels (S100B-residualised):", round(ari_s100b, 3),
    " [old: 0.791]\n")

## ---- S2: differentiation quality by subtype ----
cat("\n===== S2: diff quality ~ subtype =====\n")
aov1 <- aov(diff_quality ~ cluster_final, data = meta)
print(summary(aov1)); print(TukeyHSD(aov1))
ss <- summary(aov1)[[1]][["Sum Sq"]]
eta_sq <- ss[1] / sum(ss)
cat("Eta-squared:", round(eta_sq, 3), " [old: 0.346]\n")
print(meta %>% group_by(cluster_final) %>%
        summarise(n = n(), dq_mean = mean(diff_quality),
                  dq_sd = sd(diff_quality)))

## ---- S3: broad glial + diff-quality residualisation ----
cat("\n===== S3: broad panel residualisation =====\n")
r2 <- c()
res_broad <- sapply(AX, function(ax) {
  fit <- lm(meta[[ax]] ~ astro_score + micro_score + oligo_score +
              diff_quality, data = meta)
  r2[[ax]] <<- summary(fit)$adj.r.squared
  cat(" ", ax, "adjusted R^2:", round(summary(fit)$adj.r.squared, 3), "\n")
  residuals(fit)
})
cat("  [old: 0.257 metabolic / 0.485 cell cycle / 0.659 inflammation]\n")
m3 <- t(scale(res_broad)); colnames(m3) <- meta$sample_id
ari_broad <- recluster_ari(m3, "ccp_S3_broad")
cat("ARI vs primary labels (broad panel):", round(ari_broad, 3),
    " [old: 0.244]\n")

## ---- S4: CCT2 S-ALS-3 vs S-ALS-2 under increasing adjustment ----
cat("\n===== S4: CCT2, S-ALS-2 vs S-ALS-3 =====\n")
m23 <- meta %>% filter(cluster_final %in% c("S-ALS-2", "S-ALS-3")) %>%
  mutate(cluster_final = droplevels(cluster_final))
f1 <- lm(CCT2 ~ cluster_final, data = m23)
f2 <- lm(CCT2 ~ cluster_final + prolif_score + maturity_score +
           diff_quality, data = m23)
f3 <- lm(CCT2 ~ cluster_final + prolif_score + maturity_score +
           diff_quality + astro_score + micro_score + oligo_score,
         data = m23)
for (f in list(f1, f2, f3)) print(summary(f)$coefficients["cluster_finalS-ALS-3", ])

## ---- Save ----
sensitivity_summary <- data.frame(
  metric = c("ari_s100b", "eta_sq_diff_quality", "ari_broad",
             paste0("r2_", AX)),
  value  = c(ari_s100b, eta_sq, ari_broad, unlist(r2)),
  old_value = c(0.791, 0.346, 0.244, 0.257, 0.485, 0.659)
)
print(sensitivity_summary)
write.csv(sensitivity_summary, "sensitivity_summary.csv", row.names = FALSE)
write.csv(meta, "diff_quality_sensitivity_results.csv", row.names = FALSE)
saveRDS(list(summary = sensitivity_summary, meta = meta,
             cct2_models = list(unadjusted = summary(f1)$coefficients,
                                composition = summary(f2)$coefficients,
                                full = summary(f3)$coefficients),
             panels_verified = "identical to SensitivityAnalysis.R",
             ccp_params = "maxK3 reps1000 pItem0.8 km euclidean seed42"),
        "diff_quality_sensitivity_FINAL.rds")
cat("\nSaved diff_quality_sensitivity_FINAL.rds\n")
sessionInfo()
