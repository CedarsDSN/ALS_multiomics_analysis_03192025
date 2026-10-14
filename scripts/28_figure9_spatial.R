# =====================================================================
# FIGURE 9  Spatial validation, Target ALS postmortem cortex
#   a  CCT2 in ALS vs control, FRONTAL CORTEX ONLY
#   b  within-ALS coupling: CCT2 vs the cell cycle axis
#   c  the same with CCT2 removed from the axis gene set
#
# Output: Figure9_spatial.pdf, figure9_stats.rds
#
# WHAT WAS WRONG WITH THE PREVIOUS VERSION
#   1. It used ALL 43 samples (30 ALS, 13 control), giving Wilcoxon
#      p = 0.099 and rho = 0.55. Motor cortex contributes only 4 samples
#      (2 ALS, 2 control), too few for a group comparison, and mixing
#      the two regions is not the analysis Methods describes. Restricting
#      to frontal cortex gives 28 ALS, 11 control, p = 0.024, rho = 0.62,
#      which are the verified numbers.
#   2. Panel c was missing because the CCT2-excluded axis is stored as a
#      SEPARATE ELEMENT of the pinned object (axis_cc_noCCT2, a named
#      vector over sample_id), not as a column of results_df. The script
#      only looked for a column.
#
# Paths are resolved here rather than through p_().
# =====================================================================

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
if (file.exists("00_figure_common.R")) source("00_figure_common.R")
while (!is.null(grDevices::dev.list())) grDevices::dev.off()

SPATIAL <- local({
  cands <- path.expand(c(
    "~/Documents/DSN/ALS/targetALS/spatial_results_FINAL.rds",
    "~/Documents/targetALS/spatial_results_FINAL.rds",
    "targetALS/spatial_results_FINAL.rds",
    "spatial_results_FINAL.rds"))
  hit <- cands[file.exists(cands)]
  if (!length(hit)) stop("spatial_results_FINAL.rds not found. Tried:\n  ",
                         paste(cands, collapse = "\n  "))
  hit[1]
})
cat("spatial:", SPATIAL, "\n")

if (!exists("thm")) {
  thm <- theme_classic(base_size = 8) +
    theme(axis.text = element_text(size = 7, colour = "black"),
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.key.size = unit(3, "mm"),
          plot.title = element_text(size = 9, hjust = 0),
          plot.subtitle = element_text(size = 7),
          plot.margin = margin(3, 3, 3, 3))
}
if (!exists("W")) W <- 180 / 25.4

sp <- readRDS(SPATIAL)
cat("elements:", paste(names(sp), collapse = ", "), "\n")
stopifnot(!is.null(sp$results_df))
dd <- as.data.frame(sp$results_df)
cat("all samples:", nrow(dd), "\n")
print(table(dd$tissue, dd$clinical_group))

## ---- attach the CCT2-excluded axis, stored as a separate element -----
if (!is.null(sp$axis_cc_noCCT2)) {
  v <- sp$axis_cc_noCCT2
  dd$AXIS_CC_noCCT2 <- as.numeric(v[match(dd$sample_id, names(v))])
  cat("attached axis_cc_noCCT2 for", sum(!is.na(dd$AXIS_CC_noCCT2)),
      "of", nrow(dd), "samples\n")
} else {
  dd$AXIS_CC_noCCT2 <- NA_real_
  cat("*** axis_cc_noCCT2 not present; panel c will be skipped\n")
}

## ---- FRONTAL CORTEX ONLY. This is the analysis Methods describes. ----
dd <- dd[dd$tissue == "Frontal Cortex", ]
dd$clinical_group <- factor(dd$clinical_group,
                            levels = c("Non-Neurological Control", "ALS"))
nn <- table(dd$clinical_group)
cat("\nfrontal cortex only:", nrow(dd), "samples\n"); print(nn)
stopifnot(nrow(dd) == 39)

## ---- a: case-control ---------------------------------------------------
wt <- wilcox.test(CCT2 ~ clinical_group, data = dd)
cat(sprintf("CCT2 ALS vs control: Wilcoxon p = %.4f   (expect 0.024)\n", wt$p.value))

pA <- ggplot(dd, aes(clinical_group, CCT2, fill = clinical_group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.55) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.7, colour = "grey25") +
  scale_fill_manual(values = c("Non-Neurological Control" = "#7F7F7F",
                               "ALS" = "#5B6C8F")) +   # CHANGED: neutral, not the S-ALS-1 red
  scale_x_discrete(labels = c(paste0("Control\nn=", nn[1]),
                              paste0("ALS\nn=", nn[2]))) +
  labs(title = "a   CCT2 in postmortem frontal cortex",          # CHANGED: spelling
       subtitle = sprintf("Wilcoxon P = %.3f", wt$p.value),
       x = NULL, y = "CCT2 expression (VST)") +                    # CHANGED: units
  thm + theme(legend.position = "none")

## ---- b, c: within-ALS --------------------------------------------------
ad <- dd[dd$clinical_group == "ALS", ]
cat("within-ALS n =", nrow(ad), "\n")

## CHANGED: z-score both cell cycle axes so panels b and c share one scale.
## Spearman rho and P are computed on the raw columns below and are unchanged.
ad$cc_z    <- as.numeric(scale(ad$AXIS_CELL_CYCLE))
ad$cc_ex_z <- if (any(!is.na(ad$AXIS_CC_noCCT2))) as.numeric(scale(ad$AXIS_CC_noCCT2)) else NA_real_

ct_cc  <- cor.test(ad$CCT2, ad$AXIS_CELL_CYCLE,   method = "spearman")
ct_met <- cor.test(ad$CCT2, ad$AXIS_METABOLIC,    method = "spearman")
ct_inf <- cor.test(ad$CCT2, ad$AXIS_INFLAMMATION, method = "spearman")
cat(sprintf("cell cycle  : rho = %.3f, p = %.2g   (expect 0.62, 5.2e-4)\n",
            ct_cc$estimate, ct_cc$p.value))
cat(sprintf("metabolic   : rho = %.3f, p = %.2g\n", ct_met$estimate, ct_met$p.value))
cat(sprintf("inflammation: rho = %.3f, p = %.2g\n", ct_inf$estimate, ct_inf$p.value))

mk_scatter <- function(ycol, ttl, sub) {
  ggplot(ad, aes(CCT2, .data[[ycol]])) +
    geom_point(size = 1.3, alpha = 0.8, colour = "grey20") +   # CHANGED: neutral, not the S-ALS-3 colour
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "grey25", linewidth = 0.5) +
    labs(title = ttl, subtitle = sub,
         x = "CCT2 expression (VST)", y = "Cell cycle axis score (z)") + thm   # CHANGED
}
pB <- mk_scatter("cc_z", "b   CCT2 and cell cycle coupling",       # CHANGED: z column, no hyphen
                 sprintf("Within ALS: rho = %.2f, P = %.1g",
                         ct_cc$estimate, ct_cc$p.value))

pC <- NULL; ct_ex <- NULL
if (any(!is.na(ad$AXIS_CC_noCCT2))) {
  ct_ex <- cor.test(ad$CCT2, ad$AXIS_CC_noCCT2, method = "spearman")
  cat(sprintf("cell cycle, CCT2 excluded: rho = %.4f, p = %.2g\n",
              ct_ex$estimate, ct_ex$p.value))
  identical_ranks <- isTRUE(all.equal(unname(ct_cc$estimate),
                                      unname(ct_ex$estimate), tolerance = 1e-6))
  if (identical_ranks)
    cat("NOTE: identical to the full-axis result. Removing one gene did not\n",
        "     change any sample's rank, so CCT2 contributes negligibly to\n",
        "     the composite. Describe it that way, not as two independent\n",
        "     estimates that agree.\n")
  pC <- mk_scatter("cc_ex_z", "c   CCT2 excluded from the axis",   # CHANGED: z column
                   sprintf("rho = %.2f, P = %.1g%s", ct_ex$estimate, ct_ex$p.value,
                           if (identical_ranks) " (ranks unchanged)" else ""))
}

fig <- if (is.null(pC)) (pA | pB) else (pA | pB | pC)
## CHANGED: ggsave() does not embed fonts on this Mac. Use save_fig() from
## 00_figure_common.R, or the native quartz device if it is not loaded.
if (exists("save_fig")) {
  save_fig(fig, "Figure9_spatial.pdf", height = 2.7)
} else {
  grDevices::quartz(file = "Figure9_spatial.pdf", type = "pdf", width = W, height = 2.7)
  print(fig); grDevices::dev.off()
}
cat("\nwrote Figure9_spatial.pdf\n")
saveRDS(list(frontal = dd, casecontrol = wt,
             cc = ct_cc, metabolic = ct_met, inflammation = ct_inf,
             cc_excl = ct_ex, source = SPATIAL), "figure9_stats.rds")
cat("saved figure9_stats.rds\n")