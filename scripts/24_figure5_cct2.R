# =====================================================================
# FIGURE 5  CCT2 as a single-gene classifier for S-ALS-3 vs S-ALS-2
#   a CCT2 by subtype   b ROC + bootstrap CI
#   c bootstrap AUC     d empirical CDF (uniform shift, not outliers)
# Input:  patient_level_FINAL.rds   (earliest timepoint, 391 patients)
# Output: Figure5_CCT2_classifier.pdf, figure5_stats.rds
#
# AUDIT NOTE: the old version aggregated CCT2 by MEAN across a patient's
# samples and dropped patients whose samples spanned subtypes. The
# canonical dedup is earliest timepoint, already applied in the pinned
# patient table, so this script simply consumes it.
# =====================================================================
source("00_figure_common.R"); library(patchwork); fresh_device()

pat <- get_patients()
stopifnot("CCT2" %in% names(pat))
cat("patients:", nrow(pat), "\n"); print(table(pat$subtype))

sub <- pat[pat$subtype %in% c("S-ALS-2", "S-ALS-3"), ]
sub$subtype <- droplevels(sub$subtype)
pos <- sub$subtype == "S-ALS-3"

auc <- auc_fn(sub$CCT2, pos)
B <- 1000; i3 <- which(pos); i2 <- which(!pos)
boot <- replicate(B, {
  s <- c(sample(i3, length(i3), TRUE), sample(i2, length(i2), TRUE))
  auc_fn(sub$CCT2[s], pos[s]) })
ci <- quantile(boot, c(0.025, 0.975))
cat(sprintf("\nAUC = %.4f  95%% CI %.4f to %.4f   (n = %d S-ALS-3 vs %d S-ALS-2)\n",
            auc, ci[1], ci[2], sum(pos), sum(!pos)))

tt <- t.test(CCT2 ~ subtype, data = sub); wt <- wilcox.test(CCT2 ~ subtype, data = sub)
cat(sprintf("t = %.2f, p = %s | Wilcoxon p = %s\n", tt$statistic,
            format.pval(tt$p.value, digits = 3), format.pval(wt$p.value, digits = 3)))

## robustness: uniform shift rather than outlier-driven
max2 <- max(sub$CCT2[!pos]); pct <- 100 * mean(sub$CCT2[pos] > max2)
top  <- order(sub$CCT2, decreasing = TRUE)[1:15]
cat(sprintf("S-ALS-3 above the S-ALS-2 maximum: %.1f%%\n", pct))
cat(sprintf("AUC after dropping the 15 highest CCT2 patients: %.4f\n",
            auc_fn(sub$CCT2[-top], pos[-top])))
dip_p <- NA
if (requireNamespace("diptest", quietly = TRUE)) {
  dip_p <- diptest::dip.test(sub$CCT2[pos])$p.value
  cat(sprintf("Hartigan dip test on CCT2 (S-ALS-3), p = %.3f\n", dip_p))
}

roc_pts <- function(score, pos) {
  o <- order(score, decreasing = TRUE); l <- pos[o]
  data.frame(fpr = c(0, cumsum(!l) / sum(!l)), tpr = c(0, cumsum(l) / sum(l)))
}

pA <- ggplot(pat, aes(subtype, CCT2, fill = subtype)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.6) +
  geom_jitter(width = 0.16, size = 0.35, alpha = 0.4, colour = "grey25") +
  scale_fill_manual(values = SCOL) +
  scale_x_discrete(labels = paste0(levels(pat$subtype), "\nn=",
                                   as.integer(table(pat$subtype)))) +
  labs(title = "a   CCT2 by subtype", x = NULL, y = "CCT2 (VST)") +
  thm + theme(legend.position = "none")

pB <- ggplot(roc_pts(sub$CCT2, pos), aes(fpr, tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey75", linewidth = 0.3) +
  geom_path(colour = SCOL[["S-ALS-3"]], linewidth = 0.7) +
  annotate("text", x = 0.98, y = 0.06, size = 2.5, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\n95%% CI %.3f to %.3f", auc, ci[1], ci[2])) +
  coord_fixed() +
  labs(title = "b   S-ALS-3 vs S-ALS-2",
       x = "False positive rate", y = "True positive rate") + thm

pC <- ggplot(data.frame(auc = boot), aes(auc)) +
  geom_histogram(bins = 40, fill = "grey75", colour = NA) +
  geom_vline(xintercept = auc, colour = "#E64B35", linewidth = 0.5) +
  geom_vline(xintercept = ci, colour = "grey35", linetype = "dashed", linewidth = 0.35) +
  labs(title = sprintf("c   Bootstrap AUC (%d resamples)", B), x = "AUC", y = "Count") + thm

pD <- ggplot(sub, aes(CCT2, colour = subtype)) +
  stat_ecdf(linewidth = 0.6) + scale_colour_manual(values = SCOL, name = NULL) +
  labs(title = "d   Cumulative distribution",
       subtitle = sprintf("%.0f%% of S-ALS-3 exceed the S-ALS-2 maximum", pct),
       x = "CCT2 (VST)", y = "Cumulative proportion") +
  thm + theme(legend.position = c(0.18, 0.82))

## =====================================================================
## e, f: PROTEIN VALIDATION
## CCT2 protein in the subtyped patients with proteomics, plus the
## RNA-protein correlation. This is the same claim on an orthogonal
## platform, which is why it belongs beside the classifier rather than
## in a figure of its own.
## Landmarks: n = 265 (59/144/62), KW p = 0.005,
##            S3 vs S2 Wilcoxon p = 0.033, RNA-protein rho = 0.21 p = 8e-4,
##            TRiC composite ns (p = 0.13) -> the signal is CCT2-specific.
## =====================================================================
prot_file <- p_("cct2_protein_by_subtype.csv")
pE <- pF <- NULL; prot <- NULL
if (file.exists(prot_file)) {
  prot <- utils::read.csv(prot_file, stringsAsFactors = FALSE)
  cat("\nprotein table:", nrow(prot), "rows; cols:",
      paste(names(prot), collapse = ", "), "\n")
  pc  <- grep("protein|CCT2_prot|abundance", names(prot), value = TRUE)[1]
  psc <- grep("cluster|subtype", names(prot), value = TRUE)[1]
  prc <- grep("^CCT2$|rna|expr|vst", names(prot), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(pc) || is.na(psc))
    stop("could not identify protein / subtype columns in ", basename(prot_file))
  prot$subtype <- factor(as.character(prot[[psc]]), levels = names(SCOL))
  prot <- prot[!is.na(prot$subtype) & !is.na(prot[[pc]]), ]
  nn2 <- table(prot$subtype)
  kwp <- kruskal.test(prot[[pc]] ~ prot$subtype)$p.value
  s23 <- prot[prot$subtype %in% c("S-ALS-2", "S-ALS-3"), ]
  s23$subtype <- droplevels(s23$subtype)
  wp  <- wilcox.test(s23[[pc]] ~ s23$subtype)$p.value
  cat(sprintf("protein: n = %d (%s), KW p = %.4f, S3 vs S2 Wilcoxon p = %.4f\n",
              nrow(prot), paste(as.integer(nn2), collapse = "/"), kwp, wp))

  pE <- ggplot(prot, aes(subtype, .data[[pc]], fill = subtype)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.6) +
    geom_jitter(width = 0.16, size = 0.35, alpha = 0.4, colour = "grey25") +
    scale_fill_manual(values = SCOL) +
    scale_x_discrete(labels = paste0(names(nn2), "\nn=", as.integer(nn2))) +
    scale_y_continuous(labels = function(x) x / 1e5) +
    labs(title = "e   CCT2 protein by subtype",
         subtitle = sprintf("Kruskal-Wallis P = %.3f; S-ALS-3 vs S-ALS-2 P = %.3f", kwp, wp),
         x = NULL, y = expression(CCT2~protein~abundance~(x10^5))) +
    thm + theme(legend.position = "none")

  if (!is.na(prc) && prc != pc) {
    ok <- complete.cases(prot[[prc]], prot[[pc]])
    ct <- cor.test(prot[[prc]][ok], prot[[pc]][ok], method = "spearman")
    cat(sprintf("RNA-protein: rho = %.3f, p = %.2g (n = %d)\n",
                ct$estimate, ct$p.value, sum(ok)))
    pF <- ggplot(prot[ok, ], aes(.data[[prc]], .data[[pc]], colour = subtype)) +
      geom_point(size = 0.8, alpha = 0.7) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                  colour = "grey25", linewidth = 0.45) +
      scale_colour_manual(values = SCOL, name = NULL) +
      scale_y_continuous(labels = function(x) x / 1e5) +
      labs(title = "f   RNA versus protein",
           subtitle = sprintf("Spearman rho = %.2f, P = %.1g", ct$estimate, ct$p.value),
           x = "CCT2 RNA (VST)", y = expression(CCT2~protein~(x10^5))) +
      thm + theme(legend.position = "bottom")
  } else {
    cat("*** no RNA column alongside protein; panel f skipped.\n")
  }
} else {
  cat("\n*** cct2_protein_by_subtype.csv not found in PIN_DIRS.\n")
  cat("*** Panels e and f skipped. Add its directory to PIN_DIRS.\n")
}

if (!is.null(pE) && !is.null(pF)) {
  fig <- (pA | pB) / (pC | pD) / (pE | pF); h <- 8.0
} else if (!is.null(pE)) {
  fig <- (pA | pB) / (pC | pD) / pE; h <- 7.6
} else {
  fig <- (pA | pB) / (pC | pD); h <- 5.4
}
save_fig(fig, "Figure5_CCT2_classifier.pdf", height = h)
saveRDS(list(patients = pat, pair = sub, auc = auc, ci = ci, boot = boot,
             ttest = tt, wilcox = wt, pct_above_max = pct, dip_p = dip_p,
             protein = prot),
        "figure5_stats.rds")
cat("saved figure5_stats.rds\n")
