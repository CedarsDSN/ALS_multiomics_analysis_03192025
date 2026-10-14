# =====================================================================
# FIGURE 8  Clinical screen on the final labels
#   a subtype effect across all endpoints (broad null)
#   b CCT2 effect across the same endpoints (decline-rate hits)
#   c the ALSFRS-R result under alternative specifications
#
# Input:  clinical_screen_combined_FINAL.rds, patient_level_FINAL.rds
#         downstream_analyses_FINAL.rds (for the survival rows, if present)
# Output: Figure8_clinical_screen.pdf, figure8_stats.rds
#
# AUDIT NOTE: the old panel c read res_SURV / res_SURV2 from the SESSION
# and silently skipped if absent. Survival now comes from the pinned
# downstream object; if it is not there, panel c is replaced by the
# robustness panel rather than dropped without trace.
# =====================================================================
source("00_figure_common.R")
suppressPackageStartupMessages({ library(data.table); library(patchwork) })
fresh_device()

sc <- as.data.table(load_pin("clinical_screen_combined_FINAL.rds"))
sc <- sc[!is.na(p_S3)]
sc <- unique(sc, by = "endpoint")
cat("endpoints plotted:", nrow(sc), "\n")
sc[, lab := gsub(" \\(.*\\)$", "", endpoint)][, lab := gsub("PORTAL", "portal", lab)]

sa <- copy(sc)[order(p_S3)]
sa[, lab := factor(lab, levels = rev(lab))][, hit := fifelse(p_S3 < 0.05, "p < 0.05", "ns")]
pA <- ggplot(sa, aes(-log10(p_S3), lab, colour = hit)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = -log10(p_S3), yend = lab), colour = "grey80", linewidth = 0.3) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = c("p < 0.05" = "#E64B35", "ns" = "grey55"), name = NULL) +
  labs(title = "a   Subtype effect (S-ALS-3 vs S-ALS-2)",
       subtitle = sprintf("%d endpoints, none survives FDR (min q = %.2f)",
                          nrow(sa), min(sa$p_S3_FDR_all, na.rm = TRUE)),
       x = expression(-log[10](p)), y = NULL) +
  thm + theme(legend.position = "none", axis.text.y = element_text(size = 6))

sb <- copy(sc)[!is.na(p_CCT2_uni)][order(p_CCT2_uni)]
sb[, lab := factor(lab, levels = rev(lab))][, hit := fifelse(p_CCT2_uni < 0.05, "p < 0.05", "ns")]
pB <- ggplot(sb, aes(-log10(p_CCT2_uni), lab, colour = hit)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = -log10(p_CCT2_uni), yend = lab), colour = "grey80", linewidth = 0.3) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = c("p < 0.05" = "#00A087", "ns" = "grey55"), name = NULL) +
  labs(title = "b   CCT2 effect (continuous, per SD)",
       subtitle = sprintf("%d nominal hits, all decline rates; none survive FDR",
                          sum(sb$p_CCT2_uni < 0.05, na.rm = TRUE)),
       x = expression(-log[10](p)), y = NULL) +
  thm + theme(legend.position = "none", axis.text.y = element_text(size = 6))

cat("\nCCT2 nominal hits:\n")
print(sb[p_CCT2_uni < 0.05, .(endpoint, n, beta = round(beta_CCT2_uni, 3),
                              p = signif(p_CCT2_uni, 3))])

## ---- c: ALSFRS-R robustness, canonical covariates (age, sex, 6 SVs) ---
ds  <- load_pin("downstream_analyses_FINAL.rds")
d   <- as.data.frame(ds$df); COV <- ds$covariates
d$CCT2_sd <- as.numeric(scale(d$CCT2_expr))
d$subtype <- relevel(factor(d$subtype), ref = "S-ALS-2")
yv <- "ALSFRS_R_PROGRESSION_SLOPE"
C3 <- "+ prolif_score + maturity_score + diff_quality"
C6 <- paste(C3, "+ astro_score + micro_score + oligo_score")
rob1 <- function(data, lab, extra = "") {
  m  <- lm(as.formula(paste(yv, "~ CCT2_sd +", COV, extra)), data = data)
  ci <- confint(m)["CCT2_sd", ]
  data.frame(spec = lab, beta = coef(m)[["CCT2_sd"]], lo = ci[[1]], hi = ci[[2]],
             p = coef(summary(m))["CCT2_sd", 4], n = nobs(m),
             comp = grepl("composition", lab), stringsAsFactors = FALSE)
}
y <- d[[yv]]; q <- quantile(y, c(.01, .99), na.rm = TRUE)
rb <- rbind(
  rob1(d, "all patients"),
  rob1(d[is.na(y) | y >= -2, ], "excluding slope < -2"),
  rob1(d[is.na(y) | (y >= q[1] & y <= q[2]), ], "trimmed 1st/99th pct"),
  rob1(d, "+ subtype", "+ subtype"),
  rob1(d, "+ subtype + composition", paste("+ subtype", C3)),
  rob1(d, "+ subtype + composition + glia", paste("+ subtype", C6)))
cat("\nALSFRS-R slope, CCT2 effect by specification:\n")
print(transform(rb, beta = round(beta, 4), p = signif(p, 3)), row.names = FALSE)
stopifnot(rb$n[1] == 244, abs(rb$beta[1] + 0.1374) < 0.001, abs(rb$beta[4] + 0.1218) < 0.001)

rb$spec <- factor(rb$spec, levels = rev(rb$spec))
pC <- ggplot(rb, aes(beta, spec, colour = comp)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.35) +
  geom_point(size = 2) +
  geom_text(aes(label = sprintf("P = %.3f, n = %d", p, n)), vjust = -1.0, size = 2.1,
            colour = "grey20") +
  scale_colour_manual(values = c(`FALSE` = "grey20", `TRUE` = "grey60"), guide = "none") +
  labs(title = "c   ALSFRS-R decline: CCT2 effect by model specification",
       subtitle = "All models adjust for age at onset, sex and six surrogate variables",
       x = "Beta per SD of CCT2 (ALSFRS-R points/month)", y = NULL) + thm

save_fig((pA | pB) / pC + plot_layout(heights = c(2.4, 1.4)),
         "Figure8_clinical_screen.pdf", height = 7.0)
saveRDS(list(screen = sc, cct2_hits = sb[p_CCT2_uni < 0.05],
             robustness = if (exists("rb")) rb else NULL), "figure8_stats.rds")
cat("saved figure8_stats.rds\n")
