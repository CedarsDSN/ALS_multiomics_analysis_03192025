# =====================================================================
# FIGURE 7  Subtypes are post-genomic states
#   a  CCT2 and TRiC rare damaging variant carriage by subtype
#   b  CCT2 cis-eQTL, quantile-quantile
#   c  Genome-wide rare damaging burden by subtype
#   d  Genomic clustering vs transcriptomic subtype
#
# Output: Figure7_genomic.pdf, figure7_stats.rds
#
# WHAT WAS WRONG WITH THE PREVIOUS VERSION
#   1. Panel b read a stale cct2_cis_eqtl_results.csv and rendered
#      "3781 cis variants, n = 369, minimum FDR q = 0.82". Those are the
#      PSEUDOREPLICATED numbers: nine patients with more than one RNA
#      sample entered the regressions twice with identical genotypes. The
#      corrected run is 3,782 variants, n = 360, minimum q = 0.712. This
#      version STOPS if it sees the old file rather than plotting it.
#   2. Panel c reported only Kruskal-Wallis (P = 0.05) while the text
#      reported ANOVA (P = 0.86). Both are now computed and both are
#      shown, because a panel reading P = 0.05 beside a claim of no
#      difference is exactly what a reviewer will question.
#
# DROPPED, and not to be reinstated without a rebuild: the pathway and
# genome-wide polygenic risk panels. The van Rheenen summary statistics
# used to build those scores were on GRCh37 while the AnswerALS genotypes
# are GRCh38, so under 2% of chromosome 12 variants matched and the
# scores were assembled from coincidental positional collisions.
# =====================================================================

suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr) })
if (file.exists("00_figure_common.R")) source("00_figure_common.R")
while (!is.null(grDevices::dev.list())) grDevices::dev.off()

find1 <- function(f, dirs = c(".", "~/Documents/DSN/ALS", "~/Documents",
                              "~/Documents/DSN/ALS/cct2_eqtl",
                              "~/Documents/DSN/ALS/tric_burden_results")) {
  for (d in path.expand(dirs)) {
    p <- file.path(d, f); if (file.exists(p)) return(p)
    if (dir.exists(d)) {
      g <- list.files(d, pattern = paste0("^", f, "$"), recursive = TRUE, full.names = TRUE)
      if (length(g)) return(g[1])
    }
  }
  stop("not found: ", f)
}
if (!exists("SCOL"))
  SCOL <- c("S-ALS-1" = "#E64B35", "S-ALS-2" = "#4DBBD5", "S-ALS-3" = "#00A087")
if (!exists("thm")) {
  thm <- theme_classic(base_size = 8) +
    theme(axis.text = element_text(size = 7, colour = "black"),
          axis.title = element_text(size = 8), legend.text = element_text(size = 7),
          legend.key.size = unit(3, "mm"),
          plot.title = element_text(size = 9, hjust = 0),
          plot.subtitle = element_text(size = 7), plot.margin = margin(3, 3, 3, 3))
}
if (!exists("W")) W <- 180 / 25.4

gb <- readRDS(find1("genomic_burden_FINAL.rds"))
cat("genomic landmarks:\n"); str(gb$landmarks)

## ---- a: CCT2 and TRiC carriage ----------------------------------------
tr <- utils::read.csv(find1("tric_burden_full.csv"), stringsAsFactors = FALSE)
sc <- grep("cluster_clean|subtype", names(tr), value = TRUE)[1]
tr$subtype <- factor(as.character(tr[[sc]]), levels = names(SCOL))
nsub <- table(tr$subtype)
cat("\nTRiC table:", nrow(tr), "patients\n"); print(nsub)
cc_car   <- tapply(tr$CCT2_carrier, tr$subtype, sum)
tric_pct <- tapply(tr$TRiC_carrier, tr$subtype, function(x) 100 * mean(x))
cat("CCT2 carriers by subtype:", paste(cc_car, collapse = "/"),
    " (must be 0/0/0)\n")

ad <- data.frame(
  subtype = rep(names(nsub), 2),
  gene    = rep(c("CCT2", "TRiC complex"), each = length(nsub)),
  pct     = c(100 * as.numeric(cc_car) / as.numeric(nsub), as.numeric(tric_pct)))
ad$subtype <- factor(ad$subtype, levels = names(SCOL))
ad$gene    <- factor(ad$gene, levels = c("CCT2", "TRiC complex"))

pA <- ggplot(ad, aes(subtype, pct, fill = subtype)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 2.2) +
  facet_wrap(~ gene, nrow = 1) +
  scale_fill_manual(values = SCOL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(title = "a   Rare damaging variant carriage",
       subtitle = sprintf("No CCT2 carrier in any of %d patients", nrow(tr)),
       x = NULL, y = "Patients carrying (%)") +
  thm + theme(legend.position = "none",
              axis.text.x = element_text(angle = 25, hjust = 1))

## ---- b: eQTL QQ. Refuse the stale file. --------------------------------
## Prefer the pinned _FINAL file. The plain cct2_cis_eqtl_results.csv in
## that directory is the April 2026 PRE-CORRECTION run (n = 369, q = 0.82)
## and must not be used: nine patients with more than one RNA sample
## entered those regressions twice with identical genotypes.
.eqtl_file <- tryCatch(find1("cct2_cis_eqtl_results_FINAL.csv"),
                       error = function(e) find1("cct2_cis_eqtl_results.csv"))
cat("eQTL file:", .eqtl_file, "\n")
eq <- utils::read.csv(.eqtl_file, check.names = FALSE)
n_max <- max(eq$n, na.rm = TRUE); q_min <- min(eq$fdr_q, na.rm = TRUE)
cat(sprintf("\neQTL: %d variants, max n = %d, min q = %.3f\n",
            nrow(eq), n_max, q_min))
if (n_max > 361)
  stop("this eQTL file holds the PSEUDOREPLICATED run (max n = ", n_max,
       ", expected 360): ", .eqtl_file,
       "\n  Use cct2_cis_eqtl_results_FINAL.csv, or rerun ",
       "12_cct2_eqtl_AUDITED.R.")

p  <- sort(eq$p[!is.na(eq$p) & eq$p > 0]); nvar <- length(p)
pB <- ggplot(data.frame(exp = -log10(ppoints(nvar)), obs = -log10(p)),
             aes(exp, obs)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70",
              linetype = "dashed", linewidth = 0.35) +
  geom_point(size = 0.5, colour = "#2166AC", alpha = 0.6) +
  labs(title = "b   CCT2 cis-eQTL",
       subtitle = sprintf("%s cis variants, n = %d, minimum FDR q = %.2f",
                          format(length(p), big.mark = ","), max(eq$n, na.rm = TRUE),
                          min(eq$fdr_q, na.rm = TRUE)),
       x = expression(Expected~-log[10](p)),
       y = expression(Observed~-log[10](p))) + thm

## ---- c: rare damaging burden. Report BOTH tests. -----------------------
rb <- as.data.frame(gb$rare_burden_sub)
sc2 <- grep("cluster_clean|subtype", names(rb), value = TRUE)[1]
rb$subtype <- factor(as.character(rb[[sc2]]), levels = names(SCOL))
bcol <- grep("^rare_damaging$", names(rb), value = TRUE)[1]

aov_p <- summary(aov(rb[[bcol]] ~ rb$subtype))[[1]]$`Pr(>F)`[1]
kw_p  <- kruskal.test(rb[[bcol]] ~ rb$subtype)$p.value
cat(sprintf("\nrare damaging burden: ANOVA P = %.3f | Kruskal-Wallis P = %.3f\n",
            aov_p, kw_p))
cat("Both are shown on the panel. The means are nearly identical\n",
    "(", paste(round(tapply(rb[[bcol]], rb$subtype, mean)), collapse = ", "),
    ") so the rank test is picking up a distributional difference, not a\n",
    "shift in burden. Report both in the text.\n", sep = "")

nn2 <- table(rb$subtype)
kw <- kruskal.test(rb[[bcol]] ~ rb$subtype)
av <- summary(aov(rb[[bcol]] ~ rb$subtype))[[1]]$`Pr(>F)`[1]
pC <- ggplot(rb, aes(subtype, .data[[bcol]], fill = subtype)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.6) +
  geom_jitter(width = 0.16, size = 0.35, alpha = 0.35, colour = "grey25") +
  scale_fill_manual(values = SCOL) +
  scale_x_discrete(labels = paste0(names(nn2), "\nn=", as.integer(nn2))) +
  labs(title = "c   Genome-wide rare damaging burden",
       subtitle = sprintf("ANOVA P = %.2f; Kruskal-Wallis P = %.2f", av, kw$p.value),
       x = NULL, y = "Rare damaging variants") +
  thm + theme(legend.position = "none")

## ---- d: clustering agreement -------------------------------------------
ad2 <- data.frame(
  what = factor(c("k-means on genomic PCs", "blind consensus, k = 2"),
                levels = c("blind consensus, k = 2", "k-means on genomic PCs")),
  ARI  = c(gb$landmarks$ari_kmeans, gb$landmarks$ari_blind_k2))
pD <- ggplot(ad2, aes(ARI, what)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.3) +
  geom_col(fill = "grey65", width = 0.5) +
  geom_text(aes(label = sprintf("ARI = %+.3f", ARI)),
            hjust = ifelse(ad2$ARI < 0, 1.1, -0.1), size = 2.4) +
  scale_x_continuous(limits = c(-0.15, 0.15)) +
  labs(title = "d   Genomic vs transcriptomic partitions",
       subtitle = "No agreement; the blind k = 2 split isolates a single outlier",
       x = "Adjusted Rand index", y = NULL) + thm


save_fig((pA | pB) / (pC | pD), 
         "Figure7_genomic.pdf", height = 4.6)
saveRDS(list(carriage = ad, eqtl_file = .eqtl_file, eqtl_n = n_max, eqtl_q = q_min, eqtl_nvar = nvar,
             burden = rb, anova_p = aov_p, kw_p = kw_p,
             landmarks = gb$landmarks), "figure7_stats.rds")
cat("saved figure7_stats.rds\n")
