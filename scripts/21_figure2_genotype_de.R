# =====================================================================
# FIGURE 2  Familial-sporadic convergence at gene and pathway level
#   a-c  pairwise log2FC scatters      d  three-level Spearman summary
#
# Output: Figure2.pdf, figure2_stats.rds
#
# PATH HANDLING: this script resolves its three inputs ITSELF, with
# explicit candidate paths and a validity check. It does not rely on
# p_() from 00_figure_common.R, because that loader has repeatedly
# resolved convergence_FINAL.rds to the CONTAMINATED directory.
#
# ~/Documents/convergence_results/          script-12, CONTAMINATED
#                                           RNA row 0.761 / -0.077 / 0.104
# ~/Documents/convergence_results_corrected/ CORRECTED
#                                           RNA row 0.402 / -0.070 / 0.045
#
# VALUES (order: C9-Sporadic, C9-SOD1, SOD1-Sporadic)
#   gene-level logFC   0.326 / 0.026 / 0.263
#   RNA pathway NES    0.402 / -0.070 / 0.045
#   protein NES        0.319 / -0.282 / 0.108
#
# 0.804 / +0.170 / +0.348 is RETIRED. It has no surviving artefact and
# no construction on disk reproduces it (see RNA_NES_resolution.md).
# =====================================================================

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
if (file.exists("00_figure_common.R")) source("00_figure_common.R")
while (!is.null(grDevices::dev.list())) grDevices::dev.off()

## ---- resolve inputs, explicitly --------------------------------------
first_existing <- function(cands, label) {
  cands <- path.expand(cands)
  hit <- cands[file.exists(cands)]
  if (!length(hit))
    stop(label, " not found. Tried:\n  ", paste(cands, collapse = "\n  "))
  hit[1]
}

DE_FILE <- first_existing(c(
  "~/Documents/convergence_results/genotype_de_corrected_FINAL.rds",
  "~/Documents/DSN/ALS/genotype_de_corrected_FINAL.rds",
  "genotype_de_corrected_FINAL.rds"), "genotype_de_corrected_FINAL.rds")

## CORRECTED ONLY. The contaminated folder is deliberately not a candidate.
CONV_FILE <- first_existing(c(
  "~/Documents/convergence_results_corrected/convergence_FINAL.rds",
  "~/Documents/DSN/ALS/convergence_results_corrected/convergence_FINAL.rds",
  "convergence_results_corrected/convergence_FINAL.rds"), "corrected convergence_FINAL.rds")

PROT_FILE <- first_existing(c(
  "~/Documents/DSN/ALS/proteomics/proteomics_inputs_FINAL.rds",
  "proteomics_inputs_FINAL.rds"), "proteomics_inputs_FINAL.rds")

if (grepl("convergence_results/", CONV_FILE, fixed = TRUE))
  stop("CONV_FILE resolved to the CONTAMINATED directory: ", CONV_FILE)
cat("DE  :", DE_FILE, "\nCONV:", CONV_FILE, "\nPROT:", PROT_FILE, "\n\n")

## ---- theme, if the common file was not available ---------------------
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

## ---- panels a-c ------------------------------------------------------
de <- readRDS(DE_FILE)
pick <- function(p) {
  i <- grep(p, names(de), ignore.case = TRUE)[1]
  if (is.na(i)) stop("no '", p, "' in the DE file; have: ",
                     paste(names(de), collapse = ", "))
  as.data.frame(de[[i]])
}
res_C9 <- pick("C9"); res_SOD1 <- pick("SOD1"); res_Sp <- pick("Spor")
cat("genes:", nrow(res_C9), "\n")

rho_lfc <- function(a, b) {
  sh <- intersect(rownames(a), rownames(b))
  suppressWarnings(cor(a[sh, "log2FoldChange"], b[sh, "log2FoldChange"],
                       method = "spearman", use = "complete.obs"))
}
rho <- c(rho_lfc(res_C9, res_Sp), rho_lfc(res_C9, res_SOD1), rho_lfc(res_SOD1, res_Sp))
cat("gene-level (expect 0.326 / 0.026 / 0.263):",
    paste(sprintf("%+.3f", rho), collapse = " / "), "\n")

mk <- function(a, b, xl, yl, r, ttl) {
  sh <- intersect(rownames(a), rownames(b))
  ggplot(data.frame(x = a[sh, "log2FoldChange"], y = b[sh, "log2FoldChange"]),
         aes(x, y)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(alpha = 0.25, size = 0.5, colour = "#0f6f6b") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "#173f3d", linewidth = 0.5) +
    labs(title = ttl,
         subtitle = sprintf("Spearman rho = %.3f (n = %d genes)", r, length(sh)),
         x = xl, y = yl) + thm
}
pA <- mk(res_C9,  res_Sp,   "C9orf72 log2FC", "Sporadic log2FC", rho[1], "a   C9orf72 vs sporadic")
pB <- mk(res_C9,  res_SOD1, "C9orf72 log2FC", "SOD1 log2FC",     rho[2], "b   C9orf72 vs SOD1")
pC <- mk(res_SOD1,res_Sp,   "SOD1 log2FC",    "Sporadic log2FC", rho[3], "c   SOD1 vs sporadic")

## ---- panel d ---------------------------------------------------------
CONTR <- c("C9orf72-Sporadic", "C9orf72-SOD1", "SOD1-Sporadic")

## RNA: pattern_table is what script 16 computes and prints. Do NOT
## recompute from cv$RNA_fgsea -- in some saved files those objects are
## stale leftovers from the pre-correction run.
cv <- readRDS(CONV_FILE)
pt <- as.data.frame(cv$pattern_table)
pt <- pt[pt$omic == "RNA", ]
gg <- function(p) pt$spearman[grep(p, pt$pair, ignore.case = TRUE)][1]
rna <- c(gg("C9 vs Spor"), gg("C9 vs SOD1"), gg("SOD1 vs Spor"))
cat("RNA pathway (expect 0.402 / -0.070 / 0.045):",
    paste(sprintf("%+.3f", rna), collapse = " / "), "\n")
if (any(is.na(rna))) stop("could not read the RNA row from pattern_table")
if (rna[1] > 0.6)
  stop("RNA row reads ", round(rna[1], 3), " -- the CONTAMINATED value. ",
       "CONV_FILE is: ", CONV_FILE)

## Protein: from the pinned fgsea tables, all 50 pathways, the same
## construction as the RNA row.
pin  <- readRDS(PROT_FILE)
nesv <- function(x) { x <- as.data.frame(x); setNames(x$NES, x$pathway) }
a <- nesv(pin$fg_c9); b <- nesv(pin$fg_sod1); d <- nesv(pin$fg_spor)
sh <- Reduce(intersect, list(names(a), names(b), names(d)))
prot <- c(cor(a[sh], d[sh], method = "spearman"),
          cor(a[sh], b[sh], method = "spearman"),
          cor(b[sh], d[sh], method = "spearman"))
cat("protein NES (expect 0.319 / -0.282 / 0.108):",
    paste(sprintf("%+.3f", prot), collapse = " / "),
    " on", length(sh), "pathways\n")

sum_df <- rbind(
  data.frame(contrast = CONTR, level = "Gene-level (logFC)",    rho = rho),
  data.frame(contrast = CONTR, level = "RNA pathway (NES)",     rho = rna),
  data.frame(contrast = CONTR, level = "Protein pathway (NES)", rho = prot))
sum_df$contrast <- factor(sum_df$contrast, levels = CONTR)
sum_df$level <- factor(sum_df$level,
                       levels = c("Gene-level (logFC)", "RNA pathway (NES)", "Protein pathway (NES)"))

## y range must include the negatives or the protein C9-SOD1 bar at
## -0.28 is clipped off the bottom.
ylim <- range(c(sum_df$rho, 0)) + c(-0.08, 0.10)

pD <- ggplot(sum_df, aes(contrast, rho, fill = level)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", rho)),
            position = position_dodge(width = 0.8),
            vjust = ifelse(sum_df$rho < 0, 1.4, -0.4), size = 2.0) +
  scale_fill_manual(values = c("Gene-level (logFC)"    = "#0f6f6b",
                               "RNA pathway (NES)"     = "#7a4fa3",
                               "Protein pathway (NES)" = "#c9822f"), name = NULL) +
  scale_y_continuous(limits = ylim) +
  labs(title = "d   Genotype similarity across levels",
       subtitle = "C9orf72 shares more with sporadic at the pathway level",
       x = NULL, y = "Spearman rho") +
  thm + theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 20, hjust = 1), plot.margin = margin(3, 12, 3, 3))

#ggsave("Figure2.pdf", (pA | pB) / (pC | pD), width = W, height = 5.0, units = "in")

#cat("\nwrote Figure2.pdf\n")
save_fig((pA | pB) / (pC | pD),
        "Figure2_new.pdf", height = 4.6)
saveRDS(list(gene = rho, rna = rna, protein = prot, summary = sum_df,
             sources = c(DE_FILE, CONV_FILE, PROT_FILE)), "figure2_stats.rds")
cat("saved figure2_stats.rds\n")