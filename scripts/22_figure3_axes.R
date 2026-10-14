# =====================================================================
# FIGURE 3  The three axes: activity, correlation, dimensionality
#   a axis activity by genotype   b pairwise correlation
#   c dimensionality (PCA)        d subtype separation per axis
# Input:  backbone_FINAL.rds, patient_level_FINAL.rds
# Output: Figure3.pdf, Supplementary_Table_axes.csv, figure3_stats.rds
#
# AUDIT NOTE: gene-set sizes were hardcoded (512/256/610) and the old
# caption claimed all cross-axis rho < 0.5. Both are wrong. Sizes now
# come from backbone$axis_sets; the corrected correlations are
# 0.24 / 0.61 / 0.60 and the subtitle states them as computed.
# =====================================================================
source("00_figure_common.R"); library(patchwork); fresh_device()

b <- get_backbone(); meta <- as.data.frame(b$meta); pat <- get_patients()

## ---- supplementary table: axis composition, sizes from the backbone --
SETS <- list(
  Metabolic    = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                   "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_GLYCOLYSIS"),
  Inflammation = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                   "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_COMPLEMENT"),
  `Cell cycle` = c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
                   "HALLMARK_MYC_TARGETS_V1", "HALLMARK_DNA_REPAIR"))
key <- c(Metabolic = "AXIS_METABOLIC", Inflammation = "AXIS_INFLAMMATION",
         `Cell cycle` = "AXIS_CELL_CYCLE")
if (is.null(b$axis_sets))
  stop("backbone_FINAL.rds has no axis_sets; the audit embedded it there")
supp <- do.call(rbind, lapply(names(SETS), function(a)
  data.frame(Axis = a, Hallmark_gene_set = SETS[[a]],
             Genes_in_axis = length(b$axis_sets[[key[[a]]]]),
             stringsAsFactors = FALSE)))
write.csv(supp, "Supplementary_Table_axes.csv", row.names = FALSE)
cat("axis sizes from backbone:\n")
print(unique(supp[, c("Axis", "Genes_in_axis")]), row.names = FALSE)

## ---- a: axis activity by genotype -----------------------------------
sc <- as.data.frame(scale(t(b$scores_r1))); sc$sample_id <- rownames(sc)
meta$group <- get_group(meta)
sc$group <- meta$group[match(sc$sample_id, meta$sample_id)]
sc <- sc[!is.na(sc$group), ]
cat("\nANOVA, axis score by genotype:\n")
aov_p <- sapply(AX, function(a)
  summary(aov(sc[[a]] ~ sc$group))[[1]]$`Pr(>F)`[1])
print(signif(aov_p, 3))

lg <- do.call(rbind, lapply(AX, function(a)
  data.frame(axis = ALAB[[a]], group = sc$group, z = sc[[a]], stringsAsFactors = FALSE)))
lg$axis <- factor(lg$axis, levels = unname(ALAB))
pA <- ggplot(lg, aes(group, z, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, width = 0.65) +
  facet_wrap(~ axis, nrow = 1) +
  scale_fill_manual(values = GCOL) +
  labs(title = "a   Axis activity by genotype",
       subtitle = "The axes describe variation within the cohort, not between genotypes",
       x = NULL, y = "Axis score (z)") +
  thm + theme(axis.text.x = element_text(angle = 40, hjust = 1),
              legend.position = "none")

## ---- b: correlation structure, computed ------------------------------
rho <- cor(pat[, AX], method = "spearman")
cat("\npairwise Spearman (corrected): met~inf", round(rho[1,2], 3),
    "| met~cc", round(rho[1,3], 3), "| inf~cc", round(rho[2,3], 3), "\n")
prs <- list(c(1,2), c(1,3), c(2,3))
scat <- do.call(rbind, lapply(prs, function(k)
  data.frame(pair = sprintf("%s vs %s   rho = %.2f", ALAB[[AX[k[1]]]],
                            ALAB[[AX[k[2]]]], rho[k[1], k[2]]),
             x = pat[[AX[k[1]]]], y = pat[[AX[k[2]]]],
             subtype = pat$subtype, stringsAsFactors = FALSE)))
pB <- ggplot(scat, aes(x, y, colour = subtype)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey20", linewidth = 0.4) +
  facet_wrap(~ pair, nrow = 1, scales = "free") +
  scale_colour_manual(values = SCOL, name = NULL) +
  labs(title = "b   Axis correlation structure",
       subtitle = sprintf("Correlated but not redundant (max rho = %.2f)",
                          max(rho[upper.tri(rho)])),
       x = "Axis score (z)", y = "Axis score (z)") +
  thm + theme(legend.position = "bottom")

## ---- c: dimensionality ----------------------------------------------
pc <- prcomp(pat[, AX], scale. = TRUE); ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
cat("variance explained:", paste0("PC", 1:3, " ", round(ve, 1), "%", collapse = ", "), "\n")
print(round(pc$rotation, 3))
pC <- ggplot(data.frame(PC = paste0("PC", 1:3), ve = ve), aes(PC, ve)) +
  geom_col(fill = "grey65", width = 0.6) +
  geom_text(aes(label = sprintf("%.0f%%", ve)), vjust = -0.5, size = 2.5) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "c   Dimensionality",
       subtitle = sprintf("PC1 explains %.0f%%, so the axes do not collapse to one", ve[1]),
       x = NULL, y = "Variance explained (%)") + thm

## ---- d: separation per axis ------------------------------------------
eta <- sapply(AX, function(a) {
  s <- summary(aov(pat[[a]] ~ pat$subtype))[[1]]$`Sum Sq`; s[1] / sum(s) })
cat("eta-squared:", paste(ALAB[AX], round(eta, 3), sep = "=", collapse = "  "), "\n")
ed <- data.frame(axis = unname(ALAB[AX]), eta = as.numeric(eta))
ed$axis <- factor(ed$axis, levels = ed$axis[order(ed$eta)])
pD <- ggplot(ed, aes(eta, axis)) +
  geom_col(fill = "grey65", width = 0.55) +
  geom_text(aes(label = sprintf("%.2f", eta)), hjust = -0.25, size = 2.5) +
  scale_x_continuous(limits = c(0, max(ed$eta) * 1.3), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "d   Subtype separation per axis",
       subtitle = "Variance attributable to subtype\n(descriptive; subtypes defined on these axes)",
       x = expression(eta^2), y = NULL) + thm

save_fig(pA / pB / (pC | pD) + plot_layout(heights = c(1, 1.25, 0.95)),
         "Figure3.pdf", height = 6.6)
saveRDS(list(rho = rho, pca = pc, var_explained = ve, eta = eta,
             aov_genotype_p = aov_p, axis_sizes = supp), "figure3_stats.rds")
cat("saved figure3_stats.rds\n")
