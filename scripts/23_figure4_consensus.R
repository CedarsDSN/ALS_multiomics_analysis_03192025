# =====================================================================
# FIGURE 4  Consensus clustering of sporadic ALS
#   a consensus matrix (k=3) + subtype ribbon   b PAC
#   c axis heatmap + subtype ribbon             d distributions
# Input:  backbone_FINAL.rds
# Output: Figure4_consensus_clustering.pdf, figure4_stats.rds
# =====================================================================
source("00_figure_common.R"); library(patchwork); fresh_device()

b  <- get_backbone(); d <- b$ax_round1; cc <- b$ccp_r1
AXO <- c("AXIS_METABOLIC", "AXIS_INFLAMMATION", "AXIS_CELL_CYCLE")

## ---- a: consensus matrix --------------------------------------------
cm <- cc[[3]]$consensusMatrix
dimnames(cm) <- list(d$sample_id, d$sample_id)
ord <- cc[[3]]$consensusTree$order; cm <- cm[ord, ord]
ids <- rownames(cm); n <- length(ids)
pA <- ggplot(data.frame(x = rep(1:n, times = n), y = rep(1:n, each = n),
                        v = as.vector(cm)), aes(x, y, fill = v)) +
  geom_raster() +
  scale_fill_gradient(low = "white", high = "#2166AC", name = "Consensus", limits = c(0, 1)) +
  coord_cartesian(expand = FALSE) +
  labs(title = "a   Consensus matrix, k = 3", x = NULL, y = NULL) +
  thm + theme(axis.text = element_blank(), axis.ticks = element_blank(),
              axis.line = element_blank(), plot.margin = margin(2, 2, 0, 2))
ribA <- ggplot(data.frame(x = 1:n, y = 1,
                          subtype = d$subtype[match(ids, d$sample_id)]),
               aes(x, y, fill = subtype)) +
  geom_raster() + scale_fill_manual(values = SCOL) +
  coord_cartesian(expand = FALSE) + theme_void(base_size = 8) +
  theme(legend.position = "none", plot.margin = margin(0, 2, 2, 2))

## ---- b: PAC ----------------------------------------------------------
pacdf <- data.frame(k = seq_along(b$pac_r1) + 1, PAC = as.numeric(b$pac_r1))
cat("PAC:", paste(sprintf("k%d=%.3f", pacdf$k, pacdf$PAC), collapse = "  "), "\n")
pB <- ggplot(pacdf, aes(k, PAC)) +
  geom_line(colour = "grey40") + geom_point(size = 2.2, colour = "grey20") +
  geom_point(data = pacdf[which.min(pacdf$PAC), ], size = 3.2, colour = "#E64B35") +
  geom_text(aes(label = sprintf("%.3f", PAC)), vjust = -1.1, size = 2.3) +
  scale_x_continuous(breaks = 2:5, expand = expansion(add = 0.35)) + 
  scale_y_continuous(limits = c(0, max(pacdf$PAC) * 1.35)) +
  labs(title = "b   Cluster stability", x = "Number of clusters (k)",
       y = "Proportion of ambiguous\nclustering (PAC)") + thm

## ---- c: axis heatmap -------------------------------------------------
dC <- d[order(d$subtype, d$AXIS_METABOLIC), ]
hm <- do.call(rbind, lapply(AXO, function(a)
  data.frame(axis = ALAB[[a]], x = seq_len(nrow(dC)),
             z = pmax(pmin(dC[[a]], 2.5), -2.5), stringsAsFactors = FALSE)))
hm$axis <- factor(hm$axis, levels = rev(unname(ALAB[AXO])))
pC <- ggplot(hm, aes(x, axis, fill = z)) +
  geom_raster() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-2.5, 2.5), name = "z-score") +
  coord_cartesian(expand = FALSE) +
  labs(title = "c   Axis scores by subtype", x = NULL, y = NULL) +
  thm + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.line = element_blank(), plot.margin = margin(2, 2, 0, 2))
tb   <- table(dC$subtype); bnd <- cumsum(as.integer(tb)); mids <- bnd - as.integer(tb) / 2
ribC <- ggplot(data.frame(x = seq_len(nrow(dC)), y = 1, subtype = dC$subtype),
               aes(x, y, fill = subtype)) +
  geom_raster() + scale_fill_manual(values = SCOL) +
  annotate("text", x = mids, y = 1,
           label = paste0(levels(dC$subtype), "\nn=", as.integer(tb)),
           size = 2.2, colour = "white", lineheight = 0.9) +
  coord_cartesian(expand = FALSE) +
  labs(x = sprintf("Sporadic ALS samples (n = %d)", nrow(dC))) +
  theme_void(base_size = 8) +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 8, margin = margin(t = 2)),
        plot.margin = margin(0, 2, 2, 2))

## ---- d: distributions -------------------------------------------------
long <- do.call(rbind, lapply(AXO, function(a)
  data.frame(axis = ALAB[[a]], subtype = d$subtype, z = d[[a]], stringsAsFactors = FALSE)))
long$axis <- factor(long$axis, levels = unname(ALAB[AXO]))
pD <- ggplot(long, aes(subtype, z, fill = subtype)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, width = 0.65) +
  facet_wrap(~ axis, nrow = 1) + scale_fill_manual(values = SCOL, name = "Subtype") +
  labs(title = "d   Axis score distributions", x = NULL, y = "Axis score (z)") +
  thm + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              legend.position = "right")

colA <- pA / ribA + plot_layout(heights = c(26, 1))
colC <- pC / ribC + plot_layout(heights = c(6, 1))
save_fig((colA | pB) / colC / pD + plot_layout(heights = c(3.1, 1.5, 1.7)),
         "Figure4_consensus_clustering.pdf", height = 8.2)

cat("\n--- legend numbers ---\n")
for (s in levels(d$subtype)) {
  m <- colMeans(d[d$subtype == s, AXO])
  cat(sprintf("%s: n = %d  Met %+.2f, Infl %+.2f, CC %+.2f\n",
              s, sum(d$subtype == s), m[1], m[2], m[3]))
}
saveRDS(list(pac = pacdf, sizes = tb), "figure4_stats.rds")
cat("saved figure4_stats.rds\n")
