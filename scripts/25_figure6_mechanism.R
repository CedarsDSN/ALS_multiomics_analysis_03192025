# =====================================================================
# FIGURE 6  Mechanism: mTORC1 to LARP1 to translation, with MYC
#   a upstream regulators   b canonical pathways
#   c the mTOR/MYC/translation axis   d schematic, z-scores from IPA
# Input:  the two IPA export CSVs (noSexChr analysis)
# Output: Figure6_mechanism.pdf, Figure6_panelD_mechanism.pdf, figure6_stats.rds
#
# AUDIT NOTE: the old scripts used two different filename patterns for
# the same export (space vs underscore) and zof() returned NA silently,
# so a missed match produced a blank schematic instead of an error.
# Both are now resolved by one tolerant matcher that STOPS if nothing
# is found, and the translation z-scores (previously hardcoded 6.56 /
# 6.40) are read from the export.
# =====================================================================
source("00_figure_common.R"); library(patchwork); fresh_device()

find_ipa <- function(kind) {
  all_csv <- unlist(lapply(PIN_DIRS[dir.exists(PIN_DIRS)], function(d)
    list.files(d, pattern = paste0(kind, ".*\\.csv$"),
               ignore.case = TRUE, full.names = TRUE, recursive = TRUE)))
  f <- grep("nosexchr|noXY|no_sex", all_csv, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(f)) {
    f <- all_csv[1]
    if (!is.na(f))
      warning("using ", basename(f), " -- NOT the noSexChr export; ",
              "confirm this is the intended IPA analysis")
  }
  if (is.na(f)) stop("no IPA export matching '", kind, "' found in PIN_DIRS")
  f
}
UR <- find_ipa("Upstream.?Regulators"); CP <- find_ipa("Canonical.?Pathways")
cat("upstream: ", basename(UR), "\ncanonical:", basename(CP), "\n")

rd <- function(f) {
  x <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  names(x)[1:3] <- c("name", "p", "z")
  x$p <- as.numeric(x$p); x$z <- suppressWarnings(as.numeric(x$z))
  x$logp <- -log10(x$p); x
}
ur <- rd(UR); cp <- rd(CP)
cat("regulators:", nrow(ur), " pathways:", nrow(cp), "\n")

zfill <- scale_fill_gradient2(low = "#2166AC", mid = "grey92", high = "#D6604D",
                              midpoint = 0, limits = c(-6.6, 6.6),
                              name = "Activation\nz-score", na.value = "grey80")

a <- head(ur[order(ur$p), ], 15); a$name <- factor(a$name, levels = rev(a$name))
pA <- ggplot(a, aes(logp, name, fill = z)) +
  geom_col(width = 0.7, colour = "grey40", linewidth = 0.2) + zfill +
  labs(title = "a   Upstream regulators", subtitle = "Top 15 by P value",
       x = expression(-log[10](p)), y = NULL) +
  thm + theme(legend.position = "none")

bq <- head(cp[order(cp$p), ], 12)
bq$name <- gsub(" \\(.*\\)$", "", bq$name)
bq$name <- factor(bq$name, levels = rev(bq$name))
pB <- ggplot(bq, aes(logp, name, fill = z)) +
  geom_col(width = 0.7, colour = "grey40", linewidth = 0.2) + zfill +
  scale_y_discrete(labels = function(x) vapply(strwrap(x, 34, simplify = FALSE),
                                               paste, "", collapse = "\n")) +
  labs(title = "b   Canonical pathways", subtitle = "Translation dominates",
       x = expression(-log[10](p)), y = NULL) + thm

want <- c("RPTOR","MTOR","RICTOR","torin1","sirolimus","LARP1","MYC","MYCN","EIF4E","CIP2A")
lab  <- c(RPTOR = "RPTOR (mTORC1)", MTOR = "MTOR", RICTOR = "RICTOR (mTORC2)",
          torin1 = "torin1 (mTOR inh.)", sirolimus = "sirolimus (mTOR inh.)",
          LARP1 = "LARP1", MYC = "MYC", MYCN = "MYCN", EIF4E = "EIF4E",
          CIP2A = "CIP2A (MYC stabilizer)")
cc <- ur[ur$name %in% want & !is.na(ur$z), ]
cat("\nmechanism axis members found:\n")
print(data.frame(name = cc$name, p = signif(cc$p, 3), z = round(cc$z, 2)), row.names = FALSE)
missing <- setdiff(want, ur$name)
if (length(missing)) cat("NOT in this export:", paste(missing, collapse = ", "), "\n")

pC <- NULL
if (nrow(cc)) {
  cc$lab <- lab[cc$name]; cc <- cc[order(cc$z), ]
  cc$lab <- factor(cc$lab, levels = cc$lab)
  pC <- ggplot(cc, aes(z, lab, fill = z)) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.3) +
    geom_col(width = 0.65, colour = "grey40", linewidth = 0.2) +
    geom_text(aes(label = sprintf("P = %.0e", p), hjust = ifelse(z < 0, 1.1, -0.1)), size = 2.1) +
    zfill + scale_x_continuous(expand = expansion(mult = 0.25)) +
    labs(title = "c   mTOR, MYC and translation control",
         subtitle = "mTORC1 active, LARP1 inhibited, MYC activated",
         x = "Activation z-score", y = NULL) +
    thm + theme(legend.position = "none")
}

save_fig(if (is.null(pC)) (pA | pB) else (pA | pB) / pC + plot_layout(heights = c(1.5, 1)),
         "Figure6_mechanism.pdf", height = if (is.null(pC)) 3.4 else 6.0)

## ---- panel d: schematic, every z read from the export ----------------
zof <- function(n) { i <- match(n, ur$name); if (is.na(i)) NA_real_ else ur$z[i] }
z_transl <- { i <- grep("Translation|Elongation", cp$name, ignore.case = TRUE)
              if (length(i)) cp$z[i[which.max(cp$logp[i])]] else NA_real_ }
z_gcn2   <- zof("EIF2AK4")
if (is.na(z_gcn2)) z_gcn2 <- zof("GCN2")
cat(sprintf("\nschematic z: RPTOR %.2f | LARP1 %.2f | MYC %.2f | translation %.2f | GCN2 %.2f\n",
            zof("RPTOR"), zof("LARP1"), zof("MYC"), z_transl, z_gcn2))

nodes <- data.frame(
  id = c("mtorc1","larp1","top","myc","transl","gcn2","cct2","state"),
  label = c("mTORC1 active", "LARP1 de-repressed", "5'TOP mRNA release",
            "MYC / MYCN activated", "Translation elongation\nand termination up",
            "GCN2 amino-acid\nstress response", "CCT2 / TRiC\nchaperonin demand",
            "S-ALS-3\nfaster ALSFRS-R decline"),
  z = c(zof("RPTOR"), zof("LARP1"), NA, zof("MYC"), z_transl, z_gcn2, NA, NA),
  y = 8:1, stringsAsFactors = FALSE)
nodes$note  <- ifelse(is.na(nodes$z), "", sprintf("z = %+.2f", nodes$z))
nodes$label <- factor(nodes$label, levels = rev(nodes$label))

p_mech <- ggplot(nodes) +
  geom_segment(data = nodes[nodes$y > 1, ],
               aes(x = 0, xend = 0, y = y - 0.32, yend = y - 0.68),
               arrow = arrow(length = unit(1.6, "mm"), type = "closed"),
               linewidth = 0.35, colour = "grey45") +
  geom_tile(aes(x = 0, y = y, fill = z), width = 1.7, height = 0.62,
            colour = "grey55", linewidth = 0.3) +
  geom_text(aes(x = 0, y = y, label = label), size = 2.4, lineheight = 0.95) +
  geom_text(aes(x = 0.92, y = y, label = note), size = 2.1, hjust = 0, colour = "grey30") +
  scale_fill_gradient2(low = "#2166AC", mid = "grey94", high = "#D6604D",
                       midpoint = 0, limits = c(-6.6, 6.6),
                       na.value = "grey94", guide = "none") +
  coord_cartesian(xlim = c(-1.0, 1.9), ylim = c(0.5, 8.6)) +
  labs(title = "d   Proposed mechanism",
       subtitle = "Shading and z-scores from IPA activation predictions") +
  void_thm
save_fig(p_mech, "Figure6_panelD_mechanism.pdf", height = 120/25.4, width = 62/25.4)

saveRDS(list(upstream = ur, canonical = cp, axis_members = cc,
             ipa_files = c(UR, CP)), "figure6_stats.rds")
cat("saved figure6_stats.rds\n")
