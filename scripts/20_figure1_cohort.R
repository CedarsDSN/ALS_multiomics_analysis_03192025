# =====================================================================
# FIGURE 1  Study design and cohort
#   a design schematic  b cohort sizes  c sex  d age at onset
# Input:  backbone_FINAL.rds, patient_level_FINAL.rds
# Output: Figure1.pdf, figure1_stats.rds
# All counts computed from the pinned objects; none hardcoded.
# =====================================================================
source("00_figure_common.R"); library(patchwork); fresh_device()

b    <- get_backbone(); meta <- as.data.frame(b$meta); pat <- get_patients()
meta$group <- get_group(meta)

ppl <- meta[!duplicated(meta$Participant_ID) & !is.na(meta$group), ]
cat("participant-level counts:\n"); print(table(ppl$group))
cat("samples:", nrow(meta), " sporadic patients (dedup):", nrow(pat), "\n")

n_all <- nrow(ppl); n_sp <- sum(ppl$group == "Sporadic"); n_samp <- nrow(meta)
n_sub <- table(pat$subtype)

## ---- a: design schematic (counts injected from data) ----------------
n_samp_sp <- sum(meta$group == "Sporadic", na.rm = TRUE)    # sporadic RNA-seq samples (expect 400)

steps <- data.frame(
  i   = 1:6,
  top = c("AnswerALS", "iPSC-derived", "Variance", "ssGSEA", "Consensus", "Single-gene"),
  mid = c(paste0(n_all, " participants"), "motor neurons", "stabilized counts",
          "10 Hallmark sets", "clustering, k = 3", "classifier"),
  bot = c(paste0(n_samp, " RNA-seq samples"), "bulk RNA-seq",
          "sex + SV residualized", "three axes",
          paste0(n_samp_sp, " sporadic samples\n", n_sp, " patients"), "CCT2"),
  stringsAsFactors = FALSE)
steps$x <- steps$i

pA <- ggplot(steps) +
  geom_segment(data = steps[steps$i < 6, ],
               aes(x = x + 0.44, xend = x + 0.56, y = 0, yend = 0),
               arrow = arrow(length = unit(1.3, "mm"), type = "closed"),
               linewidth = 0.3, colour = "grey45") +
  geom_tile(aes(x = x, y = 0), width = 0.86, height = 2.0,
            fill = "grey96", colour = "grey65", linewidth = 0.3) +
  geom_text(aes(x = x, y =  0.62, label = top), size = 2.2, fontface = "bold") +
  geom_text(aes(x = x, y =  0.20, label = mid), size = 2.0) +
  geom_text(aes(x = x, y = -0.22, label = bot), size = 2.0, colour = "grey35",
            vjust = 1, lineheight = 0.9) +
  coord_cartesian(xlim = c(0.5, 6.5), ylim = c(-1.1, 1.1)) +
  labs(title = "a   Study design") + void_thm

## ---- b: cohort ------------------------------------------------------
cs <- as.data.frame(table(ppl$group)); names(cs) <- c("group", "n")
pB <- ggplot(cs, aes(group, n, fill = group)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), vjust = -0.4, size = 2.4) +
  scale_fill_manual(values = GCOL) +
  scale_y_continuous(limits = c(0, max(cs$n) * 1.15),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "b   Cohort", x = NULL, y = "Participants") +
  thm + theme(legend.position = "none",
              axis.text.x = element_text(angle = 40, hjust = 1))

## ---- c: sex ---------------------------------------------------------
sx <- ppl[!is.na(ppl$SEX) & ppl$SEX != "", ]; sx$SEX <- factor(as.character(sx$SEX))
sdf <- as.data.frame(prop.table(table(sx$group, sx$SEX), margin = 1))
names(sdf) <- c("group", "sex", "prop")
sexcol <- setNames(c("#B8860B", "#31688E")[seq_along(levels(sx$SEX))], levels(sx$SEX))
pC <- ggplot(sdf, aes(group, prop, fill = sex)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = sexcol, name = NULL) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "c   Sex", x = NULL, y = "Proportion") +
  thm + theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 40, hjust = 1))

## ---- d: age at onset ------------------------------------------------
ao <- ppl[ppl$group != "Control" & !is.na(ppl$AGE_AT_SYMPTOM_ONSET), ]
ao$group <- droplevels(ao$group)
kw  <- kruskal.test(AGE_AT_SYMPTOM_ONSET ~ group, data = ao)
nn  <- table(ao$group)
cat("\nmedian age at onset:\n")
print(round(tapply(ao$AGE_AT_SYMPTOM_ONSET, ao$group, median, na.rm = TRUE), 1))
cat("Kruskal-Wallis p =", signif(kw$p.value, 4), "\n")

pD <- ggplot(ao, aes(group, AGE_AT_SYMPTOM_ONSET, fill = group)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.35, colour = "grey25") +
  scale_fill_manual(values = GCOL[levels(ao$group)]) +
  scale_x_discrete(labels = paste0(names(nn), "\nn=", as.integer(nn))) +
  labs(title = "d   Age at symptom onset",
       subtitle = sprintf("Kruskal-Wallis p = %.3f", kw$p.value),
       x = NULL, y = "Years") +
  thm + theme(legend.position = "none")

save_fig(pA / (pB | pC | pD) + plot_layout(heights = c(1.25, 2.1)),
         "Figure1_new.pdf", height = 4.6)
saveRDS(list(counts = cs, sex = sdf, onset_kw = kw, n_subtype = n_sub),
        "figure1_stats.rds")
cat("saved figure1_stats.rds\n")
