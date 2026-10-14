# =====================================================================
# 31_pathway_endpoint_screen.R   (corrected mohamadFeedback1.R)
#
# Run 30_pathway_preamble.R first.
#
#   PHASE 1  NfL added to the ALSFRS-R model
#   PHASE 2  Pathway x clinical-endpoint association matrix (10 Hallmarks)
#   PHASE 6  Familial proximity of the sporadic subtypes
#   PHASE 7  Top-500 directional agreement on a FIXED gene set
#
# CHANGES FROM mohamadFeedback1.R
#   - Phases 3 and 5 (GO-BP drill-down, rrvgo) REMOVED: they need
#     res_clean_23, the S-ALS-3 vs S-ALS-2 DE object, which is not in
#     any pin. They belong with 32_ anyway, which does the GO work.
#   - Phase 4 REMOVED as a duplicate of 32_ Step 1.
#   - The stray console block after Phase 3 ("run this, not the above")
#     is GONE. It referenced r, sig2, universe_ens and m2e(), none of
#     which exist in that script, so it always errored.
#   - Phase 1's filter(group == "Sporadic") now works: the preamble sets
#     that column. The frame is sporadic-only by construction anyway.
#   - Phase 2 scores ssGSEA on expr_resid, NOT expr_symbol. The axes were
#     built on residualised expression; scoring the pathways on raw
#     expression compared two different things.
#   - Phase 2 runs at PATIENT level. The original joined sample-level
#     scores onto a patient-level clinical frame.
#   - Phase 6 works: the preamble builds axis_merged from scores_r1.
#
# EVERY MODEL HERE IS AN ASSOCIATION. Write "associated with", never
# "drives".
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(broom); library(ggplot2)
})
set.seed(42)

need <- c("expr_resid_mat","hallmark_list","axis_sets","axis_vars",
          "ax_clean_dedup","clin_frame","universe","sv_vars","endpoints")
miss <- need[!sapply(need, exists)]
if (length(miss)) stop("run 30_pathway_preamble.R first. Missing: ",
                       paste(miss, collapse = ", "))
cat("SV covariates:", paste(sv_vars, collapse = ", "), "\n")
cov_str <- paste(c("AGE_AT_SYMPTOM_ONSET", "SEX", sv_vars), collapse = " + ")

out_csv <- function(x, f) { write_csv(x, f); cat("wrote", f, "\n") }

# =====================================================================
# PHASE 1  NfL added to the ALSFRS-R model
# Does the metabolic axis, or the S-ALS-3 label, predict decline rate
# beyond the axonal injury NfL already captures?
# =====================================================================
cat("\n\n=========== PHASE 1: NfL in the ALSFRS-R model ===========\n")

p1 <- clin_frame %>%
  filter(group == "Sporadic",
         !is.na(ALSFRS_R_PROGRESSION_SLOPE), !is.na(logNfL)) %>%
  droplevels()
cat("n (slope & NfL):", nrow(p1), "\n")

fA <- lm(as.formula(paste("ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC +", cov_str)), data = p1)
fB <- lm(as.formula(paste("ALSFRS_R_PROGRESSION_SLOPE ~ AXIS_METABOLIC + logNfL +", cov_str)), data = p1)
cat("\n-- metabolic axis, no NfL --\n"); print(tidy(fA) %>% filter(term == "AXIS_METABOLIC"))
cat("\n-- metabolic axis + NfL --\n");   print(tidy(fB) %>% filter(term %in% c("AXIS_METABOLIC","logNfL")))
cat(sprintf("\nadj R2  A = %.3f   B = %.3f\n",
            summary(fA)$adj.r.squared, summary(fB)$adj.r.squared))
cat("ANOVA A vs B:\n"); print(anova(fA, fB))

## the subtype version of the same question
p1b <- clin_frame %>%
  filter(cluster_final %in% c("S-ALS-2","S-ALS-3"),
         !is.na(ALSFRS_R_PROGRESSION_SLOPE), !is.na(logNfL)) %>%
  droplevels() %>%
  mutate(is_S3 = as.integer(cluster_final == "S-ALS-3"))
cat("\nn (S-ALS-2 vs S-ALS-3, slope & NfL):", nrow(p1b), "\n")
f0 <- lm(as.formula(paste("ALSFRS_R_PROGRESSION_SLOPE ~ is_S3 +", cov_str)), data = p1b)
f1 <- lm(as.formula(paste("ALSFRS_R_PROGRESSION_SLOPE ~ is_S3 + logNfL +", cov_str)), data = p1b)
cat("\n-- S-ALS-3, no NfL --\n"); print(tidy(f0) %>% filter(term == "is_S3"))
cat("\n-- S-ALS-3 + NfL --\n");   print(tidy(f1) %>% filter(term %in% c("is_S3","logNfL")))
cat(sprintf("\nadj R2  without NfL = %.3f, with NfL = %.3f\n",
            summary(f0)$adj.r.squared, summary(f1)$adj.r.squared))

## CCT2 version, if the column is present
if ("CCT2z" %in% names(clin_frame)) {
  p1c <- clin_frame %>% filter(!is.na(ALSFRS_R_PROGRESSION_SLOPE), !is.na(logNfL))
  fc0 <- lm(as.formula(paste("ALSFRS_R_PROGRESSION_SLOPE ~ CCT2z +", cov_str)), data = p1c)
  fc1 <- lm(as.formula(paste("ALSFRS_R_PROGRESSION_SLOPE ~ CCT2z + logNfL +", cov_str)), data = p1c)
  cat("\n-- CCT2 without / with NfL (n =", nrow(p1c), ") --\n")
  print(bind_rows(tidy(fc0) %>% filter(term == "CCT2z") %>% mutate(model = "no NfL"),
                  tidy(fc1) %>% filter(term == "CCT2z") %>% mutate(model = "with NfL")))
}

out_csv(bind_rows(tidy(fA) %>% mutate(model = "A_axis_only"),
                  tidy(fB) %>% mutate(model = "B_axis_plus_NfL"),
                  tidy(f0) %>% mutate(model = "C_S3_only"),
                  tidy(f1) %>% mutate(model = "D_S3_plus_NfL")),
        "phase1_alsfrs_nfl_models.csv")

cat("\nINTERPRET: if the predictor survives with NfL in the model it\n")
cat("carries progression information beyond axonal damage. If it\n")
cat("attenuates, part of its signal runs through NfL. Report the n;\n")
cat("a power drop is not a true null.\n")

# =====================================================================
# PHASE 2  Pathway x clinical-endpoint association matrix
# The 10 Hallmark pathways behind the three axes, tested individually.
# =====================================================================
cat("\n\n=========== PHASE 2: pathway x endpoint matrix ===========\n")
suppressPackageStartupMessages(library(GSVA))

pathways_10 <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_GLYCOLYSIS",
  "HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MYC_TARGETS_V1","HALLMARK_DNA_REPAIR",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_COMPLEMENT")
pathway_sets <- hallmark_list[pathways_10]
stopifnot(all(lengths(pathway_sets) > 0))

## Scored on expr_resid, the same matrix the axes were built on.
pw_scores <- gsva(ssgseaParam(expr_resid_mat, pathway_sets, normalize = TRUE))
pw_df <- as.data.frame(scale(t(pw_scores)))
pw_df$sample_id <- rownames(pw_df)

## Sample level -> patient level, earliest timepoint, matching the
## clinical frame. Joining sample scores to a patient frame is what the
## original did and it pseudoreplicates.
pw_pat <- ax_clean_dedup %>%
  dplyr::select(Participant_ID, sample_id) %>%
  inner_join(pw_df, by = "sample_id") %>%
  dplyr::select(-sample_id)
cat("pathway scores at patient level:", nrow(pw_pat), "\n")

pw_clin <- clin_frame %>%
  dplyr::select(any_of(c("Participant_ID","group","SEX","AGE_AT_SYMPTOM_ONSET",
                         endpoints, sv_vars))) %>%
  inner_join(pw_pat, by = "Participant_ID")
cat("analysis frame:", nrow(pw_clin), "patients x", length(pathways_10), "pathways\n")

assoc <- list()
for (ep in endpoints) for (pw in pathways_10) {
  d <- pw_clin %>% dplyr::select(all_of(c(ep, pw, "AGE_AT_SYMPTOM_ONSET","SEX", sv_vars))) %>%
    tidyr::drop_na()
  if (nrow(d) < 30 || length(unique(d[[ep]])) < 5) next
  fml <- as.formula(paste0("`", ep, "` ~ `", pw, "` + ", cov_str))
  fit <- lm(fml, data = d)
  co  <- tidy(fit, conf.int = TRUE)
  row <- co[gsub("`", "", co$term) == pw, ]
  if (nrow(row)) assoc[[paste(ep, pw, sep = "__")]] <- data.frame(
    endpoint = ep, pathway = pw, beta = row$estimate, se = row$std.error,
    ci_low = row$conf.low, ci_high = row$conf.high, p = row$p.value,
    n = nobs(fit), adj_r2 = summary(fit)$adj.r.squared, stringsAsFactors = FALSE)
}
assoc_df <- bind_rows(assoc)
if (!nrow(assoc_df)) stop("no pathway-endpoint models were estimable")
assoc_df$FDR <- p.adjust(assoc_df$p, "BH")          # across the whole matrix
assoc_df <- assoc_df %>% arrange(FDR, p)
cat("\nTop pathway x endpoint associations:\n"); print(head(assoc_df, 15))
cat("\ncells at FDR < 0.05:", sum(assoc_df$FDR < 0.05), "of", nrow(assoc_df), "\n")
out_csv(assoc_df, "phase2_pathway_endpoint_associations.csv")

## joint model per endpoint: which pathways survive when they compete
joint <- list()
for (ep in endpoints) {
  d <- pw_clin %>% dplyr::select(all_of(c(ep, pathways_10, "AGE_AT_SYMPTOM_ONSET","SEX", sv_vars))) %>%
    tidyr::drop_na()
  if (nrow(d) < 40) next
  fml <- as.formula(paste0("`", ep, "` ~ ", paste0("`", pathways_10, "`", collapse = " + "),
                           " + ", cov_str))
  fit <- lm(fml, data = d)
  joint[[ep]] <- tidy(fit, conf.int = TRUE) %>%
    filter(gsub("`", "", term) %in% pathways_10) %>%
    mutate(endpoint = ep, n = nobs(fit))
}
joint_df <- bind_rows(joint)
if (nrow(joint_df)) {
  joint_df$FDR <- p.adjust(joint_df$p.value, "BH")
  out_csv(joint_df, "phase2_pathway_endpoint_joint.csv")
}

hm <- assoc_df %>%
  mutate(star = ifelse(FDR < 0.001, "***", ifelse(FDR < 0.01, "**",
                ifelse(FDR < 0.05, "*", ""))),
         pathway = factor(pathway, levels = rev(pathways_10)))
p_hm <- ggplot(hm, aes(endpoint, pathway, fill = beta)) +
  geom_tile(colour = "white") + geom_text(aes(label = star), size = 5) +
  scale_fill_gradient2(low = "#3E6E9E", mid = "white", high = "#C0504D", midpoint = 0) +
  labs(title = "Pathway x clinical-endpoint associations (sporadic ALS)",
       subtitle = "beta from endpoint ~ pathway + age + sex + SVs; * FDR<0.05 ** <0.01 *** <0.001",
       x = NULL, y = NULL, fill = "beta") +
  theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("phase2_pathway_endpoint_heatmap.pdf", p_hm, width = 7, height = 6)
cat("wrote phase2_pathway_endpoint_heatmap.pdf\n")

# =====================================================================
# PHASE 6  Familial proximity of the sporadic subtypes
# How close does each subtype sit to C9 and to SOD1 in 3-axis space?
# =====================================================================
cat("\n\n=========== PHASE 6: subtype proximity to familial biology ===========\n")
if (is.null(axis_merged)) {
  cat("axis_merged not available; skipped.\n")
} else {
  fam <- axis_merged
  fam$cluster_final <- ax_clean_df$cluster_clean[match(fam$sample_id, ax_clean_df$sample_id)]
  cen <- function(g) colMeans(fam[which(fam$group == g), axis_vars], na.rm = TRUE)
  cen_C9 <- cen("C9orf72"); cen_SOD1 <- cen("SOD1"); cen_Ctrl <- cen("Control")
  cat("centroids:\n")
  print(round(rbind(C9 = cen_C9, SOD1 = cen_SOD1, Control = cen_Ctrl), 3))

  sub_cen <- fam %>% filter(!is.na(cluster_final)) %>%
    group_by(cluster_final) %>%
    summarise(across(all_of(axis_vars), mean), n = n(), .groups = "drop")
  d3 <- function(v, c) sqrt(sum((v - c)^2))
  prox <- sub_cen %>% rowwise() %>%
    mutate(dist_C9   = d3(c_across(all_of(axis_vars)), cen_C9),
           dist_SOD1 = d3(c_across(all_of(axis_vars)), cen_SOD1),
           dist_Ctrl = d3(c_across(all_of(axis_vars)), cen_Ctrl),
           closest_familial = ifelse(dist_C9 < dist_SOD1, "C9-like", "SOD1-like")) %>%
    ungroup()
  cat("\ndistance to familial centroids:\n"); print(as.data.frame(prox))
  out_csv(prox, "phase6_subtype_familial_proximity.csv")
  cat("\nNOTE: centroid distance in 3-axis space is a DESCRIPTIVE summary.\n")
  cat("It has no significance test attached and should not be reported as\n")
  cat("one. The convergence claim rests on the DE and NES correlations.\n")
}

# =====================================================================
# PHASE 7  Top-500 directional agreement, FIXED reference set
# One pre-specified set of 500 genes (most significant in C9), evaluated
# in sporadic and in SOD1, so the gene set is identical on both sides.
# =====================================================================
cat("\n\n=========== PHASE 7: top-500 directional agreement ===========\n")
c9 <- res_C9; sp <- res_Sporadic; so <- res_SOD1
c9$gene <- rownames(c9); sp$gene <- rownames(sp); so$gene <- rownames(so)
ref500 <- c9 %>% filter(!is.na(padj)) %>% arrange(padj) %>%
  slice_head(n = 500) %>% pull(gene)
cat("reference: top 500 by C9 adjusted p-value\n")

eval_ref <- function(other, lab) {
  m <- inner_join(c9 %>% filter(gene %in% ref500) %>%
                    dplyr::select(gene, lfc_C9 = log2FoldChange),
                  other %>% dplyr::select(gene, lfc_other = log2FoldChange),
                  by = "gene") %>% filter(is.finite(lfc_C9), is.finite(lfc_other))
  data.frame(comparison = lab, n = nrow(m),
             dir_agree = mean(sign(m$lfc_C9) == sign(m$lfc_other)),
             spearman  = cor(m$lfc_C9, m$lfc_other, method = "spearman"),
             stringsAsFactors = FALSE)
}
p7 <- bind_rows(eval_ref(sp, "C9-significant genes in Sporadic"),
                eval_ref(so, "C9-significant genes in SOD1"))
print(p7)
out_csv(p7, "phase7_top500_fixed_directional.csv")
cat("\nReport as 'top 500 by C9 significance, evaluated in the other\n")
cat("groups', so the gene set is identical on both sides.\n")

cat("\n\nPHASES 1, 2, 6, 7 COMPLETE.\n")
cat("GO-BP work is in 32_database_independence.R.\n")
