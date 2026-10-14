# 01_rebuild_backbone.R   (AUDITED 2026-08-22)
# Canonical backbone: symbol matrix from VST (v4 recipe) + final
# configuration from backbone_FINAL.R (original 742-sample SVs, all six,
# unfiltered; the only configuration with PAC minimum at k=3, 0.299).
# Merges rebuild_backbone_v4.R + backbone_FINAL.R; finish_backbone.R is
# v3 debris and goes to archive.
#
# AUDIT ADDITIONS:
#   - All inputs from cohort_inputs_FINAL.rds; vst is recomputed from the
#     pinned corrected counts (blind VST), replacing the session object
#     vst_corrected.
#   - axis_sets is now DEFINED IN THIS FILE. It was previously a session
#     object in every version of this script, meaning the axis gene-set
#     definition itself was never under version control.
#   - Non-destructive verification: if backbone_FINAL.rds already exists,
#     the run compares its new labels against the saved ones and only
#     overwrites on exact agreement; otherwise it saves *_AUDITCHECK.rds
#     and warns.
#
# CANONICAL LABELS ARE ROUND 1. Round 2 (expression-level PCT marker
# residualization) was run on 2026-08-13, produced a degenerate solution
# (clusters 269/130/1, singleton profile +0.86/+1.81/-1.19, PAC >= 0.46
# at every k, ARI vs round 1 = 0.134), and was dropped by decision:
# cluster_final := cluster_r1. Composition sensitivity is handled at the
# axis level in 03_differentiation_sensitivity.R instead. Round 2 is
# retained below behind RUN_ROUND2_SENSITIVITY for the record.
#
# Landmarks:
#   sporadic samples 400 | round 1 clusters 89/204/107 (S-ALS-1/2/3)
#   profiles pan +0.60/+1.34/+1.17, quiescent -0.76/-0.37/-0.71,
#   metabolic +0.95/-0.41/+0.39 | round-1 PAC min at k=3 ~0.299 |
#   final labels must match saved backbone ax_final exactly (agreement 1.0);
#   patient-level agreement ~0.987 (multi-sample patients, earliest-sample rule)

suppressPackageStartupMessages({
  library(GSVA); library(ConsensusClusterPlus); library(DESeq2)
  library(AnnotationDbi); library(org.Hs.eg.db); library(msigdbr); library(dplyr)
})
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::mutate,
                               dplyr::slice, base::setdiff, base::intersect,
                               base::union, base::unname, .quiet = TRUE)
}
base::set.seed(42)
setwd("~/Documents/DSN/ALS")
AX <- c("AXIS_METABOLIC", "AXIS_INFLAMMATION", "AXIS_CELL_CYCLE")

## ---- 0. Pinned inputs ----
ci <- readRDS("~/Documents/cohort_inputs_FINAL.rds")
counts_c <- ci$counts_filtered_corrected
annot_sva_filtered <- as.data.frame(ci$annot_sva_filtered)

## ---- 1. axis_sets, now in version control ----
## Hallmark unions per the manuscript's pre-specified three-axis design.
## VERIFY ONCE against the original AxisScoring.R; the label-agreement
## check at the end will catch any divergence immediately.
h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>% distinct()
hl <- split(h$gene_symbol, h$gs_name)
axis_sets <- list(
  AXIS_METABOLIC = unique(c(hl[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],
                            hl[["HALLMARK_FATTY_ACID_METABOLISM"]],
                            hl[["HALLMARK_GLYCOLYSIS"]])),
  AXIS_CELL_CYCLE = unique(c(hl[["HALLMARK_E2F_TARGETS"]],
                             hl[["HALLMARK_G2M_CHECKPOINT"]],
                             hl[["HALLMARK_MYC_TARGETS_V1"]],
                             hl[["HALLMARK_DNA_REPAIR"]])),
  AXIS_INFLAMMATION = unique(c(hl[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
                               hl[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
                               hl[["HALLMARK_COMPLEMENT"]]))
)
cat("axis set sizes:", sapply(axis_sets, length), "\n")

## ---- 2. VST + symbol matrix (v4 recipe, from pinned counts) ----
cat("[1/5] VST + symbol matrix ...\n")
dds <- DESeqDataSetFromMatrix(counts_c, ci$annot_corrected, design = ~ 1)
expr <- assay(vst(dds, blind = TRUE)); rm(dds); gc()
ens_clean <- sub("\\..*$", "", rownames(expr))
map <- AnnotationDbi::select(org.Hs.eg.db, keys = ens_clean,
                             keytype = "ENSEMBL", columns = "SYMBOL")
map <- map[!is.na(map$SYMBOL) & !duplicated(map$ENSEMBL), ]
em <- expr[match(map$ENSEMBL, ens_clean), , drop = FALSE]
rownames(em) <- map$SYMBOL
sym <- rownames(em)
if (anyDuplicated(sym)) {
  em <- rowsum(em, group = sym, reorder = FALSE) /
        as.vector(table(sym)[unique(sym)])
  em <- em[unique(sym), , drop = FALSE]
}
E <- as.matrix(em); storage.mode(E) <- "double"
cat("   symbol matrix:", paste(dim(E), collapse = " x "), "\n")

## ---- 3. Metadata: original 742-sample SVs, corrected samples ----
meta <- annot_sva_filtered
if (is.null(meta$sample_id)) meta$sample_id <- meta$Sample_ID
meta <- meta[match(colnames(E), meta$sample_id), , drop = FALSE]
stopifnot(identical(meta$sample_id, colnames(E)))
svs <- paste0("SV", 1:6)
stopifnot(all(svs %in% colnames(meta)))
sp <- meta$sample_id[meta$group == "Sporadic"]
cat("samples", ncol(E), "| sporadic", length(sp), "(want 400) | patients",
    length(unique(meta$Participant_ID)), "\n")

## ---- helpers (verbatim from backbone_FINAL.R) ----
resid_out <- function(X, cov) {
  P <- cov %*% solve(t(cov) %*% cov) %*% t(cov)
  R <- X - X %*% P; dimnames(R) <- dimnames(X); R
}
run_ccp <- function(m) {
  grDevices::pdf(file.path(tempdir(), paste0(basename(tempfile()), ".pdf")),
                 width = 12, height = 12)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  ConsensusClusterPlus(m, maxK = 5, reps = 1000, pItem = 0.8, pFeature = 1,
                       clusterAlg = "hc", distance = "euclidean",
                       seed = 42, plot = NULL)
}
axis_df_from <- function(X) {
  s <- gsva(ssgseaParam(X, axis_sets, normalize = TRUE), verbose = FALSE)
  d <- as.data.frame(scale(t(s)))
  d <- as.data.frame(scale(d[match(sp, rownames(d)), ]))
  d$sample_id <- rownames(d); list(scores = s, df = d)
}
profile_of <- function(d, cl) {
  p <- aggregate(d[, AX], by = list(cluster = cl), FUN = mean)
  data.frame(cluster = as.character(p$cluster),
             n = as.integer(table(cl)[as.character(p$cluster)]),
             Met = round(p$AXIS_METABOLIC, 3),
             Infl = round(p$AXIS_INFLAMMATION, 3),
             CC = round(p$AXIS_CELL_CYCLE, 3))
}
label_clusters <- function(prof) {
  tot <- prof$Met + prof$Infl + prof$CC
  lab <- rep(NA_character_, nrow(prof))
  lab[which.max(tot)] <- "S-ALS-1"; lab[which.min(tot)] <- "S-ALS-2"
  lab[is.na(lab)] <- "S-ALS-3"
  setNames(lab, prof$cluster)
}
adj_rand <- function(a, b) {
  tb <- table(a, b); f <- function(x) x * (x - 1) / 2
  si <- sum(f(tb)); ra <- sum(f(rowSums(tb))); rb <- sum(f(colSums(tb)))
  n <- f(sum(tb)); (si - ra * rb / n) / ((ra + rb) / 2 - ra * rb / n)
}

## ---- 4. Round 1: SEX + SV1..SV6, all samples ----
cat("\n[2/5] round 1 ...\n")
cov1 <- model.matrix(as.formula(paste("~ SEX +", paste(svs, collapse = " + "))),
                     data = meta)
r1 <- axis_df_from(resid_out(E, cov1)); d1 <- r1$df
cr <- cor(d1[, AX], method = "spearman")
cat("cross-axis:", round(cr[1,2],3), "/", round(cr[1,3],3), "/", round(cr[2,3],3), "\n")
cc1 <- run_ccp(t(as.matrix(d1[, AX])))
pac1 <- sapply(2:5, function(k) { m <- cc1[[k]]$consensusMatrix
  sum(m > 0.1 & m < 0.9) / length(m) })
cat("PAC:", paste0("k", 2:5, "=", round(pac1, 3), collapse = "  "),
    "(want min at k=3 ~0.299)\n")
cl1 <- cc1[[3]]$consensusClass[d1$sample_id]
prof1 <- profile_of(d1, cl1)
cat("\nROUND 1 (expect 107 / 89 / 204):\n"); print(prof1, row.names = FALSE)
d1$cluster_r1 <- factor(label_clusters(prof1)[as.character(cl1)],
                        levels = paste0("S-ALS-", 1:3))
print(table(d1$cluster_r1))

## ---- 5. Round 2 sensitivity (historically degenerate; OFF by default) ----
RUN_ROUND2_SENSITIVITY <- FALSE
if (RUN_ROUND2_SENSITIVITY) {
cat("\n[3/5] round 2 (sensitivity only; expect degenerate 269/130/1) ...\n")
mk <- c("PCT_SMI32","PCT_ISL1","PCT_NKX61","PCT_TUJ1","PCT_S100b","PCT_Nestin")
mk <- mk[mk %in% colnames(meta)]
spm <- meta[match(sp, meta$sample_id), , drop = FALSE]
for (m in mk) spm[[m]] <- as.numeric(spm[[m]])
usable <- mk[sapply(mk, function(m) sum(!is.na(spm[[m]])) > 0.8 * nrow(spm))]
for (m in usable) spm[[m]][is.na(spm[[m]])] <- mean(spm[[m]], na.rm = TRUE)
cat("markers:", paste(usable, collapse = ", "), "\n")
cov2 <- model.matrix(as.formula(
  paste("~ SEX +", paste(c(svs, usable), collapse = " + "))), data = spm)
r2 <- axis_df_from(resid_out(E[, sp, drop = FALSE], cov2)); d2 <- r2$df
cc2 <- run_ccp(t(as.matrix(d2[, AX])))
pac2 <- sapply(2:5, function(k) { m <- cc2[[k]]$consensusMatrix
  sum(m > 0.1 & m < 0.9) / length(m) })
cat("PAC:", paste0("k", 2:5, "=", round(pac2, 3), collapse = "  "), "\n")
cl2 <- cc2[[3]]$consensusClass[d2$sample_id]
prof2 <- profile_of(d2, cl2)
cat("\nROUND 2:\n"); print(prof2, row.names = FALSE)
d2$cluster_final <- factor(label_clusters(prof2)[as.character(cl2)],
                           levels = paste0("S-ALS-", 1:3))
print(table(d2$cluster_final))
ari <- adj_rand(as.integer(d1$cluster_r1), as.integer(d2$cluster_final))
cat("\nARI round1 vs round2 =", round(ari, 3), "(historical: 0.134)\n")
print(table(round1 = d1$cluster_r1, final = d2$cluster_final))
} else {
  d2 <- d1; r2 <- NULL; cc2 <- NULL; pac2 <- NULL; prof2 <- NULL
  ari <- NA_real_; usable <- character(0)
  cat("\n[3/5] round 2 skipped (dropped by design; see header).\n")
}

## ---- 6. Final labels = ROUND 1; assemble, verify, save ----
d2$cluster_final <- d1$cluster_r1[match(d2$sample_id, d1$sample_id)]
d2$Participant_ID <- meta$Participant_ID[match(d2$sample_id, meta$sample_id)]
d2$cluster_r1 <- d1$cluster_r1[match(d2$sample_id, d1$sample_id)]
if ("CCT2" %in% rownames(E)) {
  d2$CCT2 <- E["CCT2", d2$sample_id]
  ct <- aggregate(CCT2 ~ cluster_final, data = d2, FUN = mean)
  cat("\nCCT2 by subtype (expect S-ALS-3 > S-ALS-2):\n")
  print(transform(ct, CCT2 = round(CCT2, 3)), row.names = FALSE)
}

backbone <- list(ax_final = d2, ax_round1 = d1, scores_r1 = r1$scores,
  scores_r2 = if (is.null(r2)) NULL else r2$scores, ccp_r1 = cc1,
  ccp_r2 = cc2, pac_r1 = pac1, pac_r2 = pac2, profile_r1 = prof1,
  profile_r2 = prof2, ari = ari,
  meta = meta, svs = svs, sv_source = "annot_sva_filtered (742-sample SVA)",
  markers = usable, expr_symbol = E, axis_sets = axis_sets,
  built = Sys.time(),
  note = paste0("AUDITED build. Corrected cohort. Symbol matrix from blind ",
    "VST of pinned corrected counts (ENSEMBL strip, org.Hs.eg.db, mean ",
    "collapse). SVs: original 742-sample SVA, all six unfiltered. ",
    "CANONICAL LABELS = ROUND 1 (SEX+SV1..6, all samples, CCP k=3 ",
    "hc/euclidean 1000 reps seed 42, labels by axis profile). Round 2 ",
    "(expression-level PCT marker residualization) was degenerate ",
    "(269/130/1, PAC>=0.46, ARI 0.134) and dropped 2026-08-13; ",
    "composition sensitivity lives at the axis level in script 03."))

cat("\n[4/5] verifying against saved backbone_FINAL.rds ...\n")
if (file.exists("backbone_FINAL.rds")) {
  old <- readRDS("backbone_FINAL.rds")
  cmp <- merge(d2[, c("sample_id", "cluster_final")],
               old$ax_final[, c("sample_id", "cluster_final")],
               by = "sample_id", suffixes = c("_new", "_old"))
  agree <- mean(as.character(cmp$cluster_final_new) ==
                as.character(cmp$cluster_final_old))
  cat("label agreement with saved backbone:", round(agree, 4), "\n")
  if (agree == 1) {
    saveRDS(backbone, "backbone_FINAL.rds")
    write.csv(d2[, c("sample_id","Participant_ID",AX,"cluster_r1","cluster_final")],
              "subtype_assignments_FINAL.csv", row.names = FALSE)
    cat("[5/5] EXACT MATCH: backbone_FINAL.rds refreshed. Script VERIFIED.\n")
  } else {
    saveRDS(backbone, "backbone_AUDITCHECK.rds")
    cat("[5/5] MISMATCH: saved as backbone_AUDITCHECK.rds; the existing",
        "backbone_FINAL.rds is untouched. Print table(cmp$cluster_final_new,",
        "cmp$cluster_final_old) and report back before anything downstream",
        "reruns.\n")
  }
} else {
  saveRDS(backbone, "backbone_FINAL.rds")
  write.csv(d2[, c("sample_id","Participant_ID",AX,"cluster_r1","cluster_final")],
            "subtype_assignments_FINAL.csv", row.names = FALSE)
  cat("[5/5] no prior backbone found; saved backbone_FINAL.rds.\n")
}
sessionInfo()
