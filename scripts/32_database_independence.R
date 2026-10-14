# =====================================================================
# 32_database_independence.R   (corrected mohamadFeedback2.R)
#
# Run 30_pathway_preamble.R first.
#
#   STEP 1  Rebuild the three axes from GO-BP and from Reactome,
#           recluster, and compare to the Hallmark partition.
#           This is the database-independence check Methods describes.
#   STEP 2  Data-driven themes from the familial-vs-control GO-BP
#           enrichments (RAW only -- see the warning below).
#   STEP 3  ARI comparison across all clusterings.
#
# CHANGES FROM mohamadFeedback2.R
#   - counts_filtered / universe now come from the preamble, so Step 2
#     no longer dies on a missing object.
#   - ConsensusClusterPlus is called with plot = "png" and a title after
#     graphics.off(). With plot = NULL it inherits whatever device is
#     current and throws "invalid graphics state" once anything has
#     written a pdf earlier in the session.
#   - rrvgo / GOSemSim are OPTIONAL. The original library()'d them at the
#     top, so the whole script failed if they were not installed, even
#     though Step 1 -- the part Methods actually describes -- does not
#     need them.
#   - STEP 2 CLEANED IS DISABLED BY DEFAULT. See below.
#
# ============ ABOUT STEP 2 "CLEANED" ============
# Residualising expression on the PCT_* differentiation markers and then
# re-clustering is ROUND 2, which was run on 2026-08-13 and abandoned:
# it produced a degenerate 269 / 130 / 1 solution. Its ARI against the
# Hallmark partition is near zero because one cluster has a single
# member, and that number means nothing.
#
# It is retained here behind RUN_CLEANED for disclosure, but if it is
# run the output must be reported as the known-degenerate result, never
# as a fresh robustness check that failed.
# ================================================
#
# LANDMARKS (manifest): axis correlations vs Hallmark 0.79-0.97;
# ARI 0.49 (GO-BP) and 0.59 (Reactome) for the discrete partition.
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr)
  library(GSVA); library(msigdbr); library(ConsensusClusterPlus); library(aricode)
  library(clusterProfiler); library(org.Hs.eg.db)
})
set.seed(42)
RUN_CLEANED <- FALSE      # see the warning above before setting TRUE

need <- c("expr_resid_mat","ax_clean_df","marker_data","universe_ens",
          "res_C9","res_SOD1","res_Sporadic","hallmark_list","axis_sets")
miss <- need[!sapply(need, exists)]
if (length(miss)) stop("run 30_pathway_preamble.R first. Missing: ",
                       paste(miss, collapse = ", "))

sporadic_ids <- ax_clean_df$sample_id[!is.na(ax_clean_df$cluster_clean)]
hallmark_labels <- setNames(as.character(ax_clean_df$cluster_clean),
                            ax_clean_df$sample_id)[sporadic_ids]
cat("sporadic samples with labels:", length(sporadic_ids), "\n")
print(table(hallmark_labels))

.msig <- function(...) {
  a <- list(...)
  tryCatch(do.call(msigdbr, c(list(species = "Homo sapiens"), a)),
           error = function(e) {
             names(a)[names(a) == "category"]    <- "collection"
             names(a)[names(a) == "subcategory"] <- "subcollection"
             do.call(msigdbr, c(list(species = "Homo sapiens"), a))
           })
}

score_axes <- function(gene_sets, expr, ids) {
  m <- as.matrix(expr); storage.mode(m) <- "double"
  sc <- gsva(ssgseaParam(m, gene_sets, normalize = TRUE))
  d  <- as.data.frame(scale(t(sc)))
  d$sample_id <- rownames(d)
  d[d$sample_id %in% ids, , drop = FALSE]
}

cluster_axes <- function(axis_df, k = 3, title = "cc_tmp") {
  ax <- axis_df[, setdiff(colnames(axis_df), "sample_id"), drop = FALSE]
  rownames(ax) <- axis_df$sample_id
  m <- t(as.matrix(as.data.frame(scale(ax))))
  ## AUDIT FIX: plot = NULL draws to the current device and errors with
  ## "invalid graphics state" after any earlier pdf(). Write to files.
  graphics.off()
  cc <- ConsensusClusterPlus(m, maxK = 6, reps = 1000, pItem = 0.8,
                             pFeature = 1.0, clusterAlg = "hc",
                             distance = "euclidean", seed = 42,
                             plot = "png", title = title)
  list(cc = cc, labels = setNames(cc[[k]]$consensusClass, colnames(m)))
}
pac <- function(cc, ks = 2:6) sapply(ks, function(k) {
  m <- cc[[k]]$consensusMatrix; sum(m > 0.1 & m < 0.9) / length(m) })

# =====================================================================
# STEP 1  GO-BP and Reactome replication
# =====================================================================
cat("\n\n===== STEP 1: GO-BP and Reactome replication =====\n")
c5   <- .msig(category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol) %>% distinct()
reac <- .msig(category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol) %>% distinct()
c5_list   <- split(c5$gene_symbol, c5$gs_name)
reac_list <- split(reac$gene_symbol, reac$gs_name)

pick <- function(lst, pattern) unique(unlist(lst[grep(pattern, names(lst))]))
go_axis_sets <- list(
  AXIS_METABOLIC    = pick(c5_list, "OXIDATIVE_PHOSPHORYLATION|FATTY_ACID_(BETA_)?OXIDATION|GLYCOLY|CHOLESTEROL_BIOSYNTH|STEROL_BIOSYNTH|TRICARBOXYLIC_ACID|ATP_SYNTHESIS"),
  AXIS_CELL_CYCLE   = pick(c5_list, "DNA_REPLICATION|MITOTIC_CELL_CYCLE|CELL_CYCLE_DNA_REPLICATION|DNA_REPAIR|CHROMOSOME_SEGREGATION"),
  AXIS_INFLAMMATION = pick(c5_list, "INFLAMMATORY_RESPONSE|INTERFERON_GAMMA|COMPLEMENT_ACTIVATION|RESPONSE_TO_TYPE_II_INTERFERON|TUMOR_NECROSIS_FACTOR"))
reac_axis_sets <- list(
  AXIS_METABOLIC    = pick(reac_list, "RESPIRATORY_ELECTRON_TRANSPORT|GLYCOLYSIS|FATTY_ACID|CHOLESTEROL_BIOSYNTHESIS|CITRIC_ACID|TCA"),
  AXIS_CELL_CYCLE   = pick(reac_list, "CELL_CYCLE|DNA_REPLICATION|DNA_REPAIR|MITOTIC|G2_M"),
  AXIS_INFLAMMATION = pick(reac_list, "INTERFERON|COMPLEMENT|CYTOKINE_SIGNALING|INNATE_IMMUNE"))
cat("GO-BP theme sizes:   "); print(lengths(go_axis_sets))
cat("Reactome theme sizes:"); print(lengths(reac_axis_sets))
stopifnot(all(lengths(go_axis_sets) > 10), all(lengths(reac_axis_sets) > 10))

go_df   <- score_axes(go_axis_sets,   expr_resid_mat, sporadic_ids)
reac_df <- score_axes(reac_axis_sets, expr_resid_mat, sporadic_ids)

hall <- ax_clean_df[match(sporadic_ids, ax_clean_df$sample_id),
                    c("AXIS_METABOLIC","AXIS_CELL_CYCLE","AXIS_INFLAMMATION")]
corr_tbl <- function(d, lab) data.frame(
  source = lab,
  metabolic    = cor(hall$AXIS_METABOLIC,    d$AXIS_METABOLIC[match(sporadic_ids, d$sample_id)],    use = "complete.obs"),
  cell_cycle   = cor(hall$AXIS_CELL_CYCLE,   d$AXIS_CELL_CYCLE[match(sporadic_ids, d$sample_id)],   use = "complete.obs"),
  inflammation = cor(hall$AXIS_INFLAMMATION, d$AXIS_INFLAMMATION[match(sporadic_ids, d$sample_id)], use = "complete.obs"))
axis_cor <- bind_rows(corr_tbl(go_df, "GO-BP"), corr_tbl(reac_df, "Reactome"))
cat("\naxis correlations vs Hallmark (manifest: 0.79 to 0.97):\n"); print(axis_cor)
write_csv(axis_cor, "step1_axis_correlations.csv")

go_cl   <- cluster_axes(go_df,   k = 3, title = "cc_go")
reac_cl <- cluster_axes(reac_df, k = 3, title = "cc_reac")
ord <- sporadic_ids
ari_go   <- ARI(as.integer(factor(hallmark_labels[ord])), as.integer(go_cl$labels[ord]))
ari_reac <- ARI(as.integer(factor(hallmark_labels[ord])), as.integer(reac_cl$labels[ord]))
cat(sprintf("\nARI Hallmark vs GO-BP    : %.3f   (manifest 0.49)\n", ari_go))
cat(sprintf("ARI Hallmark vs Reactome : %.3f   (manifest 0.59)\n", ari_reac))
cat("PAC GO-BP    k2-k6:", round(pac(go_cl$cc), 3), "\n")
cat("PAC Reactome k2-k6:", round(pac(reac_cl$cc), 3), "\n")
cat("\nCluster sizes (check for degenerate splits before reading any ARI):\n")
print(table(go_cl$labels));  print(table(reac_cl$labels))

# =====================================================================
# STEP 2  Data-driven GO-BP themes
# =====================================================================
cat("\n\n===== STEP 2: data-driven GO themes =====\n")
have_rrvgo <- requireNamespace("rrvgo", quietly = TRUE) &&
              requireNamespace("GOSemSim", quietly = TRUE)
if (!have_rrvgo)
  cat("rrvgo / GOSemSim absent -- semantic reduction skipped.\n",
      "  BiocManager::install(c('rrvgo','GOSemSim')) to enable Step 2.\n")

step2_raw <- NULL
if (have_rrvgo) {
  m2e <- function(ens) na.omit(AnnotationDbi::mapIds(org.Hs.eg.db, ens, "ENTREZID", "ENSEMBL"))
  de_sig <- function(res, lfc = 0.3) {
    d <- as.data.frame(res); d$ens <- sub("\\..*$", "", rownames(d))
    d$ens[!is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange) > lfc]
  }
  uni_entrez <- m2e(universe_ens)
  contrasts  <- list(C9 = res_C9, SOD1 = res_SOD1, Sporadic = res_Sporadic)
  go_each <- lapply(names(contrasts), function(nm) {
    g <- m2e(de_sig(contrasts[[nm]]))
    if (!length(g)) return(NULL)
    e <- enrichGO(gene = g, universe = uni_entrez, OrgDb = org.Hs.eg.db,
                  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                  qvalueCutoff = 0.10, readable = TRUE)
    if (is.null(e) || !nrow(e@result)) return(NULL)
    as.data.frame(e@result) %>% filter(p.adjust < 0.05) %>% mutate(contrast = nm)
  })
  go_all <- bind_rows(go_each)
  write_csv(go_all, "step2_GOBP_per_contrast_RAW.csv")
  recurring <- go_all %>% count(ID, Description) %>% filter(n >= 2) %>% pull(ID)
  cat("GO-BP terms significant in >= 2 contrasts:", length(recurring), "\n")

  if (length(recurring) >= 3) {
    term_score <- go_all %>% filter(ID %in% recurring) %>%
      group_by(ID) %>% summarise(s = max(-log10(p.adjust)), .groups = "drop")
    sc_named <- setNames(term_score$s, term_score$ID)
    sim <- rrvgo::calculateSimMatrix(recurring, orgdb = "org.Hs.eg.db",
                                     ont = "BP", method = "Rel")
    red <- rrvgo::reduceSimMatrix(sim, sc_named[recurring], threshold = 0.7,
                                  orgdb = "org.Hs.eg.db")
    write_csv(red, "step2_themes_RAW.csv")
    cat("parent themes:\n"); print(table(red$parentTerm))

    gs <- list()
    for (pt in unique(red$parentTerm)) {
      ids <- red$go[red$parentTerm == pt]
      genes <- go_all %>% filter(ID %in% ids) %>% pull(geneID) %>%
        strsplit("/") %>% unlist() %>% unique()
      if (length(genes) >= 10) gs[[make.names(pt)]] <- genes
    }
    cat("data-derived theme sizes:\n"); print(lengths(gs))
    if (length(gs) >= 2) {
      td <- score_axes(gs, expr_resid_mat, sporadic_ids)
      cl <- cluster_axes(td, k = 3, title = "cc_step2_RAW")
      pacs <- pac(cl$cc); names(pacs) <- 2:6
      k_star <- as.integer(names(pacs)[which.min(pacs)])
      labs <- setNames(cl$cc[[k_star]]$consensusClass, td$sample_id)
      o2  <- intersect(sporadic_ids, names(labs))
      ari <- ARI(as.integer(factor(hallmark_labels[o2])), as.integer(labs[o2]))
      cat("PAC by k:\n"); print(round(pacs, 3))
      cat(sprintf("most stable k = %d; ARI vs Hallmark = %.3f\n", k_star, ari))
      cat("cluster sizes:\n"); print(table(labs))
      cat("A near-zero ARI with one tiny cluster is DEGENERATE, not a\n")
      cat("failed replication. Check the sizes before interpreting.\n")
      step2_raw <- list(themes = gs, labels = labs, k = k_star, ari = ari)
    }
  } else cat("too few recurring terms for semantic reduction.\n")
}

if (RUN_CLEANED) {
  cat("\n*** STEP 2 CLEANED: this is the ABANDONED round-2 PCT marker\n")
  cat("*** residualisation (2026-08-13, degenerate 269/130/1). Any ARI\n")
  cat("*** it produces must be reported as the known-degenerate result.\n")
  mk <- intersect(c("PCT_SMI32","PCT_ISL1","PCT_NKX61","PCT_TUJ1",
                    "PCT_S100b","PCT_Nestin"), colnames(marker_data))
  mrk <- marker_data[match(colnames(expr_symbol), marker_data$sample_id), mk, drop = FALSE]
  ok  <- which(complete.cases(mrk))
  cat("samples with complete markers:", length(ok), "of", ncol(expr_symbol), "\n")
  Xm  <- model.matrix(~ ., data = as.data.frame(mrk[ok, , drop = FALSE]))
  emr <- expr_symbol[, ok, drop = FALSE]
  emr <- t(apply(emr, 1, function(y) resid(lm(y ~ Xm - 1))))
  colnames(emr) <- colnames(expr_symbol)[ok]
  cat("marker-residualised matrix:", paste(dim(emr), collapse = " x "), "\n")
} else {
  cat("\nSTEP 2 CLEANED skipped (RUN_CLEANED = FALSE).\n")
  cat("It is the abandoned round-2 residualisation; see the header.\n")
}

# =====================================================================
# STEP 3  ARI comparison matrix
# =====================================================================
cat("\n\n===== STEP 3: ARI comparison =====\n")
labs_list <- list(
  Hallmark      = factor(hallmark_labels[sporadic_ids]),
  GO_BP_axes    = factor(go_cl$labels[sporadic_ids]),
  Reactome_axes = factor(reac_cl$labels[sporadic_ids]))
if (!is.null(step2_raw))
  labs_list$DataDriven_RAW <- factor(
    step2_raw$labels[intersect(sporadic_ids, names(step2_raw$labels))])

nm <- names(labs_list)
A <- matrix(NA_real_, length(nm), length(nm), dimnames = list(nm, nm))
for (i in nm) for (j in nm) {
  common <- intersect(names(labs_list[[i]]), names(labs_list[[j]]))
  if (length(common) > 10)
    A[i, j] <- ARI(as.integer(labs_list[[i]][common]),
                   as.integer(labs_list[[j]][common]))
}
cat("\nARI matrix (1 = identical, 0 = random):\n"); print(round(A, 3))
write_csv(as.data.frame(A) %>% rownames_to_column("clustering"),
          "step3_ARI_matrix.csv")

saveRDS(list(axis_cor = axis_cor, ari_go = ari_go, ari_reac = ari_reac,
             pac_go = pac(go_cl$cc), pac_reac = pac(reac_cl$cc),
             ari_matrix = A, step2_raw = step2_raw),
        "database_independence_FINAL.rds")
cat("\nsaved database_independence_FINAL.rds\n")

cat("\nINTERPRET\n")
cat(" STEP 1 is the reviewer-facing result: if the GO-BP and Reactome\n")
cat("        axes correlate highly with Hallmark and the partitions\n")
cat("        agree at ARI ~0.5-0.6, the framework is not an artefact of\n")
cat("        one database. Check cluster SIZES first: a high correlation\n")
cat("        with a degenerate partition is not a replication.\n")
cat(" STEP 2 RAW is exploratory. A differentiation or ECM theme\n")
cat("        appearing is the known confound, not a new finding.\n")
