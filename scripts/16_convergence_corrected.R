# 16_convergence_corrected.R
# RNA vs proteomics convergence, rerun on the CORRECTED genotype DE
# layer (res_C9 / res_SOD1 / res_Sporadic from
# genotype_de_corrected_FINAL.rds, script 15).
#
# Supersedes the script-12 run, which consumed res_*_filtered -- the
# December 2025 PRE-correction fit whose control group contains the
# replicated reference line. Every RNA-side number from that run is
# retired along with the objects. Protein inputs are unchanged.
#
# Outputs go to ./convergence_results_corrected/ so nothing from the
# contaminated run can be picked up by mistake.
#
# Gene-level cross-genotype correlations (Section 10) now use logFC as
# the primary metric, matching the manuscript convention (Spearman of
# log2 fold changes across all shared genes); the Wald-stat version is
# kept as a secondary column, clearly labelled, and must not appear in
# the paper or on slides.

#Run this first
#setwd("~/Documents/DSN/ALS")

#source("16a_pin_proteomics_inputs.R")

#de <- readRDS(path.expand("~/Documents/convergence_results/genotype_de_corrected_FINAL.rds"))
#res_C9 <- de$res_C9; res_SOD1 <- de$res_SOD1; res_Sporadic <- de$res_Sporadic

#source("14_tals_rebuild.R")          # optional, unblocks Section 11
#source("16_convergence_corrected.R") # only now you can generate figure 2 

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

set.seed(42)

if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename,
                               dplyr::mutate, dplyr::summarise, dplyr::count,
                               dplyr::slice, base::setdiff, base::intersect,
                               base::union, base::unname, .quiet = TRUE)
}

out_dir <- "./convergence_results_corrected"
dir.create(out_dir, showWarnings = FALSE)

## ---- 0. LOAD OBJECTS ----
## Corrected genotype DE layer, saved by 15_genotype_de_rebuild.R.
## Only load if not already in the session (do not clobber a session
## where 15 just ran).
de_file <- "convergence_results/genotype_de_corrected_FINAL.rds"
if (!all(vapply(c("res_C9", "res_SOD1", "res_Sporadic"), exists, logical(1))) &&
    file.exists(de_file)) {
  de <- readRDS(de_file)
  res_C9 <- de$res_C9; res_SOD1 <- de$res_SOD1; res_Sporadic <- de$res_Sporadic
  cat("Loaded corrected genotype DE (", de$design, ",", de$n_genes, "genes )\n")
}
## Protein objects: auto-load the pinned proteomics inputs if present.
prot_file <- path.expand("~/Documents/DSN/ALS/proteomics/proteomics_inputs_FINAL.rds")
prot_req <- c("res_c9", "res_sod1", "res_spor", "fg_c9", "fg_sod1", "fg_spor")
if (!all(vapply(prot_req, exists, logical(1))) && file.exists(prot_file)) {
  pi <- readRDS(prot_file)
  list2env(pi, envir = .GlobalEnv)
  cat("Loaded pinned proteomics inputs (", paste(names(pi), collapse = ", "), ")\n")
}

required <- c(
  "res_C9", "res_SOD1", "res_Sporadic",                             # corrected RNA DESeq2
  "res_c9", "res_sod1", "res_spor",                                 # protein limma
  "fg_c9", "fg_sod1", "fg_spor"                                     # protein pathway fgsea tables
)
optional <- c("C9_tab", "SOD1_tab", "SPORADIC_tab",                 # RNA clusterProfiler GSEA
              "res_C9_tals", "res_SOD1_tals", "res_Spor_tals")      # TargetALS bulk

missing_req <- required[!vapply(required, exists, logical(1))]
if (length(missing_req) > 0) {
  stop("Load these objects first: ", paste(missing_req, collapse = ", "))
}
missing_opt <- optional[!vapply(optional, exists, logical(1))]
if (length(missing_opt) > 0) {
  cat("Optional objects not loaded (those sections will be skipped or",
      "use fallbacks):", paste(missing_opt, collapse = ", "), "\n")
}

## ---- 1. Helpers ----

strip_ens_version <- function(x) sub("\\..*$", "", x)

split_core <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  if (is.list(x)) return(unique(unlist(x)))
  if (is.na(x) || x == "") return(character(0))
  unlist(strsplit(x, "/"))
}

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

## DESeq2 stat-ranked vector in version-stripped ENSEMBL space,
## deduplicated by strongest |stat|.
make_rank_from_deseq <- function(res) {
  df <- as.data.frame(res)
  r <- df$stat
  names(r) <- strip_ens_version(rownames(df))
  r <- r[is.finite(r) & !is.na(names(r)) & names(r) != ""]
  tmp <- tibble(g = names(r), v = as.numeric(r)) %>%
    group_by(g) %>%
    slice_max(order_by = abs(v), n = 1, with_ties = FALSE) %>%
    ungroup()
  out <- tmp$v
  names(out) <- tmp$g
  sort(out, decreasing = TRUE)
}

ens_to_symbol <- function(ens_vec) {
  ens_vec <- strip_ens_version(ens_vec)
  map <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(ens_vec),
                               keytype = "ENSEMBL", columns = "SYMBOL") %>%
    filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    distinct(ENSEMBL, .keep_all = TRUE)
  sym <- map$SYMBOL
  names(sym) <- map$ENSEMBL
  sym
}

run_fgsea <- function(rank_vec, pathways_list, tag) {
  stopifnot(is.numeric(rank_vec), !is.null(names(rank_vec)))
  r <- rank_vec[is.finite(rank_vec) & !is.na(names(rank_vec)) & names(rank_vec) != ""]
  r <- sort(r, decreasing = TRUE)
  fg <- fgsea(pathways = pathways_list, stats = r,
              minSize = 10, maxSize = 500, eps = 0)
  fg %>%
    as.data.frame() %>%
    transmute(pathway, NES, pval, FDR = padj, leadingEdge, contrast = tag)
}

## ---- 2. RNA gene tables (ENSEMBL and SYMBOL space) ----

rna_tbl <- function(res, label) {
  df <- as.data.frame(res)
  df$gene <- strip_ens_version(rownames(df))
  df %>%
    transmute(gene, logFC = log2FoldChange, tstat = stat,
              pval = pvalue, FDR = padj, contrast = label) %>%
    filter(!is.na(gene) & gene != "")
}

rna_c9  <- rna_tbl(res_C9,       "C9_vs_CTRL")
rna_sod <- rna_tbl(res_SOD1,     "SOD1_vs_CTRL")
rna_sp  <- rna_tbl(res_Sporadic, "SPORADIC_vs_CTRL")

to_symbol_space <- function(rna_df) {
  ens_map <- ens_to_symbol(rna_df$gene)
  rna_df %>%
    mutate(symbol = unname(ens_map[gene])) %>%
    filter(!is.na(symbol) & symbol != "") %>%
    group_by(symbol) %>%
    slice_max(order_by = abs(tstat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(gene = symbol) %>%
    dplyr::select(-symbol)
}

rna_c9_sym  <- to_symbol_space(rna_c9)
rna_sod_sym <- to_symbol_space(rna_sod)
rna_sp_sym  <- to_symbol_space(rna_sp)

## ---- 3. Protein gene tables ----

prot_gene_table <- function(res, contrast_label) {
  df <- as.data.frame(res)
  stopifnot("gene" %in% names(df))
  df %>%
    transmute(gene, logFC, tstat = t, pval = P.Value,
              FDR = adj.P.Val, contrast = contrast_label) %>%
    filter(!is.na(gene) & gene != "")
}

prot_c9_full  <- prot_gene_table(res_c9,   "C9_vs_CTRL")
prot_sod_full <- prot_gene_table(res_sod1, "SOD1_vs_CTRL")
prot_sp_full  <- prot_gene_table(res_spor, "SPORADIC_vs_CTRL")

## ---- 4. Gene-level overlap: strict FDR and nominal-protein versions ----

gene_overlap_report <- function(rna_df, prot_df, rna_fdr = 0.10,
                                prot_fdr = NULL, prot_p = NULL) {
  common <- intersect(rna_df$gene, prot_df$gene)
  rna <- rna_df %>% filter(gene %in% common)
  pro <- prot_df %>% filter(gene %in% common)
  
  rna_sig <- rna %>% filter(!is.na(FDR) & FDR < rna_fdr)
  pro_sig <- if (!is.null(prot_fdr)) {
    pro %>% filter(!is.na(FDR) & FDR < prot_fdr)
  } else {
    pro %>% filter(!is.na(pval) & pval < prot_p)
  }
  
  overlap <- intersect(rna_sig$gene, pro_sig$gene)
  a <- length(overlap)
  b <- nrow(rna_sig) - a
  c <- nrow(pro_sig) - a
  d <- length(common) - a - b - c
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE),
                    alternative = "greater")
  
  dir_agree <- NA_real_
  if (a > 0) {
    joined <- inner_join(
      rna %>% filter(gene %in% overlap) %>% dplyr::select(gene, rna_logFC = logFC),
      pro  %>% filter(gene %in% overlap) %>% dplyr::select(gene, prot_logFC = logFC),
      by = "gene"
    ) %>% mutate(agree = sign(rna_logFC) == sign(prot_logFC))
    dir_agree <- mean(joined$agree, na.rm = TRUE)
  }
  
  merged_all <- inner_join(rna %>% dplyr::select(gene, r = tstat),
                           pro %>% dplyr::select(gene, p = tstat), by = "gene")
  
  tibble(
    universe_n = length(common),
    rna_sig_n = nrow(rna_sig),
    prot_sig_n = nrow(pro_sig),
    overlap_n = a,
    fisher_or = unname(ft$estimate),
    fisher_p = ft$p.value,
    dir_agree_rate = dir_agree,
    spearman_tstat_all = suppressWarnings(
      cor(merged_all$r, merged_all$p, method = "spearman",
          use = "pairwise.complete.obs"))
  )
}

gene_overlap_all <- bind_rows(
  gene_overlap_report(rna_c9_sym,  prot_c9_full,  prot_fdr = 0.10) %>% mutate(contrast = "C9",       protein_cut = "FDR<0.10"),
  gene_overlap_report(rna_sod_sym, prot_sod_full, prot_fdr = 0.10) %>% mutate(contrast = "SOD1",     protein_cut = "FDR<0.10"),
  gene_overlap_report(rna_sp_sym,  prot_sp_full,  prot_fdr = 0.10) %>% mutate(contrast = "Sporadic", protein_cut = "FDR<0.10"),
  gene_overlap_report(rna_c9_sym,  prot_c9_full,  prot_p = 0.05)   %>% mutate(contrast = "C9",       protein_cut = "p<0.05"),
  gene_overlap_report(rna_sod_sym, prot_sod_full, prot_p = 0.05)   %>% mutate(contrast = "SOD1",     protein_cut = "p<0.05"),
  gene_overlap_report(rna_sp_sym,  prot_sp_full,  prot_p = 0.05)   %>% mutate(contrast = "Sporadic", protein_cut = "p<0.05")
)
print(gene_overlap_all)
write.csv(gene_overlap_all, file.path(out_dir, "gene_level_overlap.csv"), row.names = FALSE)

## ---- 5. RNA-significant sets tested on the protein rank ----

rna_sets_on_prot <- function(rna_sym_df, prot_df, tag, rna_fdr = 0.10) {
  sig  <- rna_sym_df %>% filter(FDR < rna_fdr)
  up   <- sig %>% filter(logFC > 0) %>% pull(gene) %>% unique()
  down <- sig %>% filter(logFC < 0) %>% pull(gene) %>% unique()
  
  ranks <- prot_df$tstat
  names(ranks) <- prot_df$gene
  ranks <- ranks[is.finite(ranks) & !is.na(names(ranks)) & names(ranks) != ""]
  ranks <- sort(ranks, decreasing = TRUE)
  
  sets <- list(RNA_UP = intersect(up, names(ranks)),
               RNA_DOWN = intersect(down, names(ranks)))
  sets <- sets[lengths(sets) >= 10]
  if (length(sets) == 0) return(NULL)
  fgsea(pathways = sets, stats = ranks, eps = 0) %>%
    as.data.frame() %>% mutate(contrast = tag)
}

fg_sets <- bind_rows(
  rna_sets_on_prot(rna_c9_sym,  prot_c9_full,  "C9"),
  rna_sets_on_prot(rna_sod_sym, prot_sod_full, "SOD1"),
  rna_sets_on_prot(rna_sp_sym,  prot_sp_full,  "Sporadic")
)
print(fg_sets %>% dplyr::select(-leadingEdge))
saveRDS(fg_sets, file.path(out_dir, "rna_sets_on_protein_fgsea.rds"))

## ---- 6. Hallmark pathway level: RNA (fgsea) ----

pathways_rna <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, ensembl_gene) %>%
  distinct() %>%
  filter(!is.na(ensembl_gene) & ensembl_gene != "") %>%
  { split(.$ensembl_gene, .$gs_name) }

rank_C9   <- make_rank_from_deseq(res_C9)
rank_SOD1 <- make_rank_from_deseq(res_SOD1)
rank_SP   <- make_rank_from_deseq(res_Sporadic)

RNA_C9_fg   <- run_fgsea(rank_C9,   pathways_rna, "C9_vs_CTRL")
RNA_SOD1_fg <- run_fgsea(rank_SOD1, pathways_rna, "SOD1_vs_CTRL")
RNA_SP_fg   <- run_fgsea(rank_SP,   pathways_rna, "SPORADIC_vs_CTRL")

## RNA pathway tables for concordance: use the clusterProfiler tables if
## loaded (they carry core_enrichment); otherwise fall back to fgsea.
prep_rna_path <- function(tab, contrast) {
  df <- as.data.frame(tab)
  stopifnot(all(c("ID", "NES", "p.adjust") %in% names(df)))
  df %>%
    transmute(pathway = ID, NES, FDR = p.adjust,
              core = if ("core_enrichment" %in% names(df)) core_enrichment else NA_character_,
              omic = "RNA", contrast = contrast)
}
fg_as_path <- function(fg, contrast) {
  fg %>% transmute(pathway, NES, FDR,
                   core = vapply(leadingEdge, paste, collapse = "/", FUN.VALUE = ""),
                   omic = "RNA", contrast = contrast)
}

rna_path_c9   <- if (exists("C9_tab"))       prep_rna_path(C9_tab, "C9_vs_CTRL")             else fg_as_path(RNA_C9_fg, "C9_vs_CTRL")
rna_path_sod1 <- if (exists("SOD1_tab"))     prep_rna_path(SOD1_tab, "SOD1_vs_CTRL")         else fg_as_path(RNA_SOD1_fg, "SOD1_vs_CTRL")
rna_path_sp   <- if (exists("SPORADIC_tab")) prep_rna_path(SPORADIC_tab, "SPORADIC_vs_CTRL") else fg_as_path(RNA_SP_fg, "SPORADIC_vs_CTRL")

## ---- 7. Hallmark pathway level: protein ----
## fg_c9 / fg_sod1 / fg_spor are fgsea outputs (pathway, NES, padj,
## leadingEdge). The helper also accepts wide NES_*/padj_* tables.

prep_prot_path <- function(tab, contrast, wide_nes = NULL, wide_fdr = NULL) {
  df <- as.data.frame(tab)
  stopifnot("pathway" %in% names(df))
  if (!is.null(wide_nes) && all(c(wide_nes, wide_fdr) %in% names(df))) {
    df %>% transmute(pathway, NES = .data[[wide_nes]], FDR = .data[[wide_fdr]],
                     core = NA_character_, omic = "Protein", contrast = contrast)
  } else if (all(c("NES", "padj") %in% names(df))) {
    df %>% transmute(pathway, NES, FDR = padj,
                     core = if ("leadingEdge" %in% names(df))
                       vapply(leadingEdge, function(z) paste(unlist(z), collapse = "/"), "")
                     else NA_character_,
                     omic = "Protein", contrast = contrast)
  } else {
    stop("Unrecognised protein pathway table for ", contrast,
         ". Columns present: ", paste(names(df), collapse = ", "))
  }
}

prot_path_c9   <- prep_prot_path(fg_c9,   "C9_vs_CTRL")
prot_path_sod1 <- prep_prot_path(fg_sod1, "SOD1_vs_CTRL")
prot_path_sp   <- prep_prot_path(fg_spor, "SPORADIC_vs_CTRL")

## ---- 8. Concordance and overlap summaries ----

make_nes_wide <- function(df_long) {
  df_long %>% dplyr::select(contrast, pathway, NES) %>% distinct() %>%
    pivot_wider(names_from = contrast, values_from = NES)
}
pair_stats <- function(wide_df, a, b) {
  x <- wide_df[[a]]; y <- wide_df[[b]]
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  tibble(n = length(x),
         spearman = suppressWarnings(cor(x, y, method = "spearman")),
         sign_agree = mean(sign(x) == sign(y)))
}

rna_nes  <- make_nes_wide(bind_rows(rna_path_c9, rna_path_sod1, rna_path_sp))
prot_nes <- make_nes_wide(bind_rows(prot_path_c9, prot_path_sod1, prot_path_sp))

pattern_table <- bind_rows(
  tibble(omic = "RNA",     pair = "C9 vs Sporadic", pair_stats(rna_nes,  "C9_vs_CTRL", "SPORADIC_vs_CTRL")),
  tibble(omic = "RNA",     pair = "C9 vs SOD1",     pair_stats(rna_nes,  "C9_vs_CTRL", "SOD1_vs_CTRL")),
  tibble(omic = "RNA",     pair = "SOD1 vs Sporadic", pair_stats(rna_nes, "SOD1_vs_CTRL", "SPORADIC_vs_CTRL")),
  tibble(omic = "Protein", pair = "C9 vs Sporadic", pair_stats(prot_nes, "C9_vs_CTRL", "SPORADIC_vs_CTRL")),
  tibble(omic = "Protein", pair = "C9 vs SOD1",     pair_stats(prot_nes, "C9_vs_CTRL", "SOD1_vs_CTRL")),
  tibble(omic = "Protein", pair = "SOD1 vs Sporadic", pair_stats(prot_nes, "SOD1_vs_CTRL", "SPORADIC_vs_CTRL"))
)
print(pattern_table)
write.csv(pattern_table, file.path(out_dir, "pattern_replication_table.csv"), row.names = FALSE)

nes_concordance <- function(rna_df, prot_df, fdr_thresh = 0.10) {
  merged <- inner_join(
    rna_df  %>% dplyr::select(pathway, NES_RNA = NES, FDR_RNA = FDR),
    prot_df %>% dplyr::select(pathway, NES_Prot = NES, FDR_Prot = FDR),
    by = "pathway")
  sig_both <- merged %>% filter(FDR_RNA < fdr_thresh, FDR_Prot < fdr_thresh)
  tibble(
    n_overlap = nrow(merged),
    spearman_all = suppressWarnings(cor(merged$NES_RNA, merged$NES_Prot, method = "spearman")),
    sign_agree_all = mean(sign(merged$NES_RNA) == sign(merged$NES_Prot)),
    n_sig_both = nrow(sig_both),
    spearman_sig_both = if (nrow(sig_both) >= 2) suppressWarnings(cor(sig_both$NES_RNA, sig_both$NES_Prot, method = "spearman")) else NA_real_,
    sign_agree_sig_both = if (nrow(sig_both) >= 1) mean(sign(sig_both$NES_RNA) == sign(sig_both$NES_Prot)) else NA_real_
  )
}

conc_table <- bind_rows(
  nes_concordance(rna_path_c9,   prot_path_c9)   %>% mutate(contrast = "C9", .before = 1),
  nes_concordance(rna_path_sod1, prot_path_sod1) %>% mutate(contrast = "SOD1", .before = 1),
  nes_concordance(rna_path_sp,   prot_path_sp)   %>% mutate(contrast = "Sporadic", .before = 1)
)
print(conc_table)
write.csv(conc_table, file.path(out_dir, "rna_protein_nes_concordance.csv"), row.names = FALSE)

sig_paths <- function(df, fdr = 0.10) df %>% filter(FDR < fdr) %>% pull(pathway) %>% unique()
overlap_summary <- function(A, B, label) {
  tibble(pair = label, A = length(A), B = length(B),
         overlap = length(intersect(A, B)),
         jaccard = if (length(union(A, B)) > 0) length(intersect(A, B)) / length(union(A, B)) else NA_real_)
}

sig_overlaps <- bind_rows(
  overlap_summary(sig_paths(rna_path_c9),  sig_paths(rna_path_sod1), "RNA C9 vs SOD1"),
  overlap_summary(sig_paths(rna_path_c9),  sig_paths(rna_path_sp),   "RNA C9 vs Sporadic"),
  overlap_summary(sig_paths(prot_path_c9), sig_paths(prot_path_sod1), "Protein C9 vs SOD1"),
  overlap_summary(sig_paths(prot_path_c9), sig_paths(prot_path_sp),   "Protein C9 vs Sporadic")
)
print(sig_overlaps)
write.csv(sig_overlaps, file.path(out_dir, "significant_pathway_overlaps.csv"), row.names = FALSE)

## Direction agreement of C9-significant pathways in the other groups
dir_agree_on_c9 <- function(base_df, other_df, label) {
  base_sig <- base_df %>% filter(FDR < 0.10) %>% dplyr::select(pathway, NES_base = NES)
  m <- inner_join(base_sig, other_df %>% dplyr::select(pathway, NES_other = NES), by = "pathway")
  tibble(comparison = label, n_c9_sig = nrow(m),
         dir_agree = mean(sign(m$NES_base) == sign(m$NES_other)))
}
dir_table <- bind_rows(
  dir_agree_on_c9(rna_path_c9,  rna_path_sp,   "RNA: Sporadic on C9-sig"),
  dir_agree_on_c9(rna_path_c9,  rna_path_sod1, "RNA: SOD1 on C9-sig"),
  dir_agree_on_c9(prot_path_c9, prot_path_sp,  "Protein: Sporadic on C9-sig"),
  dir_agree_on_c9(prot_path_c9, prot_path_sod1,"Protein: SOD1 on C9-sig")
)
print(dir_table)
write.csv(dir_table, file.path(out_dir, "c9_sig_direction_agreement.csv"), row.names = FALSE)

## Leading-edge jaccard where core genes are available
le_jaccard_tbl <- function(dfA, dfB, paths, label) {
  if (length(paths) == 0) return(tibble())
  bind_rows(lapply(paths, function(p) {
    leA <- split_core(dfA$core[dfA$pathway == p][1])
    leB <- split_core(dfB$core[dfB$pathway == p][1])
    tibble(comparison = label, pathway = p, n_A = length(leA),
           n_B = length(leB), overlap = length(intersect(leA, leB)),
           jaccard = jaccard(leA, leB))
  })) %>% arrange(desc(jaccard))
}
le_rna_c9_sod1 <- le_jaccard_tbl(rna_path_c9, rna_path_sod1,
                                 intersect(sig_paths(rna_path_c9), sig_paths(rna_path_sod1)),
                                 "RNA C9 vs SOD1")
print(le_rna_c9_sod1)
write.csv(le_rna_c9_sod1, file.path(out_dir, "leading_edge_jaccard.csv"), row.names = FALSE)

## ---- 9. Plots ----

c9_sig_prot_tbl <- prot_path_c9 %>%
  filter(FDR < 0.10) %>%
  dplyr::select(pathway, NES_C9 = NES, FDR_C9 = FDR) %>%
  left_join(prot_path_sp   %>% dplyr::select(pathway, NES_Spor = NES), by = "pathway") %>%
  left_join(prot_path_sod1 %>% dplyr::select(pathway, NES_SOD1 = NES), by = "pathway") %>%
  arrange(FDR_C9)
write.csv(c9_sig_prot_tbl, file.path(out_dir, "c9_sig_protein_pathways.csv"), row.names = FALSE)

plot_df <- c9_sig_prot_tbl %>%
  mutate(pathway_short = gsub("^HALLMARK_", "", pathway)) %>%
  dplyr::select(pathway_short, NES_C9, NES_Spor, NES_SOD1) %>%
  pivot_longer(starts_with("NES_"), names_to = "contrast", values_to = "NES") %>%
  mutate(contrast = recode(contrast, NES_C9 = "C9", NES_Spor = "Sporadic", NES_SOD1 = "SOD1"),
         pathway_short = factor(pathway_short,
                                levels = rev(gsub("^HALLMARK_", "", c9_sig_prot_tbl$pathway))))

p_dot <- ggplot(plot_df, aes(x = NES, y = pathway_short)) +
  geom_vline(xintercept = 0) +
  geom_point(size = 3) +
  facet_grid(. ~ contrast) +
  labs(title = "Proteomics: C9-significant Hallmark pathways across groups",
       x = "NES", y = NULL) +
  theme_bw()
ggsave(file.path(out_dir, "prot_hallmarks_dotplot.png"), p_dot, width = 12, height = 8, dpi = 300)

scatter_hallmarks <- function(wide_df, xcol, ycol, title) {
  df <- wide_df %>%
    dplyr::select(pathway, x = all_of(xcol), y = all_of(ycol)) %>%
    filter(is.finite(x) & is.finite(y))
  sp <- suppressWarnings(cor(df$x, df$y, method = "spearman"))
  ggplot(df, aes(x = x, y = y)) +
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    geom_point(size = 2) +
    labs(title = title,
         subtitle = paste0("Spearman = ", round(sp, 3),
                           " | sign agree = ", round(mean(sign(df$x) == sign(df$y)), 2),
                           " | n = ", nrow(df)),
         x = xcol, y = ycol) +
    theme_bw()
}

ggsave(file.path(out_dir, "RNA_C9_vs_Sporadic_scatter.png"),
       scatter_hallmarks(rna_nes, "C9_vs_CTRL", "SPORADIC_vs_CTRL", "RNA Hallmarks: C9 vs Sporadic"),
       width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "RNA_C9_vs_SOD1_scatter.png"),
       scatter_hallmarks(rna_nes, "C9_vs_CTRL", "SOD1_vs_CTRL", "RNA Hallmarks: C9 vs SOD1"),
       width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "PROT_C9_vs_Sporadic_scatter.png"),
       scatter_hallmarks(prot_nes, "C9_vs_CTRL", "SPORADIC_vs_CTRL", "Protein Hallmarks: C9 vs Sporadic"),
       width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "PROT_C9_vs_SOD1_scatter.png"),
       scatter_hallmarks(prot_nes, "C9_vs_CTRL", "SOD1_vs_CTRL", "Protein Hallmarks: C9 vs SOD1"),
       width = 6, height = 5, dpi = 300)

## ---- 10. AnswerALS RNA gene-level cross-genotype correlations ----
## PRIMARY metric: Spearman of log2 fold changes across all shared
## genes (manuscript convention). Stat-based version is secondary and
## never goes in the paper or on slides.

pair_cor <- function(df1, df2, col = "logFC") {
  m <- inner_join(df1 %>% dplyr::select(gene, s1 = all_of(col)),
                  df2 %>% dplyr::select(gene, s2 = all_of(col)), by = "gene")
  cor(m$s1, m$s2, method = "spearman", use = "pairwise.complete.obs")
}
topK_dir_agree <- function(df1, df2, K = 500) {
  d1 <- df1 %>% filter(is.finite(tstat)) %>% arrange(desc(abs(tstat))) %>%
    slice(1:K) %>% transmute(gene, s1 = sign(tstat))
  d2 <- df2 %>% transmute(gene, s2 = sign(tstat))
  m <- inner_join(d1, d2, by = "gene")
  mean(m$s1 == m$s2)
}

aals_gene_cor <- tibble(
  pair = c("C9 vs SOD1", "C9 vs Sporadic", "SOD1 vs Sporadic"),
  spearman_logFC = c(pair_cor(rna_c9, rna_sod), pair_cor(rna_c9, rna_sp), pair_cor(rna_sod, rna_sp)),
  spearman_stat_secondary = c(pair_cor(rna_c9, rna_sod, "tstat"),
                              pair_cor(rna_c9, rna_sp, "tstat"),
                              pair_cor(rna_sod, rna_sp, "tstat")),
  top500_dir_agree = c(topK_dir_agree(rna_c9, rna_sod), topK_dir_agree(rna_c9, rna_sp), topK_dir_agree(rna_sod, rna_sp))
)
print(aals_gene_cor)   # logFC column must read 0.013 / 0.289 / 0.260
write.csv(aals_gene_cor, file.path(out_dir, "aals_gene_level_correlations.csv"), row.names = FALSE)

## ---- 11. TargetALS bulk gene-level ----

tals_file <- "./convergence_results/tals_de_results.rds"
if (!all(vapply(c("res_C9_tals", "res_SOD1_tals", "res_Spor_tals"), exists, logical(1))) &&
    file.exists(tals_file)) {
  tals <- readRDS(tals_file)
  res_C9_tals <- tals$res_C9_tals; res_SOD1_tals <- tals$res_SOD1_tals
  res_Spor_tals <- tals$res_Spor_tals
  cat("Loaded saved TALS DE results.\n")
}
if (all(vapply(c("res_C9_tals", "res_SOD1_tals", "res_Spor_tals"), exists, logical(1)))) {
  tals_c9   <- rna_tbl(res_C9_tals,   "C9_vs_CTRL")   %>% filter(is.finite(tstat))
  tals_sod1 <- rna_tbl(res_SOD1_tals, "SOD1_vs_CTRL") %>% filter(is.finite(tstat))
  tals_sp   <- rna_tbl(res_Spor_tals, "SPORADIC_vs_CTRL") %>% filter(is.finite(tstat))
  
  tals_gene_cor <- tibble(
    pair = c("C9 vs SOD1", "C9 vs Sporadic", "SOD1 vs Sporadic"),
    spearman_logFC = c(pair_cor(tals_c9, tals_sod1), pair_cor(tals_c9, tals_sp), pair_cor(tals_sod1, tals_sp)),
    top500_dir_agree = c(topK_dir_agree(tals_c9, tals_sod1), topK_dir_agree(tals_c9, tals_sp), topK_dir_agree(tals_sod1, tals_sp))
  )
  print(tals_gene_cor)
  write.csv(tals_gene_cor, file.path(out_dir, "tals_gene_level_correlations.csv"), row.names = FALSE)
} else {
  cat("TALS objects not loaded; Section 11 skipped.\n")
}

saveRDS(list(gene_overlap_all = gene_overlap_all, fg_sets = fg_sets,
             pattern_table = pattern_table, conc_table = conc_table,
             sig_overlaps = sig_overlaps, dir_table = dir_table,
             c9_sig_prot_tbl = c9_sig_prot_tbl,
             RNA_fgsea = list(C9 = RNA_C9_fg, SOD1 = RNA_SOD1_fg, SP = RNA_SP_fg)),
        file.path(out_dir, "convergence_FINAL.rds"))

sessionInfo()

## ---- inputs for the CCT2 protein block ------------------------------
## These come from proteomicsDE.R, NOT from the convergence work above.
## Script 16's `required` check does not list them, so without this
## guard the script passes its own check and then dies here.
if (!exists("prot_expr_collapsed") || !exists("prot_gene"))
  stop("The CCT2 protein block needs prot_expr_collapsed and prot_gene, ",
       "which only proteomicsDE.R builds.\n",
       "  Either source('proteomicsDE.R') first, or skip this block: the ",
       "results are already pinned in proteomics_cct2_FINAL.rds.")

## patient table with subtype labels
if (!exists("pat")) {
  pat <- as.data.frame(readRDS(
    path.expand("~/Documents/DSN/ALS/patient_level_FINAL.rds"))$pat)
  cat("loaded pat from patient_level_FINAL.rds:", nrow(pat), "patients\n")
}

## Label and expression column names vary by object (subtype on the
## patient table, cluster_clean on ax_clean_df, cluster_final on
## ax_clean_dedup). Resolve rather than hardcode.
sub_col <- intersect(c("subtype", "cluster_clean", "cluster_final"), names(pat))[1]
rna_col <- intersect(c("CCT2", "CCT2_rna", "CCT2_expr"), names(pat))[1]
if (is.na(sub_col))
  stop("no subtype column in pat. Columns: ", paste(names(pat), collapse = ", "))
cat("sub_col =", sub_col, " | rna_col =", rna_col, "\n")

## Rescue CCT2 from the pre-mapping matrix
cct2_row <- "1/sp|P78371|TCPB_HUMAN"
stopifnot(cct2_row %in% rownames(prot_expr_collapsed))
cct2_vec <- as.numeric(prot_expr_collapsed[cct2_row, ])
names(cct2_vec) <- colnames(prot_expr_collapsed)
cat("CCT2 quantified in", sum(is.finite(cct2_vec)), "of", length(cct2_vec), "participants\n")

## Rebuild the test frame around the rescued vector
df <- tibble::tibble(Participant_ID = names(cct2_vec),
                     CCT2_prot = cct2_vec) %>%
  dplyr::filter(is.finite(CCT2_prot)) %>%
  inner_join(pat %>% dplyr::distinct(Participant_ID, .keep_all = TRUE) %>%
               dplyr::select(Participant_ID, subtype = all_of(sub_col),
                             CCT2_rna = any_of(rna_col)),
             by = "Participant_ID")
cat("n =", nrow(df), "\n")
print(table(df$subtype))

## Primary: S-ALS-3 vs S-ALS-2 at the protein level
d23 <- df %>% dplyr::filter(grepl("2$|3$", as.character(subtype)))
print(wilcox.test(CCT2_prot ~ droplevels(factor(subtype)), data = d23))

## Context: all three subtypes, medians
print(kruskal.test(CCT2_prot ~ subtype, data = df))
print(df %>% dplyr::group_by(subtype) %>%
        dplyr::summarise(n = dplyr::n(), median = median(CCT2_prot), IQR = IQR(CCT2_prot)))

## CCT2 transcript-protein coupling
if ("CCT2_rna" %in% names(df))
  print(cor.test(df$CCT2_rna, df$CCT2_prot, method = "spearman"))

## Secondary: TRiC complex score (all 8 subunits, z-scored then averaged)
tric_syms <- c("TCP1", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8")
stopifnot(all(tric_syms %in% rownames(prot_gene)))

common_ids <- intersect(colnames(prot_gene), names(cct2_vec))
cat("TRiC score computable for", length(common_ids), "participants\n")

tric_mat <- rbind(prot_gene[tric_syms, common_ids, drop = FALSE],
                  CCT2 = cct2_vec[common_ids])
tric_z <- t(scale(t(tric_mat)))
tric_score <- colMeans(tric_z, na.rm = TRUE)

df$TRiC_score  <- tric_score[df$Participant_ID]
d23 <- df %>% dplyr::filter(grepl("2$|3$", as.character(subtype)),
                            is.finite(TRiC_score))
print(wilcox.test(TRiC_score ~ droplevels(factor(subtype)), data = d23))
print(df %>% dplyr::filter(is.finite(TRiC_score)) %>%
        dplyr::group_by(subtype) %>%
        dplyr::summarise(n = dplyr::n(), median_TRiC = median(TRiC_score)))

#readr::write_csv(df, "~/Documents/DSN/ALS/proteomics/cct2_protein_by_subtype.csv")

# =====================================================================
# 16b_pin_proteomics.R
#
# Pins the CCT2 protein validation from 16_convergence_corrected.R.
# Run it directly after that script's proteomics block, with df,
# cct2_vec, prot_gene and tric_score still in session.
#
# Why this exists: the block wrote cct2_protein_by_subtype.csv but never
# pinned the test objects, so the manifest row for the proteomics layer
# had no reproducible artefact behind it. The CCT2 rescue from
# "1/sp|P78371|TCPB_HUMAN" is a manual step and must be recorded.
#
# Landmarks: n = 265 (59/144/62); Kruskal-Wallis p = 0.005;
#            S-ALS-3 vs S-ALS-2 Wilcoxon p = 0.033;
#            RNA-protein Spearman rho = 0.21, p = 8e-4;
#            TRiC composite not significant (p = 0.13) -> CCT2-specific.
# =====================================================================

stopifnot(exists("df"), "CCT2_prot" %in% names(df))
CCT2_ROW <- "1/sp|P78371|TCPB_HUMAN"   # provenance: rescued pre-mapping row

df$subtype <- factor(as.character(df$subtype),
                     levels = c("S-ALS-1", "S-ALS-2", "S-ALS-3"))

kw   <- kruskal.test(CCT2_prot ~ subtype, data = df)
d23  <- droplevels(subset(df, subtype %in% c("S-ALS-2", "S-ALS-3")))
wil  <- wilcox.test(CCT2_prot ~ subtype, data = d23)
rp   <- if ("CCT2_rna" %in% names(df))
  cor.test(df$CCT2_rna, df$CCT2_prot, method = "spearman") else NULL
tri  <- if ("TRiC_score" %in% names(df)) {
  dt <- droplevels(subset(df, subtype %in% c("S-ALS-2", "S-ALS-3") &
                            is.finite(TRiC_score)))
  wilcox.test(TRiC_score ~ subtype, data = dt)
} else NULL

cat("\n=== LANDMARKS (compare to manifest) ===\n")
cat("n =", nrow(df), " sizes:", paste(as.integer(table(df$subtype)), collapse = "/"), "\n")
cat(sprintf("Kruskal-Wallis p = %.4f   (expect 0.005)\n", kw$p.value))
cat(sprintf("S-ALS-3 vs S-ALS-2 Wilcoxon p = %.4f   (expect 0.033)\n", wil$p.value))
if (!is.null(rp))
  cat(sprintf("RNA-protein rho = %.3f, p = %.2g   (expect 0.21, 8e-4)\n",
              rp$estimate, rp$p.value))
if (!is.null(tri))
  cat(sprintf("TRiC composite p = %.3f   (expect 0.13, not significant)\n",
              tri$p.value))

payload <- list(
  df = df, cct2_row = CCT2_ROW,
  n = nrow(df), sizes = table(df$subtype),
  kruskal = kw, wilcox_S3_vs_S2 = wil,
  rna_protein = rp, tric_composite = tri,
  note = paste("CCT2 rescued from the pre-mapping proteomics matrix row",
               CCT2_ROW, "-- it is absent from the gene-symbol collapsed",
               "matrix, so this manual step is required to reproduce."),
  run_date = Sys.Date())

FINAL <- "~/Documents/DSN/ALS/proteomics/proteomics_cct2_FINAL.rds"
FINAL <- path.expand(FINAL)
dir.create(dirname(FINAL), showWarnings = FALSE, recursive = TRUE)
if (!file.exists(FINAL)) {
  saveRDS(payload, FINAL); cat("\nSaved", FINAL, "\n")
} else {
  saveRDS(payload, sub("_FINAL", "_AUDITCHECK", FINAL))
  old <- readRDS(FINAL)
  cat("\nFINAL exists -- wrote the AUDITCHECK copy instead.\n")
  cat("  n pinned vs rerun:", old$n, "vs", nrow(df), "\n")
  cat("  Wilcoxon p pinned vs rerun:", signif(old$wilcox_S3_vs_S2$p.value, 3),
      "vs", signif(wil$p.value, 3), "\n")
}

## Figure 5 panels e/f read the csv; keep both in step.
readr::write_csv(df, "~/Documents/DSN/ALS/proteomics/cct2_protein_by_subtype.csv")
cat("csv refreshed alongside the pin\n")
