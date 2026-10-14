# =====================================================================
# 16a_pin_proteomics_inputs.R
#
# Builds ~/Documents/DSN/ALS/proteomics/proteomics_inputs_FINAL.rds,
# the file 16_convergence_corrected.R auto-loads at its line 54.
# It never existed because proteomicsDE.R writes TSVs but no pin.
#
# WHICH OBJECTS ARE WHICH (this cost an hour, so it is written down):
#   res_C9, res_SOD1, res_Sporadic   UPPERCASE = RNA, DESeq2,
#                                    log2FoldChange / padj
#   res_c9, res_sod1, res_spor       lowercase = PROTEIN, limma,
#                                    gene / logFC / t / P.Value / adj.P.Val
#   fg_c9, fg_sod1, fg_spor          PROTEIN Hallmark fgsea, ranked on the
#                                    limma t statistic
# The two DE sets differ only by letter case. Do not infer from the name.
#
# TWO WAYS TO RUN
#   A. proteomicsDE.R just finished and its objects are in the session
#      -> this script pins them directly.
#   B. fresh session -> it rebuilds from the TSVs proteomicsDE.R wrote
#      into outputs/proteomics_guardrail/ and recomputes fgsea from the
#      limma t statistic exactly as do_fgsea() did (nperm = 10000).
# =====================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr) })
set.seed(42)

PROT_DIR <- path.expand("~/Documents/DSN/ALS/proteomics")
FINAL    <- file.path(PROT_DIR, "proteomics_inputs_FINAL.rds")
TSV_DIRS <- c("outputs/proteomics_guardrail",
              path.expand("~/Documents/DSN/ALS/outputs/proteomics_guardrail"),
              path.expand("~/Documents/DSN/ALS/proteomics"),
              ".")

have_limma <- all(sapply(c("res_c9","res_sod1","res_spor"), exists))
have_fg    <- all(sapply(c("fg_c9","fg_sod1","fg_spor"), exists))

## The expression MATRICES cannot be rebuilt from the TSVs. Sections 12+
## of script 16 need them:
##   prot_expr_collapsed  proteins x participants, PRE gene-symbol mapping.
##                        Required for the CCT2 rescue, because CCT2 is
##                        absent from the symbol-collapsed matrix and has
##                        to be pulled from row "1/sp|P78371|TCPB_HUMAN".
##   prot_gene            HGNC symbol x participants, used for the TRiC
##                        composite (TCP1, CCT3-CCT8).
## If they are not in the session and not already pinned, proteomicsDE.R
## must be sourced once to build them.
have_mats <- all(sapply(c("prot_expr_collapsed","prot_gene"), exists))
if (!have_mats && file.exists(FINAL)) {
  .old <- readRDS(FINAL)
  if (!is.null(.old$prot_expr_collapsed)) {
    prot_expr_collapsed <- .old$prot_expr_collapsed
    prot_gene           <- .old$prot_gene
    if (!is.null(.old$meta_prot2) && !exists("meta_prot2"))
      meta_prot2 <- .old$meta_prot2
    have_mats <- TRUE
    cat("expression matrices restored from the existing pin\n")
  }
}

## ---- route B: rebuild the limma tables from the TSVs -----------------
if (!have_limma) {
  find_tsv <- function(f) {
    for (d in TSV_DIRS) { p <- file.path(d, f); if (file.exists(p)) return(p) }
    hits <- unlist(lapply(TSV_DIRS[dir.exists(TSV_DIRS)], function(d)
      list.files(d, pattern = paste0("^", f, "$"), recursive = TRUE, full.names = TRUE)))
    if (length(hits)) return(hits[1])
    stop("cannot find ", f, " -- searched: ", paste(TSV_DIRS, collapse = ", "))
  }
  res_c9   <- as.data.frame(read_tsv(find_tsv("limma_C9orf72_vs_Control.tsv"),  show_col_types = FALSE))
  res_sod1 <- as.data.frame(read_tsv(find_tsv("limma_SOD1_vs_Control.tsv"),     show_col_types = FALSE))
  res_spor <- as.data.frame(read_tsv(find_tsv("limma_Sporadic_vs_Control.tsv"), show_col_types = FALSE))
  cat("rebuilt limma tables from TSV\n")
}
for (nm in c("res_c9","res_sod1","res_spor")) {
  d <- get(nm)
  need <- c("gene","logFC","t","P.Value","adj.P.Val")
  miss <- setdiff(need, names(d))
  if (length(miss))
    stop(nm, " is missing limma columns: ", paste(miss, collapse = ", "),
         "\n  columns present: ", paste(names(d), collapse = ", "),
         "\n  (if these look like log2FoldChange/padj you have loaded the RNA",
         " tables into the protein slots -- see the header)")
}
cat("protein limma tables:",
    paste(sapply(list(res_c9, res_sod1, res_spor), nrow), collapse = " / "), "genes\n")

## ---- route B: recompute fgsea exactly as proteomicsDE.R did ----------
if (!have_fg) {
  suppressPackageStartupMessages({ library(fgsea); library(msigdbr) })
  hall <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol)
  pathways <- split(hall$gene_symbol, hall$gs_name)
  do_fgsea <- function(res_tbl, label) {          # same as proteomicsDE.R
    ranks <- res_tbl$t; names(ranks) <- res_tbl$gene
    ranks <- sort(ranks, decreasing = TRUE)
    fg <- fgsea(pathways = pathways, stats = ranks, nperm = 10000)
    fg <- fg[order(fg$padj), ]; fg$comparison <- label; fg
  }
  fg_c9   <- do_fgsea(res_c9,   "C9orf72_vs_Control")
  fg_sod1 <- do_fgsea(res_sod1, "SOD1_vs_Control")
  fg_spor <- do_fgsea(res_spor, "Sporadic_vs_Control")
  cat("recomputed protein fgsea from the limma t statistic\n")
}

## ---- landmark check: the PROTEIN NES correlations ---------------------
## Verified 21 Aug run: 0.319 / -0.282 / 0.108
## (order: C9-Sporadic, C9-SOD1, SOD1-Sporadic)
nes <- function(f) { f <- as.data.frame(f); setNames(f$NES, f$pathway) }
a <- nes(fg_c9); b <- nes(fg_sod1); d <- nes(fg_spor)
sh <- Reduce(intersect, list(names(a), names(b), names(d)))
got <- c(cor(a[sh], d[sh], method = "spearman"),
         cor(a[sh], b[sh], method = "spearman"),
         cor(b[sh], d[sh], method = "spearman"))
cat(sprintf("\nprotein NES rho on %d pathways: %s\n", length(sh),
            paste(sprintf("%+.3f", got), collapse = " / ")))
cat("verified landmarks:                  +0.319 / -0.282 / +0.108\n")
if (!isTRUE(all.equal(got, c(0.319, -0.282, 0.108), tolerance = 0.03)))
  warning("protein NES correlations do not reproduce the verified run. ",
          "fgsea uses permutations; if the gap is small this is nperm noise, ",
          "if it is large the input tables differ.")

## ---- pin ---------------------------------------------------------------
if (!have_mats) {
  cat("\n*** prot_expr_collapsed / prot_gene are NOT available.\n")
  cat("*** They cannot be rebuilt from the TSVs -- source proteomicsDE.R once,\n")
  cat("*** then rerun this script to pin them alongside the DE tables.\n")
  cat("*** Script 16 will run up to section 12 and fail at the CCT2 rescue.\n")
} else {
  cat("\nexpression matrices present: prot_expr_collapsed ",
      paste(dim(prot_expr_collapsed), collapse = " x "),
      " | prot_gene ", paste(dim(prot_gene), collapse = " x "), "\n", sep = "")
  CCT2_ROW <- "1/sp|P78371|TCPB_HUMAN"
  cat("CCT2 rescue row present:", CCT2_ROW %in% rownames(prot_expr_collapsed), "\n")
  tric <- c("TCP1","CCT3","CCT4","CCT5","CCT6A","CCT7","CCT8")
  cat("TRiC subunits in prot_gene:", sum(tric %in% rownames(prot_gene)), "of 7\n")
}

payload <- list(res_c9 = res_c9, res_sod1 = res_sod1, res_spor = res_spor,
                fg_c9 = fg_c9, fg_sod1 = fg_sod1, fg_spor = fg_spor,
                prot_expr_collapsed = if (have_mats) prot_expr_collapsed else NULL,
                prot_gene           = if (have_mats) prot_gene else NULL,
                meta_prot2          = if (exists("meta_prot2")) meta_prot2 else NULL,
                nes_rho_C9_Spor_SOD1 = got,
                source = if (have_limma) "session (proteomicsDE.R)" else "TSV rebuild",
                has_matrices = have_mats,
                run_date = Sys.Date())
dir.create(PROT_DIR, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(FINAL)) {
  saveRDS(payload, FINAL); cat("\nSaved", FINAL, "\n")
} else {
  saveRDS(payload, sub("_FINAL", "_AUDITCHECK", FINAL))
  old <- readRDS(FINAL)
  cat("\nFINAL exists -- wrote the AUDITCHECK copy instead.\n")
  cat("  pinned NES rho:", paste(sprintf("%+.3f", old$nes_rho_C9_Spor_SOD1),
                                 collapse = " / "), "\n")
}
cat("\nNext: load the RNA side explicitly, then source 16.\n")
cat('  de <- readRDS("~/Documents/convergence_results/genotype_de_corrected_FINAL.rds")\n')
cat('  res_C9 <- de$res_C9; res_SOD1 <- de$res_SOD1; res_Sporadic <- de$res_Sporadic\n')