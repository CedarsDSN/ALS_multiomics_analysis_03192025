# =====================================================================
# 16a_preamble.R
#
# Supplies the six objects 16_convergence_corrected.R checks for:
#   res_c9, res_sod1, res_spor      <- genotype_de_corrected_FINAL.rds
#   fg_c9,  fg_sod1,  fg_spor       <- recomputed, or reused from the pin
#
# source() this, then source 16_convergence_corrected.R.
#
# The fgsea tables can come from either direction:
#   REUSE_FGSEA = TRUE   read them back from convergence_FINAL.rds
#                        (fast, and guarantees the figures match the pin,
#                        but is circular if you are rebuilding that file)
#   REUSE_FGSEA = FALSE  recompute from the DE tables with fgsea
#                        (the honest path when rerunning from scratch)
# =====================================================================

suppressPackageStartupMessages({ library(dplyr) })
set.seed(42)
REUSE_FGSEA <- TRUE   # default TRUE: the pinned tables are the VERIFIED ones
# (21 Aug run). Recomputation below does NOT reproduce
# them -- see the note at the landmark check.

find_file <- function(f, dirs = c(".", "~/Documents/DSN/ALS",
                                  "~/Documents/convergence_results",
                                  "~/Documents")) {
  for (d in path.expand(dirs)) {
    p <- file.path(d, f)
    if (file.exists(p)) return(p)
    if (dir.exists(d)) {
      g <- list.files(d, pattern = paste0("^", f, "$"),
                      recursive = TRUE, full.names = TRUE)
      if (length(g)) return(g[1])
    }
  }
  stop("not found: ", f)
}

## ---- the three DE tables ---------------------------------------------
de <- readRDS(find_file("genotype_de_corrected_FINAL.rds"))
cat("genotype_de_corrected_FINAL.rds elements:",
    paste(names(de), collapse = ", "), "\n")
pick <- function(patt) {
  i <- grep(patt, names(de), ignore.case = TRUE)[1]
  if (is.na(i)) stop("no element matching '", patt, "'; available: ",
                     paste(names(de), collapse = ", "))
  as.data.frame(de[[i]])
}
res_c9   <- pick("^res_C9$|C9")
res_sod1 <- pick("SOD1")
res_spor <- pick("Spor")

## ---- add the columns 16 expects on each DE frame ---------------------
## 16's prot_gene_table() requires a `gene` column; the pinned tables carry
## versioned ENSEMBL rownames only. Add `ens` (stripped) and `gene` (SYMBOL).
suppressPackageStartupMessages({ library(AnnotationDbi); library(org.Hs.eg.db) })
add_gene <- function(res) {
  res <- as.data.frame(res)
  res$ens <- sub("\\..*$", "", rownames(res))
  map <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(res$ens),
                               keytype = "ENSEMBL", columns = "SYMBOL")
  map <- map[!is.na(map$SYMBOL) & !duplicated(map$ENSEMBL), ]
  res$gene <- map$SYMBOL[match(res$ens, map$ENSEMBL)]
  res
}
res_c9 <- add_gene(res_c9); res_sod1 <- add_gene(res_sod1); res_spor <- add_gene(res_spor)
cat("gene symbols mapped:", sum(!is.na(res_c9$gene)), "of", nrow(res_c9), "\n")
cat("res_c9 / res_sod1 / res_spor:",
    paste(sapply(list(res_c9, res_sod1, res_spor), nrow), collapse = " / "),
    "genes\n")

## some versions of 16 want the *_corrected / *_filtered aliases too
res_C9_corrected <- res_C9_filtered <- res_c9
res_SOD1_corrected <- res_SOD1_filtered <- res_sod1
res_Sporadic_corrected <- res_Sporadic_filtered <- res_spor

## ---- the three fgsea tables -------------------------------------------
if (REUSE_FGSEA) {
  cv <- readRDS(find_file("convergence_FINAL.rds"))
  stopifnot(!is.null(cv$RNA_fgsea))
  fg_c9   <- as.data.frame(cv$RNA_fgsea$C9)
  fg_sod1 <- as.data.frame(cv$RNA_fgsea$SOD1)
  fg_spor <- as.data.frame(cv$RNA_fgsea$SP)
  cat("fgsea tables reused from convergence_FINAL.rds\n")
} else {
  suppressPackageStartupMessages({
    library(fgsea); library(msigdbr)
    library(AnnotationDbi); library(org.Hs.eg.db)
  })
  h <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>% distinct()
  hallmark_list <- split(h$gene_symbol, h$gs_name)
  cat("Hallmark sets:", length(hallmark_list), "\n")
  
  ## DESeq2 rownames are versioned ENSEMBL; rank on the Wald stat, one
  ## value per symbol (mean across duplicated mappings).
  rank_from_res <- function(res) {
    ens <- sub("\\..*$", "", rownames(res))
    map <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(ens),
                                 keytype = "ENSEMBL", columns = "SYMBOL")
    map <- map[!is.na(map$SYMBOL) & !duplicated(map$ENSEMBL), ]
    d <- data.frame(ens = ens, stat = res$stat, stringsAsFactors = FALSE)
    d <- merge(d, map, by.x = "ens", by.y = "ENSEMBL")
    d <- d[is.finite(d$stat), ]
    v <- tapply(d$stat, d$SYMBOL, mean)
    sort(v, decreasing = TRUE)
  }
  run_fg <- function(res, lab) {
    r <- rank_from_res(res)
    cat("  ", lab, ": ", length(r), " ranked symbols\n", sep = "")
    as.data.frame(fgsea(pathways = hallmark_list, stats = r))
  }
  cat("recomputing fgsea from the corrected DE tables:\n")
  fg_c9   <- run_fg(res_c9,   "C9")
  fg_sod1 <- run_fg(res_sod1, "SOD1")
  fg_spor <- run_fg(res_spor, "Sporadic")
}

## ---- landmark check: recomputed vs pinned ----------------------------
## RESOLVED 2026-08-26. The AUTHORITATIVE pathway NES correlations come
## from the verified 21 Aug run of script 16 on the corrected fit-3 layer:
##
##   RNA (iPSC-MN)   C9-Spor  0.804 | C9-SOD1 +0.170 | SOD1-Spor 0.348
##   Protein (DIA)   C9-Spor  0.319 | C9-SOD1 -0.282 | SOD1-Spor 0.108
##
## Recomputing fgsea here from the Wald stat over all 50 Hallmark sets
## gives 0.431 / -0.071 / 0.025, which does NOT reproduce them. Note the
## C9-SOD1 sign: -0.071 is the CONTAMINATED-layer sign, so the ranking
## used below is closer to the old layer than the corrected one. The
## ranking construction differs; the gene sets do not (50 matched).
##
## => REUSE_FGSEA is TRUE by default. Use the pinned tables. The
##    recomputation path is retained only as a documented mismatch and
##    must not be used to generate a reportable number.
##
## SWEEP NOTE for the writing pass: the pathway NES values above are a
## DIFFERENT quantity from the gene-level logFC correlations in Fig 2D
## (0.326 / 0.263 / 0.026). Same three pairs, different granularity, both
## changed from pre-correction. Anywhere the C9-SOD1 pathway comparison is
## quoted, the sign flipped from -0.08 to +0.17 and needs updating.
cmp_pinned <- function() {
  f <- try(find_file("convergence_FINAL.rds"), silent = TRUE)
  if (inherits(f, "try-error")) return(invisible(NULL))
  cv <- readRDS(f)
  if (is.null(cv$RNA_fgsea)) return(invisible(NULL))
  g <- function(x) { x <- as.data.frame(x); setNames(x$NES, x$pathway) }
  a <- g(cv$RNA_fgsea$C9); b <- g(cv$RNA_fgsea$SOD1); d <- g(cv$RNA_fgsea$SP)
  sh <- Reduce(intersect, list(names(a), names(b), names(d)))
  cat("\nPINNED fgsea (", length(sh), " pathways):\n", sep = "")
  cat(sprintf("  C9-Sporadic %.3f | C9-SOD1 %.3f | SOD1-Sporadic %.3f\n",
              cor(a[sh], d[sh], method = "spearman"),
              cor(a[sh], b[sh], method = "spearman"),
              cor(b[sh], d[sh], method = "spearman")))
  cat("  pinned table columns:", paste(names(as.data.frame(cv$RNA_fgsea$C9)),
                                       collapse = ", "), "\n")
  cat("  pinned n pathways per contrast:",
      paste(sapply(cv$RNA_fgsea, nrow), collapse = " / "), "\n")
}

## ---- landmark check: the NES convergence rho --------------------------
nes <- function(f) setNames(f$NES, f$pathway)
v9 <- nes(fg_c9); vs <- nes(fg_sod1); vp <- nes(fg_spor)
sh <- Reduce(intersect, list(names(v9), names(vs), names(vp)))
cat("\nNES rho on", length(sh), "shared pathways (manifest 0.804 / 0.170 / 0.348):\n")
cat(sprintf("  C9-Sporadic %.3f | C9-SOD1 %.3f | SOD1-Sporadic %.3f\n",
            cor(v9[sh], vp[sh], method = "spearman"),
            cor(v9[sh], vs[sh], method = "spearman"),
            cor(vs[sh], vp[sh], method = "spearman")))

cmp_pinned()
cat("\nIf the two blocks above disagree, the manifest value needs adjudicating\n")
cat("before Figure 2 panel d is finalised. Set REUSE_FGSEA <- TRUE to run 16\n")
cat("against the pinned tables in the meantime.\n")
cat("\nReady. Now source 16_convergence_corrected.R\n")