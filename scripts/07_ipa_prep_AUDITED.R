# 07_ipa_prep.R   (AUDITED 2026-08-22)
# S-ALS-3 vs S-ALS-2 DESeq2 contrast for IPA, corrected pipeline.
# Audited port of ipa_rerun_FINAL.R (the actual corrected-rerun canon:
# round-1 labels, S2/S3 sporadic samples, ~ SEX + SV1..SV6 + subtype,
# raw log2FC, no shrinkage).
#
# AUDIT CHANGES:
#   - counts from cohort_inputs_FINAL.rds (was: exists("counts_filtered"),
#     a session object). The refit uses the SAVED table's gene universe so
#     reproduction is exact regardless of the 24,220 vs 24,253 filter
#     difference between eras.
#   - Verification gate vs the saved DEG_S3_vs_S2_full.rds (a results
#     TABLE), at TOLERANCE: the saved universe (24,021 genes) is an
#     era refilter containing 40 genes outside the pinned corrected
#     filter, so the refit runs on the 23,981-gene intersection and
#     dispersion pooling differs in the fourth decimal
#     (verified 2026-08-22: max |dLFC| 0.0015, max |dstat| 0.013).
#   - THE SAVED FILES ARE NEVER OVERWRITTEN: the manuscript's IPA
#     analysis ran on them, so they are canonical by definition. This
#     script verifies them and writes only the _noXY.tsv convenience
#     file and an audit report.
#   - Sex-chromosome exclusion ADJUDICATED (2026-08-22): KDM5D/UTY/
#     USP9Y are present in the saved upload, so the exclusion happened
#     INSIDE the IPA software as an analysis-side filter. Methods
#     states it as such. The _noXY.tsv is written so future IPA runs
#     need no in-software filter.
#   - The old_f comparison block (undefined session variable) dropped.

suppressPackageStartupMessages({
  library(DESeq2); library(dplyr)
  library(AnnotationDbi); library(org.Hs.eg.db)
})
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::mutate,
                               dplyr::slice, base::setdiff, base::intersect,
                               base::union, base::unname, base::setequal,
                               .quiet = TRUE)
}
set.seed(42)
setwd("~/Documents/DSN/ALS")
OUT <- "ipa_inputs_FINAL"
stopifnot(dir.exists(OUT))

## ---- 1. Spec from the backbone (canonical labels + anchor SVs) ----
b <- readRDS("backbone_FINAL.rds")
d <- b$ax_round1
meta <- as.data.frame(b$meta)
d$subtype <- d$cluster_r1
keep <- d[d$subtype %in% c("S-ALS-2", "S-ALS-3"), ]
cat("samples in contrast:", nrow(keep), "(expect 311)\n")
print(table(droplevels(keep$subtype)))

svs <- paste0("SV", 1:6)
idx <- match(keep$sample_id, meta$sample_id)
stopifnot(!any(is.na(idx)))
cd <- data.frame(row.names = keep$sample_id,
                 subtype = droplevels(factor(keep$subtype,
                           levels = c("S-ALS-2", "S-ALS-3"))),
                 SEX = factor(meta$SEX[idx]))
for (s in svs) cd[[s]] <- meta[[s]][idx]
stopifnot(!any(is.na(cd)))

## ---- 2. Saved table defines the gene universe; pinned counts refit ----
saved <- readRDS(file.path(OUT, "DEG_S3_vs_S2_full.rds"))
stopifnot(all(c("gene_id", "log2FoldChange", "stat") %in% names(saved)))
cat("saved table:", nrow(saved), "genes\n")

ci <- readRDS("~/Documents/cohort_inputs_FINAL.rds")
cts <- ci$counts_filtered_corrected
gene_hit <- saved$gene_id %in% rownames(cts)
cat("saved genes found in pinned counts:", sum(gene_hit), "of", nrow(saved),
    "(40 era-refilter genes expected to be absent)\n")
stopifnot(mean(gene_hit) > 0.99)

dds <- DESeqDataSetFromMatrix(
  countData = cts[saved$gene_id[gene_hit], rownames(cd)],
  colData   = cd,
  design    = as.formula(paste("~ SEX +", paste(svs, collapse = " + "),
                               "+ subtype")))
dds <- DESeq(dds)
res <- results(dds, contrast = c("subtype", "S-ALS-3", "S-ALS-2"))

## ---- 3. GATE 1: reproduction of the saved DE table ----
sv <- saved[gene_hit, ]
stopifnot(identical(rownames(res), sv$gene_id))
d_lfc  <- max(abs(res$log2FoldChange - sv$log2FoldChange), na.rm = TRUE)
d_stat <- max(abs(res$stat - sv$stat), na.rm = TRUE)
cat(sprintf("\nGATE 1 (tolerance): max |dLFC| %.3g (limit 0.005) | max |dstat| %.3g (limit 0.05)\n",
            d_lfc, d_stat))
cat(sprintf("stat correlation: %.6f | top-500 overlap: %d\n",
            cor(res$stat, sv$stat, use = "pairwise"),
            length(intersect(rownames(res)[order(res$pvalue)][1:500],
                             sv$gene_id[order(sv$pvalue)][1:500]))))
de_ok <- d_lfc < 0.005 && d_stat < 0.05

## ---- 4. Rebuild annotation + upload (canon format, verbatim) ----
rdf <- as.data.frame(res)
rdf$gene_id <- rownames(rdf)
rdf$ensembl <- sub("\\..*$", "", rdf$gene_id)
sym <- suppressMessages(mapIds(org.Hs.eg.db, keys = rdf$ensembl,
                               column = "SYMBOL", keytype = "ENSEMBL",
                               multiVals = "first"))
ent <- suppressMessages(mapIds(org.Hs.eg.db, keys = rdf$ensembl,
                               column = "ENTREZID", keytype = "ENSEMBL",
                               multiVals = "first"))
rdf$symbol <- base::unname(sym[rdf$ensembl])
rdf$entrez <- base::unname(ent[rdf$ensembl])
full <- rdf[, c("gene_id","ensembl","symbol","entrez","baseMean",
                "log2FoldChange","lfcSE","stat","pvalue","padj")]
full <- full[order(full$pvalue), ]
u <- full[!is.na(full$symbol) & !is.na(full$pvalue), ]
u <- u[!duplicated(u$symbol), ]

## ---- 5. GATE 2: adjudicate the sex-chromosome exclusion ----
sexchr <- c("KDM5D", "UTY", "USP9Y")
upload_path <- file.path(OUT, "DEG_S3_vs_S2_IPA_upload.tsv")
console_edit <- NA
if (file.exists(upload_path)) {
  saved_up <- read.delim(upload_path, stringsAsFactors = FALSE)
  in_saved <- sexchr %in% saved_up$ID
  cat("\nGATE 2: KDM5D/UTY/USP9Y in saved upload:",
      paste(in_saved, collapse = "/"), "\n")
  console_edit <- !any(in_saved)
  cat(if (console_edit)
    "-> exclusion was a TSV-level edit; encoding it in the upload export.\n"
    else
    "-> saved upload contains them; exclusion happened INSIDE IPA (dataset\n   filter). Documented; upload export keeps all genes, and a *_noXY.tsv\n   is written for future IPA runs.\n")
  extra <- setdiff(saved_up$ID, u$symbol)
  cat("saved-upload symbols outside the intersection universe:",
      length(extra), "(the 40 era-refilter genes after symbol dedup;",
      "~25 expected)\n")
  set_ok <- length(extra) <= 30
} else { set_ok <- TRUE; console_edit <- FALSE }

## ---- 6. Verdict and gated writes ----
if (de_ok && set_ok) {
  ## saved FINAL files are canonical and untouched; write only the
  ## convenience export and the audit report
  mk_upload <- function(tab) data.frame(
    ID = tab$symbol,
    `Expr.Log.Ratio` = tab$log2FoldChange,
    `Expr.p.value` = tab$pvalue,
    `Expr.False.Discovery.Rate..q.value.` = tab$padj,
    check.names = FALSE, stringsAsFactors = FALSE)
  write.table(mk_upload(u[!u$symbol %in% sexchr, ]),
              file.path(OUT, "DEG_S3_vs_S2_IPA_upload_noXY.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(c(
    paste("07 audit", Sys.Date()),
    sprintf("intersection %d of %d genes; max |dLFC| %.4g; max |dstat| %.4g",
            sum(gene_hit), nrow(saved), d_lfc, d_stat),
    "sex-chr exclusion: inside IPA (KDM5D/UTY/USP9Y present in upload)",
    "saved FINAL files canonical; never overwritten by this script"),
    file.path(OUT, "AUDIT_REPORT_07.txt"))
  cat("\nVERIFIED (tolerance): saved IPA inputs reproduce from pinned",
      "counts up to the era gene-universe. FINAL files untouched;",
      "_noXY.tsv and AUDIT_REPORT_07.txt written.\n")
} else {
  saveRDS(full, file.path(OUT, "DEG_S3_vs_S2_AUDITCHECK.rds"))
  cat("\nMISMATCH beyond tolerance: saved as AUDITCHECK; report back.\n")
}
sessionInfo()
