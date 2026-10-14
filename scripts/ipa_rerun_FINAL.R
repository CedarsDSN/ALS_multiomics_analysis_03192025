################################################################################
## ipa_rerun_FINAL.R
##
## Regenerates the S-ALS-3 vs S-ALS-2 DESeq2 contrast on the FINAL subtype
## labels and writes a fresh IPA upload file.
##
## The previous upload was built from cluster_final under the old labels,
## so every IPA result derived from it (including MYC activation at
## p = 1.20E-51) is stale.
##
## Design matches the original: ~ SEX + SV1..SV6 + subtype, sporadic only,
## S-ALS-2 as reference.
##
## Output: ipa_inputs_FINAL/DEG_S3_vs_S2_IPA_upload.tsv
##         ipa_inputs_FINAL/DEG_S3_vs_S2_full.tsv
##         ipa_inputs_FINAL/DEG_S3_vs_S2_full.rds
##
## Runtime: 10-20 min for the DESeq2 fit.
################################################################################

suppressPackageStartupMessages({
  library(DESeq2); library(dplyr)
  library(AnnotationDbi); library(org.Hs.eg.db)
})
set.seed(42)

OUT <- "ipa_inputs_FINAL"
dir.create(OUT, showWarnings = FALSE)

b  <- readRDS("backbone_FINAL.rds")
d  <- b$ax_round1
meta <- as.data.frame(b$meta)
if (is.null(d$Participant_ID))
  d$Participant_ID <- meta$Participant_ID[match(d$sample_id, meta$sample_id)]
d$subtype <- d$cluster_r1

## ---------------------------------------------------------------------------
## 1. sample set: sporadic, S-ALS-2 and S-ALS-3 only
## ---------------------------------------------------------------------------
keep <- d[d$subtype %in% c("S-ALS-2", "S-ALS-3"), ]
cat("samples in contrast:", nrow(keep), "\n")
print(table(keep$subtype))

stopifnot(exists("counts_filtered"))
cat("counts_filtered:", paste(dim(counts_filtered), collapse = " x "), "\n")

ids <- intersect(keep$sample_id, colnames(counts_filtered))
cat("matched to counts:", length(ids), "of", nrow(keep), "\n")
if (length(ids) < nrow(keep))
  warning("some samples missing from counts_filtered; check it is the corrected object")

cd <- meta[match(ids, meta$sample_id), , drop = FALSE]
cd$subtype <- droplevels(factor(keep$subtype[match(ids, keep$sample_id)],
                                levels = c("S-ALS-2", "S-ALS-3")))
cd$SEX <- factor(cd$SEX)
rownames(cd) <- ids
stopifnot(identical(rownames(cd), ids))

svs <- paste0("SV", 1:6)
stopifnot(all(svs %in% colnames(cd)))
stopifnot(all(is.finite(cd$SV1)))

## ---------------------------------------------------------------------------
## 2. DESeq2
## ---------------------------------------------------------------------------
design_f <- as.formula(paste("~ SEX +", paste(svs, collapse = " + "), "+ subtype"))
cat("\ndesign:", deparse(design_f), "\n")

dds <- DESeqDataSetFromMatrix(countData = counts_filtered[, ids, drop = FALSE],
                              colData = cd, design = design_f)
dds <- dds[rowSums(counts(dds) >= 10) >= round(0.1 * ncol(dds)), ]
cat("genes tested:", nrow(dds), "\n")

dds <- DESeq(dds)
nc <- sum(!mcols(dds)$betaConv, na.rm = TRUE)
if (nc > 0) {
  cat("non-converged:", nc, "- refitting with maxit 500\n")
  dds <- nbinomWaldTest(dds, maxit = 500)
  cat("still non-converged:", sum(!mcols(dds)$betaConv, na.rm = TRUE), "\n")
}

res <- results(dds, contrast = c("subtype", "S-ALS-3", "S-ALS-2"))
cat("\nsummary:\n"); summary(res)
cat("padj < 0.05:", sum(res$padj < 0.05, na.rm = TRUE), "\n")
cat("padj < 0.05 and |LFC| > 0.5:",
    sum(res$padj < 0.05 & abs(res$log2FoldChange) > 0.5, na.rm = TRUE), "\n")

## ---------------------------------------------------------------------------
## 3. annotate and write
## ---------------------------------------------------------------------------
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
cat("\nmapped to symbol:", sum(!is.na(rdf$symbol)), "of", nrow(rdf), "\n")

full <- rdf[, c("gene_id","ensembl","symbol","entrez","baseMean",
                "log2FoldChange","lfcSE","stat","pvalue","padj")]
full <- full[order(full$pvalue), ]

saveRDS(full, file.path(OUT, "DEG_S3_vs_S2_full.rds"))
write.table(full, file.path(OUT, "DEG_S3_vs_S2_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## IPA upload: one row per gene, identifier plus fold change and p-value
full <- readRDS(file.path(OUT, "DEG_S3_vs_S2_full.rds"))
u <- full[!is.na(full$symbol) & !is.na(full$pvalue), ]
u <- u[!duplicated(u$symbol), ]

up_ipa <- data.frame(
  ID                                    = u$symbol,
  `Expr.Log.Ratio`                      = u$log2FoldChange,
  `Expr.p.value`                        = u$pvalue,
  `Expr.False.Discovery.Rate..q.value.` = u$padj,
  check.names = FALSE, stringsAsFactors = FALSE)

write.table(up_ipa, file.path(OUT, "DEG_S3_vs_S2_IPA_upload.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", nrow(up_ipa), "genes\n")
print(head(up_ipa, 3))

## comparison against the previous upload
old <- read.delim(old_f, stringsAsFactors = FALSE)
m <- merge(up_ipa[, c("ID", "Expr.Log.Ratio")],
           old[,    c("ID", "Expr.Log.Ratio")],
           by = "ID", suffixes = c("_new", "_old"))
cat("shared genes:", nrow(m), "\n")
cat("Pearson :", round(cor(m$Expr.Log.Ratio_new, m$Expr.Log.Ratio_old,
                           use = "complete.obs"), 3), "\n")
cat("Spearman:", round(cor(m$Expr.Log.Ratio_new, m$Expr.Log.Ratio_old,
                           method = "spearman", use = "complete.obs"), 3), "\n")





library(DESeq2); library(clusterProfiler); library(org.Hs.eg.db)
res <- as.data.frame(results(dds, contrast = c("subtype", "S-ALS-3", "S-ALS-2")))
res <- res[!is.na(res$padj), ]
ids <- sub("\\..*$", "", rownames(res))                 # strip Ensembl versions if present
up  <- ids[res$padj < 0.001 & res$log2FoldChange >=  0.3]
dn  <- ids[res$padj < 0.001 & res$log2FoldChange <= -0.3]
cat("up:", length(up), " down:", length(dn), "\n")      # should be close to IPA's 771 / 2,094
go <- function(g) simplify(enrichGO(g, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
                                    universe = ids, pAdjustMethod = "BH"))
go_up <- go(up); go_dn <- go(dn)
head(as.data.frame(go_up)[, c("Description", "GeneRatio", "p.adjust")], 10)
head(as.data.frame(go_dn)[, c("Description", "GeneRatio", "p.adjust")], 10)
packageVersion("clusterProfiler")
