# ============================================================
# Refactored entrypoint for original script: proteomicsDE.R
# Use --config, --base-dir, --outdir, and --objects. See scripts/README.md.
# ============================================================
# ---- ALS project bootstrap ----
.get_script_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- "--file="
  hit <- grep(file_arg, args, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub(file_arg, "", hit[[1]]), mustWork = FALSE)))
  getwd()
}
.find_project_root <- function(start) {
  path <- normalizePath(start, mustWork = FALSE)
  for (i in seq_len(8)) {
    if (file.exists(file.path(path, "R", "project_io.R"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }
  stop("Cannot find R/project_io.R. Run this script from the repository root.", call. = FALSE)
}
ALSPROJ_ROOT <- .find_project_root(.get_script_dir())
source(file.path(ALSPROJ_ROOT, "R", "project_io.R"))
setwd(ALSPROJ_ROOT)
ALS_ARGS <- als_parse_args()
if (isTRUE(ALS_ARGS$help)) { als_print_help("02_proteomics_de"); quit(save = "no", status = 0) }
CFG <- als_load_config(ALS_ARGS$config)
ALSPROJ_BASE_DIR <- als_choose(ALS_ARGS$base_dir, CFG, "project.base_dir", getwd())
ALSPROJ_OUTDIR <- als_choose(ALS_ARGS$outdir, CFG, "outputs.root", file.path(getwd(), "results", "02_proteomics_de"))
ALSPROJ_SCRIPT_ID <- "02_proteomics_de"
als_init(ALSPROJ_SCRIPT_ID, ALSPROJ_BASE_DIR, ALSPROJ_OUTDIR, ALS_ARGS$redirect_outputs)
als_register_output_wrappers()
als_load_objects(ALS_ARGS$objects)
on.exit(als_write_session_info(ALSPROJ_SCRIPT_ID), add = TRUE)
# ---- End bootstrap ----

library(arrow)
library(NMF)
library(DESeq2)
library(dplyr)
library(biomaRt)
library(data.table)
library(ggplot2)
# library(uwot)
# library(Rtsne)
library(BiocParallel)
register(MulticoreParam(20))

mainDir <- file.path(ALSPROJ_BASE_DIR, "proteomics/4_matrix")
setwd(mainDir)
prot_filename <- read.table("AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt", sep = "\t", header = TRUE)
prot_matrix <- as.data.frame(prot_filename)
sample_mapping <- read.csv("Sample_Mapping_Information/Sample_Mapping_File_Feb_2024.csv")

library(dplyr)
library(stringr)
library(tibble)

sample_cols <- setdiff(colnames(prot_matrix), c("Protein","nFragment","nPeptide"))

prot_colmap <- tibble(sample_col = sample_cols) %>%
  mutate(
    Participant.ID = str_replace(sample_col, "\\..*$", ""),  # before first dot
    Participant.ID = str_replace_all(Participant.ID, "_", "-"),
    case_ctrl_from_col = case_when(
      str_starts(sample_col, "CASE_") ~ "CASE",
      str_starts(sample_col, "CTRL_") ~ "CTRL",
      TRUE ~ NA_character_
    )
  )

# join to mapping at participant level
prot_colmap2 <- prot_colmap %>%
  left_join(sample_mapping, by = "Participant.ID")

prot_colmap2 %>% summarise(
  n_cols = n(),
  n_participant_mapped = sum(!is.na(Sample.ID)),
  n_participant_unmapped = sum(is.na(Sample.ID)),
  n_unique_participants = n_distinct(Participant.ID)
)

prot_expr <- prot_matrix %>%
  dplyr::select(Protein, all_of(sample_cols)) %>%
  column_to_rownames("Protein") %>%
  as.matrix()

# build participant groups (same order as columns)
participant_ids <- prot_colmap$Participant.ID
names(participant_ids) <- prot_colmap$sample_col  # map column -> participant

# collapse columns by participant
unique_participants <- unique(participant_ids)

prot_expr_collapsed <- sapply(unique_participants, function(pid) {
  cols <- names(participant_ids)[participant_ids == pid]
  mat <- prot_expr[, cols, drop = FALSE]
  if (ncol(mat) == 1) return(mat[,1])
  apply(mat, 1, median, na.rm = TRUE)
})

prot_expr_collapsed <- as.matrix(prot_expr_collapsed)  # proteins × participants

#!/usr/bin/env Rscript

# ============================================================
# AnswerALS Proteomics (Level 4, 436-P corrected+imputed)
# End-to-end script:
#   1) Load proteomics protein matrix (wide)
#   2) Collapse repeated proteomics runs to Participant_ID (median)
#   3) Build ALS/Control subtype metadata (Mutation_Subtype) from your existing pipeline inputs
#   4) Normalize proteomics (log2 + per-sample median-centering)
#   5) Map IDs:
#        - Guardrail Ensembl -> HGNC symbol via biomaRt
#        - Proteomics UniProt -> HGNC symbol via biomaRt
#   6) Guardrail score (UP–DOWN) at protein level
#   7) Differential abundance (limma) for:
#        - C9orf72 vs Control
#        - SOD1 vs Control
#        - Sporadic vs Control
#   8) Hallmark GSEA (fgsea) on limma t-stat ranks
#   9) Save outputs
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(biomaRt)
  library(limma)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
})

# ----------------------------
# User parameters (EDIT PATHS)
# ----------------------------

# Proteomics protein matrix (corrected+imputed)
prot_matrix_file <- "AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt"

# Your RNA metadata inputs (as in your snippet)
repeat_df_file <- file.path(ALSPROJ_BASE_DIR, "metadata/clinical/aals_dataportal_datatable.csv")
mutation_subtypes_file <- file.path(ALSPROJ_BASE_DIR, "ALS_Mutation_Subtypes_with_ALS_Control.csv")

# Your transcriptomics sample_map_tx object is assumed to come from a file; set this path or replace with your own load
# This should contain at least: Participant_ID, Sample_ID, SUBJECT_GROUP (or SUBJECT_GROUP.y after merges)
sample_map_tx <- read.csv(file.path(ALSPROJ_BASE_DIR, "transcriptomics/4_matrix/Sample_Mapping_Information/Sample_Mapping_File_Feb_2024.csv"))
colnames(sample_map_tx) <- gsub("\\.", "_", colnames(sample_map_tx))
sample_map_tx <- sample_map_tx %>%
  filter(grepl("-T$", Sample_ID)) %>%
  mutate(Sample_ID = trimws(Sample_ID),
         Participant_ID = trimws(Participant_ID))

# Optional: your RNA counts matrix (for the exact filtering step you used)
# If you already have "counts" in memory from another script, you can comment this and the filter will skip.
#counts_rds_file <- "data/transcriptomics/counts.rds"
counts <- read.csv(
  file.path(ALSPROJ_BASE_DIR, "transcriptomics/4_matrix/AnswerALS-651-T-v1-release6_raw-counts.csv"),
  check.names = FALSE
)
stopifnot("Geneid" %in% colnames(counts))
counts <- column_to_rownames(counts, var = "Geneid")
# keep only transcriptomics samples
counts <- counts[, grepl("-T$", colnames(counts)), drop=FALSE]
colnames(counts) <- trimws(colnames(counts))

# Guardrail gene lists (Ensembl IDs). Provide as plain text files, one ID per line.
guardrail_table_file <- "familial_signature_strict.csv"

# Output directory
outdir <- file.path(outdir, "proteomics_guardrail")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Helpers
# ----------------------------

strip_ens_ver <- function(x) sub("\\.\\d+$", "", x)

safe_factor2 <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- NA_character_
  factor(x)
}

write_tsv_safe <- function(df, file) {
  readr::write_tsv(df, file)
  message("Wrote: ", file)
}

get_participant_core <- function(x) {
  # 1) remove the run/suffix chunk like:
  #    -9014-P_D3, -5730-P-A01, .5730.P_A01, etc.
  core <- sub("([._-]\\d{3,5}[._-]?P[._-].*)$", "", x)
  
  # 2) standardize CASE_/CTRL_ (and CASE./CTRL.) to CASE-/CTRL-
  core <- sub("^(CASE|CTRL)[_.]", "\\1-", core)
  
  core
}

# ----------------------------
# 1) Load proteomics matrix (wide)
# ----------------------------

message("[1/9] Loading proteomics protein matrix...")
prot_matrix <- readr::read_tsv(prot_matrix_file, show_col_types = FALSE)

stopifnot("Protein" %in% colnames(prot_matrix))
prot_annot_cols <- intersect(c("nFragment", "nPeptide"), colnames(prot_matrix))

sample_cols <- setdiff(colnames(prot_matrix), c("Protein", prot_annot_cols))
stopifnot(length(sample_cols) > 10)

# Matrix: proteins x proteomics-run-columns
prot_expr <- prot_matrix %>%
  dplyr::select(Protein, all_of(sample_cols)) %>%
  column_to_rownames("Protein") %>%
  as.matrix()

storage.mode(prot_expr) <- "numeric"

# ----------------------------
# 2) Collapse proteomics runs -> Participant_ID (median)
# ----------------------------

message("[2/9] Collapsing proteomics columns to Participant_ID (median across repeats)...")

# Derive participant core from column name:
#   CASE_NEUAB000NKC.5730.P_A01  -> CASE-NEUAB000NKC
prot_colmap <- data.frame(
  sample_col = sample_cols,
  Participant_ID = get_participant_core(sample_cols),
  stringsAsFactors = FALSE
)

# sanity checks
head(prot_colmap, 10)
sum(is.na(prot_colmap$Participant_ID) | prot_colmap$Participant_ID == "")
head(prot_colmap$Participant_ID[!prot_colmap$Participant_ID %in% meta_prot$Participant_ID], 20)


if (any(is.na(prot_colmap$Participant_ID))) {
  bad <- prot_colmap %>% filter(is.na(Participant_ID)) %>% slice_head(n = 20)
  stop("Could not derive Participant_ID for some proteomics columns. Examples:\n",
       paste(bad$sample_col, collapse = "\n"))
}

participant_ids <- prot_colmap$Participant_ID
names(participant_ids) <- prot_colmap$sample_col
unique_participants <- unique(participant_ids)

prot_expr_collapsed <- sapply(unique_participants, function(pid) {
  cols <- names(participant_ids)[participant_ids == pid]
  mat <- prot_expr[, cols, drop = FALSE]
  if (ncol(mat) == 1) return(mat[, 1])
  apply(mat, 1, median, na.rm = TRUE)
})
prot_expr_collapsed <- as.matrix(prot_expr_collapsed) # proteins x participant cores


message("Proteomics collapsed dims: ", nrow(prot_expr_collapsed), " proteins x ", ncol(prot_expr_collapsed), " participants")

# ----------------------------
# 3) Build subtype metadata (Mutation_Subtype) using your pipeline inputs
# ----------------------------

message("[3/9] Building subtype metadata (Mutation_Subtype) from your existing RNA pipeline inputs...")

# Load sample_map_tx and repeat_df, then join + subtype join
# sample_map_tx <- readr::read_csv(sample_map_tx_file, show_col_types = FALSE) %>%
#   mutate(Participant_ID = trimws(as.character(Participant_ID)))

repeat_df <- readr::read_csv(repeat_df_file, show_col_types = FALSE) %>%
  mutate(Participant_ID = trimws(as.character(Participant_ID)))

annot <- dplyr::left_join(sample_map_tx, repeat_df, by = "Participant_ID")

# Optional: retain only samples with RNA counts present (mirrors your snippet)
# counts <- NULL
# if (file.exists(counts_rds_file)) {
#   message("Loading counts from: ", counts_rds_file)
#   counts <- readRDS(counts_rds_file)
#   stopifnot(is.matrix(counts) || is.data.frame(counts))
#   counts <- as.matrix(counts)
#   annot <- annot %>% filter(Sample_ID %in% colnames(counts))
# } else {
#   message("Counts file not found (", counts_rds_file, "). Skipping Sample_ID-in-counts filter.")
# }

counts <- as.matrix(counts)
annot <- annot %>% filter(Sample_ID %in% colnames(counts))
mutation_subtypes_with_ALS_control <- read.csv(mutation_subtypes_file, stringsAsFactors = FALSE)
annot_als_test <- annot %>%
  left_join(mutation_subtypes_with_ALS_control, by = "Participant_ID")

stopifnot(nrow(annot_als_test) > 0)

# Your earlier script used SUBJECT_GROUP.y after merges; handle both possibilities
subj_col <- if ("SUBJECT_GROUP.y" %in% colnames(annot_als_test)) "SUBJECT_GROUP.y" else if ("SUBJECT_GROUP" %in% colnames(annot_als_test)) "SUBJECT_GROUP" else NA_character_
stopifnot(!is.na(subj_col))

annot_als_filtered <- annot_als_test %>%
  filter(.data[[subj_col]] %in% c("ALS", "Healthy Control")) %>%
  filter(Mutation_Subtype %in% c("C9orf72", "SOD1", "Sporadic", "Control")) %>%
  droplevels()

# For proteomics, we want participant-level metadata
meta_prot <- annot_als_filtered %>%
  mutate(
    Participant_ID = as.character(Participant_ID),
    Mutation_Subtype = as.character(Mutation_Subtype),
    SEX = if ("SEX" %in% colnames(.)) as.character(SEX) else NA_character_
  ) %>%
  distinct(Participant_ID, .keep_all = TRUE)

# Restrict to participants present in proteomics
meta_prot2 <- meta_prot %>%
  filter(Participant_ID %in% colnames(prot_expr_collapsed)) %>%
  mutate(
    Mutation_Subtype = factor(Mutation_Subtype, levels = c("Control", "C9orf72", "SOD1", "Sporadic")),
    SEX = safe_factor2(SEX)
  )

message("Proteomics participants with Mutation_Subtype: ")
print(table(meta_prot2$Mutation_Subtype, useNA = "ifany"))

# Align proteomics matrix columns to meta order
prot_expr_collapsed2 <- prot_expr_collapsed[, meta_prot2$Participant_ID, drop = FALSE]
stopifnot(identical(colnames(prot_expr_collapsed2), meta_prot2$Participant_ID))

# ----------------------------
# 4) Normalize proteomics (log2 + sample median-centering) + filter
# ----------------------------

message("[4/9] Normalizing proteomics (log2 + per-sample median-centering) and filtering...")

min_pos <- min(prot_expr_collapsed2[prot_expr_collapsed2 > 0], na.rm = TRUE)
stopifnot(is.finite(min_pos))
offset <- min_pos / 2

prot_log2 <- log2(prot_expr_collapsed2 + offset)

# median-center each sample
prot_norm <- sweep(prot_log2, 2, apply(prot_log2, 2, median, na.rm = TRUE), FUN = "-")

# filter low-information proteins
miss_rate <- rowMeans(is.na(prot_norm))
row_sd <- apply(prot_norm, 1, sd, na.rm = TRUE)
keep <- (miss_rate <= 0.2) & (row_sd > 0.1)

prot_filt <- prot_norm[keep, , drop = FALSE]
message("Filtered proteomics dims: ", nrow(prot_filt), " proteins x ", ncol(prot_filt), " participants")

# ----------------------------
# 5) Map IDs via biomaRt:
#    - Guardrail Ensembl -> HGNC symbol
#    - UniProt -> HGNC symbol (proteomics)
# ----------------------------

message("[5/9] Mapping IDs using biomaRt (Ensembl + UniProt -> HGNC symbols)...")

# Read guardrail Ensembl lists
# Read guardrail signature table (single file with direction column)
guardrail_tbl <- readr::read_csv(guardrail_table_file, show_col_types = FALSE)

# Expect columns like:
# ensembl_id, direction (UP/DOWN), external_gene_name, gene_biotype, etc.
stopifnot(all(c("ensembl_id", "direction") %in% colnames(guardrail_tbl)))

guardrail_tbl <- guardrail_tbl %>%
  mutate(
    ensembl_id = strip_ens_ver(as.character(ensembl_id)),
    direction = toupper(trimws(as.character(direction)))
  ) %>%
  filter(!is.na(ensembl_id), ensembl_id != "", direction %in% c("UP", "DOWN"))

# Optional: keep protein-coding only (recommended if that was your original definition)
if ("gene_biotype" %in% colnames(guardrail_tbl)) {
  guardrail_tbl <- guardrail_tbl %>%
    filter(tolower(gene_biotype) == "protein_coding")
}

up_ens <- guardrail_tbl %>%
  filter(direction == "UP") %>%
  pull(ensembl_id) %>%
  unique()

down_ens <- guardrail_tbl %>%
  filter(direction == "DOWN") %>%
  pull(ensembl_id) %>%
  unique()

stopifnot(length(up_ens) > 0, length(down_ens) > 0)

message("Guardrail table loaded: UP Ensembl=", length(up_ens), " DOWN Ensembl=", length(down_ens))


mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Ensembl gene -> HGNC
ens_all <- unique(c(up_ens, down_ens))
ens_map <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ens_all,
  mart = mart
) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

up_symbols <- ens_map %>% filter(ensembl_gene_id %in% up_ens) %>% pull(hgnc_symbol) %>% unique()
down_symbols <- ens_map %>% filter(ensembl_gene_id %in% down_ens) %>% pull(hgnc_symbol) %>% unique()

message("Guardrail Ensembl mapped to HGNC: UP=", length(up_symbols), " DOWN=", length(down_symbols))

# UniProt -> HGNC for proteomics features
# UniProt -> HGNC for proteomics features (your IDs look like: 1/sp|A0AVT1|UBA6_HUMAN)

uniprot_raw <- rownames(prot_filt)

extract_uniprot <- function(x) {
  x <- as.character(x)
  
  # remove leading "1/" or anything like "<digit>/"
  x <- sub("^\\d+/", "", x)
  
  # if "sp|ACC|NAME" or "tr|ACC|NAME", extract ACC
  x <- ifelse(grepl("\\|", x), sub("^.*\\|([^|]+)\\|.*$", "\\1", x), x)
  
  # keep first if multiple separated by ; or ,
  x <- sub("[;,].*$", "", x)
  
  # drop isoform suffix -2, -3, ...
  x <- sub("-\\d+$", "", x)
  
  trimws(x)
}

uniprot_ids_clean <- extract_uniprot(uniprot_raw)

# quick sanity
print(head(data.frame(raw = uniprot_raw, clean = uniprot_ids_clean), 10))
cat("Empty cleaned IDs:", sum(is.na(uniprot_ids_clean) | uniprot_ids_clean == ""), "\n")

uniprot_map <- biomaRt::getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "uniprotswissprot",
  values = unique(uniprot_ids_clean),
  mart = mart
) %>%
  dplyr::filter(!is.na(hgnc_symbol), hgnc_symbol != "",
                !is.na(uniprotswissprot), uniprotswissprot != "") %>%
  dplyr::distinct(uniprotswissprot, .keep_all = TRUE)

idx <- match(uniprot_ids_clean, uniprot_map$uniprotswissprot)
keep_map <- which(!is.na(idx))

cat("Mapped proteins:", length(keep_map), "of", length(uniprot_ids_clean), "\n")
stopifnot(length(keep_map) > 0)

prot_gene <- prot_filt[keep_map, , drop = FALSE]
rownames(prot_gene) <- uniprot_map$hgnc_symbol[idx[keep_map]]

# Collapse duplicated gene symbols by median (robust)
# prot_gene currently: rows = (HGNC symbols, possibly duplicated), cols = participants
stopifnot(!is.null(colnames(prot_gene)))
cn <- colnames(prot_gene)

genes <- unique(rownames(prot_gene))

prot_gene_collapsed <- matrix(NA_real_, nrow = length(genes), ncol = length(cn),
                              dimnames = list(genes, cn))

for (g in genes) {
  mat <- prot_gene[rownames(prot_gene) == g, , drop = FALSE]
  if (nrow(mat) == 1) {
    prot_gene_collapsed[g, ] <- as.numeric(mat[1, ])
  } else {
    prot_gene_collapsed[g, ] <- apply(mat, 2, median, na.rm = TRUE)
  }
}

prot_gene <- prot_gene_collapsed

rm(prot_gene_collapsed)


cat("Gene-symbol matrix dims:", nrow(prot_gene), "x", ncol(prot_gene), "\n")


message("Gene-symbol proteomics dims: ", nrow(prot_gene), " genes x ", ncol(prot_gene), " participants")

# ----------------------------
# 6) Guardrail score (protein) - don't do this there is no overlap between proteomic genes and guardrail genes
# ----------------------------

message("[6/9] Computing proteomics guardrail score...") - 

up_in <- intersect(up_symbols, rownames(prot_gene))
down_in <- intersect(down_symbols, rownames(prot_gene))
message("Guardrail genes present in proteomics: UP=", length(up_in), " DOWN=", length(down_in))

stopifnot(length(up_in) >= 5, length(down_in) >= 5)

score_prot <- colMeans(prot_gene[up_in, , drop = FALSE], na.rm = TRUE) -
  colMeans(prot_gene[down_in, , drop = FALSE], na.rm = TRUE)

score_df <- tibble(
  Participant_ID = names(score_prot),
  score_protein = as.numeric(score_prot)
) %>%
  left_join(meta_prot2 %>% select(Participant_ID, Mutation_Subtype, SEX), by = "Participant_ID")

write_tsv_safe(score_df, file.path(outdir, "proteomics_guardrail_scores.tsv"))

# Pairwise Wilcoxon tests vs Control
pairwise_test <- function(df, case_label) {
  sub <- df %>% filter(Mutation_Subtype %in% c("Control", case_label))
  if (nrow(sub) < 3) return(NULL)
  wilcox.test(score_protein ~ Mutation_Subtype, data = sub)
}

wt_c9 <- pairwise_test(score_df, "C9orf72")
wt_sod1 <- pairwise_test(score_df, "SOD1")
wt_spor <- pairwise_test(score_df, "Sporadic")

message("Wilcoxon C9orf72 vs Control:")
print(wt_c9)
message("Wilcoxon SOD1 vs Control:")
print(wt_sod1)
message("Wilcoxon Sporadic vs Control:")
print(wt_spor)

# Plot guardrail scores
p1 <- ggplot(score_df, aes(x = Mutation_Subtype, y = score_protein)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
  theme_bw() +
  labs(title = "Proteomics guardrail score (UP - DOWN)", x = NULL, y = "Score (log2-centered)")

ggsave(file.path(outdir, "guardrail_score_boxplot.png"), p1, width = 7, height = 4, dpi = 200)

# ----------------------------
# 7) Differential abundance (limma) per subtype vs Control
# ----------------------------

message("[7/9] Running limma differential abundance (subtype vs Control)...")

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(readr)
})

stopifnot(identical(colnames(prot_gene), meta_prot2$Participant_ID))

run_limma_subtype <- function(expr_gene, meta_tbl, case_label, ctrl_label = "Control") {
  idx <- meta_tbl$Mutation_Subtype %in% c(ctrl_label, case_label)
  meta_sub <- meta_tbl[idx, , drop = FALSE]
  expr_sub <- expr_gene[, meta_sub$Participant_ID, drop = FALSE]
  
  meta_sub$Mutation_Subtype <- factor(as.character(meta_sub$Mutation_Subtype),
                                      levels = c(ctrl_label, case_label))
  
  use_sex <- ("SEX" %in% names(meta_sub)) && (nlevels(droplevels(meta_sub$SEX)) >= 2)
  
  design <- if (use_sex) {
    model.matrix(~ 0 + Mutation_Subtype + SEX, data = meta_sub)
  } else {
    model.matrix(~ 0 + Mutation_Subtype, data = meta_sub)
  }
  colnames(design) <- make.names(colnames(design))
  
  fit <- lmFit(expr_sub, design)
  
  coef_case <- make.names(paste0("Mutation_Subtype", case_label))
  coef_ctrl <- make.names(paste0("Mutation_Subtype", ctrl_label))
  
  # makeContrasts wants a string expression using column names
  contrast_str <- paste0(coef_case, " - ", coef_ctrl)
  cont <- makeContrasts(contrasts = contrast_str, levels = design)
  
  fit2 <- eBayes(contrasts.fit(fit, cont))
  
  res <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")
  res$gene <- rownames(res)
  res
}


res_c9   <- run_limma_subtype(prot_gene, meta_prot2, "C9orf72", "Control")
res_sod1 <- run_limma_subtype(prot_gene, meta_prot2, "SOD1", "Control")
res_spor <- run_limma_subtype(prot_gene, meta_prot2, "Sporadic", "Control")

write_tsv(res_c9,   file.path(outdir, "limma_C9orf72_vs_Control.tsv"))
write_tsv(res_sod1, file.path(outdir, "limma_SOD1_vs_Control.tsv"))
write_tsv(res_spor, file.path(outdir, "limma_Sporadic_vs_Control.tsv"))

count_sig <- function(res, alpha = 0.05) sum(res$adj.P.Val <= alpha, na.rm = TRUE)
cat("Significant genes (FDR<=0.05): C9 =", count_sig(res_c9),
    "SOD1 =", count_sig(res_sod1),
    "Sporadic =", count_sig(res_spor), "\n")

head(res_c9[order(res_c9$adj.P.Val), c("gene","logFC","adj.P.Val","t")], 20)

# ----------------------------
# 8) Hallmark GSEA (fgsea) per comparison
# ----------------------------

message("[8/9] Running Hallmark GSEA (fgsea) using limma t-stat ranks...")

hall <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
pathways <- split(hall$gene_symbol, hall$gs_name)

do_fgsea <- function(res_tbl, label) {
  if (is.null(res_tbl) || nrow(res_tbl) == 0) return(NULL)
  ranks <- res_tbl$t
  names(ranks) <- res_tbl$gene
  ranks <- sort(ranks, decreasing = TRUE)
  
  fg <- fgsea(pathways = pathways, stats = ranks, nperm = 10000)
  fg <- fg[order(fg$padj), ]
  fg$comparison <- label
  fg
}

fg_c9 <- do_fgsea(res_c9, "C9orf72_vs_Control")
fg_sod1 <- do_fgsea(res_sod1, "SOD1_vs_Control")
fg_spor <- do_fgsea(res_spor, "Sporadic_vs_Control")

fg_all <- bind_rows(fg_c9, fg_sod1, fg_spor)
if (nrow(fg_all) > 0) {
  write_tsv_safe(fg_all, file.path(outdir, "fgsea_hallmark_all.tsv"))
}

# Save top 20 per comparison
if (nrow(fg_all) > 0) {
  top20 <- fg_all %>%
    group_by(comparison) %>%
    arrange(padj, .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    ungroup()
  write_tsv_safe(top20, file.path(outdir, "fgsea_hallmark_top20_per_comparison.tsv"))
}

# ----------------------------
# 9) Quick PCA QC (optional but recommended)
# ----------------------------

message("[9/9] PCA QC on normalized proteomics (gene-symbol matrix)...")

# PCA on filtered gene-symbol matrix (transposed: samples x genes)
pca_mat <- t(prot_gene)
# Remove zero-variance columns if any
sd_cols <- apply(pca_mat, 2, sd, na.rm = TRUE)
pca_mat2 <- pca_mat[, sd_cols > 0, drop = FALSE]

pca <- prcomp(pca_mat2, center = TRUE, scale. = TRUE)
pcs <- as.data.frame(pca$x[, 1:5, drop = FALSE])
pcs$Participant_ID <- rownames(pcs)

pcs <- pcs %>%
  left_join(meta_prot2 %>% select(Participant_ID, Mutation_Subtype, SEX), by = "Participant_ID")

write_tsv_safe(pcs, file.path(outdir, "proteomics_pcs.tsv"))

p_pca <- ggplot(pcs, aes(x = PC1, y = PC2, color = Mutation_Subtype)) +
  geom_point(alpha = 0.8) +
  labs(title = "Proteomics PCA (normalized, filtered)", x = "PC1", y = "PC2")

ggsave(file.path(outdir, "proteomics_PCA_PC1_PC2.png"), p_pca, width = 6, height = 4, dpi = 200)

message("Done. Outputs are in: ", outdir)

p_pca <- ggplot(pcs, aes(x = PC1, y = PC2, color = Batch.x)) +
  geom_point(alpha = 0.8) +
  labs(title = "Proteomics PCA (normalized, filtered)", x = "PC1", y = "PC2")

ggsave(file.path(outdir, "proteomics_PCA_PC1_PC2_colored_by_batch.png"), p_pca, width = 6, height = 4, dpi = 200)


pca <- prcomp(t(prot_gene), center = TRUE, scale. = TRUE)
pcs <- as.data.frame(pca$x[, 1:10, drop = FALSE])
pcs$Participant_ID <- rownames(pca$x)

pcs <- dplyr::left_join(
  pcs,
  meta_prot2,
  by = "Participant_ID"
)


# continuous vars you might have
cont_vars <- c("AGE_AT_SYMPTOM_ONSET", "Age_at_First_PBMC_Collection",
               "ALSFRS_R_Baseline_Value", "ALSFRS_R_PROGRESSION_SLOPE")

# categorical vars you might have
cat_vars <- c("SEX", "Site", "Batch.x", "Shipment", "External_Control")

pc_names <- paste0("PC", 1:5)

pc_assoc <- list()

for (pc in pc_names) {
  # continuous
  for (v in cont_vars) {
    if (!v %in% colnames(pcs)) next
    x <- pcs[[v]]
    y <- pcs[[pc]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) next
    ct <- cor.test(x[ok], y[ok], method = "spearman")
    pc_assoc[[length(pc_assoc)+1]] <- data.frame(
      PC = pc, var = v, type = "continuous",
      stat = unname(ct$estimate), p = ct$p.value
    )
  }
  
  # categorical
  for (v in cat_vars) {
    if (!v %in% colnames(pcs)) next
    f <- droplevels(as.factor(pcs[[v]]))
    if (nlevels(f) < 2) next
    df <- data.frame(y = pcs[[pc]], f = f)
    a <- anova(lm(y ~ f, data = df))
    pc_assoc[[length(pc_assoc)+1]] <- data.frame(
      PC = pc, var = v, type = "categorical",
      stat = NA_real_, p = a$`Pr(>F)`[1]
    )
  }
}

pc_assoc_df <- dplyr::bind_rows(pc_assoc) %>%
  dplyr::mutate(p_adj = p.adjust(p, method = "BH")) %>%
  dplyr::arrange(p_adj)

pc_assoc_df %>% head(30)

ggplot(pcs, aes(PC1, PC2, color = Mutation_Subtype)) +
  geom_point(alpha = 0.8) +
  stat_ellipse(aes(group = Mutation_Subtype), type = "t", linetype = 2)


pairs_df <- expand.grid(
  x = paste0("PC", 1:5),
  y = paste0("PC", 1:5),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(x < y)


for (i in seq_len(nrow(pairs_df))) {
  x <- pairs_df$x[i]; y <- pairs_df$y[i]
  print(
    ggplot(pcs, aes(.data[[x]], .data[[y]], shape = Mutation_Subtype)) +
      geom_point(alpha = 0.7) +
      theme_bw() +
      labs(title = paste(x, "vs", y))
  )
}

library(limma)
library(dplyr)
library(ggplot2)

stopifnot(identical(colnames(prot_gene), meta_prot2$Participant_ID))

# regress out SEX only (since that's what PCs flagged)
design_sex <- model.matrix(~ SEX, data = meta_prot2)
prot_resid <- removeBatchEffect(prot_gene, covariates = design_sex[, -1, drop = FALSE])

pca_res <- prcomp(t(prot_resid), center = TRUE, scale. = TRUE)
pcs_res <- as.data.frame(pca_res$x[, 1:5, drop = FALSE])
pcs_res$Participant_ID <- rownames(pca_res$x)

pcs_res <- left_join(pcs_res, meta_prot2 %>% dplyr::select(Participant_ID, Mutation_Subtype), by = "Participant_ID")

p_pca_sex <- ggplot(pcs_res, aes(PC1, PC2, color = Mutation_Subtype)) +
  geom_point(alpha = 0.8) +
  labs(title = "Proteomics PCA after regressing SEX")
ggsave(file.path(outdir, "proteomics_PCA_PC1_PC2_after_reg_sex.png"), p_pca_sex, width = 6, height = 4, dpi = 200)

summary(lm(PC1 ~ Mutation_Subtype, data = pcs_res))
summary(lm(PC2 ~ Mutation_Subtype, data = pcs_res))


pcs_res <- pcs_res %>%
  mutate(dist12 = sqrt(PC1^2 + PC2^2)) %>%
  arrange(desc(dist12))

head(pcs_res %>% dplyr::select(Participant_ID, Mutation_Subtype, PC1, PC2, dist12), 10)

table(meta_prot2$Mutation_Subtype, meta_prot2$SEX)
if ("Site" %in% colnames(meta_prot2)) table(meta_prot2$Mutation_Subtype, meta_prot2$Site)
if ("Batch.x" %in% colnames(meta_prot2)) table(meta_prot2$Mutation_Subtype, meta_prot2$Batch.x)
if ("Shipment" %in% colnames(meta_prot2)) table(meta_prot2$Mutation_Subtype, meta_prot2$Shipment)


run_limma_subtype <- function(expr_gene, meta_tbl, case_label, ctrl_label = "Control") {
  idx <- meta_tbl$Mutation_Subtype %in% c(ctrl_label, case_label)
  meta_sub <- meta_tbl[idx, , drop = FALSE]
  expr_sub <- expr_gene[, meta_sub$Participant_ID, drop = FALSE]
  
  meta_sub$Mutation_Subtype <- factor(as.character(meta_sub$Mutation_Subtype),
                                      levels = c(ctrl_label, case_label))
  
  use_sex <- ("SEX" %in% names(meta_sub)) && (nlevels(droplevels(meta_sub$SEX)) >= 2)
  
  design <- {
    terms <- c("0 + Mutation_Subtype")
    
    if ("SEX" %in% names(meta_sub) && nlevels(droplevels(meta_sub$SEX)) >= 2)
      terms <- c(terms, "SEX")
    
    if ("Shipment" %in% names(meta_sub) && nlevels(droplevels(as.factor(meta_sub$Shipment))) >= 2)
      terms <- c(terms, "as.factor(Shipment)")
    
    model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), data = meta_sub)
  }
  colnames(design) <- make.names(colnames(design))
  
  fit <- lmFit(expr_sub, design)
  
  coef_case <- make.names(paste0("Mutation_Subtype", case_label))
  coef_ctrl <- make.names(paste0("Mutation_Subtype", ctrl_label))
  
  # makeContrasts wants a string expression using column names
  contrast_str <- paste0(coef_case, " - ", coef_ctrl)
  cont <- makeContrasts(contrasts = contrast_str, levels = design)
  
  fit2 <- eBayes(contrasts.fit(fit, cont))
  
  res <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")
  res$gene <- rownames(res)
  res
}

res_c9_ship   <- run_limma_subtype(prot_gene, meta_prot2, "C9orf72", "Control")
res_sod1_ship <- run_limma_subtype(prot_gene, meta_prot2, "SOD1", "Control")
res_spor_ship <- run_limma_subtype(prot_gene, meta_prot2, "Sporadic", "Control")

min(res_c9_ship$adj.P.Val, na.rm=TRUE)
head(res_c9_ship[order(res_c9_ship$P.Value), c("gene","logFC","P.Value","adj.P.Val","t")], 10)
count_sig <- function(res, alpha = 0.05) sum(res$adj.P.Val <= alpha, na.rm = TRUE)
cat("Significant genes (FDR<=0.05): C9 =", count_sig(res_c9_ship),
    "SOD1 =", count_sig(res_sod1_ship),
    "Sporadic =", count_sig(res_spor_ship), "\n")


# using your pcs_res with dist12 already computed
keep_ids <- pcs_res %>%
  mutate(cut = quantile(dist12, 0.95, na.rm=TRUE)) %>%
  filter(dist12 <= cut) %>%
  pull(Participant_ID)

prot_gene_no <- prot_gene[, keep_ids, drop=FALSE]
meta_no <- meta_prot2 %>% filter(Participant_ID %in% keep_ids)
meta_no <- meta_no[match(colnames(prot_gene_no), meta_no$Participant_ID), ]
stopifnot(identical(colnames(prot_gene_no), meta_no$Participant_ID))

res_c9_no <- run_limma_subtype(prot_gene_no, meta_no, "C9orf72", "Control")
min(res_c9_no$adj.P.Val, na.rm=TRUE)


prot_C9_sig <- fg_c9 %>%
  transmute(
    ID = pathway,
    NES_Prot = NES,
    FDR_Prot = padj,
    sig_Prot = padj < 0.05
  )

C9_join <- C9_sig %>% inner_join(prot_C9_sig, by = "ID")

cor.test(C9_join$NES_C9, C9_join$NES_Prot, method = "spearman")
cor.test(C9_join$NES_C9, C9_join$NES_Prot, method = "pearson")

C9_join %>%
  filter(sig_C9 | sig_Prot) %>%
  arrange(FDR_C9, FDR_Prot) %>%
  dplyr::select(ID, NES_C9, FDR_C9, sig_C9, NES_Prot, FDR_Prot, sig_Prot)

meta_c9_ctrl <- meta_prot2 %>%
  dplyr::filter(Mutation_Subtype %in% c("Control", "C9orf72")) %>%
  droplevels()

gene <- "CBX1"

expr_df <- data.frame(
  Participant_ID = colnames(prot_gene),
  expr = as.numeric(prot_gene[gene, ]),
  stringsAsFactors = FALSE
)

expr_df <- expr_df %>%
  dplyr::inner_join(meta_c9_ctrl, by = "Participant_ID")
pdf(file.path(outdir, "CBX1_protein_abundance_boxplot_C9_vs_Control.pdf"), width = 5, height = 4)
boxplot(
  expr ~ Mutation_Subtype,
  data = expr_df,
  main = "CBX1 protein abundance",
  xlab = "Mutation Subtype",
  ylab = "Normalized protein expression",
  col = c("grey80", "grey60"),
  outline = TRUE
)
dev.off()

wilcox.test(expr ~ Mutation_Subtype, data = expr_df)

tapply(df$expr[df$Mutation_Subtype %in% c("Control","C9orf72")],
       df$Mutation_Subtype[df$Mutation_Subtype %in% c("Control","C9orf72")], mean, na.rm=TRUE)



library(dplyr)

# pick your threshold(s)
alpha <- 0.10   # use 0.10 given low power; also try 0.05

c9_sig <- fg_c9 %>%
  filter(!is.na(padj), padj < alpha) %>%
  transmute(pathway, NES_C9 = NES, FDR_C9 = padj,
            dir_C9 = ifelse(NES > 0, "Up_in_C9", "Down_in_C9"))

sod1_sig <- fg_sod1 %>%
  filter(!is.na(padj), padj < alpha) %>%
  transmute(pathway, NES_SOD1 = NES, FDR_SOD1 = padj,
            dir_SOD1 = ifelse(NES > 0, "Up_in_SOD1", "Down_in_SOD1"))

overlap_sig <- inner_join(c9_sig, sod1_sig, by = "pathway") %>%
  mutate(direction_agree = (sign(NES_C9) == sign(NES_SOD1))) %>%
  arrange(pmax(FDR_C9, FDR_SOD1), desc(abs(NES_C9) + abs(NES_SOD1)))

# Summary counts
n_c9 <- nrow(c9_sig)
n_sod1 <- nrow(sod1_sig)
n_overlap <- nrow(overlap_sig)
jaccard <- ifelse((n_c9 + n_sod1 - n_overlap) > 0,
                  n_overlap / (n_c9 + n_sod1 - n_overlap), NA_real_)

cat("alpha =", alpha, "\n")
cat("C9 sig:", n_c9, " | SOD1 sig:", n_sod1, " | overlap:", n_overlap, "\n")
cat("Jaccard overlap:", round(jaccard, 3), "\n")
cat("Direction agree in overlap:", sum(overlap_sig$direction_agree), "/", n_overlap, "\n")

# Top overlapped pathways
print(head(overlap_sig, 20))

# Join all pathways by name, then correlate NES across contrasts
c9_all <- fg_c9 %>% transmute(pathway, NES_C9 = NES, padj_C9 = padj)
sod1_all <- fg_sod1 %>% transmute(pathway, NES_SOD1 = NES, padj_SOD1 = padj)

nes_join <- full_join(c9_all, sod1_all, by = "pathway") %>%
  filter(is.finite(NES_C9), is.finite(NES_SOD1)) %>%
  mutate(direction_agree = (sign(NES_C9) == sign(NES_SOD1)))

# Correlations
spearman_nes <- cor(nes_join$NES_C9, nes_join$NES_SOD1, method = "spearman")
pearson_nes  <- cor(nes_join$NES_C9, nes_join$NES_SOD1, method = "pearson")

cat("NES correlation (Spearman):", round(spearman_nes, 3), "\n")
cat("NES correlation (Pearson): ", round(pearson_nes, 3), "\n")
cat("Direction agreement:", round(mean(nes_join$direction_agree), 3), "\n")

# Top concordant (large |NES| in both, same direction)
top_concordant <- nes_join %>%
  filter(direction_agree) %>%
  mutate(score = abs(NES_C9) + abs(NES_SOD1)) %>%
  arrange(desc(score))

print(head(top_concordant, 20))

# Top discordant (opposite direction)
top_discordant <- nes_join %>%
  filter(!direction_agree) %>%
  mutate(score = abs(NES_C9) + abs(NES_SOD1)) %>%
  arrange(desc(score))

print(head(top_discordant, 20))

topK <- 15

c9_top <- fg_c9 %>% arrange(padj, desc(abs(NES))) %>% slice(1:topK) %>% dplyr::select(pathway, NES_C9 = NES)
sod1_top <- fg_sod1 %>% arrange(padj, desc(abs(NES))) %>% slice(1:topK) %>% dplyr::select(pathway, NES_SOD1 = NES)

topK_overlap <- inner_join(c9_top, sod1_top, by = "pathway") %>%
  mutate(direction_agree = (sign(NES_C9) == sign(NES_SOD1))) %>%
  arrange(desc(abs(NES_C9) + abs(NES_SOD1)))

cat("Top", topK, "overlap n =", nrow(topK_overlap), "\n")
print(topK_overlap)


alpha <- 0.10

c9_sig <- fg_c9 %>%
  dplyr::filter(!is.na(padj), padj < alpha) %>%
  dplyr::transmute(pathway, NES_C9 = NES, FDR_C9 = padj,
                   dir_C9 = ifelse(NES > 0, "Up_in_C9", "Down_in_C9"))

spor_sig <- fg_spor %>%
  dplyr::filter(!is.na(padj), padj < alpha) %>%
  dplyr::transmute(pathway, NES_Spor = NES, FDR_Spor = padj,
                   dir_Spor = ifelse(NES > 0, "Up_in_Spor", "Down_in_Spor"))

overlap_sig_c9_spor <- dplyr::inner_join(c9_sig, spor_sig, by = "pathway") %>%
  dplyr::mutate(direction_agree = (sign(NES_C9) == sign(NES_Spor))) %>%
  dplyr::arrange(pmax(FDR_C9, FDR_Spor), dplyr::desc(abs(NES_C9) + abs(NES_Spor)))

n_c9   <- nrow(c9_sig)
n_spor <- nrow(spor_sig)
n_ov   <- nrow(overlap_sig_c9_spor)
jacc   <- ifelse((n_c9 + n_spor - n_ov) > 0, n_ov / (n_c9 + n_spor - n_ov), NA_real_)

cat("alpha =", alpha, "\n")
cat("C9 sig:", n_c9, " | Sporadic sig:", n_spor, " | overlap:", n_ov, "\n")
cat("Jaccard overlap:", round(jacc, 3), "\n")
cat("Direction agree in overlap:", sum(overlap_sig_c9_spor$direction_agree), "/", n_ov, "\n")

print(head(overlap_sig_c9_spor, 20))

c9_all <- fg_c9 %>% dplyr::transmute(pathway, NES_C9 = NES, padj_C9 = padj)
sp_all <- fg_spor %>% dplyr::transmute(pathway, NES_Spor = NES, padj_Spor = padj)

nes_join_c9_spor <- dplyr::full_join(c9_all, sp_all, by = "pathway") %>%
  dplyr::filter(is.finite(NES_C9), is.finite(NES_Spor)) %>%
  dplyr::mutate(direction_agree = (sign(NES_C9) == sign(NES_Spor)))

spearman_nes <- cor(nes_join_c9_spor$NES_C9, nes_join_c9_spor$NES_Spor, method = "spearman")
pearson_nes  <- cor(nes_join_c9_spor$NES_C9, nes_join_c9_spor$NES_Spor, method = "pearson")

cat("NES correlation (Spearman):", round(spearman_nes, 3), "\n")
cat("NES correlation (Pearson): ", round(pearson_nes, 3), "\n")
cat("Direction agreement:", round(mean(nes_join_c9_spor$direction_agree), 3), "\n")

top_concordant <- nes_join_c9_spor %>%
  dplyr::filter(direction_agree) %>%
  dplyr::mutate(score = abs(NES_C9) + abs(NES_Spor)) %>%
  dplyr::arrange(dplyr::desc(score))

top_discordant <- nes_join_c9_spor %>%
  dplyr::filter(!direction_agree) %>%
  dplyr::mutate(score = abs(NES_C9) + abs(NES_Spor)) %>%
  dplyr::arrange(dplyr::desc(score))

cat("\nTop concordant pathways:\n")
print(head(top_concordant, 15))

cat("\nTop discordant pathways:\n")
print(head(top_discordant, 15))

topK <- 15

c9_top <- fg_c9 %>%
  dplyr::arrange(padj, dplyr::desc(abs(NES))) %>%
  dplyr::slice(1:topK) %>%
  dplyr::select(pathway, NES_C9 = NES)

sp_top <- fg_spor %>%
  dplyr::arrange(padj, dplyr::desc(abs(NES))) %>%
  dplyr::slice(1:topK) %>%
  dplyr::select(pathway, NES_Spor = NES)

topK_overlap <- dplyr::inner_join(c9_top, sp_top, by = "pathway") %>%
  dplyr::mutate(direction_agree = (sign(NES_C9) == sign(NES_Spor))) %>%
  dplyr::arrange(dplyr::desc(abs(NES_C9) + abs(NES_Spor)))

cat("Top", topK, "overlap n =", nrow(topK_overlap), "\n")
print(topK_overlap)
cat("Direction agree in TopK overlap:", sum(topK_overlap$direction_agree), "/", nrow(topK_overlap), "\n")

