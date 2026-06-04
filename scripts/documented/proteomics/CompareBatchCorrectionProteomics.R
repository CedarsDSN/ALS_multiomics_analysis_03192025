# ============================================================
# Script: CompareBatchCorrectionProteomics.R
# Project: ALS multi-omics analysis
# Source status: documented/sanitized for GitHub handoff
# ============================================================
# Documentation summary from scripts/README.md
# 
# ## `CompareBatchCorrectionProteomics.R`
# 
# **Category:** Proteomics / batch correction sensitivity
# 
# **Purpose**
# 
# Compares C9orf72-vs-control proteomics results across pre-batch and post-batch correction strategies. Maps protein IDs to genes, collapses to gene-level values, runs limma variants with/without shipment adjustment or ComBat, compares t-statistic rank correlations, and runs Hallmark camera pathway analysis.
# 
# **Inputs**
# 
# - `AnswerALS-447-P_proteomics-protein-matrix_prebatchcorrected.txt`.
# - `AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt`.
# - `metadata/clinical/aals_dataportal_datatable.csv`.
# - `ALS_Mutation_Subtypes_with_ALS_Control.csv`.
# - MSigDB Hallmark sets via `msigdbr`.
# 
# **Outputs**
# 
# - Folder `C9vsCtrl_proteoDA_genelevel_outputs`: `protein_to_gene_mapping_prebatch.csv`, `protein_to_gene_mapping_postbatch.csv`, `limma_genelevel_<variant>.csv`, `summary_genelevel_across_variants.csv`, `tstat_spearman_corr_genelevel.csv`, `tstat_spearman_corr_genelevel_heatmap.png`, `camera_hallmark_genelevel_across_variants.csv`, `camera_hallmark_genelevel_top20_per_variant.csv`.
# 
# **Upstream dependencies**
# 
# - No required upstream script if all data files exist. It reuses the same metadata conventions as proteomicsDE.
# 
# **How to run**
# 
# Legacy documented version:
# 
# ```bash
# Rscript scripts/documented/proteomics/CompareBatchCorrectionProteomics.R
# ```
# 
# Refactored command-line version:
# 
# ```bash
# Rscript scripts/refactored/proteomics/03_compare_batch_correction.R --config config/config.yml --base-dir /path/to/private/ALS --outdir analysis/batch_05_proteomics/results/batch_correction
# ```
# 
# **Refactor notes**
# 
# - Make shipment variable and C9/control group definitions explicit.
# - Remove dead/commented RNA count code unless it is truly needed.
# - Save the aligned participant metadata used for each variant.
# ============================================================

#!/usr/bin/env Rscript

## ============================================================
## AnswerALS Proteomics: UniProt-like Protein IDs -> Gene Symbols
## + Collapse to gene-level (median across proteins per gene)
## + C9 vs Control batch-strategy bake-off
## + Hallmark pathway testing (camera) at gene-level
##
## Works with row IDs like:
##   "1/sp|A0AVT1|UBA6_HUMAN"
##   "1/sp|P21964|COMT_HUMAN"
##   "1/iRT_protein"  (will be dropped)
##
## IMPORTANT:
## - This script supports two mapping modes:
##   (A) FAST guess from entry name (UBA6_HUMAN -> UBA6)
##   (B) FORMAL UniProt accession -> HGNC mapping via biomaRt (recommended, needs internet)
##
## - camera/Hallmark requires gene symbols, so we build a gene-level matrix.
##
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(limma)
})

has_sva <- requireNamespace("sva", quietly = TRUE)
has_biomart <- requireNamespace("biomaRt", quietly = TRUE)
has_msigdbr <- requireNamespace("msigdbr", quietly = TRUE)
mainDir <- file.path(ALSPROJ_BASE_DIR, "proteomics/4_matrix")
setwd(mainDir)
## -----------------------------
## USER INPUTS (EDIT THESE)
## -----------------------------

# REQUIRED objects you must load/define before running:
#   pre_mat   : pre-batch protein matrix (rows=Protein IDs, cols=sample IDs)
#   post_mat  : post-batch protein matrix (rows=Protein IDs, cols=sample IDs)
#   meta_prot2: metadata with sample IDs + Mutation_Subtype + Shipment (+ SEX optional)
#
# Example loading:
pre_mat  <- read.table("AnswerALS-447-P_proteomics-protein-matrix_prebatchcorrected.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
post_mat <- read.table("AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
repeat_df_file <- file.path(ALSPROJ_BASE_DIR, "metadata/clinical/aals_dataportal_datatable.csv")
mutation_subtypes_file <- file.path(ALSPROJ_BASE_DIR, "ALS_Mutation_Subtypes_with_ALS_Control.csv")
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
# meta_prot2 <- read.csv("meta_prot2.csv", stringsAsFactors = FALSE)

# Metadata columns (edit to match your metadata)
sample_col   <- "Participant_ID"
subtype_col  <- "Mutation_Subtype"      # "Control", "C9orf72", "Sporadic", "SOD1", ...
shipment_col <- "Shipment"
sex_col      <- "SEX"                   # set to NULL if not present

# Contrast
case_label <- "C9orf72"
ctrl_label <- "Control"

# Missingness filtering (protein-level before gene collapse)
max_na_prop_default <- 0.20

# Balanced shipments subset (optional variant)
min_ctrl_per_ship <- 5
min_case_per_ship <- 1

# Mapping mode:
#   "guess"  = gene symbol guessed from entry name (UBA6_HUMAN -> UBA6)
#   "biomart"= map UniProt accession -> HGNC symbol (recommended; requires biomaRt + internet)
mapping_mode <- "guess"

# Output directory
out_dir <- "C9vsCtrl_proteoDA_genelevel_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## -----------------------------
## SANITY CHECKS
## -----------------------------

stop_if_missing_object <- function(nm) {
  if (!exists(nm, inherits = TRUE)) stop("Missing required object: ", nm)
}
stop_if_missing_object("pre_mat")
stop_if_missing_object("post_mat")
stop_if_missing_object("meta_prot2")

cat("Pre-batch dims:  ", paste(dim(pre_mat), collapse = " x "), "\n")
cat("Post-batch dims: ", paste(dim(post_mat), collapse = " x "), "\n")

## -----------------------------
## PARSERS: Protein ID -> UniProt accession / entry name / gene guess
## -----------------------------

extract_uniprot_accession <- function(x) {
  # Extract ACC from sp|ACC|NAME pattern
  acc <- ifelse(grepl("\\|", x), sub("^.*\\|([^|]+)\\|.*$", "\\1", x), NA_character_)
  acc
}

extract_entry_name <- function(x) {
  # Extract NAME from sp|ACC|NAME (e.g., UBA6_HUMAN)
  name <- ifelse(grepl("\\|", x), sub("^.*\\|[^|]+\\|([^|]+)$", "\\1", x), NA_character_)
  name
}

guess_gene_from_entry <- function(entry_name) {
  # UBA6_HUMAN -> UBA6; handles NA safely
  gene <- entry_name
  gene[is.na(gene)] <- NA_character_
  gene <- sub("_HUMAN$", "", gene)
  # if entry doesn't end with _HUMAN, keep as-is (may be imperfect)
  gene
}

is_irt <- function(x) grepl("iRT_protein", x, fixed = TRUE)

## -----------------------------
## FILTER + IMPUTE HELPERS
## -----------------------------

filter_by_missing <- function(x, max_na_prop = 0.2) {
  keep <- rowMeans(is.na(x)) <= max_na_prop
  x[keep, , drop = FALSE]
}

impute_row_median <- function(x) {
  x2 <- x
  row_med <- apply(x2, 1, median, na.rm = TRUE)
  for (i in seq_len(nrow(x2))) {
    idx <- is.na(x2[i, ])
    if (any(idx)) x2[i, idx] <- row_med[i]
  }
  x2
}

## -----------------------------
## ALIGN MATRIX WITH METADATA
## -----------------------------

align_matrix_and_meta <- function(mat, meta, sample_col) {
  stopifnot(sample_col %in% colnames(meta))
  ids <- intersect(colnames(mat), meta[[sample_col]])
  if (length(ids) < 10) stop("Too few overlapping samples between matrix and metadata.")
  
  meta2 <- meta %>%
    filter(.data[[sample_col]] %in% ids) %>%
    mutate(.sample_id = .data[[sample_col]]) %>%
    arrange(match(.sample_id, ids))
  
  mat2 <- mat[, ids, drop = FALSE]
  stopifnot(all(colnames(mat2) == meta2$.sample_id))
  list(mat = mat2, meta = meta2)
}

prep_c9_vs_ctrl_meta <- function(meta, subtype_col, shipment_col, sex_col,
                                 case_label, ctrl_label) {
  meta2 <- meta %>%
    filter(.data[[subtype_col]] %in% c(ctrl_label, case_label)) %>%
    mutate(
      group = factor(ifelse(.data[[subtype_col]] == case_label, "C9", "Control"),
                     levels = c("Control", "C9")),
      Shipment = factor(.data[[shipment_col]])
    )
  if (!is.null(sex_col) && sex_col %in% colnames(meta2)) {
    meta2 <- meta2 %>% mutate(SEX = factor(.data[[sex_col]]))
  }
  meta2
}

make_balanced_shipment_subset <- function(meta, min_ctrl = 5, min_case = 1) {
  tab <- table(meta$group, meta$Shipment)
  ok_ship <- colnames(tab)[(tab["Control", ] >= min_ctrl) & (tab["C9", ] >= min_case)]
  meta %>% filter(Shipment %in% ok_ship)
}

## -----------------------------
## MAPPING: Protein IDs -> gene symbols
## -----------------------------

map_proteins_to_genes_guess <- function(protein_ids) {
  entry <- extract_entry_name(protein_ids)
  gene_guess <- guess_gene_from_entry(entry)
  tibble(
    Protein = protein_ids,
    uniprot_acc = extract_uniprot_accession(protein_ids),
    entry_name = entry,
    gene_symbol = gene_guess,
    is_irt = is_irt(protein_ids)
  )
}

map_proteins_to_genes_biomart <- function(protein_ids) {
  if (!has_biomart) stop("biomaRt not installed. Install with install.packages('biomaRt') or set mapping_mode='guess'.")
  suppressPackageStartupMessages(library(biomaRt))
  
  acc <- extract_uniprot_accession(protein_ids)
  acc_u <- unique(acc[!is.na(acc)])
  # drop iRT pseudo-id if it sneaks in
  acc_u <- acc_u[acc_u != "iRT_protein"]
  
  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  map_swiss <- biomaRt::getBM(
    attributes = c("uniprotswissprot", "hgnc_symbol"),
    filters = "uniprotswissprot",
    values = acc_u,
    mart = mart
  ) %>% rename(uniprot_acc = uniprotswissprot)
  
  map_trembl <- biomaRt::getBM(
    attributes = c("uniprotsptrembl", "hgnc_symbol"),
    filters = "uniprotsptrembl",
    values = acc_u,
    mart = mart
  ) %>% rename(uniprot_acc = uniprotsptrembl)
  
  map_all <- bind_rows(map_swiss, map_trembl) %>%
    filter(hgnc_symbol != "") %>%
    distinct(uniprot_acc, hgnc_symbol)
  
  entry <- extract_entry_name(protein_ids)
  
  tibble(
    Protein = protein_ids,
    uniprot_acc = acc,
    entry_name = entry,
    is_irt = is_irt(protein_ids)
  ) %>%
    left_join(map_all, by = "uniprot_acc") %>%
    mutate(gene_symbol = hgnc_symbol) %>%
    select(-hgnc_symbol)
}

## -----------------------------
## COLLAPSE TO GENE-LEVEL (median across proteins per gene)
## -----------------------------

collapse_protein_to_gene <- function(expr_mat, protein_ids, mapping_tbl) {
  stopifnot(nrow(expr_mat) == length(protein_ids))
  stopifnot(all(mapping_tbl$Protein %in% protein_ids))
  
  # attach mapping in the same order as expr_mat rows
  map_ord <- mapping_tbl[match(protein_ids, mapping_tbl$Protein), ]
  stopifnot(all(map_ord$Protein == protein_ids))
  
  # drop iRT and rows without a gene symbol
  keep <- (!map_ord$is_irt) & (!is.na(map_ord$gene_symbol)) & (map_ord$gene_symbol != "")
  expr2 <- expr_mat[keep, , drop = FALSE]
  genes <- map_ord$gene_symbol[keep]
  
  # collapse: for each gene, take median across proteins (row-wise median per sample)
  # Implemented by splitting indices per gene for speed/clarity.
  uniq_genes <- unique(genes)
  out <- matrix(NA_real_, nrow = length(uniq_genes), ncol = ncol(expr2),
                dimnames = list(uniq_genes, colnames(expr2)))
  
  for (g in uniq_genes) {
    idx <- which(genes == g)
    if (length(idx) == 1) {
      out[g, ] <- expr2[idx, ]
    } else {
      out[g, ] <- apply(expr2[idx, , drop = FALSE], 2, median, na.rm = TRUE)
    }
  }
  out
}

## -----------------------------
## LIMMA + CAMERA
## -----------------------------
major_ship <- c("14","20","4","15","16","9")

pre_meta0 <- pre_meta0 %>%
  mutate(
    ShipmentMajor = as.character(Shipment),
    ShipmentMajor = ifelse(ShipmentMajor %in% major_ship, ShipmentMajor, "Other"),
    ShipmentMajor = factor(ShipmentMajor)
  )

run_limma_C9vsCtrl <- function(expr, meta,
                               include_shipment = TRUE,
                               include_sex = TRUE,
                               robust = TRUE,
                               trend = TRUE,
                               id_col = "Participant_ID") {
  
  al <- align_expr_meta_by_id(expr, meta, id_col = id_col)
  expr <- al$expr
  meta2 <- al$meta
  
  # Reference: Control
  meta2$group <- relevel(factor(meta2$group), ref = "Control")
  
  f_terms <- c("group")
  
  if (include_sex && "SEX" %in% colnames(meta2) &&
      nlevels(droplevels(meta2$SEX)) >= 2) {
    f_terms <- c(f_terms, "SEX")
  }
  
  if (include_shipment && "ShipmentMajor" %in% colnames(meta2) &&
      nlevels(droplevels(meta2$ShipmentMajor)) >= 2) {
    meta2$ShipmentMajor <- droplevels(meta2$ShipmentMajor)
    f_terms <- c(f_terms, "ShipmentMajor")
  }
  
  fmla <- as.formula(paste("~", paste(f_terms, collapse = " + ")))
  design <- model.matrix(fmla, data = meta2)
  colnames(design) <- make.names(colnames(design))
  
  fit <- lmFit(expr, design)
  fit <- eBayes(fit, robust = robust, trend = trend)
  
  if (!("groupC9" %in% colnames(design))) {
    stop("Design missing groupC9. Check group coding.")
  }
  
  tt <- topTable(fit, coef = "groupC9", number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  
  list(top = tt, fit = fit, design = design, meta_aligned = meta2, expr_aligned = expr)
}


combat_protect_C9 <- function(expr, meta, batch_col = "Shipment") {
  if (!has_sva) stop("Package 'sva' is required for ComBat. Install it: install.packages('sva')")
  batch <- meta[[batch_col]]
  
  # Protect group effect (and SEX if present)
  if ("SEX" %in% colnames(meta) && nlevels(droplevels(meta$SEX)) >= 2) {
    mod <- model.matrix(~ group + SEX, data = meta)
  } else {
    mod <- model.matrix(~ group, data = meta)
  }
  
  expr_imp <- impute_row_median(expr)
  sva::ComBat(dat = expr_imp, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
}

get_hallmark_sets <- function() {
  if (!has_msigdbr) stop("msigdbr not installed. Install with install.packages('msigdbr')")
  msig <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  split(msig$gene_symbol, msig$gs_name)
}

run_camera_hallmark <- function(expr, design, fit, hallmark_sets) {
  idx <- limma::ids2indices(hallmark_sets, rownames(expr), remove.empty = TRUE)
  # camera can take the fitted design/contrast, but easiest is direct contrast vector name:
  # We used makeContrasts(C9vsCtrl = groupC9 - groupControl), so use contrast = "groupC9 - groupControl"
  # But camera takes either a contrast vector or contrast specification.
  # We'll compute a contrast vector explicitly from design columns:
  cn <- colnames(design)
  if (!all(c("groupC9", "groupControl") %in% cn)) stop("Design missing group columns for camera.")
  contrast_vec <- rep(0, length(cn)); names(contrast_vec) <- cn
  contrast_vec["groupC9"] <- 1
  contrast_vec["groupControl"] <- -1
  
  cam <- limma::camera(expr, index = idx, design = design, contrast = contrast_vec)
  cam$pathway <- rownames(cam)
  cam
}

## ------------------------------------------------------------
## Align protein matrix to metadata by "core" Participant ID
## ------------------------------------------------------------

extract_core_id2 <- function(x) {
  # Normalize first: convert "CASE_NEU..." or "CTRL_NEU..." to "CASE-NEU..." / "CTRL-NEU..."
  x2 <- sub("^(CASE|CTRL)_", "\\1-", x)
  
  # Extract the core participant ID: "CASE-NEU<alnum>" or "CTRL-NEU<alnum>"
  core <- sub("^((CASE|CTRL)-NEU[0-9A-Za-z]+).*$", "\\1", x2)
  
  # If it didn't match, return NA (so we can detect issues)
  core[!grepl("^(CASE|CTRL)-NEU[0-9A-Za-z]+$", core)] <- NA_character_
  core
}

align_matrix_and_meta_by_core2 <- function(mat, meta, sample_col) {
  stopifnot(sample_col %in% colnames(meta))
  
  mat_cols_raw <- colnames(mat)
  mat_core <- extract_core_id2(mat_cols_raw)
  
  meta_core <- extract_core_id2(meta[[sample_col]])
  
  map_mat <- data.frame(
    mat_col = mat_cols_raw,
    core_id = mat_core,
    stringsAsFactors = FALSE
  )
  
  meta2 <- meta %>% mutate(core_id = meta_core)
  
  # Drop rows where core_id couldn't be parsed
  map_mat <- map_mat %>% filter(!is.na(core_id))
  meta2   <- meta2 %>% filter(!is.na(core_id))
  
  common_core <- intersect(unique(map_mat$core_id), unique(meta2$core_id))
  
  if (length(common_core) < 10) {
    message("Example matrix core IDs: ", paste(head(unique(map_mat$core_id), 6), collapse = ", "))
    message("Example meta core IDs:   ", paste(head(unique(meta2$core_id), 6), collapse = ", "))
    stop("Too few overlapping *core* sample IDs between matrix and metadata.")
  }
  
  # Keep matrix columns whose core IDs exist in metadata
  map_keep <- map_mat %>% filter(core_id %in% common_core)
  
  # If multiple matrix columns map to the same participant, keep the first
  # (Alternative: average replicates — tell me if you want that.)
  map_keep <- map_keep %>%
    group_by(core_id) %>%
    slice(1) %>%
    ungroup()
  
  mat2 <- mat[, map_keep$mat_col, drop = FALSE]
  colnames(mat2) <- map_keep$core_id
  
  meta3 <- meta2 %>%
    filter(core_id %in% map_keep$core_id) %>%
    distinct(core_id, .keep_all = TRUE) %>%
    arrange(match(core_id, map_keep$core_id))
  
  stopifnot(all(colnames(mat2) == meta3$core_id))
  
  list(mat = mat2, meta = meta3)
}

align_expr_meta_by_id <- function(expr, meta, id_col = "Participant_ID") {
  expr <- as.matrix(expr)
  stopifnot(id_col %in% colnames(meta))
  if (is.null(colnames(expr))) stop("Expression matrix has NULL colnames.")
  
  ids_expr <- colnames(expr)
  ids_meta <- meta[[id_col]]
  
  common <- intersect(ids_expr, ids_meta)
  if (length(common) < 10) {
    stop("Too few overlapping sample IDs between expr colnames and meta[[", id_col, "]]. Overlap=", length(common))
  }
  
  # subset + reorder meta to match expr (or reorder expr to match meta; pick one and stick to it)
  # I recommend: reorder expr to match meta, because meta carries grouping variables
  meta2 <- meta[match(common, ids_meta), , drop = FALSE]  # in meta order
  expr2 <- expr[, meta2[[id_col]], drop = FALSE]
  
  stopifnot(all(colnames(expr2) == meta2[[id_col]]))
  list(expr = expr2, meta = meta2)
}

## -----------------------------
## MAIN: prepare meta + align matrices
## -----------------------------

meta_c9 <- prep_c9_vs_ctrl_meta(meta_prot2, subtype_col, shipment_col, sex_col, case_label, ctrl_label)

post_al <- align_matrix_and_meta_by_core2(post_mat, meta_c9, sample_col)
pre_al  <- align_matrix_and_meta_by_core2(pre_mat,  meta_c9, sample_col)

dim(post_al$mat)
head(colnames(post_al$mat))
head(post_al$meta$core_id)



pre_expr0 <- as.matrix(pre_al$mat); mode(pre_expr0) <- "numeric"
post_expr0 <- as.matrix(post_al$mat); mode(post_expr0) <- "numeric"

pre_meta0 <- pre_al$meta
post_meta0 <- post_al$meta

# Apply missingness filter at protein level before collapsing
pre_expr_f <- filter_by_missing(pre_expr0, max_na_prop = max_na_prop_default)
post_expr_f <- filter_by_missing(post_expr0, max_na_prop = max_na_prop_default)

cat("After protein-level missingness filter (<= ", max_na_prop_default*100, "%):\n", sep = "")
cat("  Pre-batch proteins kept:  ", nrow(pre_expr_f), "\n")
cat("  Post-batch proteins kept: ", nrow(post_expr_f), "\n")

# Map protein IDs -> gene symbols
prot_ids_pre <- rownames(pre_expr_f)
prot_ids_post <- rownames(post_expr_f)

cat("Mapping mode: ", mapping_mode, "\n", sep = "")

if (mapping_mode == "guess") {
  map_pre <- map_proteins_to_genes_guess(prot_ids_pre)
  map_post <- map_proteins_to_genes_guess(prot_ids_post)
} else if (mapping_mode == "biomart") {
  map_pre <- map_proteins_to_genes_biomart(prot_ids_pre)
  map_post <- map_proteins_to_genes_biomart(prot_ids_post)
} else {
  stop("Unknown mapping_mode. Use 'guess' or 'biomart'.")
}

# Save mapping tables for transparency
write.csv(map_pre,  file.path(out_dir, "protein_to_gene_mapping_prebatch.csv"), row.names = FALSE)
write.csv(map_post, file.path(out_dir, "protein_to_gene_mapping_postbatch.csv"), row.names = FALSE)

# Collapse to gene-level matrices (median across proteins per gene)
pre_gene <- collapse_protein_to_gene(pre_expr_f, prot_ids_pre, map_pre)
post_gene <- collapse_protein_to_gene(post_expr_f, prot_ids_post, map_post)

cat("Gene-level dims:\n")
cat("  Pre-batch genes:  ", nrow(pre_gene), " x ", ncol(pre_gene), "\n", sep = "")
cat("  Post-batch genes: ", nrow(post_gene), " x ", ncol(post_gene), "\n", sep = "")

# Impute remaining NAs at gene level (for modeling stability)
pre_gene_imp  <- impute_row_median(pre_gene)
post_gene_imp <- impute_row_median(post_gene)

## -----------------------------
## VARIANTS: gene-level analysis
## -----------------------------
major_ship <- c("14","20","4","15","16","9")

pre_meta0 <- pre_meta0 %>%
  mutate(
    ShipmentMajor = as.character(Shipment),
    ShipmentMajor = ifelse(ShipmentMajor %in% major_ship, ShipmentMajor, "Other"),
    ShipmentMajor = factor(ShipmentMajor)
  )


variants <- list()

# V1: Post-batch gene matrix + limma (SEX + group)
variants[["V1_post_gene_SEX_group"]] <- list(
  expr = post_gene_imp,
  meta = post_meta0,
  include_shipment = FALSE
)

# V1b: Post-batch gene matrix + limma (SEX + Shipment + group)
variants[["V1b_post_gene_SEX_Ship_group"]] <- list(
  expr = post_gene_imp,
  meta = post_meta0,
  include_shipment = TRUE
)

# V2: Pre-batch gene matrix + limma (SEX + Shipment + group)
variants[["V2_pre_gene_SEX_Ship_group"]] <- list(
  expr = pre_gene_imp,
  meta = pre_meta0,
  include_shipment = TRUE
)
has_sva <- requireNamespace("sva", quietly = TRUE)
if (!has_sva) {
  stop("Package 'sva' is required for ComBat. Install it with: install.packages('sva')")
}

# V3: Pre-batch gene matrix -> ComBat protecting C9 vs Ctrl -> limma (SEX + group)
expr_cb <- combat_protect_C9(pre_gene_imp, pre_meta0, batch_col = "Shipment")
variants[["V3_pre_gene_ComBatProtectC9_SEX_group"]] <- list(
  expr = expr_cb,
  meta = pre_meta0,
  include_shipment = FALSE
)

# V4 (optional): Balanced shipments only on pre-batch gene matrix + limma (SEX + Shipment + group)
# IMPORTANT: meta_bal must still have Participant_ID, group, Shipment, (optional SEX)
# expr stays the full matrix

meta_bal <- make_balanced_shipment_subset(pre_meta0,
                                          min_ctrl = min_ctrl_per_ship,
                                          min_case = min_case_per_ship)

variants[["V4_pre_gene_balancedShip_SEX_Ship_group"]] <- list(
  expr = pre_gene_imp,   # full matrix (has correct colnames)
  meta = meta_bal,       # balanced subset
  include_shipment = TRUE
)



## -----------------------------
## RUN LIMMA + SAVE RESULTS
## -----------------------------

all_top <- list()
summary_tbl <- list()

summarize_top <- function(tt, run_id) {
  tibble(
    run_id = run_id,
    n_genes_tested = nrow(tt),
    best_p = min(tt$P.Value, na.rm = TRUE),
    best_FDR = min(tt$adj.P.Val, na.rm = TRUE),
    n_FDR_0.05 = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
    n_FDR_0.10 = sum(tt$adj.P.Val < 0.10, na.rm = TRUE)
  )
}

for (nm in names(variants)) {
  v <- variants[[nm]]
  cat("Running ", nm, " ...\n", sep = "")
  
  res <- run_limma_C9vsCtrl(
    expr = v$expr,
    meta = v$meta,
    include_shipment = v$include_shipment,
    include_sex = TRUE,
    robust = TRUE,
    trend = TRUE
  )
  
  tt <- res$top
  all_top[[nm]] <- tt
  
  summary_tbl[[nm]] <- summarize_top(tt, nm)
  
  write.csv(tt, file.path(out_dir, paste0("limma_genelevel_", nm, ".csv")), row.names = FALSE)
}

summary_df <- bind_rows(summary_tbl) %>% arrange(best_FDR)
print(summary_df)
write.csv(summary_df, file.path(out_dir, "summary_genelevel_across_variants.csv"), row.names = FALSE)

## -----------------------------
## RANK STABILITY: t-stat Spearman correlations
## -----------------------------

t_wide <- bind_rows(lapply(names(all_top), function(nm) {
  all_top[[nm]] %>% dplyr::select(gene, t) %>% mutate(run_id = nm)
})) %>%
  distinct() %>%
  pivot_wider(names_from = run_id, values_from = t)

run_ids <- setdiff(colnames(t_wide), "gene")
corr_mat <- matrix(NA_real_, nrow = length(run_ids), ncol = length(run_ids),
                   dimnames = list(run_ids, run_ids))

for (a in run_ids) {
  for (b in run_ids) {
    xa <- t_wide[[a]]
    xb <- t_wide[[b]]
    ok <- is.finite(xa) & is.finite(xb)
    if (sum(ok) > 50) corr_mat[a, b] <- cor(xa[ok], xb[ok], method = "spearman")
  }
}
write.csv(round(corr_mat, 6), file.path(out_dir, "tstat_spearman_corr_genelevel.csv"), row.names = TRUE)

corr_df <- as.data.frame(as.table(corr_mat)) %>% rename(x = Var1, y = Var2, spearman = Freq)
p_corr <- ggplot(corr_df, aes(x = x, y = y, fill = spearman)) +
  geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene-level Spearman correlation of limma t-stats across variants", x = "", y = "")
ggsave(file.path(out_dir, "tstat_spearman_corr_genelevel_heatmap.png"),
       plot = p_corr, width = 12, height = 8, dpi = 150)

## -----------------------------
## PATHWAY ANALYSIS: Hallmark via camera (gene symbols required)
## -----------------------------

if (has_msigdbr) {
  hallmark_sets <- get_hallmark_sets()
  
  camera_all <- list()
  for (nm in names(variants)) {
    v <- variants[[nm]]
    # rebuild design exactly like limma run to keep consistent
    res <- run_limma_C9vsCtrl(
      expr = v$expr,
      meta = v$meta,
      include_shipment = v$include_shipment,
      include_sex = TRUE,
      robust = TRUE,
      trend = TRUE
    )
    
    cam <- run_camera_hallmark(v$expr, res$design, res$fit, hallmark_sets)
    cam$variant <- nm
    camera_all[[nm]] <- cam
  }
  
  camera_df <- bind_rows(camera_all) %>%
    select(variant, pathway, NGenes, Direction, PValue, FDR) %>%
    arrange(variant, FDR, PValue)
  
  write.csv(camera_df, file.path(out_dir, "camera_hallmark_genelevel_across_variants.csv"), row.names = FALSE)
  
  # Quick convenience: top 20 pathways per variant
  top_cam <- camera_df %>%
    group_by(variant) %>%
    slice_min(order_by = FDR, n = 20, with_ties = FALSE) %>%
    ungroup()
  write.csv(top_cam, file.path(out_dir, "camera_hallmark_genelevel_top20_per_variant.csv"), row.names = FALSE)
  
} else {
  cat("Skipping Hallmark camera: msigdbr not installed.\n")
}

cat("Done. Outputs written to: ", out_dir, "\n", sep = "")
