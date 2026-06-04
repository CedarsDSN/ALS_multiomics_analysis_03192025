# ============================================================
# Refactored entrypoint for original script: genomic_stratification.R
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
if (isTRUE(ALS_ARGS$help)) { als_print_help("03_genomic_stratification"); quit(save = "no", status = 0) }
CFG <- als_load_config(ALS_ARGS$config)
ALSPROJ_BASE_DIR <- als_choose(ALS_ARGS$base_dir, CFG, "project.base_dir", getwd())
ALSPROJ_OUTDIR <- als_choose(ALS_ARGS$outdir, CFG, "outputs.root", file.path(getwd(), "results", "03_genomic_stratification"))
ALSPROJ_SCRIPT_ID <- "03_genomic_stratification"
als_init(ALSPROJ_SCRIPT_ID, ALSPROJ_BASE_DIR, ALSPROJ_OUTDIR, ALS_ARGS$redirect_outputs)
als_register_output_wrappers()
als_load_objects(ALS_ARGS$objects)
on.exit(als_write_session_info(ALSPROJ_SCRIPT_ID), add = TRUE)
# ---- End bootstrap ----

# ============================================================
# STEP 1: VARIANT SELECTION FROM ANSWERAS WGS DATA
# Using axis-defining genes + ALS panel (not DisGeNET)
# ============================================================

library(dplyr)
library(data.table)
library(msigdbr)
library(tidyr)

# ============================================================
# PART A: DEFINE YOUR GENE SETS
# Three tiers in priority order
# ============================================================

# --- Tier 1: Axis-defining genes (from your ssGSEA Hallmark sets)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")

metabolic_genes <- hallmark |>
  filter(gs_name %in% c(
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_FATTY_ACID_METABOLISM",
    "HALLMARK_ADIPOGENESIS",
    "HALLMARK_GLYCOLYSIS"
  )) |> pull(gene_symbol) |> unique()

inflam_genes <- hallmark |>
  filter(gs_name %in% c(
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
  )) |> pull(gene_symbol) |> unique()

cellcycle_genes <- hallmark |>
  filter(gs_name %in% c(
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_DNA_REPAIR"
  )) |> pull(gene_symbol) |> unique()

# --- Tier 2: Your 34-gene ALS panel (already validated in your analysis)
als_panel_genes <- c(
  "SOD1", "TARDBP", "FUS", "C9orf72", "TBK1", "NEK1",
  "OPTN", "SQSTM1", "ANXA11", "CHCHD10", "DCTN1", "SETX",
  "UBQLN2", "VCP", "ATXN2", "TIA1", "MATR3", "ALS2",
  "SIGMAR1", "VAPB", "TRPM7", "GRN", "CHMP2B", "ERBB4",
  "DAO", "EWSR1", "HNRNPA1", "HNRNPA2B1", "MAPT",
  "PFN1", "TUBA4A", "KIF5A", "CCNF", "DNAJC7"
  # NEFH excluded — hypervariable repeat artifact
)

# --- Combined gene universe for this analysis
target_genes <- unique(c(
  metabolic_genes,
  inflam_genes,
  cellcycle_genes,
  als_panel_genes
))

cat("=== GENE SET SUMMARY ===\n")
cat("Metabolic axis genes:   ", length(metabolic_genes), "\n")
cat("Inflammation axis genes:", length(inflam_genes), "\n")
cat("Cell cycle axis genes:  ", length(cellcycle_genes), "\n")
cat("ALS panel genes:        ", length(als_panel_genes), "\n")
cat("Total unique genes:     ", length(target_genes), "\n")

# Save gene set membership for later profiling of results
gene_tier_df <- data.frame(
  gene = target_genes,
  in_metabolic   = target_genes %in% metabolic_genes,
  in_inflam      = target_genes %in% inflam_genes,
  in_cellcycle   = target_genes %in% cellcycle_genes,
  in_als_panel   = target_genes %in% als_panel_genes
)

# Assign primary axis label (for later edge profiling)
gene_tier_df <- gene_tier_df |>
  mutate(primary_axis = case_when(
    in_als_panel                                   ~ "ALS_panel",
    in_metabolic & !in_inflam & !in_cellcycle      ~ "metabolic",
    in_inflam    & !in_metabolic & !in_cellcycle   ~ "inflammatory",
    in_cellcycle & !in_metabolic & !in_inflam      ~ "cell_cycle",
    TRUE                                            ~ "multi_axis"
  ))

table(gene_tier_df$primary_axis)
write.csv(gene_tier_df, "gene_tier_membership.csv", row.names=FALSE)

# ============================================================
# PART B: EXTRACT VARIANTS FROM ANNOVAR FILES
# Filter to: rare (MAF < 0.01), damaging, in target genes
# ============================================================

# --- Define damaging variant classes (same as your burden analysis)
damaging_classes <- c(
  "nonsynonymous SNV",
  "stopgain",
  "stoploss",
  "frameshift insertion",
  "frameshift deletion"
)

# --- Function to process one chromosome annotation file
process_chromosome <- function(chr_file, target_genes, 
                               damaging_classes, maf_threshold = 0.01) {
  
  cat("Processing:", chr_file, "\n")
  
  # Load ANNOVAR file — adjust column names to match your actual files
  # Standard ANNOVAR columns include: Chr, Start, End, Ref, Alt,
  # Func.refGene, Gene.refGene, ExonicFunc.refGene,
  # Max_Freq_of_3 (or similar MAF column)
  
  df <- fread(chr_file, sep="\t", header=TRUE, 
              showProgress=FALSE, fill=TRUE)
  
  # Standardise column names (adjust to match your actual headers)
  # Check your file: names(fread(chr_file, nrows=0))
  colnames(df) <- make.names(colnames(df))
  
  # Identify MAF column — your files use "Max.Freq.of.3" or similar
  # From your prior analysis you used "Max Freq of 3"
  maf_col <- grep("Max.*Freq|AF|MAF|freq", colnames(df), 
                  value=TRUE, ignore.case=TRUE)[1]
  cat("  MAF column detected:", maf_col, "\n")
  
  # Identify gene column
  gene_col <- grep("Gene.refGene|Gene_refGene|Gene", 
                   colnames(df), value=TRUE)[1]
  
  # Identify functional class column
  func_col <- grep("ExonicFunc|exonic_func", 
                   colnames(df), value=TRUE)[1]
  
  cat("  Gene column:", gene_col, "\n")
  cat("  Function column:", func_col, "\n")
  cat("  Rows before filtering:", nrow(df), "\n")
  
  # Step 1: Filter to damaging variants
  df_damaging <- df |>
    filter(.data[[func_col]] %in% damaging_classes)
  
  cat("  After damaging filter:", nrow(df_damaging), "\n")
  
  # Step 2: Filter to rare variants (MAF < threshold)
  # Handle missing/dot entries as rare (conservative)
  df_rare <- df_damaging |>
    mutate(maf_numeric = suppressWarnings(
      as.numeric(.data[[maf_col]]))) |>
    filter(is.na(maf_numeric) | maf_numeric < maf_threshold)
  
  cat("  After MAF < ", maf_threshold, " filter:", nrow(df_rare), "\n")
  
  # Step 3: Filter to target genes
  # Handle cases where gene column contains multiple genes (e.g. "GENE1;GENE2")
  df_target <- df_rare |>
    mutate(gene_clean = gsub("[,;].*", "", .data[[gene_col]])) |>
    filter(gene_clean %in% target_genes)
  
  cat("  After target gene filter:", nrow(df_target), "\n\n")
  
  return(df_target)
}

# --- Process all chromosomes
# Adjust path pattern to match your HPC file structure
chr_files <- list.files(
  path = file.path(ALSPROJ_BASE_DIR),   # UPDATE THIS
  pattern = "chr[0-9]+.*\\.txt$",
  full.names = TRUE
)

cat("Found", length(chr_files), "chromosome files\n\n")

# Process all chromosomes and combine
# Note: this will take a few minutes on HPC
all_variants <- lapply(chr_files, function(f) {
  tryCatch(
    process_chromosome(f, target_genes, damaging_classes),
    error = function(e) {
      cat("ERROR processing", f, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
}) |>
  bind_rows()

cat("=== TOTAL VARIANTS AFTER FILTERING ===\n")
cat("Total rare damaging variants in target genes:", nrow(all_variants), "\n")
cat("Unique genes hit:", length(unique(all_variants$gene_clean)), "\n")

# Check which target genes have no variants (important for PPI later)
genes_with_variants <- unique(all_variants$gene_clean)
genes_without_variants <- setdiff(target_genes, genes_with_variants)
cat("Target genes with zero rare damaging variants:", 
    length(genes_without_variants), "\n")
cat("(These will be excluded from PPI subnetwork)\n")

# ============================================================
# PART C: IDENTIFY SAMPLE IDs IN VARIANT DATA
# Your ANNOVAR files contain sample columns
# ============================================================

# Check column structure — sample genotypes are in final columns
# Standard ANNOVAR output has sample columns after annotation columns
# Column names are typically the sample/participant IDs

# ============================================================
# CORRECTED STEP 1 — BASED ON YOUR ACTUAL COLUMN STRUCTURE
# ============================================================

library(dplyr)
library(data.table)
library(tibble)

# Enforce dplyr namespace explicitly
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise

# ============================================================
# CONFIRMED COLUMN MAPPING FROM YOUR FILES
# ============================================================
# Col 1:  Gene                        <- gene name
# Col 5:  ExonicFunc                  <- variant class
# Col 34: "Max Freq of 3"             <- MAF column
# Col 36+: CASE-NEU* / CTRL-NEU*     <- sample genotypes

# Annotation columns = everything before sample columns
ANNOTATION_COLS <- c(
  "Gene", "Func.knownGene", "Func.refGene", "Func.ensGene",
  "ExonicFunc", "AAChange.knownGene", "AAChange.refGene",
  "AAChange.ensGene", "chrom", "position", "FILTER", "REF", "ALT",
  "rsID", "GTEx (spinal_cord/cortex/frontal_cortex)",
  "Up regulated Protein Atlas Gene", "OMIM disease",
  "ClinVar sig (CLINSIG)", "ClinVar trait (CLNDBN)",
  "CGD condition", "CGD manifestation", "CGD refs (PubMed Ids)",
  "ALSoD/Lander/DisGeNet", "Harms Combining Criteria",
  "Harms Meta", "Harms Leaning Towards", "Known ALS Variant",
  "ACMG variant", "ACMG gene", "in silico Predictions",
  "RVIS perc", "GDI Phred / Level", "LoFtool perc",
  "Max Freq of 3", "Freq:1KG/ESP/ExAC"
)

# Damaging variant classes in your ExonicFunc column
DAMAGING_CLASSES <- c(
  "nonsynonymous SNV",
  "stopgain",
  "stoploss",
  "frameshift insertion",
  "frameshift deletion",
  "frameshift substitution"  # include in case present
)

# ============================================================
# PART B (CORRECTED): PROCESS CHROMOSOME FILES
# ============================================================

process_chromosome_corrected <- function(chr_file, target_genes,
                                         damaging_classes,
                                         maf_threshold = 0.01) {
  cat("Processing:", basename(chr_file), "... ")
  
  df <- data.table::fread(
    chr_file, sep = "\t", header = TRUE,
    showProgress = FALSE, fill = TRUE,
    quote = ""
  )
  
  # Identify sample columns by pattern (CASE- or CTRL-)
  sample_cols <- grep("^(CASE|CTRL)-", colnames(df), value = TRUE)
  cat(length(sample_cols), "sample cols | ")
  
  # --- Filter 1: damaging variants
  df_damaging <- df[df$ExonicFunc %in% damaging_classes, ]
  cat("damaging:", nrow(df_damaging), "| ")
  
  if (nrow(df_damaging) == 0) {
    cat("SKIP (no damaging variants)\n")
    return(NULL)
  }
  
  # --- Filter 2: rare variants (MAF < threshold)
  # "Max Freq of 3" column — dots and NA treated as rare
  df_damaging$maf_num <- suppressWarnings(
    as.numeric(df_damaging[["Max Freq of 3"]])
  )
  df_rare <- df_damaging[is.na(df_damaging$maf_num) |
                           df_damaging$maf_num < maf_threshold, ]
  cat("rare:", nrow(df_rare), "| ")
  
  if (nrow(df_rare) == 0) {
    cat("SKIP (no rare variants)\n")
    return(NULL)
  }
  
  # --- Filter 3: target genes
  # Gene column may contain multiple genes separated by ";" or ","
  # Keep primary gene (first listed)
  df_rare$gene_clean <- gsub("[,;(].*", "", df_rare$Gene)
  df_rare$gene_clean <- trimws(df_rare$gene_clean)
  
  df_target <- df_rare[df_rare$gene_clean %in% target_genes, ]
  cat("in target genes:", nrow(df_target), "\n")
  
  if (nrow(df_target) == 0) return(NULL)
  
  # Keep only gene column + sample columns
  keep_cols <- c("gene_clean", sample_cols)
  return(df_target[, keep_cols, with = FALSE])
}

# --- Run on all chromosome files
chr_files <- list.files(
  path = file.path(ALSPROJ_BASE_DIR),  # UPDATE THIS PATH
  pattern = "\\_exon.txt$",
  full.names = TRUE
)

# If files are per-chromosome, confirm
cat("Chromosome files found:", length(chr_files), "\n")
cat("First file:", basename(chr_files[1]), "\n\n")

all_variants_list <- lapply(chr_files, function(f) {
  tryCatch(
    process_chromosome_corrected(f, target_genes, DAMAGING_CLASSES),
    error = function(e) {
      cat("ERROR in", basename(f), ":", conditionMessage(e), "\n")
      NULL
    }
  )
})

# Remove NULLs and combine
all_variants_list <- Filter(Negate(is.null), all_variants_list)
cat("\nChromosomes with target variants:", length(all_variants_list), "\n")

all_variants <- data.table::rbindlist(all_variants_list, fill = TRUE)
cat("Total rare damaging variants in target genes:", nrow(all_variants), "\n")
cat("Unique target genes hit:", length(unique(all_variants$gene_clean)), "\n")

# ============================================================
# PART C (CORRECTED): IDENTIFY SAMPLE COLUMNS
# ============================================================

# All sample columns in your data
all_sample_cols <- grep("^(CASE|CTRL)-", colnames(all_variants),
                        value = TRUE)

# Separate CASE (ALS) from CTRL (healthy) columns
case_cols <- grep("^CASE-", all_sample_cols, value = TRUE)
ctrl_cols  <- grep("^CTRL-", all_sample_cols, value = TRUE)

cat("\nTotal sample columns:", length(all_sample_cols), "\n")
cat("CASE (ALS) columns:", length(case_cols), "\n")
cat("CTRL columns:", length(ctrl_cols), "\n")

# For this analysis: use CASE samples only
# (we are subtyping sporadic ALS patients)
# CTRL included later as reference if needed

# ============================================================
# PART D (CORRECTED): BUILD PATIENT × GENE VARIANT COUNT MATRIX
# ============================================================

# Genotype values in your files are typically:
# "0/0" = homozygous reference (no variant)
# "0/1" or "1/0" = heterozygous (carries variant)
# "1/1" = homozygous alt (carries variant)
# "./." or "." = missing

is_variant_gt <- function(gt) {
  # Strip leading single quote (Excel artifact in your files)
  gt_clean <- sub("^'", "", gt)
  # Strip trailing whitespace (note '././ ' has trailing space)
  gt_clean <- trimws(gt_clean)
  
  ifelse(
    is.na(gt_clean) |
      gt_clean %in% c(".", "./.", "0/0", "0|0", ""),
    0L, 1L
  )
}

# ============================================================
# VERIFY THE FIX on your test data before rerunning
# ============================================================

# Test on actual values from your files
test_values <- c("'0/0", "'0/1", "'1/1", "'./. ", NA, "'./.")
cat("=== GENOTYPE PARSING VERIFICATION ===\n")
result <- data.frame(
  raw      = test_values,
  parsed   = is_variant_gt(test_values),
  expected = c(0, 1, 1, 0, 0, 0)
)
print(result)
cat("All correct:", all(result$parsed == result$expected), "\n\n")
# Convert genotype strings to 0/1 for CASE samples only
# Then count variants per gene per patient
# Use data.table for speed


variants_dt <- as.data.table(all_variants)

case_cols_present <- intersect(case_cols, colnames(variants_dt))

# Apply corrected variant indicator
variants_dt[,
            (case_cols_present) := lapply(.SD, is_variant_gt),
            .SDcols = case_cols_present
]

# Sum per gene per patient
gene_patient_counts <- variants_dt[,
                                   lapply(.SD, sum, na.rm = TRUE),
                                   by = gene_clean,
                                   .SDcols = case_cols_present
]

cat("Gene x patient matrix:", nrow(gene_patient_counts),
    "genes x", ncol(gene_patient_counts)-1, "patients\n")

# Transpose
patient_gene_matrix <- as.data.frame(t(gene_patient_counts[, -1]))
colnames(patient_gene_matrix) <- gene_patient_counts$gene_clean
rownames(patient_gene_matrix) <- case_cols_present

# ============================================================
# PART E (CORRECTED): MAP SAMPLE IDs TO PARTICIPANT IDs
# ============================================================

# Your sample IDs are in CASE-NEU format
# Your transcriptomic data uses Participant_ID format
# Apply your existing ID crosswalk

# The CASE-NEU IDs in the VCF need to be mapped to 
# CGND-HDA identifiers and then to Participant_ID
# Using your existing NYGC_CGND_ID crosswalk

# Load your crosswalk file (adjust path)
# crosswalk <- read.csv("NYGC_CGND_crosswalk.csv")

# For now, strip to get base ID format
# CASE-NEUMN012EVP -> check what format your sporadic_scores uses

# Check what format sporadic_scores Participant_IDs are in
cat("\nSample of sporadic_scores Participant_IDs:\n")
print(head(sporadic_final_check$Participant_ID, 5))

cat("\nSample of CASE column names from VCF:\n")
print(head(case_cols_present, 5))

# If they match directly:
# patient_ids_vcf <- case_cols_present

# If they need stripping (e.g. CASE-NEUMN012EVP -> CASE-NEUMN012EVP):
# They appear to already be in CASE-NEU format
# Check if your sporadic_scores uses this format or CGND format

# Transpose to patient x gene
patient_gene_matrix <- as.data.frame(
  t(gene_patient_counts[, -1])  # remove gene_clean column
)
colnames(patient_gene_matrix) <- gene_patient_counts$gene_clean
rownames(patient_gene_matrix) <- case_cols_present

cat("\nTransposed matrix:", nrow(patient_gene_matrix),
    "patients x", ncol(patient_gene_matrix), "genes\n")

# ============================================================
# PART F: MATCH TO SPORADIC SUBTYPE LABELS
# ============================================================

# First check ID format overlap
cat("\n=== ID FORMAT CHECK ===\n")
cat("VCF IDs (first 3):\n"); print(head(rownames(patient_gene_matrix), 3))
cat("sporadic_scores IDs (first 3):\n")
print(head(sporadic_final_check$Participant_ID, 3))

# Direct match attempt
direct_overlap <- intersect(
  rownames(patient_gene_matrix),
  sporadic_final_check$Participant_ID
)
cat("Direct ID overlap:", length(direct_overlap), "patients\n")

# If direct overlap is low, try ID conversion
# Option 1: Strip trailing suffix
rownames_stripped <- gsub("-[0-9]+-T$", "",
                          rownames(patient_gene_matrix))
stripped_overlap <- intersect(rownames_stripped,
                              sporadic_scores$Participant_ID)
cat("After suffix stripping:", length(stripped_overlap), "patients\n")

# Option 2: Use crosswalk if needed
# Uncomment and adjust if above options don't work:
# crosswalk <- read.csv("path/to/crosswalk.csv")
# id_map <- setNames(crosswalk$Participant_ID,
#                    crosswalk$CASE_NEU_ID)
# rownames(patient_gene_matrix) <- id_map[rownames(patient_gene_matrix)]

# Use whichever mapping worked
if (length(direct_overlap) >= length(stripped_overlap)) {
  cat("Using direct IDs\n")
  final_ids <- rownames(patient_gene_matrix)
} else {
  cat("Using stripped IDs\n")
  rownames(patient_gene_matrix) <- rownames_stripped
  final_ids <- rownames(patient_gene_matrix)
}

# Filter to sporadic samples with subtype labels
sporadic_ids <- sporadic_final_check$Participant_ID
overlap_final <- intersect(rownames(patient_gene_matrix), sporadic_ids)

cat("\n=== FINAL MATCHING ===\n")
cat("Sporadic samples with subtypes:", length(sporadic_ids), "\n")
cat("Matched in VCF matrix:", length(overlap_final), "\n")
cat("Unmatched sporadic samples:",
    length(setdiff(sporadic_ids, rownames(patient_gene_matrix))), "\n")

# Subset to matched sporadic samples
patient_gene_matrix_sporadic <- patient_gene_matrix[
  rownames(patient_gene_matrix) %in% overlap_final, ]

# Add zero rows for any sporadic samples not in VCF
missing_samples <- setdiff(sporadic_ids,
                           rownames(patient_gene_matrix))
if (length(missing_samples) > 0) {
  cat("Adding", length(missing_samples),
      "sporadic samples with zero target-gene variants\n")
  zero_mat <- matrix(
    0L,
    nrow = length(missing_samples),
    ncol = ncol(patient_gene_matrix_sporadic),
    dimnames = list(missing_samples,
                    colnames(patient_gene_matrix_sporadic))
  )
  patient_gene_matrix_sporadic <- rbind(
    patient_gene_matrix_sporadic, zero_mat)
}

# Reorder rows to match sporadic_scores
keep_order <- intersect(sporadic_ids,
                        rownames(patient_gene_matrix_sporadic))
patient_gene_matrix_sporadic <- patient_gene_matrix_sporadic[
  keep_order, ]

# ============================================================
# PART G: SANITY CHECKS
# ============================================================

cat("\n=== SANITY CHECKS ===\n")

# How many genes have variants in at least one patient?
genes_with_vars <- colSums(patient_gene_matrix_sporadic > 0)
cat("Genes with ≥1 variant in ≥1 patient:",
    sum(genes_with_vars > 0), "/",
    ncol(patient_gene_matrix_sporadic), "\n")

# Genes with variants in ≥5% of patients
cat("Genes hit in ≥5% of patients:",
    sum(genes_with_vars >= 0.05 *
          nrow(patient_gene_matrix_sporadic)), "\n")

# Genes with NO variants (will be dropped from PPI)
empty_genes <- names(genes_with_vars[genes_with_vars == 0])
cat("Genes with zero variants (excluded from PPI):",
    length(empty_genes), "\n")

# Remove empty genes from matrix
patient_gene_matrix_sporadic <- patient_gene_matrix_sporadic[,
                                                             genes_with_vars > 0]

cat("\nFinal matrix dimensions:",
    nrow(patient_gene_matrix_sporadic), "patients x",
    ncol(patient_gene_matrix_sporadic), "genes\n")

# Variants per patient
patient_totals <- rowSums(patient_gene_matrix_sporadic > 0)
cat("Median genes hit per patient:", median(patient_totals), "\n")
cat("Range:", min(patient_totals), "-", max(patient_totals), "\n")
cat("Patients with zero target-gene variants:",
    sum(patient_totals == 0), "\n")

# Coverage by axis
genes_hit <- colnames(patient_gene_matrix_sporadic)
cat("\nVariant coverage by gene tier:\n")
cat("  Metabolic axis genes with variants:",
    sum(genes_hit %in% metabolic_genes), "/",
    length(metabolic_genes), "\n")
cat("  Inflammation axis genes with variants:",
    sum(genes_hit %in% inflam_genes), "/",
    length(inflam_genes), "\n")
cat("  Cell cycle axis genes with variants:",
    sum(genes_hit %in% cellcycle_genes), "/",
    length(cellcycle_genes), "\n")
cat("  ALS panel genes with variants:",
    sum(genes_hit %in% als_panel_genes), "/",
    length(als_panel_genes), "\n")

# ============================================================
# SAVE
# ============================================================

saveRDS(patient_gene_matrix_sporadic,
        "patient_gene_matrix_sporadic.rds")
write.csv(patient_gene_matrix_sporadic,
          "patient_gene_matrix_sporadic.csv")

genes_with_variants_final <- colnames(patient_gene_matrix_sporadic)
saveRDS(genes_with_variants_final, "genes_with_variants.rds")

cat("\n=== STEP 1 COMPLETE ===\n")
cat("Output saved: patient_gene_matrix_sporadic.rds\n")
cat("Shape:", nrow(patient_gene_matrix_sporadic),
    "patients x", ncol(patient_gene_matrix_sporadic), "genes\n")
cat("Proceed to Step 2: build motor neuron PPI\n")





# Step 2
# ============================================================
# STEP 2: BUILD MOTOR-NEURON-RELEVANT PPI SUBNETWORK
# Using STRINGdb + your iPSC MN co-expression data
# ============================================================

library(STRINGdb)
library(igraph)
library(dplyr)
library(WGCNA)

# Load Step 1 outputs
patient_gene_matrix_sporadic <- readRDS("patient_gene_matrix_sporadic.rds")
genes_with_variants_final    <- readRDS("genes_with_variants.rds")

cat("Genes with variants (nodes to include in PPI):", 
    length(genes_with_variants_final), "\n")

# ============================================================
# PART A: BUILD STRING PPI FOR YOUR GENE SET
# ============================================================

# Initialise STRINGdb — version 11.5, human, medium confidence
# score_threshold=400 to start broad, we'll filter later
string_db <- STRINGdb$new(
  version          = "12.0",
  species          = 9606,
  score_threshold  = 400,
  network_type = "full",
  input_directory  = "stringdb_cache"  # cached locally after first download
)

cat("Mapping", length(genes_with_variants_final), 
    "genes to STRING IDs...\n")

# Map gene symbols to STRING protein IDs
gene_df <- data.frame(gene = genes_with_variants_final,
                      stringsAsFactors = FALSE)
mapped  <- string_db$mp(gene_df$gene)

# How many mapped successfully?
cat("Genes mapped to STRING:", nrow(mapped), "/",
    length(genes_with_variants_final), "\n")

# Get interactions between mapped proteins
cat("Fetching STRING interactions...\n")
interactions <- string_db$get_interactions(mapped)

cat("Raw STRING interactions found:", nrow(interactions), "\n")
cat("Score distribution:\n")
print(quantile(interactions$combined_score))

# ============================================================
# PART B: BUILD CO-EXPRESSION NETWORK FROM YOUR iPSC MN DATA
# This replaces the brain-specific PPI with ALS motor neuron data
# ============================================================

cat("\nBuilding co-expression network from iPSC MN expression...\n")

# Your SV-adjusted expression matrix: genes x samples
# Subset to genes that are both expressed AND have variants
library(AnnotationDbi)
library(org.Hs.eg.db)

# expression matrix rownames
ensembl_ids <- rownames(counts_filtered)

# remove version suffix if present, e.g. ENSG00000123456.12 -> ENSG00000123456
ensembl_ids_clean <- sub("\\..*$", "", ensembl_ids)

# map Ensembl -> gene symbol
gene_annot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_ids_clean),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

# keep non-missing symbols
gene_annot <- gene_annot[!is.na(gene_annot$SYMBOL) & gene_annot$SYMBOL != "", ]

# symbols represented in your expression matrix
expressed_symbols <- unique(gene_annot$SYMBOL)

# intersect with your variant genes
target_expressed <- intersect(genes_with_variants_final, expressed_symbols)

cat("Genes with variants AND expressed in iPSC MNs:",
    length(target_expressed), "/",
    length(genes_with_variants_final), "\n")

# keep one symbol per Ensembl row
gene_annot2 <- gene_annot[match(ensembl_ids_clean, gene_annot$ENSEMBL), ]

# rows whose mapped symbol is in target_expressed
keep <- !is.na(gene_annot2$SYMBOL) & gene_annot2$SYMBOL %in% target_expressed

counts_target <- counts_filtered[keep, , drop = FALSE]

# optional: replace rownames with gene symbols
rownames(counts_target) <- gene_annot2$SYMBOL[keep]


# keep only rows in mapping file for sporadic participants
map_sporadic <- sporadic_final_check[
  sporadic_final_check$Participant_ID %in% sporadic_sample_ids,
  c("Participant_ID", "Sample_ID")
]

# if there may be multiple Sample_IDs per participant, decide how to handle that
# for now take the first Sample_ID per participant
map_sporadic <- map_sporadic[!duplicated(map_sporadic$Participant_ID), ]

# reorder mapping to match patient_gene_matrix_sporadic row order
map_sporadic <- map_sporadic[match(sporadic_sample_ids, map_sporadic$Participant_ID), ]

# keep only mapped samples that exist in counts_target
keep <- !is.na(map_sporadic$Sample_ID) & map_sporadic$Sample_ID %in% colnames(counts_target)
map_sporadic <- map_sporadic[keep, ]

# subset counts_target in matching order
expr_sporadic_target <- counts_target[, map_sporadic$Sample_ID, drop = FALSE]

# optional: rename columns to Participant_ID so they align with mutation matrix rows
colnames(expr_sporadic_target) <- map_sporadic$Participant_ID

dim(expr_sporadic_target)
head(colnames(expr_sporadic_target))

cat("Expression matrix for co-expression:",
    nrow(expr_sporadic_target), "genes x",
    ncol(expr_sporadic_target), "samples\n")

# Check for genes with too-low variance (remove before WGCNA)
gene_vars <- apply(expr_sporadic_target, 1, var)
low_var_genes <- names(gene_vars[gene_vars < quantile(gene_vars, 0.05)])
cat("Removing", length(low_var_genes), 
    "low-variance genes before co-expression\n")

expr_for_coexpr <- expr_sporadic_target[
  !rownames(expr_sporadic_target) %in% low_var_genes, ]

# Compute pairwise Pearson correlation between genes
# Transpose: WGCNA expects samples x genes
cat("Computing gene-gene correlation matrix...\n")
cor_matrix <- cor(t(expr_for_coexpr), method = "pearson",
                  use = "pairwise.complete.obs")

cat("Correlation matrix dimensions:", dim(cor_matrix), "\n")

# Define co-expression edges: |r| > 0.3
# This is the connectivity threshold
CO_EXPR_THRESHOLD <- 0.3

coexpr_adj <- (abs(cor_matrix) > CO_EXPR_THRESHOLD) * 1
diag(coexpr_adj) <- 0

cat("Co-expression edges at |r| >", CO_EXPR_THRESHOLD, ":",
    sum(coexpr_adj) / 2, "\n")

# ============================================================
# PART C: INTERSECT STRING AND CO-EXPRESSION NETWORKS
# Keep edges supported by BOTH data sources
# ============================================================

cat("\nIntersecting STRING and co-expression networks...\n")

# Build STRING edge set for our gene subset
# Map STRING IDs back to gene symbols
# 1. map genes to STRING
mapped_df <- string_db$map(
  data.frame(gene = unique(target_expressed), stringsAsFactors = FALSE),
  "gene",
  takeFirst = TRUE,
  removeUnmappedRows = TRUE
)

# 2. get interactions
interactions <- string_db$get_interactions(mapped_df$STRING_id)

# 3. convert STRING IDs back to symbols
id_to_gene <- setNames(mapped_df$gene, mapped_df$STRING_id)
interactions$gene1 <- id_to_gene[interactions$from]
interactions$gene2 <- id_to_gene[interactions$to]

# 4. clean
interactions <- interactions[
  !is.na(interactions$gene1) & !is.na(interactions$gene2),
]

cat("Mapped genes:", nrow(mapped_df), "\n")
cat("STRING interactions:", nrow(interactions), "\n")

# Convert STRING interactions to gene-gene pairs
interactions_genes <- interactions |>
  mutate(
    gene_from = string_gene_map_rev[from],
    gene_to   = string_gene_map_rev[to]
  ) |>
  filter(!is.na(gene_from), !is.na(gene_to)) |>
  # Filter by confidence score
  filter(combined_score >= 400)



library(dplyr)
library(igraph)

# build STRING ID -> gene lookup
string_gene_map_rev <- setNames(mapped_df$gene, mapped_df$STRING_id)

# convert STRING interactions to gene-gene edges
interactions_genes <- interactions |>
  mutate(
    gene_from = string_gene_map_rev[from],
    gene_to   = string_gene_map_rev[to]
  ) |>
  filter(!is.na(gene_from), !is.na(gene_to)) |>
  filter(combined_score >= 400) |>
  mutate(
    edge_a = pmin(gene_from, gene_to),
    edge_b = pmax(gene_from, gene_to)
  ) |>
  select(edge_a, edge_b, combined_score) |>
  distinct()

cat("STRING edges (score ≥ 400):", nrow(interactions_genes), "\n")

# build coexpression edge set from 0/1 adjacency
coexpr_edges_df <- which(coexpr_adj == 1, arr.ind = TRUE) |>
  as.data.frame() |>
  mutate(
    gene_r = rownames(coexpr_adj)[row],
    gene_c = rownames(coexpr_adj)[col]
  ) |>
  filter(row < col) |>
  mutate(
    edge_a = pmin(gene_r, gene_c),
    edge_b = pmax(gene_r, gene_c)
  ) |>
  select(edge_a, edge_b) |>
  distinct()

cat("Co-expression edges:", nrow(coexpr_edges_df), "\n")

validated_edges <- inner_join(
  interactions_genes,
  coexpr_edges_df,
  by = c("edge_a", "edge_b")
)

cat("Validated edges (STRING + co-expressed):", nrow(validated_edges), "\n")

if (nrow(validated_edges) == 0) {
  cat("No validated edges found. Falling back to STRING-only network.\n")
  
  als_mn_network_lcc <- graph_from_data_frame(
    interactions_genes[, c("edge_a", "edge_b")],
    directed = FALSE
  )
  
  comp2 <- components(als_mn_network_lcc)
  lcc2  <- names(comp2$membership[comp2$membership == which.max(comp2$csize)])
  als_mn_network_lcc <- induced_subgraph(als_mn_network_lcc, lcc2)
  
} else {
  als_mn_network <- graph_from_data_frame(
    validated_edges[, c("edge_a", "edge_b")],
    directed = FALSE
  )
  
  comp <- components(als_mn_network)
  lcc_nodes <- names(comp$membership[comp$membership == which.max(comp$csize)])
  als_mn_network_lcc <- induced_subgraph(als_mn_network, lcc_nodes)
  
  if (ecount(als_mn_network_lcc) < 50) {
    cat("Validated network too sparse. Falling back to STRING-only network.\n")
    
    als_mn_network_lcc <- graph_from_data_frame(
      interactions_genes[, c("edge_a", "edge_b")],
      directed = FALSE
    )
    
    comp2 <- components(als_mn_network_lcc)
    lcc2  <- names(comp2$membership[comp2$membership == which.max(comp2$csize)])
    als_mn_network_lcc <- induced_subgraph(als_mn_network_lcc, lcc2)
  }
}

network_genes <- V(als_mn_network_lcc)$name
genes_in_both <- intersect(network_genes, colnames(patient_gene_matrix_sporadic))

patient_gene_matrix_network <- patient_gene_matrix_sporadic[, genes_in_both, drop = FALSE]

als_mn_network_final <- induced_subgraph(
  als_mn_network_lcc,
  intersect(V(als_mn_network_lcc)$name, genes_in_both)
)

cat("Final network:", vcount(als_mn_network_final), "nodes,", ecount(als_mn_network_final), "edges\n")
cat("Patient matrix:", nrow(patient_gene_matrix_network), "patients x", ncol(patient_gene_matrix_network), "genes\n")
saveRDS(mapped_df, "mapped_df.rds")
saveRDS(interactions_genes, "interactions_genes.rds")
saveRDS(coexpr_edges_df, "coexpr_edges_df.rds")
saveRDS(validated_edges, "validated_edges.rds")
saveRDS(als_mn_network_final, "als_mn_network_final.rds")
saveRDS(patient_gene_matrix_network, "patient_gene_matrix_network.rds")
# ============================================================
# PART D: ASSESS NETWORK AND DECIDE ON STRATEGY
# ============================================================

# Build the validated network
als_mn_network <- igraph::graph_from_data_frame(
  validated_edges[, c("edge_a", "edge_b")],
  directed = FALSE
)

cat("\n=== VALIDATED MOTOR NEURON PPI ===\n")
cat("Nodes:", igraph::vcount(als_mn_network), "\n")
cat("Edges:", igraph::ecount(als_mn_network), "\n")
cat("Density:", round(igraph::edge_density(als_mn_network), 4), "\n")

# Check connected components
components <- igraph::components(als_mn_network)
cat("Connected components:", components$no, "\n")
cat("Largest component size:", max(components$csize), "nodes\n")

# Keep only largest connected component for clustering
lcc_nodes <- names(components$membership[
  components$membership == which.max(components$csize)])
als_mn_network_lcc <- igraph::induced_subgraph(
  als_mn_network, lcc_nodes)

cat("\nLargest connected component:\n")
cat("  Nodes:", igraph::vcount(als_mn_network_lcc), "\n")
cat("  Edges:", igraph::ecount(als_mn_network_lcc), "\n")

# ============================================================
# PART E: FALLBACK IF VALIDATED NETWORK IS TOO SPARSE
# If validated edges < 50, use STRING-only at lower threshold
# ============================================================

if (igraph::ecount(als_mn_network_lcc) < 50) {
  
  cat("\nWARNING: Validated network too sparse for clustering\n")
  cat("Falling back to STRING-only network at score ≥ 400\n")
  
  als_mn_network_lcc <- igraph::graph_from_data_frame(
    interactions_genes[, c("edge_a", "edge_b")],
    directed = FALSE
  )
  
  # Keep largest connected component
  comp2 <- igraph::components(als_mn_network_lcc)
  lcc2  <- names(comp2$membership[
    comp2$membership == which.max(comp2$csize)])
  als_mn_network_lcc <- igraph::induced_subgraph(
    als_mn_network_lcc, lcc2)
  
  cat("STRING-only LCC:\n")
  cat("  Nodes:", igraph::vcount(als_mn_network_lcc), "\n")
  cat("  Edges:", igraph::ecount(als_mn_network_lcc), "\n")
}

# ============================================================
# PART F: ALIGN PATIENT MATRIX TO NETWORK GENES
# Only patients AND genes present in LCC go forward
# ============================================================

network_genes <- igraph::V(als_mn_network_lcc)$name
cat("\nNetwork genes available:", length(network_genes), "\n")

# Subset patient matrix to network genes only
genes_in_both <- intersect(network_genes,
                           colnames(patient_gene_matrix_sporadic))
cat("Genes in network AND patient matrix:", length(genes_in_both), "\n")

patient_gene_matrix_network <- patient_gene_matrix_sporadic[,
                                                            genes_in_both, drop = FALSE]

# Subset network to genes with patient data
als_mn_network_final <- igraph::induced_subgraph(
  als_mn_network_lcc,
  intersect(igraph::V(als_mn_network_lcc)$name, genes_in_both)
)

cat("\n=== FINAL NETWORK FOR EDGE SCORING ===\n")
cat("Nodes:", igraph::vcount(als_mn_network_final), "\n")
cat("Edges:", igraph::ecount(als_mn_network_final), "\n")

# ============================================================
# PART G: SAVE STEP 2 OUTPUTS
# ============================================================

saveRDS(als_mn_network_final,       "als_mn_network_final.rds")
saveRDS(patient_gene_matrix_network,"patient_gene_matrix_network.rds")

# Save edge list for inspection
edge_list_final <- igraph::as_data_frame(als_mn_network_final,
                                         what = "edges")
write.csv(edge_list_final, "als_mn_network_edges.csv", row.names=FALSE)

# Save node metadata: which axis does each gene belong to?
node_metadata <- data.frame(
  gene = igraph::V(als_mn_network_final)$name
) |>
  mutate(
    in_metabolic   = gene %in% metabolic_genes,
    in_inflam      = gene %in% inflam_genes,
    in_cellcycle   = gene %in% cellcycle_genes,
    in_als_panel   = gene %in% als_panel_genes,
    primary_axis   = case_when(
      in_als_panel                                      ~ "ALS_panel",
      in_metabolic & !in_inflam & !in_cellcycle        ~ "metabolic",
      in_inflam    & !in_metabolic & !in_cellcycle      ~ "inflammatory",
      in_cellcycle & !in_metabolic & !in_inflam         ~ "cell_cycle",
      TRUE                                              ~ "multi_axis"
    )
  )

write.csv(node_metadata, "als_mn_network_nodes.csv", row.names=FALSE)
cat("\nNode axis distribution:\n")
print(table(node_metadata$primary_axis))

cat("\n=== STEP 2 COMPLETE ===\n")
cat("Final network:", igraph::vcount(als_mn_network_final),
    "nodes,", igraph::ecount(als_mn_network_final), "edges\n")
cat("Patient matrix:", nrow(patient_gene_matrix_network),
    "patients x", ncol(patient_gene_matrix_network), "genes\n")
cat("Ready for Step 3: compute per-patient edge conservation scores\n")

# Quick check — do the 31 zero-variant patients cluster together?
# They will all have edge scores of exactly 1.0 everywhere
# This is fine biologically but worth knowing before interpreting clusters

# Use sporadic_final_check as the bridge
# It has both Participant_ID and Sample_ID columns
zero_subtypes <- sporadic_final_check |>
  filter(Participant_ID %in% zero_variant_patients) |>
  left_join(
    ax_clean_df |> select(sample_id, cluster_final),
    by = c("Sample_ID" = "sample_id")
  ) |>
  pull(cluster_final)

cat("Zero-variant patients subtype distribution:\n")
print(table(zero_subtypes, useNA = "always"))
cat("Total zero-variant patients:", length(zero_variant_patients), "\n")

# Overall subtype counts for comparison
overall_counts <- ax_clean_df |>
  filter(sample_id %in% sporadic_final_check$Sample_ID) |>
  count(cluster_final)

cat("Overall subtype distribution:\n")
print(overall_counts)

# Zero-variant patients as percentage of each subtype
zero_df <- data.frame(
  subtype      = c("S-ALS-1", "S-ALS-2", "S-ALS-3"),
  zero_variant = c(9, 20, 12)
) |>
  left_join(overall_counts, by = c("subtype" = "cluster_final")) |>
  mutate(pct_of_subtype = round(100 * zero_variant / n, 1))

cat("\nZero-variant patients as % of each subtype:\n")
print(zero_df)

# Step 3

library(igraph)
library(Matrix)

# ── Step 3: Per-patient edge conservation scores ──────────────────────────────

# 1. Extract edge list as gene symbol pairs
edge_list <- as_edgelist(als_mn_network_final, names = TRUE)
# Should be 171 × 2 matrix of gene symbols
stopifnot(nrow(edge_list) == 171)

# Verify all edge genes are in the variant matrix
edge_genes <- unique(c(edge_list[,1], edge_list[,2]))
missing_genes <- setdiff(edge_genes, colnames(patient_gene_matrix_network))
if (length(missing_genes) > 0) {
  cat("WARNING — genes in network but missing from variant matrix:\n")
  print(missing_genes)
} else {
  cat("All edge genes confirmed in variant matrix.\n")
}

# 2. Compute cohort-level totals per gene (denominator)
gene_totals <- colSums(patient_gene_matrix_network)
cat(sprintf("Gene totals range: %d – %d (median: %.1f)\n",
            min(gene_totals), max(gene_totals), median(gene_totals)))

# 3. Build edge names for column labels
edge_names <- paste(edge_list[,1], edge_list[,2], sep = "--")

# 4. Initialise conservation matrix
n_patients <- nrow(patient_gene_matrix_network)   # 391
n_edges    <- nrow(edge_list)                      # 171

edge_conservation_matrix <- matrix(
  NA_real_,
  nrow = n_patients,
  ncol = n_edges,
  dimnames = list(
    rownames(patient_gene_matrix_network),  # Participant_IDs
    edge_names
  )
)

# 5. Fill conservation scores
skipped_edges <- character(0)

for (e in seq_len(n_edges)) {
  gA <- edge_list[e, 1]
  gB <- edge_list[e, 2]
  
  # Cohort denominator
  denom <- gene_totals[gA] + gene_totals[gB]
  
  if (denom == 0) {
    # No cohort variants in either gene — universally conserved
    edge_conservation_matrix[, e] <- 1.0
    skipped_edges <- c(skipped_edges, edge_names[e])
    next
  }
  
  # Per-patient numerator
  numerator <- patient_gene_matrix_network[, gA] + patient_gene_matrix_network[, gB]
  
  perturbation <- numerator / denom
  edge_conservation_matrix[, e] <- 1 - perturbation
}

# 6. Sanity checks
cat(sprintf("\nEdge conservation matrix dimensions: %d × %d\n",
            nrow(edge_conservation_matrix), ncol(edge_conservation_matrix)))
cat(sprintf("NA cells: %d\n", sum(is.na(edge_conservation_matrix))))
cat(sprintf("Edges with zero cohort variants (universally conserved): %d\n",
            length(skipped_edges)))
cat(sprintf("Conservation score range: %.4f – %.4f\n",
            min(edge_conservation_matrix, na.rm=TRUE),
            max(edge_conservation_matrix, na.rm=TRUE)))
cat(sprintf("Median conservation score: %.4f\n",
            median(edge_conservation_matrix, na.rm=TRUE)))

# 7. Distribution of perturbation across all cells
perturbation_all <- 1 - edge_conservation_matrix
cat(sprintf("\nFraction of patient-edge pairs with any perturbation: %.1f%%\n",
            100 * mean(perturbation_all > 0, na.rm=TRUE)))

# 8. Per-patient total perturbation load (sum across edges)
patient_perturbation_load <- rowSums(1 - edge_conservation_matrix)
cat(sprintf("\nPer-patient total perturbation load (sum across 171 edges):\n"))
print(summary(patient_perturbation_load))

# 9. Quick check: do the 31 zero-variant patients show conservation = 1 everywhere?
# (These patients had zero variants in target genes)
# You'll need to cross-reference with your zero-variant list from earlier
zero_var_patients <- rownames(patient_gene_matrix_network)[
  rowSums(patient_gene_matrix_network) == 0
]
if (length(zero_var_patients) > 0) {
  cat(sprintf("\nZero-variant patients: %d — verifying all-1 conservation...\n",
              length(zero_var_patients)))
  all_ones <- all(edge_conservation_matrix[zero_var_patients, ] == 1.0)
  cat(sprintf("All conservation = 1.0 for zero-variant patients: %s\n", all_ones))
}

# 10. Save for Step 4
saveRDS(edge_conservation_matrix, "edge_conservation_matrix.rds")
cat("\nSaved: edge_conservation_matrix.rds\n")

# ── Step 4 Prep: Edge variance filtering ──────────────────────────────────────

# Per-edge variance in conservation scores
edge_variance <- apply(edge_conservation_matrix, 2, var)

cat(sprintf("Edge variance range: %.6f – %.6f\n", min(edge_variance), max(edge_variance)))
cat(sprintf("Edges with zero variance (constant = 1.0 across all patients): %d\n",
            sum(edge_variance == 0)))
cat(sprintf("Edges with var > 0: %d\n", sum(edge_variance > 0)))

# Distribution of variance among non-zero edges
non_zero_var <- edge_variance[edge_variance > 0]
cat(sprintf("\nAmong informative edges (var > 0):\n"))
print(summary(non_zero_var))

# How many patients are perturbed per edge?
patients_perturbed_per_edge <- colSums(edge_conservation_matrix < 1.0)
cat(sprintf("\nPatients perturbed per edge (conservation < 1.0):\n"))
print(summary(patients_perturbed_per_edge))
cat(sprintf("Edges with ≥5 patients perturbed: %d\n",  sum(patients_perturbed_per_edge >= 5)))
cat(sprintf("Edges with ≥10 patients perturbed: %d\n", sum(patients_perturbed_per_edge >= 10)))
cat(sprintf("Edges with ≥20 patients perturbed: %d\n", sum(patients_perturbed_per_edge >= 20)))

# Check: which subtype are the high perturbation-load patients?
# Merge with cluster assignments
patient_load_df <- data.frame(
  Participant_ID = names(patient_perturbation_load),
  perturbation_load = patient_perturbation_load
)
library(dplyr)

ax_clean_df <- ax_clean_df |>
  left_join(
    sporadic_final_check |>
      select(Sample_ID, Participant_ID),
    by = c("sample_id" = "Sample_ID")
  )
# Join with subtype labels from ax_clean_df
patient_load_df <- patient_load_df |>
  left_join(
    ax_clean_df |> select(Participant_ID, cluster_final)
  )

cat(sprintf("\nMean perturbation load by subtype:\n"))
patient_load_df |>
  group_by(cluster_final) |>
  summarise(
    n = n(),
    mean_load = mean(perturbation_load, na.rm=TRUE),
    median_load = median(perturbation_load, na.rm=TRUE),
    max_load = max(perturbation_load, na.rm=TRUE)
  ) |>
  print()

# ANOVA: does perturbation load differ by subtype?
aov_load <- aov(perturbation_load ~ cluster_final, data = patient_load_df)
cat(sprintf("\nANOVA perturbation load ~ subtype: p = %.4f\n",
            summary(aov_load)[[1]][["Pr(>F)"]][1]))

kruskal.test(perturbation_load ~ cluster_final, data = patient_load_df)
# ── Step 4: SNF + Spectral Clustering ─────────────────────────────────────────

# Install SNFtool if needed
# install.packages("SNFtool")
library(SNFtool)

# Filter to informative edges only (var > 0)
ecm_filtered <- edge_conservation_matrix[, edge_variance > 0]
cat(sprintf("\nFiltered conservation matrix: %d patients × %d edges\n",
            nrow(ecm_filtered), ncol(ecm_filtered)))

# Standardise columns to mean=0, sd=1 before SNF
# (SNF distance is Euclidean-based; standardisation prevents high-variance edges dominating)
ecm_scaled <- scale(ecm_filtered)

# SNF parameters
# K = neighbourhood size (~sqrt(n) or n/10; try both)
K_vals <- c(20, 39, 50)
mu_vals <- c(0.3, 0.5, 0.8)
k_clusters <- 2:5

results_grid <- list()
idx <- 1

for (K in K_vals) {
  for (mu in mu_vals) {
    # Affinity matrix from conservation score matrix
    # SNFtool::affinityMatrix expects a distance matrix
    dist_mat <- as.matrix(dist(ecm_scaled, method = "euclidean"))
    W <- affinityMatrix(dist_mat, K = K, sigma = mu)
    
    for (k in k_clusters) {
      labels <- spectralClustering(W, k)
      sil <- mean(cluster::silhouette(labels, dist_mat)[, 3])
      ari_val <- mclust::adjustedRandIndex(labels,
                                           patient_load_df$cluster_final[match(rownames(ecm_scaled),
                                                                               patient_load_df$Participant_ID)])
      
      results_grid[[idx]] <- data.frame(
        K = K, mu = mu, k = k,
        silhouette = round(sil, 4),
        ARI_vs_transcriptomic = round(ari_val, 4)
      )
      idx <- idx + 1
    }
  }
}

results_df <- do.call(rbind, results_grid)
cat(sprintf("\nSNF grid search results (%d parameter combinations):\n", nrow(results_df)))
print(results_df[order(-results_df$silhouette), ])

# Best by silhouette
best <- results_df[which.max(results_df$silhouette), ]
cat(sprintf("\nBest parameters: K=%d, mu=%.1f, k=%d\n", best$K, best$mu, best$k))
cat(sprintf("  Silhouette = %.4f\n", best$silhouette))
cat(sprintf("  ARI vs transcriptomic = %.4f\n", best$ARI_vs_transcriptomic))


# ── Diagnose the high-silhouette k=5 solution ─────────────────────────────────

# Refit the best solution explicitly
dist_mat_final <- as.matrix(dist(ecm_scaled, method = "euclidean"))
W_best <- affinityMatrix(dist_mat_final, K = 50, sigma = 0.3)
labels_k5 <- spectralClustering(W_best, 5)

# Cluster sizes
cat("Cluster sizes at best k=5 solution:\n")
print(table(labels_k5))

# How many zero-variant patients fall in each cluster?
zero_var_idx <- which(rownames(ecm_scaled) %in% zero_var_patients)
cat(sprintf("\nZero-variant patients (n=%d) by cluster:\n", length(zero_var_idx)))
print(table(labels_k5[zero_var_idx]))

# Identify the zero-variant cluster (should be one dominant cluster)
zero_cluster <- as.integer(names(which.max(table(labels_k5[zero_var_idx]))))
cat(sprintf("Dominant zero-variant cluster: Cluster %d\n", zero_cluster))
cat(sprintf("Fraction of zero-variant patients in this cluster: %.1f%%\n",
            100 * mean(labels_k5[zero_var_idx] == zero_cluster)))

# What is the within-cluster mean conservation for each cluster?
for (cl in 1:5) {
  members <- labels_k5 == cl
  cat(sprintf("Cluster %d (n=%d): mean conservation = %.4f, mean perturbation load = %.4f\n",
              cl, sum(members),
              mean(edge_conservation_matrix[members, ]),
              mean(patient_perturbation_load[members])))
}

# Cross-tabulate genetic clusters vs transcriptomic subtypes
cat("\nGenetic clusters (rows) × Transcriptomic subtypes (cols):\n")
subtype_labels <- patient_load_df$cluster_final[
  match(rownames(ecm_scaled), patient_load_df$Participant_ID)
]
print(table(Genetic = labels_k5, Transcriptomic = subtype_labels))

# ARI at k=3 (forcing match to transcriptomic k) for comparison
labels_k3 <- spectralClustering(W_best, 3)
ari_k3 <- mclust::adjustedRandIndex(labels_k3, subtype_labels)
cat(sprintf("\nARI at k=3 (K=50, mu=0.3): %.4f\n", ari_k3))

# ── Summary stats for paper ────────────────────────────────────────────────────
cat("\n--- PAPER-READY SUMMARY ---\n")
cat(sprintf("Edge conservation matrix: %d patients × %d edges\n",
            nrow(edge_conservation_matrix), ncol(edge_conservation_matrix)))
cat(sprintf("Patient-edge pairs with any perturbation: %.1f%%\n",
            100 * mean(edge_conservation_matrix < 1.0)))
cat(sprintf("Perturbation load by subtype ANOVA: p = 0.759\n"))
cat(sprintf("SNF grid search (K ∈ {20,39,50}, μ ∈ {0.3,0.5,0.8}, k ∈ {2..5}):\n"))
cat(sprintf("  ARI range across all 36 combinations: %.4f to %.4f\n",
            min(results_df$ARI_vs_transcriptomic),
            max(results_df$ARI_vs_transcriptomic)))
cat(sprintf("  Best silhouette: %.4f (k=%d) — likely zero-variant cluster artifact\n",
            max(results_df$silhouette), results_df$k[which.max(results_df$silhouette)]))
