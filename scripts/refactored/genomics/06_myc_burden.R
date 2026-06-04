# ============================================================
# Refactored entrypoint for original script: TRICBurdenMYCBurden.R
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
if (isTRUE(ALS_ARGS$help)) { als_print_help("06_myc_burden"); quit(save = "no", status = 0) }
CFG <- als_load_config(ALS_ARGS$config)
ALSPROJ_BASE_DIR <- als_choose(ALS_ARGS$base_dir, CFG, "project.base_dir", getwd())
ALSPROJ_OUTDIR <- als_choose(ALS_ARGS$outdir, CFG, "outputs.root", file.path(getwd(), "results", "06_myc_burden"))
ALSPROJ_SCRIPT_ID <- "06_myc_burden"
als_init(ALSPROJ_SCRIPT_ID, ALSPROJ_BASE_DIR, ALSPROJ_OUTDIR, ALS_ARGS$redirect_outputs)
als_register_output_wrappers()
als_load_objects(ALS_ARGS$objects)
on.exit(als_write_session_info(ALSPROJ_SCRIPT_ID), add = TRUE)
# ---- End bootstrap ----

# ============================================================
# HALLMARK_MYC_TARGETS_V1 pathway rare variant burden tests
# Drop-in script, parallel to tric_burden_tests.R
# ============================================================
# Logic: CCT2 (a direct MYC target) shows zero rare coding burden but
# clear expression variation across subtypes. The expression variation
# must therefore be trans-regulatory. This script tests whether rare
# coding burden anywhere in the MYC target network -- the upstream
# regulators and co-regulated genes that could plausibly drive CCT2
# expression -- differs by subtype.
#
# Five tests (parallel to TRiC analysis):
#   1. Pathway-wide burden by subtype (full + clean)
#   2. Pan-MYC carrier rate by subtype
#   3. Per-gene Fisher within clean subset (FDR-corrected)
#   4. Logistic: pathway burden -> P(S-ALS-3), age + sex adjusted
#   5. Pathway burden -> CCT2 expression (does upstream burden
#      predict the downstream regulatory signal we already see?)
#
# Assumes already in workspace:
#   - use_cols, col_to_id      (sample column infrastructure)
#   - rare_burden_sub          (subtype + clinical, keyed by sample)
#   - cct2_df                  (CCT2 expression by Participant_ID)
#   - axis_sets_clean OR similar object containing the Hallmark
#     gene sets (we read MYC_TARGETS_V1 from msigdbr if absent)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(broom)
})

# ---------------- Config ----------------
EXON_DIR <- file.path(ALSPROJ_BASE_DIR)
RARE_MAF <- 0.01

damaging_classes <- c("nonsynonymous_SNV", "stopgain", "stoploss",
                      "frameshift_deletion", "frameshift_insertion")
lof_classes      <- c("stopgain", "stoploss",
                      "frameshift_deletion", "frameshift_insertion")

# ---------------- Step 0: load MYC target gene set --------------
# Prefer the in-workspace object if it already contains MYC_TARGETS_V1;
# fall back to msigdbr. This avoids version drift across sessions.
get_myc_targets <- function() {
  if (exists("axis_sets_clean")) {
    nm <- names(axis_sets_clean)
    hit <- grep("MYC_TARGETS_V1", nm, value = TRUE, ignore.case = TRUE)
    if (length(hit)) return(unique(axis_sets_clean[[hit[1]]]))
  }
  if (exists("axis_sets")) {
    nm <- names(axis_sets)
    hit <- grep("MYC_TARGETS_V1", nm, value = TRUE, ignore.case = TRUE)
    if (length(hit)) return(unique(axis_sets[[hit[1]]]))
  }
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Need MYC_TARGETS_V1 -- not in axis_sets and msigdbr unavailable")
  m <- msigdbr::msigdbr(species = "Homo sapiens",
                        category = "H")
  unique(m$gene_symbol[m$gs_name == "HALLMARK_MYC_TARGETS_V1"])
}

MYC_GENES <- get_myc_targets()
cat("MYC_TARGETS_V1 genes loaded:", length(MYC_GENES), "\n")
cat("Includes CCT2 (sanity check):", "CCT2" %in% MYC_GENES, "\n")

# ---------------- Step 1: parse rare damaging MYC variants ----------
# Same parser logic as TRiC, but reads ALL chromosomes since
# MYC targets are scattered genome-wide (~200 genes).
parse_myc_chr <- function(chr,
                          genes      = MYC_GENES,
                          classes    = damaging_classes,
                          maf_thresh = RARE_MAF) {
  f <- file.path(EXON_DIR, paste0("ALS.", chr, ".anno_report_exon.txt"))
  if (!file.exists(f)) { message("Missing: ", f); return(NULL) }
  
  dt <- fread(f, sep = "\t", header = TRUE, na.strings = c(".", "NA", ""))
  
  if (!"Gene" %in% names(dt)) stop("Expected 'Gene' in ", f)
  dt <- dt[Gene %in% genes]
  if (nrow(dt) == 0) return(NULL)
  
  if (!"ExonicFunc" %in% names(dt)) stop("Expected 'ExonicFunc' in ", f)
  dt <- dt[ExonicFunc %in% classes]
  if (nrow(dt) == 0) return(NULL)
  
  maf_candidates <- c("Max Freq of 3", "gnomAD_genome_ALL", "gnomAD_exome_ALL",
                      "AF_popmax", "ExAC_ALL", "AF")
  maf_col <- intersect(maf_candidates, names(dt))[1]
  if (is.na(maf_col))
    stop("No MAF column in ", f, "\nFound: ", paste(names(dt), collapse=", "))
  
  af <- suppressWarnings(as.numeric(dt[[maf_col]]))
  af[is.na(af)] <- 0
  dt <- dt[af < maf_thresh]
  dt
}

all_chrs <- paste0("chr", c(1:22, "X"))
cat("\nReading ", length(all_chrs), " exon files...\n")
myc_variants <- rbindlist(
  lapply(all_chrs, parse_myc_chr),
  use.names = TRUE, fill = TRUE
)

cat("Rare damaging MYC-target variants:", nrow(myc_variants), "\n")
cat("Genes hit:", length(unique(myc_variants$Gene)), "of", length(MYC_GENES), "\n")
cat("Top 15 genes by variant count:\n")
print(head(sort(table(myc_variants$Gene), decreasing = TRUE), 15))
cat("\nBy variant class:\n"); print(table(myc_variants$ExonicFunc))

# ---------------- Step 2: per-sample carrier matrix --------------
stopifnot(exists("use_cols"))

# Same genotype parser as TRiC script -- partial-missing GTs treated
# conservatively as non-carriers.
is_carrier <- function(x) {
  x <- as.character(x)
  hom_ref  <- x %in% c("0", "0/0", "0|0")
  has_miss <- is.na(x) | x == "" | grepl("\\.", x)
  has_alt  <- grepl("[1-9]", x) & !has_miss
  as.integer(has_alt & !hom_ref)
}

build_burden_matrix <- function(variants, samples,
                                classes = damaging_classes,
                                gene_universe = MYC_GENES) {
  v <- variants[ExonicFunc %in% classes]
  if (nrow(v) == 0) {
    return(matrix(0L, nrow = length(gene_universe), ncol = length(samples),
                  dimnames = list(gene_universe, samples)))
  }
  carr <- vapply(samples,
                 function(s) as.integer(is_carrier(v[[s]])),
                 integer(nrow(v)))
  rownames(carr) <- v$Gene
  gb <- rowsum(carr, group = rownames(carr))
  miss <- setdiff(gene_universe, rownames(gb))
  if (length(miss)) {
    pad <- matrix(0L, nrow = length(miss), ncol = ncol(gb),
                  dimnames = list(miss, colnames(gb)))
    gb <- rbind(gb, pad)
  }
  gb[gene_universe, , drop = FALSE]
}

myc_burden_dmg <- build_burden_matrix(myc_variants, use_cols, damaging_classes)
myc_burden_lof <- build_burden_matrix(myc_variants, use_cols, lof_classes)

# Pathway-level totals per sample
myc_total_dmg <- colSums(myc_burden_dmg)
myc_total_lof <- colSums(myc_burden_lof)

cat("\nPathway-level damaging burden distribution:\n")
print(summary(myc_total_dmg))
cat("Carriers of >=1 rare damaging MYC-target variant:",
    sum(myc_total_dmg > 0), "/", length(myc_total_dmg), "\n")
cat("Carriers of >=1 rare LoF MYC-target variant:",
    sum(myc_total_lof > 0), "/", length(myc_total_lof), "\n")

# ---------------- Step 3: attach subtype + clinical ----------------
stopifnot(exists("col_to_id"), exists("rare_burden_sub"))

myc_df <- data.frame(
  sample_col      = colnames(myc_burden_dmg),
  cgnd_id         = unname(col_to_id[colnames(myc_burden_dmg)]),
  MYC_burden      = myc_total_dmg,
  MYC_burden_lof  = myc_total_lof,
  stringsAsFactors = FALSE
)
myc_df$MYC_carrier     <- as.integer(myc_df$MYC_burden > 0)
myc_df$MYC_carrier_lof <- as.integer(myc_df$MYC_burden_lof > 0)

key_candidates <- c("cgnd_id", "CGND_HDA_ID", "NYGC_CGND_ID",
                    "Participant_ID", "sample_id")
key_col <- intersect(key_candidates, names(rare_burden_sub))[1]
if (is.na(key_col))
  stop("No join key in rare_burden_sub. Cols: ",
       paste(names(rare_burden_sub), collapse = ", "))
cat("\nJoining myc_df to rare_burden_sub on:", key_col, "\n")

myc_df <- myc_df %>%
  rename(!!key_col := cgnd_id) %>%
  inner_join(
    rare_burden_sub %>%
      dplyr::select(any_of(c(key_col, "Participant_ID",
                      "cluster_clean", "cluster_final",
                      "AGE_AT_SYMPTOM_ONSET", "Age_at_Onset",
                      "SEX", "Sex"))),
    by = key_col
  )

if (!"cluster_clean" %in% names(myc_df))
  myc_df$cluster_clean <- myc_df$cluster_final
if (!"AGE" %in% names(myc_df)) {
  age_col <- intersect(c("AGE_AT_SYMPTOM_ONSET", "Age_at_Onset"),
                       names(myc_df))[1]
  if (!is.na(age_col)) myc_df$AGE <- myc_df[[age_col]]
}
if (!"SEX" %in% names(myc_df)) {
  sex_col <- intersect(c("Sex"), names(myc_df))[1]
  if (!is.na(sex_col)) myc_df$SEX <- myc_df[[sex_col]]
}

myc_df <- myc_df %>% filter(!is.na(cluster_clean))
myc_df$cluster_clean <- factor(myc_df$cluster_clean,
                               levels = c("S-ALS-1", "S-ALS-2", "S-ALS-3"))
myc_df$MYC_carrier     <- factor(myc_df$MYC_carrier,     levels = c(0, 1))
myc_df$MYC_carrier_lof <- factor(myc_df$MYC_carrier_lof, levels = c(0, 1))

cat("\nFinal analysis frame: n =", nrow(myc_df), "\n")
cat("Subtype distribution:\n"); print(table(myc_df$cluster_clean))

myc_clean <- myc_df %>%
  filter(cluster_clean %in% c("S-ALS-2", "S-ALS-3")) %>%
  mutate(cluster_clean = droplevels(cluster_clean),
         is_S3 = as.integer(cluster_clean == "S-ALS-3"))

# ============================================================
# TEST 1 -- Pathway-wide burden by subtype (count comparison)
# ============================================================
cat("\n=== TEST 1: MYC-target pathway-wide burden ===\n")

cat("\nMean burden (rare damaging count) by subtype:\n")
print(myc_df %>% group_by(cluster_clean) %>%
        summarise(mean_burden = mean(MYC_burden),
                  median_burden = median(MYC_burden),
                  sd          = sd(MYC_burden),
                  n           = n()))

cat("\nANOVA on burden count:\n")
print(summary(aov(MYC_burden ~ cluster_clean, data = myc_df)))

cat("\nKruskal-Wallis (non-parametric, robust to overdispersion):\n")
print(kruskal.test(MYC_burden ~ cluster_clean, data = myc_df))

if (requireNamespace("MASS", quietly = TRUE)) {
  cat("\nNegative binomial regression on pathway burden:\n")
  print(summary(MASS::glm.nb(MYC_burden ~ cluster_clean, data = myc_df)))
}

# Clean subset 2-group test
cat("\nClean subset (S-ALS-2 vs S-ALS-3):\n")
print(wilcox.test(MYC_burden ~ cluster_clean, data = myc_clean))
print(myc_clean %>% group_by(cluster_clean) %>%
        summarise(mean_burden = mean(MYC_burden), n = n()))

# ============================================================
# TEST 2 -- Pan-MYC carrier rate (binary, any rare variant)
# ============================================================
cat("\n=== TEST 2: Pan-MYC carrier rate by subtype ===\n")

t2_full <- table(myc_df$cluster_clean, myc_df$MYC_carrier,
                 dnn = c("subtype","MYC_carrier"))
cat("Full cohort 3x2:\n"); print(t2_full)
cat("Carrier rates:\n"); print(prop.table(t2_full, margin = 1))
print(chisq.test(t2_full, simulate.p.value = TRUE, B = 10000))

t2_clean <- table(myc_clean$cluster_clean, myc_clean$MYC_carrier,
                  dnn = c("subtype","MYC_carrier"))
cat("\nClean subset Fisher:\n"); print(t2_clean)
print(fisher.test(t2_clean))

# LoF-only is the cleaner test (missense burden noisy in 200-gene set)
t2_lof <- table(myc_df$cluster_clean, myc_df$MYC_carrier_lof,
                dnn = c("subtype","MYC_carrier_lof"))
cat("\nLoF-only carrier rate (full cohort):\n"); print(t2_lof)
cat("Carrier rates:\n"); print(prop.table(t2_lof, margin = 1))
print(chisq.test(t2_lof, simulate.p.value = TRUE, B = 10000))

# ============================================================
# TEST 3 -- Per-gene Fisher within clean subset, FDR-corrected
# This is the discovery test -- which MYC target(s), if any,
# carry differential rare burden between S-ALS-2 and S-ALS-3?
# ============================================================
cat("\n=== TEST 3: Per-gene Fisher (clean subset, BH-FDR) ===\n")

per_gene_fisher <- function(gene, mat, df) {
  carrier <- as.integer(mat[gene, df$sample_col] > 0)
  tab <- table(df$cluster_clean, carrier, dnn = c("subtype","carrier"))
  if (any(dim(tab) < 2) || !"1" %in% colnames(tab) ||
      sum(tab[,"1"]) == 0) {
    return(tibble(gene = gene, n_S2 = 0L, n_S3 = 0L,
                  rate_S2 = 0, rate_S3 = 0,
                  OR = NA_real_, p = NA_real_))
  }
  ft <- fisher.test(tab)
  n_S2 <- if ("S-ALS-2" %in% rownames(tab)) tab["S-ALS-2","1"] else 0L
  n_S3 <- if ("S-ALS-3" %in% rownames(tab)) tab["S-ALS-3","1"] else 0L
  N_S2 <- sum(df$cluster_clean == "S-ALS-2")
  N_S3 <- sum(df$cluster_clean == "S-ALS-3")
  tibble(gene = gene,
         n_S2 = n_S2, n_S3 = n_S3,
         rate_S2 = n_S2 / N_S2, rate_S3 = n_S3 / N_S3,
         OR = unname(ft$estimate),
         p  = ft$p.value)
}

# Only test genes that have at least 1 carrier across the clean subset
# (no point testing genes with zero observed variants in this cohort)
gene_has_carrier <- rowSums(
  myc_burden_dmg[, myc_clean$sample_col, drop = FALSE] > 0
) > 0
testable_genes <- MYC_GENES[gene_has_carrier[MYC_GENES]]
cat("Testable genes (>=1 carrier in clean subset):",
    length(testable_genes), "of", length(MYC_GENES), "\n")

per_gene_dmg <- bind_rows(lapply(testable_genes,
                                 per_gene_fisher,
                                 mat = myc_burden_dmg,
                                 df  = myc_clean)) %>%
  mutate(p_fdr = p.adjust(p, "BH")) %>%
  arrange(p)

cat("\nTop 20 genes (rare damaging, S-ALS-3 vs S-ALS-2):\n")
print(head(per_gene_dmg, 20))

cat("\nGenes with FDR-q < 0.10:\n")
sig_dmg <- per_gene_dmg %>% filter(p_fdr < 0.10)
if (nrow(sig_dmg) == 0) cat("  None.\n") else print(sig_dmg)

# Same scan on LoF-only burden
gene_has_lof <- rowSums(
  myc_burden_lof[, myc_clean$sample_col, drop = FALSE] > 0
) > 0
testable_lof <- MYC_GENES[gene_has_lof[MYC_GENES]]
cat("\nGenes with >=1 LoF carrier in clean subset:", length(testable_lof), "\n")

if (length(testable_lof) > 0) {
  per_gene_lof <- bind_rows(lapply(testable_lof,
                                   per_gene_fisher,
                                   mat = myc_burden_lof,
                                   df  = myc_clean)) %>%
    mutate(p_fdr = p.adjust(p, "BH")) %>%
    arrange(p)
  cat("LoF-only per-gene results:\n")
  print(per_gene_lof)
} else {
  per_gene_lof <- tibble()
  cat("(no LoF carriers in any single MYC-target gene -- expected at this N)\n")
}

# ============================================================
# TEST 4 -- Logistic: MYC pathway burden -> P(S-ALS-3)
# ============================================================
cat("\n=== TEST 4: MYC burden -> P(S-ALS-3), age + sex adjusted ===\n")

myc_clean$SEX <- factor(myc_clean$SEX)

t4a <- glm(is_S3 ~ MYC_burden, data = myc_clean, family = binomial)
t4b <- glm(is_S3 ~ MYC_burden + AGE + SEX,
           data = myc_clean, family = binomial)
t4c <- glm(is_S3 ~ MYC_burden_lof + AGE + SEX,
           data = myc_clean, family = binomial)

cat("Damaging burden, unadjusted:\n");           print(summary(t4a))
cat("\nDamaging burden, age + sex adjusted:\n"); print(summary(t4b))
cat("\nLoF-only burden, age + sex adjusted:\n"); print(summary(t4c))

# ============================================================
# TEST 5 -- MYC pathway burden -> CCT2 expression
# Does upstream rare burden predict the downstream regulatory
# signature? If MYC-target rare burden disrupts the MYC program
# globally, CCT2 expression -- as a MYC target -- should track.
# Direction tells you the mechanism:
#   beta < 0  -> burden disrupts MYC output -> CCT2 down
#   beta > 0  -> compensatory upregulation across the program
#   beta ~ 0  -> rare burden is silent at the regulatory level
# ============================================================
cat("\n=== TEST 5: MYC burden -> CCT2 expression ===\n")

stopifnot(exists("cct2_df"))
expr_col_candidates <- c("CCT2", "CCT2_expr", "expression", "expr")
expr_col <- intersect(expr_col_candidates, names(cct2_df))[1]
if (is.na(expr_col))
  stop("No CCT2 expression column in cct2_df. Looked for: ",
       paste(expr_col_candidates, collapse = ", "))

cct2_unique <- cct2_df %>%
  distinct(Participant_ID, .keep_all = TRUE) %>%
  transmute(Participant_ID, CCT2_expr = .data[[expr_col]])

if (!"Participant_ID" %in% names(myc_clean))
  stop("myc_clean missing Participant_ID -- check Step 3 join")

mech_df <- myc_clean %>%
  inner_join(cct2_unique, by = "Participant_ID") %>%
  filter(!is.na(CCT2_expr))

cat("Patients with both MYC burden and CCT2 expression:", nrow(mech_df), "\n")

# Note: include cluster_clean as covariate -- we're asking whether
# rare burden adds explanatory power BEYOND the subtype label
m1 <- lm(CCT2_expr ~ MYC_burden + AGE + SEX + cluster_clean,
         data = mech_df)
cat("\nDamaging burden, subtype-adjusted:\n"); print(summary(m1))

m2 <- lm(CCT2_expr ~ MYC_burden_lof + AGE + SEX + cluster_clean,
         data = mech_df)
cat("\nLoF burden, subtype-adjusted:\n"); print(summary(m2))

# Also test without subtype as covariate -- is there a marginal
# burden effect that gets absorbed by subtype?
m3 <- lm(CCT2_expr ~ MYC_burden + AGE + SEX, data = mech_df)
cat("\nDamaging burden, NOT subtype-adjusted (marginal):\n")
print(summary(m3))

cat("\nBurden distribution by subtype in mechanism subset:\n")
print(mech_df %>% group_by(cluster_clean) %>%
        summarise(mean_burden = mean(MYC_burden),
                  mean_lof    = mean(MYC_burden_lof),
                  n           = n()))

# ============================================================
# Save outputs
# ============================================================
out_dir <- "myc_burden_results"
dir.create(out_dir, showWarnings = FALSE)
fwrite(myc_df,        file.path(out_dir, "myc_burden_full.csv"))
fwrite(myc_clean,     file.path(out_dir, "myc_burden_clean.csv"))
fwrite(per_gene_dmg,  file.path(out_dir, "myc_per_gene_dmg.csv"))
if (nrow(per_gene_lof) > 0)
  fwrite(per_gene_lof, file.path(out_dir, "myc_per_gene_lof.csv"))
saveRDS(list(t4a = t4a, t4b = t4b, t4c = t4c,
             m1 = m1, m2 = m2, m3 = m3),
        file.path(out_dir, "myc_models.rds"))

cat("\nResults written to:", out_dir, "\n")
