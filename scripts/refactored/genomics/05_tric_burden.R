# ============================================================
# Refactored entrypoint for original script: tric_burden.R
# Use --config, --base-dir, --outdir, and --objects. See docs/SCRIPT_SPECS_FULL.md.
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
if (isTRUE(ALS_ARGS$help)) { als_print_help("05_tric_burden"); quit(save = "no", status = 0) }
CFG <- als_load_config(ALS_ARGS$config)
ALSPROJ_BASE_DIR <- als_choose(ALS_ARGS$base_dir, CFG, "project.base_dir", getwd())
ALSPROJ_OUTDIR <- als_choose(ALS_ARGS$outdir, CFG, "outputs.root", file.path(getwd(), "results", "05_tric_burden"))
ALSPROJ_SCRIPT_ID <- "05_tric_burden"
als_init(ALSPROJ_SCRIPT_ID, ALSPROJ_BASE_DIR, ALSPROJ_OUTDIR, ALS_ARGS$redirect_outputs)
als_register_output_wrappers()
als_load_objects(ALS_ARGS$objects)
on.exit(als_write_session_info(ALSPROJ_SCRIPT_ID), add = TRUE)
# ---- End bootstrap ----

# ============================================================
# CCT2 + TRiC complex rare variant burden tests
# Drop-in script, parallel to parse_gene_burden_chr / rare_burden_sub
# ============================================================
# Five tests:
#   1. CCT2 single-gene chi-square (full) + Fisher (clean subset)
#   2. TRiC complex pooled burden by subtype
#   3. Per-gene Fisher within clean subset
#   4. Logistic: TRiC burden -> P(S-ALS-3), adjusted age + sex
#   5. CCT2 burden -> CCT2 expression (LoF mechanism check)
#
# Assumes already in workspace from your existing pipeline:
#   - use_cols          (CASE-NEU sample columns matched to metadata)
#   - col_to_id         (named vector: CASE-NEU column -> CGND-HDA ID)
#   - rare_burden_sub   (carries cluster_clean + clinical, keyed by sample id)
#   - cct2_df           (CCT2 expression keyed on Participant_ID)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(broom)
})

# ---------------- Config ----------------
TRIC_GENES <- c("TCP1", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8")
EXON_DIR   <- file.path(ALSPROJ_BASE_DIR)        # path to ALS.chrN.anno_report_exon.txt files
RARE_MAF   <- 0.01

damaging_classes <- c("nonsynonymous_SNV", "stopgain", "stoploss",
                      "frameshift_deletion", "frameshift_insertion")
lof_classes      <- c("stopgain", "stoploss",
                      "frameshift_deletion", "frameshift_insertion")

# TRiC genes span 7 chromosomes -- only read the relevant exon files.
tric_chr_map <- c(TCP1 = "chr6",  CCT2 = "chr12", CCT3 = "chr1",  CCT4 = "chr2",
                  CCT5 = "chr5",  CCT6A = "chr7", CCT7 = "chr2",  CCT8 = "chr21")
tric_chrs <- unique(tric_chr_map)

# ---------------- Step 1: parse rare damaging TRiC variants -------------
# Mirrors parse_gene_burden_chr() but pre-filters to TRiC + rare MAF
# at the read step. Variant class column is left in for downstream
# LoF subsetting.
parse_tric_chr <- function(chr,
                           genes      = TRIC_GENES,
                           classes    = damaging_classes,
                           maf_thresh = RARE_MAF) {
  f <- file.path(EXON_DIR, paste0("ALS.", chr, ".anno_report_exon.txt"))
  if (!file.exists(f)) { message("Missing: ", f); return(NULL) }
  
  dt <- fread(f, sep = "\t", header = TRUE, na.strings = c(".", "NA", ""))
  
  # gene + class filter
  if (!"Gene" %in% names(dt))
    stop("Expected column 'Gene' in ", f)
  dt <- dt[Gene %in% genes]
  if (nrow(dt) == 0) return(NULL)
  
  if (!"ExonicFunc" %in% names(dt))
    stop("Expected column 'ExonicFunc' in ", f)
  dt <- dt[ExonicFunc %in% classes]
  if (nrow(dt) == 0) return(NULL)
  
  # rare filter -- "Max Freq of 3" is the max across 1KG/ESP/ExAC panels.
  # Stricter than any single panel: a variant only passes if it's rare
  # in ALL three references.
  maf_candidates <- c("Max Freq of 3", "gnomAD_genome_ALL", "gnomAD_exome_ALL",
                      "AF_popmax", "ExAC_ALL", "AF")
  maf_col <- intersect(maf_candidates, names(dt))[1]
  if (is.na(maf_col))
    stop("No recognised MAF column in ", f,
         "\nLooked for: ", paste(maf_candidates, collapse = ", "),
         "\nFound cols: ", paste(names(dt), collapse = ", "))
  
  af <- suppressWarnings(as.numeric(dt[[maf_col]]))
  af[is.na(af)] <- 0   # absent in reference panel -> treat as rare
  dt <- dt[af < maf_thresh]
  dt
}

tric_variants <- rbindlist(
  lapply(tric_chrs, parse_tric_chr),
  use.names = TRUE, fill = TRUE
)
cat("Rare damaging TRiC variants:", nrow(tric_variants), "\n")
cat("By gene:\n");  print(table(tric_variants$Gene))
cat("By class:\n"); print(table(tric_variants$ExonicFunc))

# ---------------- Step 2: per-sample x gene burden matrix -------------
# Sample columns in the ANNOVAR exon report are the genotype columns
# (everything that isn't an annotation column). use_cols already holds
# the 360 sample columns matched to metadata in your pipeline.
stopifnot(exists("use_cols"))

# Standard VCF genotype parsing.
#   0/0, 0|0, "0"           -> non-carrier
#   0/1, 1/0, 1/1, 1|2 etc. -> carrier
#   ./., 0/., ./1, ".", ""  -> missing (treated as non-carrier, conservative)
# Tighter than a simple exclusion list: anything containing "." is missing,
# so partial-missing genotypes can't silently inflate burden counts.
is_carrier <- function(x) {
  x <- as.character(x)
  hom_ref  <- x %in% c("0", "0/0", "0|0")
  has_miss <- is.na(x) | x == "" | grepl("\\.", x)
  has_alt  <- grepl("[1-9]", x) & !has_miss
  as.integer(has_alt & !hom_ref)
}

build_burden_matrix <- function(variants, samples,
                                classes = damaging_classes) {
  v <- variants[ExonicFunc %in% classes]
  if (nrow(v) == 0) {
    return(matrix(0L, nrow = length(TRIC_GENES), ncol = length(samples),
                  dimnames = list(TRIC_GENES, samples)))
  }
  # variant x sample carrier matrix (0/1)
  carr <- vapply(samples,
                 function(s) as.integer(is_carrier(v[[s]])),
                 integer(nrow(v)))
  rownames(carr) <- v$Gene
  # collapse to gene-level: count of rare damaging variants per gene per sample
  gb <- rowsum(carr, group = rownames(carr))
  # ensure every TRiC gene present (zero rows for genes with no variants)
  miss <- setdiff(TRIC_GENES, rownames(gb))
  if (length(miss)) {
    pad <- matrix(0L, nrow = length(miss), ncol = ncol(gb),
                  dimnames = list(miss, colnames(gb)))
    gb <- rbind(gb, pad)
  }
  gb[TRIC_GENES, , drop = FALSE]
}

tric_burden_dmg <- build_burden_matrix(tric_variants, use_cols, damaging_classes)
tric_burden_lof <- build_burden_matrix(tric_variants, use_cols, lof_classes)

cat("\nDamaging burden -- per-gene carrier counts (n samples =",
    ncol(tric_burden_dmg), "):\n")
print(rowSums(tric_burden_dmg > 0))
cat("LoF-only burden -- per-gene carrier counts:\n")
print(rowSums(tric_burden_lof > 0))

# ---------------- Step 3: attach subtype + clinical metadata --------
stopifnot(exists("col_to_id"), exists("rare_burden_sub"))

tric_df <- data.frame(
  sample_col      = colnames(tric_burden_dmg),
  cgnd_id         = unname(col_to_id[colnames(tric_burden_dmg)]),
  CCT2_burden     = tric_burden_dmg["CCT2", ],
  CCT2_burden_lof = tric_burden_lof["CCT2", ],
  TRiC_burden     = colSums(tric_burden_dmg),
  TRiC_burden_lof = colSums(tric_burden_lof),
  stringsAsFactors = FALSE
)
tric_df$CCT2_carrier <- as.integer(tric_df$CCT2_burden > 0)
tric_df$TRiC_carrier <- as.integer(tric_df$TRiC_burden > 0)

# Find a join key that exists in rare_burden_sub
key_candidates <- c("cgnd_id", "CGND_HDA_ID", "NYGC_CGND_ID",
                    "Participant_ID", "sample_id")
key_col <- intersect(key_candidates, names(rare_burden_sub))[1]
if (is.na(key_col))
  stop("No join key found in rare_burden_sub. Columns: ",
       paste(names(rare_burden_sub), collapse = ", "))
cat("\nJoining tric_df to rare_burden_sub on:", key_col, "\n")

tric_df <- tric_df %>%
  rename(!!key_col := cgnd_id) %>%
  inner_join(
    rare_burden_sub %>%
      dplyr::select(any_of(c(key_col, "Participant_ID",
                      "cluster_clean",
                      "AGE_AT_SYMPTOM_ONSET",
                      "SEX"))),
    by = key_col
  )

# Resolve cluster + age + sex column aliases
if (!"cluster_clean" %in% names(tric_df) && "cluster_final" %in% names(tric_df))
  tric_df$cluster_clean <- tric_df$cluster_final
if (!"AGE" %in% names(tric_df)) {
  age_col <- intersect(c("AGE_AT_SYMPTOM_ONSET", "Age_at_Onset"),
                       names(tric_df))[1]
  if (!is.na(age_col)) tric_df$AGE <- tric_df[[age_col]]
}
if (!"SEX" %in% names(tric_df)) {
  sex_col <- intersect(c("Sex"), names(tric_df))[1]
  if (!is.na(sex_col)) tric_df$SEX <- tric_df[[sex_col]]
}

# Sporadic only with subtype assignment
tric_df <- tric_df %>% filter(!is.na(cluster_clean))
tric_df$cluster_clean <- factor(tric_df$cluster_clean,
                                levels = c("S-ALS-1", "S-ALS-2", "S-ALS-3"))

cat("\nFinal analysis frame: n =", nrow(tric_df), "\n")
cat("Subtype distribution:\n"); print(table(tric_df$cluster_clean))
cat("CCT2 carriers x subtype:\n")
print(table(tric_df$cluster_clean, tric_df$CCT2_carrier,
            dnn = c("subtype","CCT2_carrier")))
cat("TRiC carriers x subtype:\n")
print(table(tric_df$cluster_clean, tric_df$TRiC_carrier,
            dnn = c("subtype","TRiC_carrier")))

# Clean subset (S-ALS-2 + S-ALS-3 -- the differentiation-clean half)
tric_clean <- tric_df %>%
  filter(cluster_clean %in% c("S-ALS-2", "S-ALS-3")) %>%
  mutate(cluster_clean = droplevels(cluster_clean),
         is_S3 = as.integer(cluster_clean == "S-ALS-3"))

# ============================================================
# TEST 1 -- CCT2 single-gene rare damaging burden by subtype
# ============================================================
cat("\n=== TEST 1: CCT2 single-gene burden ===\n")

t1_full <- table(tric_df$cluster_clean, tric_df$CCT2_carrier,
                 dnn = c("subtype","CCT2_carrier"))
cat("Full cohort 3x2 carrier table:\n"); print(t1_full)
print(chisq.test(t1_full, simulate.p.value = TRUE, B = 10000))
cat("Carrier rates by subtype:\n"); print(prop.table(t1_full, margin = 1))

t1_clean <- table(tric_clean$cluster_clean, tric_clean$CCT2_carrier,
                  dnn = c("subtype","CCT2_carrier"))
cat("\nClean subset (S-ALS-2 vs S-ALS-3) Fisher:\n"); print(t1_clean)
print(fisher.test(t1_clean))

# Poisson on count -- accounts for patients carrying multiple rare variants
cat("\nPoisson regression on rare damaging variant count (clean subset):\n")
print(summary(glm(CCT2_burden ~ cluster_clean,
                  data = tric_clean, family = poisson)))

# ============================================================
# TEST 2 -- TRiC complex pooled burden by subtype
# ============================================================
cat("\n=== TEST 2: TRiC complex pooled burden ===\n")

t2_full <- table(tric_df$cluster_clean, tric_df$TRiC_carrier,
                 dnn = c("subtype","TRiC_carrier"))
cat("Full cohort TRiC carrier table:\n"); print(t2_full)
print(chisq.test(t2_full, simulate.p.value = TRUE, B = 10000))
cat("Carrier rates:\n"); print(prop.table(t2_full, margin = 1))

cat("\nMean TRiC burden by subtype (rare damaging count):\n")
print(tric_df %>% group_by(cluster_clean) %>%
        summarise(mean = mean(TRiC_burden), sd = sd(TRiC_burden), n = n()))
print(summary(aov(TRiC_burden ~ cluster_clean, data = tric_df)))

# Negative binomial -- variant counts often overdispersed
if (requireNamespace("MASS", quietly = TRUE)) {
  cat("\nNegative binomial regression on TRiC burden:\n")
  print(summary(MASS::glm.nb(TRiC_burden ~ cluster_clean, data = tric_df)))
}

# ============================================================
# TEST 3 -- Per-gene Fisher tests (clean subset)
# Identifies which TRiC subunit(s) drive any pooled signal
# ============================================================
cat("\n=== TEST 3: Per-gene Fisher in clean subset ===\n")

per_gene_fisher <- function(gene, df = tric_clean) {
  carrier <- as.integer(tric_burden_dmg[gene, df$sample_col] > 0)
  tab <- table(df$cluster_clean, carrier, dnn = c("subtype", "carrier"))
  if (any(dim(tab) < 2) || !"1" %in% colnames(tab) ||
      sum(tab[, "1"]) == 0) {
    return(tibble(gene = gene, n_carriers_S2 = 0L, n_carriers_S3 = 0L,
                  OR = NA_real_, p = NA_real_))
  }
  ft <- fisher.test(tab)
  tibble(gene          = gene,
         n_carriers_S2 = if ("S-ALS-2" %in% rownames(tab)) tab["S-ALS-2","1"] else 0L,
         n_carriers_S3 = if ("S-ALS-3" %in% rownames(tab)) tab["S-ALS-3","1"] else 0L,
         OR            = unname(ft$estimate),
         p             = ft$p.value)
}

per_gene <- bind_rows(lapply(TRIC_GENES, per_gene_fisher)) %>%
  mutate(p_fdr = p.adjust(p, "BH")) %>%
  arrange(p)
cat("Per-gene Fisher (S-ALS-3 vs S-ALS-2, clean subset):\n")
print(per_gene)

# ============================================================
# TEST 4 -- Logistic: TRiC burden -> P(S-ALS-3), adjusted
# ============================================================
cat("\n=== TEST 4: TRiC burden -> P(S-ALS-3), adjusted age + sex ===\n")

tric_clean$SEX <- factor(tric_clean$SEX)

t4_unadj <- glm(is_S3 ~ TRiC_burden, data = tric_clean, family = binomial)
t4_adj   <- glm(is_S3 ~ TRiC_burden + AGE + SEX,
                data = tric_clean, family = binomial)

cat("Unadjusted:\n"); print(summary(t4_unadj))
cat("\nAge + sex adjusted:\n"); print(summary(t4_adj))

t4_lof <- glm(is_S3 ~ TRiC_burden_lof + AGE + SEX,
              data = tric_clean, family = binomial)
cat("\nLoF-only TRiC burden, adjusted:\n"); print(summary(t4_lof))

# ============================================================
# TEST 5 -- CCT2 burden -> CCT2 expression (LoF mechanism)
# Direction of beta tells you mechanism:
#   beta > 0  compensatory upregulation (TRiC stress response)
#   beta < 0  haploinsufficiency / direct LoF reduction
#   beta ~ 0  no transcriptional consequence (typical of missense)
# ============================================================
cat("\n=== TEST 5: CCT2 burden -> CCT2 expression ===\n")

stopifnot(exists("cct2_df"))

# Dedup CCT2 expression to one row per patient (longitudinal hazard)
expr_col_candidates <- c("CCT2", "CCT2_expr", "expression", "expr")
expr_col <- intersect(expr_col_candidates, names(cct2_df))[1]
if (is.na(expr_col))
  stop("Could not find CCT2 expression column in cct2_df. Looked for: ",
       paste(expr_col_candidates, collapse = ", "))

cct2_unique <- cct2_df %>%
  distinct(Participant_ID, .keep_all = TRUE) %>%
  transmute(Participant_ID, CCT2_expr = .data[[expr_col]])

if (!"Participant_ID" %in% names(tric_clean))
  stop("tric_clean is missing Participant_ID -- check the join in Step 3 ",
       "(rare_burden_sub may not have carried Participant_ID through)")

mech_df <- tric_clean %>%
  inner_join(cct2_unique, by = "Participant_ID") %>%
  filter(!is.na(CCT2_expr))

cat("Patients with both CCT2 burden and expression:", nrow(mech_df), "\n")
cat("Carrier counts in mechanism subset:\n")
print(table(CCT2_carrier = mech_df$CCT2_carrier))

m1 <- lm(CCT2_expr ~ CCT2_burden + AGE + SEX + cluster_clean, data = mech_df)
cat("\nDamaging variant burden -> CCT2 expression:\n"); print(summary(m1))

m2 <- lm(CCT2_expr ~ CCT2_burden_lof + AGE + SEX + cluster_clean, data = mech_df)
cat("\nLoF-only burden -> CCT2 expression:\n"); print(summary(m2))

cat("\nMean CCT2 expression by carrier status (descriptive):\n")
print(mech_df %>% group_by(CCT2_carrier) %>%
        summarise(n = n(),
                  mean_expr = mean(CCT2_expr),
                  sd        = sd(CCT2_expr)))

# ============================================================
# Save outputs
# ============================================================
out_dir <- "tric_burden_results"
dir.create(out_dir, showWarnings = FALSE)
fwrite(tric_df,    file.path(out_dir, "tric_burden_full.csv"))
fwrite(tric_clean, file.path(out_dir, "tric_burden_clean.csv"))
fwrite(per_gene,   file.path(out_dir, "tric_per_gene_fisher.csv"))
saveRDS(list(t4_unadj = t4_unadj, t4_adj = t4_adj, t4_lof = t4_lof,
             m1_dmg = m1, m2_lof = m2),
        file.path(out_dir, "tric_models.rds"))

cat("\nResults written to:", out_dir, "\n")

