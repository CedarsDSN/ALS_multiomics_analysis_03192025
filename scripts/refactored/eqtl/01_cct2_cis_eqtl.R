# ============================================================
# Refactored entrypoint for original script: ciseQTLTesting.R
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
if (isTRUE(ALS_ARGS$help)) { als_print_help("01_cct2_cis_eqtl"); quit(save = "no", status = 0) }
CFG <- als_load_config(ALS_ARGS$config)
ALSPROJ_BASE_DIR <- als_choose(ALS_ARGS$base_dir, CFG, "project.base_dir", getwd())
ALSPROJ_OUTDIR <- als_choose(ALS_ARGS$outdir, CFG, "outputs.root", file.path(getwd(), "results", "01_cct2_cis_eqtl"))
ALSPROJ_SCRIPT_ID <- "01_cct2_cis_eqtl"
als_init(ALSPROJ_SCRIPT_ID, ALSPROJ_BASE_DIR, ALSPROJ_OUTDIR, ALS_ARGS$redirect_outputs)
als_register_output_wrappers()
als_load_objects(ALS_ARGS$objects)
on.exit(als_write_session_info(ALSPROJ_SCRIPT_ID), add = TRUE)
# ---- End bootstrap ----

# ============================================================
# CCT2 cis-eQTL analysis
# Tests whether common variants in chr12:CCT2 ±500kb predict
# CCT2 expression. Then looks up top variants in van Rheenen
# ALS GWAS to ask: do CCT2 cis-eQTLs also affect ALS risk?
#
# Pipeline:
#   1. Load dosage matrix + CCT2 expression + covariates
#   2. Per-SNP linear regression with sex + 6 SVs
#   3. Bonferroni + FDR correction
#   4. Top-hit lookup in van Rheenen sumstats
#   5. Visualization (Manhattan + locus zoom)
# ============================================================

library(dplyr)
library(data.table)

eqtldir <- file.path(ALSPROJ_BASE_DIR, "cct2_eqtl")
gwasdir <- "/common/compbiomed-dsn/projects/ALS_multiomics_analysis_03192025/large_data/genomics/genomics/4_JointGenotyping/prset"

# ----------------------------------------------------------------
# Step 1: Load PLINK dosage matrix (output of --recode A)
# Format: FID IID PAT MAT SEX PHENOTYPE then one column per SNP
# Each SNP column is named SNPID_A1; values are 0/1/2 dosages of A1
# ----------------------------------------------------------------
dosage <- fread(file.path(eqtldir, "cct2_cis_dosage.raw"))
cat("Dosage matrix loaded:", nrow(dosage), "samples ×",
    ncol(dosage) - 6, "SNPs\n")

# Identify SNP columns (everything past col 6)
snp_cols <- colnames(dosage)[7:ncol(dosage)]
cat("First few SNPs:", head(snp_cols, 3), "...\n")

# ----------------------------------------------------------------
# Step 2: Build phenotype + covariate frame
# CCT2 expression keyed on Participant_ID; need to map IID
# (CGND-HDA-XXXXX-b38) → Participant_ID via sporadic_final_check
# ----------------------------------------------------------------
cct2_vec <- expr_clean_mat["CCT2", ]
pid_from_sample <- sub("-[0-9]+-T$", "", names(cct2_vec))
cct2_df <- data.frame(Participant_ID = pid_from_sample,
                      CCT2_expr = as.numeric(cct2_vec))

# ID crosswalk: VCF_ID (CGND-HDA-XXXXX-b38) ↔ Participant_ID
crosswalk_v2 <- sporadic_final_check %>%
filter(!is.na(NYGC_CGND_ID), !is.na(Participant_ID)) %>%
mutate(VCF_ID = paste0(NYGC_CGND_ID, "-b38")) %>%
distinct(VCF_ID, .keep_all = TRUE) %>%   # dedupe on VCF_ID alone
select(VCF_ID, Participant_ID, SEX, SV1, SV2, SV3, SV4, SV5, SV6)

dosage$IID <- sub("-b38$", "", dosage$IID)
crosswalk_v2$VCF_ID <- sub("-b38$", "", crosswalk_v2$VCF_ID)
# Merge expression + covariates onto dosage rows by IID
geno_pheno <- dosage %>%
  rename(VCF_ID = IID) %>%
  inner_join(crosswalk_v2, by = "VCF_ID") %>%
  inner_join(cct2_df, by = "Participant_ID")

cat("Final analysis n:", nrow(geno_pheno), "\n")
cat("Missing CCT2:", sum(is.na(geno_pheno$CCT2_expr)), "\n")
cat("Missing SVs:", sum(is.na(geno_pheno$SV1)), "\n")

# Drop any incomplete rows
geno_pheno <- geno_pheno %>%
  filter(!is.na(CCT2_expr), !is.na(SV1))
cat("After complete-case filter:", nrow(geno_pheno), "\n")

# ----------------------------------------------------------------
# Step 3: Per-SNP linear regression
# Model: CCT2 ~ dosage + SEX + 6 SVs
# SVs included to match the same correction used in the
# transcriptomic pipeline (DESeq2 with 6 SVs)
# ----------------------------------------------------------------
cat("\n=== Running cis-eQTL regressions ===\n")

results <- data.frame(
  SNP    = snp_cols,
  beta   = NA_real_,
  se     = NA_real_,
  t      = NA_real_,
  p      = NA_real_,
  n      = NA_integer_,
  maf    = NA_real_,
  stringsAsFactors = FALSE
)

# Pre-extract covariate matrix once (avoids repeated lookups)
covars <- geno_pheno[, c("SEX.x", paste0("SV", 1:6)), with = FALSE]

pb <- txtProgressBar(min = 0, max = length(snp_cols), style = 3)
for (i in seq_along(snp_cols)) {
  snp <- snp_cols[i]
  d   <- geno_pheno[[snp]]
  
  # Skip if monomorphic in this subset or too much missingness
  if (length(unique(na.omit(d))) < 2) next
  if (mean(is.na(d)) > 0.05) next
  
  af <- mean(d, na.rm = TRUE) / 2
  maf <- min(af, 1 - af)
  if (maf < 0.01) next
  
  fit <- tryCatch(
    lm(geno_pheno$CCT2_expr ~ d + as.matrix(covars)),
    error = function(e) NULL)
  if (is.null(fit)) next
  
  s <- summary(fit)$coefficients
  if (!"d" %in% rownames(s)) next
  
  results[i, "beta"] <- s["d", "Estimate"]
  results[i, "se"]   <- s["d", "Std. Error"]
  results[i, "t"]    <- s["d", "t value"]
  results[i, "p"]    <- s["d", "Pr(>|t|)"]
  results[i, "n"]    <- sum(!is.na(d) & !is.na(geno_pheno$CCT2_expr))
  results[i, "maf"]  <- maf
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Keep only successful tests
results <- results %>% filter(!is.na(p))
cat("\nSuccessful regressions:", nrow(results), "\n")

# ----------------------------------------------------------------
# Step 4: Multiple testing correction
# Bonferroni against the number of independent SNPs in the cis-window
# (conservative: actual N tests). FDR-BH for discovery.
# ----------------------------------------------------------------
results <- results %>%
  mutate(
    bonf_p = pmin(1, p * nrow(results)),
    fdr_q  = p.adjust(p, method = "BH")
  ) %>%
  arrange(p)

cat("\n=== Top cis-eQTL hits ===\n")
print(head(results, 20))

cat("\nGenome-wide significant (p < 5e-8):",
    sum(results$p < 5e-8), "\n")
cat("Cis-significant (Bonferroni p < 0.05):",
    sum(results$bonf_p < 0.05), "\n")
cat("FDR < 0.05:", sum(results$fdr_q < 0.05), "\n")
cat("Nominal (p < 0.001):", sum(results$p < 1e-3), "\n")

# Save full results
fwrite(results, file.path(eqtldir, "cct2_cis_eqtl_results.csv"))

# ----------------------------------------------------------------
# Step 5: Lookup top variants in van Rheenen 2021 GWAS
# Question: do CCT2 cis-eQTLs also affect ALS risk?
# This is the mechanistic anchor — a SNP that both regulates
# CCT2 AND affects ALS risk is a candidate causal variant
# ----------------------------------------------------------------
cat("\n=== Looking up top eQTL hits in van Rheenen ALS GWAS ===\n")

# Top 50 cis-eQTLs by p-value
top_eqtls <- head(results, 50)

# Strip _A1 suffix from PLINK SNP names → chr:pos format
top_eqtls$snp_id <- sub("_[ACGT]+$", "", top_eqtls$SNP)

# Load van Rheenen sumstats (chr:pos format from prior PRSet prep)
gwas <- read.table(file.path(ALSPROJ_BASE_DIR, "van_rheenen2021_als.tsv"), header = TRUE, stringsAsFactors = FALSE)
cat("GWAS sumstats columns:", colnames(gwas), "\n")
cat("GWAS sumstats rows:", nrow(gwas), "\n")

# Top 50 cis-eQTLs by p-value
top_eqtls <- head(results, 50)

# Strip _A1 suffix from PLINK SNP names → chr:pos format
top_eqtls$snp_id <- sub("_[ACGT]+$", "", top_eqtls$SNP)

# Build snp_id on GWAS side from chr:pos (chr12 only, much faster)
gwas_chr12 <- gwas[gwas$chromosome == 12, ]
gwas_chr12$snp_id <- paste0(gwas_chr12$chromosome, ":", gwas_chr12$base_pair_location)

cat("GWAS chr12 variants:", nrow(gwas_chr12), "\n")
cat("Top eQTL SNP IDs (first 5):", head(top_eqtls$snp_id, 5), "\n")
cat("GWAS SNP IDs (first 5):", head(gwas_chr12$snp_id, 5), "\n")

# Join eQTL results to GWAS by chr:pos
eqtl_gwas <- top_eqtls %>%
  left_join(gwas_chr12 %>%
              select(snp_id, rsid,
                     gwas_a1     = effect_allele,
                     gwas_a2     = other_allele,
                     gwas_eaf    = effect_allele_frequency,
                     gwas_beta   = beta,
                     gwas_se     = standard_error,
                     gwas_p      = p_value),
            by = "snp_id")

cat("\nTop eQTLs with GWAS lookup:\n")
print(eqtl_gwas %>%
        select(snp_id, rsid, beta, se, p, fdr_q, gwas_beta, gwas_p) %>%
        head(20))

cat("\nN top eQTLs with GWAS data:",
    sum(!is.na(eqtl_gwas$gwas_p)), "of 50\n")

# Top eQTLs nominally ALS-associated
hits <- eqtl_gwas %>% filter(!is.na(gwas_p) & gwas_p < 0.05)
cat("N nominally ALS-associated (GWAS p < 0.05):", nrow(hits), "\n")
if (nrow(hits) > 0) {
  cat("\nNominally ALS-associated CCT2 cis-eQTLs:\n")
  print(hits %>% select(snp_id, rsid, beta, p, fdr_q,
                        gwas_beta, gwas_p, gwas_a1, gwas_a2))
}

fwrite(eqtl_gwas, file.path(eqtldir, "cct2_eqtl_gwas_lookup.csv"))

# ----------------------------------------------------------------
# Step 6: Visualization — Manhattan plot of cis-window
# ----------------------------------------------------------------
library(ggplot2)

# Extract position from SNP ID (chr:pos_A1 → pos)
results$pos <- as.integer(sub(".*:([0-9]+)_.*", "\\1", results$SNP))

# CCT2 gene boundaries
gene_start <- 69716872
gene_end   <- 69750479

p_manhattan <- ggplot(results, aes(x = pos / 1e6, y = -log10(p))) +
  geom_point(alpha = 0.6, size = 1.5,
             color = ifelse(results$fdr_q < 0.05, "firebrick", "gray40")) +
  geom_hline(yintercept = -log10(0.05 / nrow(results)),
             linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(5e-8),
             linetype = "dotted", color = "red") +
  annotate("rect", xmin = gene_start / 1e6, xmax = gene_end / 1e6,
           ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "darkgreen") +
  annotate("text", x = (gene_start + gene_end) / 2e6,
           y = max(-log10(results$p)) * 0.95,
           label = "CCT2", size = 4, fontface = "italic") +
  labs(x = "Position on chr12 (Mb)",
       y = "-log10(p)",
       title = "CCT2 cis-eQTL — chr12 ±500kb",
       subtitle = paste0("n = ", max(results$n), " sporadic patients; ",
                         nrow(results), " variants tested")) +
  theme_bw() + theme(panel.grid.minor = element_blank())

ggsave(file.path(eqtldir, "cct2_cis_eqtl_manhattan.pdf"),
       p_manhattan, width = 8, height = 5)

# ----------------------------------------------------------------
# Step 7: Headline numbers for the slide
# ----------------------------------------------------------------
cat("\n========================================\n")
cat("CCT2 cis-eQTL — summary for slide\n")
cat("========================================\n")
cat("N samples:           ", max(results$n), "\n")
cat("N variants tested:   ", nrow(results), "\n")
cat("Top SNP:             ", results$SNP[1], "\n")
cat("Top SNP β (per allele):", round(results$beta[1], 4), "\n")
cat("Top SNP p:           ", format.pval(results$p[1], digits = 3), "\n")
cat("Top SNP FDR-q:       ", format.pval(results$fdr_q[1], digits = 3), "\n")
cat("Bonferroni-significant SNPs:", sum(results$bonf_p < 0.05), "\n")
cat("FDR-significant SNPs:", sum(results$fdr_q < 0.05), "\n")
if (nrow(hits) > 0) {
  cat("Top eQTL also nominally ALS-associated:", hits$snp_id[1],
      "(eQTL p =", format.pval(hits$p[1], 2),
      ", GWAS p =", format.pval(hits$gwas_p[1], 2), ")\n")
}
