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

eqtldir <- "~/Documents/DSN/ALS/cct2_eqtl"
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
gwas <- read.table("~/Documents/DSN/ALS/van_rheenen2021_als.tsv", header = TRUE, stringsAsFactors = FALSE)
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












# ============================================================
# 12_cct2_eqtl.R   (AUDITED 2026-08-24)
# ============================================================
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

eqtldir <- "~/Documents/DSN/ALS/cct2_eqtl"
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
## AUDIT REBUILD: the circulated script read expr_clean_mat and
## sporadic_final_check from a session that no longer exists. Both are
## rebuilt here from the pinned layer.
##   - CCT2 expression comes from patient_level_FINAL.rds (one VST value
##     per participant, the same dedup the rest of the pipeline uses)
##   - covariates + CGND crosswalk come from backbone_FINAL.rds meta
BACKBONE <- "backbone_FINAL.rds"
PATIENTS <- "patient_level_FINAL.rds"
for (f in c(BACKBONE, PATIENTS))
  if (!file.exists(f)) stop("not found in working directory: ", f)

b    <- readRDS(BACKBONE)
pl   <- readRDS(PATIENTS)
meta <- as.data.frame(b$meta)
pat  <- as.data.frame(pl$pat)

cct2_df <- pat %>%
  distinct(Participant_ID, .keep_all = TRUE) %>%
  transmute(Participant_ID, CCT2_expr = CCT2)
cat("CCT2 expression: one value per participant, n =", nrow(cct2_df), "\n")

## PROVENANCE NOTE: the original built cct2_df from the sample-level VST
## matrix WITHOUT deduplicating participants before the inner_join below.
## Any participant with >1 RNA sample would have entered the regressions
## twice with identical genotypes. The gate at the bottom against the saved
## cct2_cis_eqtl_results.csv adjudicates whether that happened: if the
## pinned n exceeds the participant-level n printed here, it did.

cg <- grep("CGND", colnames(meta), value = TRUE, ignore.case = TRUE)[1]
if (is.na(cg)) stop("no CGND column in backbone meta")

crosswalk_v2 <- meta %>%
  filter(!is.na(.data[[cg]]), .data[[cg]] != "",
         !is.na(Participant_ID)) %>%
  transmute(VCF_ID = sub("-b[0-9]+$", "", as.character(.data[[cg]])),
            Participant_ID, SEX,
            SV1, SV2, SV3, SV4, SV5, SV6) %>%
  distinct(VCF_ID, .keep_all = TRUE)
cat("Crosswalk:", nrow(crosswalk_v2), "participants with CGND IDs\n")

dosage$IID <- sub("-b38$", "", dosage$IID)
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
## AUDIT NOTE: the circulated script used SEX.x = the SEX column of the
## PLINK .raw file. Here sex comes from the pinned metadata. If the .fam
## sex was populated from the same metadata these are identical; the gate
## below will say so. Suffix-proof lookup:
sex_col <- intersect(c("SEX", "SEX.y", "SEX.x"), colnames(geno_pheno))[1]
cat("Sex covariate column used:", sex_col, "\n")
covars <- geno_pheno[, c(sex_col, paste0("SV", 1:6)), with = FALSE]
covars[[sex_col]] <- as.numeric(factor(covars[[sex_col]]))

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

# Save full results -- non-destructive gate against the pinned csv
RES_FINAL <- file.path(eqtldir, "cct2_cis_eqtl_results.csv")
if (!file.exists(RES_FINAL)) {
  fwrite(results, RES_FINAL)
  cat("Saved", RES_FINAL, "\n")
} else {
  fwrite(results, file.path(eqtldir, "cct2_cis_eqtl_results_AUDITCHECK.csv"))
  old <- fread(RES_FINAL)
  cat("\nPinned results csv exists -- wrote AUDITCHECK instead. GATE:\n")
  cat("  n tested (pinned vs rerun): ", nrow(old), "vs", nrow(results), "\n")
  cat("  max n per SNP:              ", max(old$n, na.rm=TRUE), "vs",
      max(results$n, na.rm=TRUE), "\n")
  cat("  top SNP:                    ", old$SNP[which.min(old$p)], "vs",
      results$SNP[1], "\n")
  cat("  min FDR-q (manifest ~0.71): ", round(min(old$fdr_q, na.rm=TRUE), 3),
      "vs", round(min(results$fdr_q, na.rm=TRUE), 3), "\n")
  shared <- intersect(old$SNP, results$SNP)
  rho <- cor(old$beta[match(shared, old$SNP)],
             results$beta[match(shared, results$SNP)],
             use = "complete.obs", method = "spearman")
  cat("  beta rank correlation on", length(shared), "shared SNPs:",
      round(rho, 4), "\n")
}

# ----------------------------------------------------------------
# Step 5: Lookup top variants in van Rheenen 2021 GWAS
# Question: do CCT2 cis-eQTLs also affect ALS risk?
# This is the mechanistic anchor — a SNP that both regulates
# CCT2 AND affects ALS risk is a candidate causal variant
# ----------------------------------------------------------------
cat("\n=== Looking up top eQTL hits in van Rheenen ALS GWAS ===\n")

## AUDIT REPLACEMENT (2026-08-26): the previous file
## (van_rheenen2021_als.tsv) was the GWAS Catalog's Build37 release, while
## the AnswerALS genotypes are GRCh38 -- a 393,780 bp offset on chr12 that
## reduced the join to chance level. Use the HARMONISED release instead:
## GCST90027164 .h.tsv.gz, whose hm_* columns are GRCh38. Coordinates were
## already GRCh38 in the submitted block, but 45% of variants had their
## effect allele swapped during harmonisation (hm_beta negated), so the
## hm_* columns must be used as a set -- never mixed with the submitted ones.
GWAS_BASE <- "vanrheenen_GRCh38_harmonised.base"   # SNP CHR BP A1 A2 BETA SE P
if (!file.exists(GWAS_BASE)) stop("harmonised base not found: ", GWAS_BASE)

gwas <- data.table::fread(GWAS_BASE)
cat("GWAS sumstats columns:", colnames(gwas), "\n")
cat("GWAS sumstats rows:", nrow(gwas), "\n")

## Harmonisation could not map every variant (hm_code 15); those rows carry
## empty CHR/BP. `gwas$CHR == 12` yields NA for them and NA row indices
## return all-NA rows rather than dropping them, which silently inflated
## chr12 from 503,184 to 518,098. Drop them explicitly.
.n_bad <- sum(is.na(gwas$CHR) | is.na(gwas$BP))
gwas   <- gwas[!is.na(gwas$CHR) & !is.na(gwas$BP), ]
cat("Dropped unmapped (hm_code 15) variants:", .n_bad, "\n")

# Top 50 cis-eQTLs by p-value
top_eqtls <- head(results, 50)

# Strip _A1 suffix from PLINK SNP names -> chr:pos format
top_eqtls$snp_id <- sub("_[ACGT]+$", "", top_eqtls$SNP)

# Build snp_id on GWAS side from chr:pos (chr12 only, much faster)
gwas_chr12 <- gwas[gwas$CHR == 12, ]
gwas_chr12$snp_id <- paste0(gwas_chr12$CHR, ":", gwas_chr12$BP)

cat("GWAS chr12 variants:", nrow(gwas_chr12), "\n")     # expect 503,184
cat("Top eQTL SNP IDs (first 5):", head(top_eqtls$snp_id, 5), "\n")
cat("GWAS SNP IDs (first 5):", head(gwas_chr12$snp_id, 5), "\n")
stopifnot(!any(is.na(gwas_chr12$CHR)))

# Join eQTL results to GWAS by chr:pos
## Column names follow the harmonised base file. Note gwas_eaf is absent:
## allele frequency was not carried through when the base was built.
eqtl_gwas <- top_eqtls %>%
  left_join(gwas_chr12 %>%
              dplyr::select(snp_id, rsid = SNP,
                            gwas_a1   = A1,
                            gwas_a2   = A2,
                            gwas_beta = BETA,
                            gwas_se   = SE,
                            gwas_p    = P),
            by = "snp_id")

cat("\nTop eQTLs with GWAS lookup:\n")
print(eqtl_gwas %>%
        dplyr::select(snp_id, rsid, beta, se, p, fdr_q, gwas_beta, gwas_p) %>%
        head(20))

cat("\nN top eQTLs with GWAS data:",
    sum(!is.na(eqtl_gwas$gwas_p)), "of 50\n")

## ---------------------------------------------------------------------------
## ALLELE ORIENTATION CHECK
## The eQTL A1 is the PLINK counted allele (the letter after the underscore
## in the SNP name); the GWAS A1 is the harmonised effect allele. Betas are
## only comparable in direction when these agree. Flag rows where they do
## not, so no direction-of-effect claim is made on mismatched alleles.
## ---------------------------------------------------------------------------
eqtl_gwas$eqtl_a1 <- sub(".*_([ACGT]+)$", "\\1", eqtl_gwas$SNP)
eqtl_gwas$allele_match <- with(eqtl_gwas,
                               ifelse(is.na(gwas_a1), NA,
                                      ifelse(toupper(eqtl_a1) == toupper(gwas_a1), "same",
                                             ifelse(toupper(eqtl_a1) == toupper(gwas_a2), "flipped", "mismatch"))))
cat("\nAllele orientation of matched SNPs:\n")
print(table(eqtl_gwas$allele_match, useNA = "ifany"))
cat("For 'flipped' rows the GWAS beta sign must be reversed before",
    "comparing direction with the eQTL beta.\n")
eqtl_gwas$gwas_beta_aligned <- with(eqtl_gwas,
                                    ifelse(allele_match == "flipped", -gwas_beta,
                                           ifelse(allele_match == "same", gwas_beta, NA)))

## ===========================================================================
## AUDIT BLOCK -- GENOME BUILD MISMATCH (DIAGNOSED 2026-08-24, FIXED 2026-08-26)
##
## RESOLVED: Step 5 now reads the GWAS Catalog HARMONISED release
## (GCST90027164 .h.tsv.gz -> vanrheenen_GRCh38_harmonised.base.gz). The
## lookup matches 42 of the top 50 cis-eQTLs and returns a clean result:
## none is nominally ALS-associated (minimum GWAS p = 0.116). The block
## below is retained as a permanent regression check.
##
## ORIGINAL DIAGNOSIS (kept for the record):
##
## FINDING: the GWAS lookup matches essentially nothing, and this is a
## COORDINATE-SYSTEM FAILURE, not a null result. Evidence from the rerun:
##
##   - 3,782 eQTL variants tested; 4,767 GWAS variants lie inside the same
##     1.03 Mb window; only 22 positions coincide. Chance alone predicts
##     ~18, so the overlap is at background level.
##   - NOT a frequency artifact: matched and unmatched eQTL variants have
##     near-identical MAF distributions (median 0.145 vs 0.144), and 2,771
##     tested variants exceed MAF 0.05.
##   - The pairwise position-difference histogram shows one sharp spike at
##     +393,780 bp, occurring 356 times against a background of 11-12.
##     That is the GRCh37-minus-GRCh38 offset for this region of chr12.
##
##   WARNING to anyone rerunning: a NEAREST-NEIGHBOUR offset check is
##   useless here and was initially misread during this audit. GWAS variant
##   density in the window is ~1 per 216 bp, so the nearest neighbour to any
##   position sits ~100 bp away by chance whatever the build. Use the
##   difference histogram below, not nearest neighbours.
##
## CONSEQUENCE: no statement of the form "no CCT2 cis-eQTL is ALS-associated"
## is supported. Redo the lookup joining on rsID (preferred) or after a
## proper liftOver. Do NOT just add 393,780 to every position: that offset
## is a property of this region between two liftover chains, not a global
## constant, and it will be wrong across any chain break inside the window.
## Use the measured offset only to sanity check a correct liftOver.
##
## KNOCK-ON: script 13 (PRSet) uses this same GWAS file against the same
## GRCh38 genotypes. If PRSice matched on position, the pathway PRS were
## built from almost no variants and the specificity matrix is
## uninterpretable. Check the PRSice log BEFORE trusting any PRS result.
## ===========================================================================
.all_pos  <- unique(as.integer(sub(".*:([0-9]+)_.*", "\\1", results$SNP)))
.gwas_pos <- unique(gwas_chr12$BP)   # AUDIT: harmonised base uses BP
.win      <- range(.all_pos)
.in_win   <- sum(.gwas_pos >= .win[1] & .gwas_pos <= .win[2])
.overlap  <- length(intersect(.all_pos, .gwas_pos))
.expected <- length(.all_pos) * .in_win / diff(.win)

cat("\n--- build/join diagnostic ---\n")
cat("eQTL variants tested:   ", length(.all_pos), "\n")
cat("GWAS variants in window:", .in_win, "\n")
cat("Positional matches:     ", .overlap,
    sprintf(" (chance expectation ~%.0f)\n", .expected))

## Name the offset rather than inferring it from nearest neighbours.
.gw  <- .gwas_pos[.gwas_pos > .win[1] - 8e5 & .gwas_pos < .win[2] + 8e5]
set.seed(1)
.sub <- sample(.all_pos, min(400, length(.all_pos)))
.d   <- as.vector(outer(.gw, .sub, "-"))
.d   <- .d[abs(.d) < 8e5]
.tt  <- sort(table(.d), decreasing = TRUE)
cat("Top position-difference modes (a spike = systematic build offset):\n")
print(head(.tt, 5))

if (length(.gwas_pos) == 0)
  stop("diagnostic has no GWAS positions -- check the coordinate column name")
if (.overlap <= 3 * .expected) {
  cat("\n*** JOIN FAILED: positional overlap is at chance level.\n")
  if (length(.tt) > 0 &&
      as.integer(.tt[1]) > 20 * stats::median(as.integer(.tt)))
    cat("*** Systematic offset:", names(.tt)[1], "bp (",
        as.integer(.tt[1]), "occurrences ) -- GENOME BUILD MISMATCH.\n")
  cat("*** Do NOT report the GWAS lookup as a null result.\n")
  cat("*** Remedy: join on rsID, or liftOver the eQTL positions properly.\n")
  cat("*** Also inspect the PRSice log before trusting script 13.\n")
} else {
  cat("Positional overlap exceeds chance; lookup can be interpreted.\n")
}
# Top eQTLs nominally ALS-associated
hits <- eqtl_gwas %>% filter(!is.na(gwas_p) & gwas_p < 0.05)
cat("\nN nominally ALS-associated (GWAS p < 0.05):", nrow(hits),
    "of", sum(!is.na(eqtl_gwas$gwas_p)), "with GWAS data\n")
cat("Note: 50 SNPs tested at alpha 0.05 -> ~2.5 expected by chance.",
    "Do not over-read a small count.\n")
if (nrow(hits) > 0) {
  cat("\nNominally ALS-associated CCT2 cis-eQTLs:\n")
  print(hits %>% dplyr::select(snp_id, rsid, beta, p, fdr_q,
                               gwas_beta, gwas_beta_aligned, gwas_p,
                               eqtl_a1, gwas_a1, gwas_a2, allele_match))
  cat("\nDirection concordance (eQTL beta vs allele-aligned GWAS beta):\n")
  print(table(sign(hits$beta) == sign(hits$gwas_beta_aligned), useNA = "ifany"))
} else {
  cat("None. With the join now working, this IS a reportable null:",
      "\nno top CCT2 cis-eQTL is nominally associated with ALS risk.\n")
}
cat("\nMinimum GWAS p among matched top eQTLs:",
    if (any(!is.na(eqtl_gwas$gwas_p)))
      signif(min(eqtl_gwas$gwas_p, na.rm = TRUE), 3) else NA, "\n")

LOOKUP_FINAL <- file.path(eqtldir, "cct2_eqtl_gwas_lookup.csv")
if (!file.exists(LOOKUP_FINAL)) {
  fwrite(eqtl_gwas, LOOKUP_FINAL)
} else {
  fwrite(eqtl_gwas, file.path(eqtldir, "cct2_eqtl_gwas_lookup_AUDITCHECK.csv"))
  cat("Pinned GWAS lookup exists -- wrote AUDITCHECK instead.\n")
}

# ----------------------------------------------------------------
# Step 6: Visualization — Manhattan plot of cis-window
# ----------------------------------------------------------------
library(ggplot2)

# Extract position from SNP ID (chr:pos_A1 → pos)
results$pos <- as.integer(sub(".*:([0-9]+)_.*", "\\1", results$SNP))

# CCT2 gene boundaries
## AUDIT: confirm which build these coordinates are on before using them to
## annotate the plot -- they must match the GENOTYPE build (GRCh38), not the
## GWAS build. Mislabelling shifts the shaded gene body by ~400 kb.
gene_start <- 69716872
gene_end   <- 69750479
if (gene_start < min(.all_pos) || gene_end > max(.all_pos))
  warning("CCT2 gene boundaries fall outside the tested cis-window -- ",
          "check the genome build of these coordinates")

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