# ============================================================
# Script: genomic_analysis.R
# Project: ALS multi-omics analysis
# Source status: documented/sanitized for GitHub handoff
# ============================================================
# Documentation summary from scripts/README.md
# 
# ## `genomic_analysis.R`
# 
# **Category:** Genomics / rare variant burden
# 
# **Purpose**
# 
# Parses ANNOVAR-style exon annotation files, creates sample-level rare damaging/LoF burden metrics, compares burden across sporadic ALS subtypes, associates LoF burden with NfL, builds gene-level burden matrices, and explores ALS-panel burden/clustering.
# 
# **Inputs**
# 
# - Exon annotation files such as `ALS.chr2.anno_report_exon.txt` and corresponding files for all chromosomes.
# - Workspace object `sporadic_meta` with sample/participant/CGND IDs and subtype labels.
# - NfL and ALSFRS fields in merged subtype/clinical metadata for downstream associations.
# 
# **Outputs**
# 
# - `genomic_burden_per_sample_all_chrom.csv`.
# - `genomic_burden_by_subtype.pdf`.
# - `nfl_vs_lof_burden_scatter.pdf`.
# - Important in-memory objects currently not saved consistently: `rare_burden_all`, `rare_burden_sub`, `gene_burden_combined`, `gene_burden_filt`, `col_to_id`, `use_cols`, `crosswalk`.
# 
# **Upstream dependencies**
# 
# - Run after `sporadic_subtyping.R` has produced sporadic subtype metadata.
# - Requires WGS/ANNOVAR preprocessing outside R.
# 
# **How to run**
# 
# Legacy documented version:
# 
# ```bash
# Rscript scripts/documented/genomics/genomic_analysis.R
# ```
# 
# Refactored command-line version:
# 
# ```bash
# Rscript scripts/refactored/genomics/02_genomic_analysis.R --config config/config.yml --base-dir /path/to/private/ALS --outdir analysis/batch_06_genomics_prs_eqtl/results/genomic_analysis
# ```
# 
# **Refactor notes**
# 
# - Separate sample-level burden, gene-level burden, ALS-panel analysis, and clinical association into separate scripts.
# - Save burden matrices and crosswalk objects as RDS for TRiC/MYC/DEG burden scripts.
# - Make variant class names and MAF column names configurable because ANNOVAR headers vary.
# ============================================================

# ============================================================
# GENOMIC MUTATION BURDEN FROM ANNOTATION FILE
# ============================================================

library(data.table)
library(dplyr)
library(ggplot2)

# ============================================================
# STEP 1: LOAD ONE CHROMOSOME ANNOTATION FILE AND INSPECT
# ============================================================

anno <- fread(
  "ALS.chr2.anno_report_exon.txt",
  header    = TRUE,
  sep       = "\t",
  fill      = TRUE,
  quote     = ""         # critical — apostrophes in genotype calls
)

cat("Dimensions:", dim(anno), "\n")
cat("First 35 column names:\n")
print(colnames(anno)[1:35])

# Identify where sample columns start
# Sample columns start with CASE- or CTRL-
sample_cols <- grep("^(CASE|CTRL)-", colnames(anno), value = TRUE)
annot_cols  <- setdiff(colnames(anno), sample_cols)

cat("\nAnnotation columns:", length(annot_cols), "\n")
cat("Sample columns:", length(sample_cols), "\n")
cat("\nFirst 5 sample IDs in annotation file:\n")
print(head(sample_cols, 5))

# Check genotype call format
cat("\nExample genotype values (first variant, first 10 samples):\n")
print(unlist(anno[1, sample_cols[1:10], with=FALSE]))

# ============================================================
# STEP 2: BUILD SAMPLE ID CROSSWALK
# annotation file: CASE-NEUMN012EVP  (expression sample ID)
# sporadic_meta:   NYGC_CGND_ID = CGND-HDA-00239
# Need: CASE-NEU... → CGND-HDA-... → subtype label
# ============================================================

# Load sample mapping (already in your environment from pipeline)
# sample_map_tx has Sample_ID (CASE-NEU...) and Participant_ID (CASE-NEU...)
# sporadic_meta has NYGC_CGND_ID

# Check what ID columns are available in sporadic_meta
cat("\nID columns in sporadic_meta:\n")
print(grep("ID|id|CGND|Sample|Participant",
           colnames(sporadic_meta), ignore.case=TRUE, value=TRUE))

# Build crosswalk: expression sample ID → CGND ID → subtype
# sporadic_meta should have both Sample_ID and NYGC_CGND_ID
crosswalk <- sporadic_meta %>%
  filter(!is.na(NYGC_CGND_ID), NYGC_CGND_ID != "") %>%
  dplyr::select(
    expr_sample_id = Participant_ID,    # CASE-NEUMN012EVP format
    cgnd_id        = NYGC_CGND_ID  # CGND-HDA-00239 format
  ) %>%
  distinct()

cat("\nCrosswalk dimensions:", dim(crosswalk), "\n")
cat("Example crosswalk rows:\n")
print(head(crosswalk, 5))

# annotation file columns are CASE-NEU... → match directly to expr_sample_id
# VCF IDs are CGND-HDA-...-b38 → strip -b38 to match cgnd_id
# We don't need the VCF directly — annotation file already has genotypes

# Which annotation sample columns have a crosswalk match?
matched_samples <- intersect(sample_cols, crosswalk$expr_sample_id)
cat("\nAnnotation sample columns matched to sporadic metadata:",
    length(matched_samples), "of", length(sample_cols), "\n")

# If 0 matches, the annotation file uses CGND IDs not CASE IDs
# Check:
if (length(matched_samples) == 0) {
  cat("\nNo direct match — checking if annotation columns are CGND IDs...\n")
  # Strip -b38 from sample_cols and try
  sample_cols_stripped <- sub("-b38$", "", sample_cols)
  matched_cgnd <- intersect(sample_cols_stripped, crosswalk$cgnd_id)
  cat("Matched after stripping -b38:", length(matched_cgnd), "\n")
  
  # If matched, rename for easy lookup
  names(sample_cols) <- sample_cols_stripped
}

# ============================================================
# STEP 3: FILTER TO DAMAGING VARIANTS
# ============================================================

# ============================================================
# GENOMIC BURDEN — CORRECTED PIPELINE
# 360 matched samples, CASE-NEU IDs, underscore ExonicFunc
# ============================================================

library(data.table)
library(dplyr)
library(aricode)
library(ggplot2)

# ============================================================
# STEP 1: DEFINE DAMAGING CLASSES (match your ExonicFunc exactly)
# ============================================================

damaging_classes <- c(
  "nonsynonymous_SNV",
  "stopgain",
  "stoploss",
  "frameshift_deletion",
  "frameshift_insertion"
)

lof_classes <- c(
  "stopgain",
  "stoploss",
  "frameshift_deletion",
  "frameshift_insertion"
)

# use_cols = matched CASE-NEU sample columns
use_cols <- matched_samples  # 360 samples

# col_to_id: maps CASE-NEU... → CGND-HDA-...
col_to_id <- setNames(
  crosswalk$cgnd_id[match(use_cols, crosswalk$expr_sample_id)],
  use_cols
)

cat("Samples to process:", length(use_cols), "\n")

# ============================================================
# STEP 2: GENOTYPE PARSER
# Format is '0/0  '0/1  '1/1 (leading apostrophe)
# ============================================================

has_alt <- function(gt_vec) {
  gt_clean <- gsub("'|\\s", "", as.character(gt_vec))
  as.integer(grepl("^(0[/|]1|1[/|]0|1[/|]1)$", gt_clean))
}

# ============================================================
# STEP 3: PER-CHROMOSOME BURDEN FUNCTION
# ============================================================

parse_burden_chr <- function(filepath, use_cols, col_to_id,
                             damaging_classes, lof_classes) {
  cat("  Reading:", basename(filepath), "... ")
  
  anno <- fread(filepath, header=TRUE, sep="\t",
                fill=TRUE, quote="", showProgress=FALSE)
  
  # Filter to damaging variants
  anno_d <- anno[ExonicFunc %in% damaging_classes]
  cat(nrow(anno_d), "damaging variants\n")
  
  if (nrow(anno_d) == 0) return(NULL)
  
  # Keep only matched sample columns present in this file
  present <- intersect(use_cols, colnames(anno_d))
  if (length(present) == 0) {
    cat("  WARNING: no matched samples in this file\n")
    return(NULL)
  }
  
  # Extract and convert genotype matrix
  gt_mat  <- as.matrix(anno_d[, present, with=FALSE])
  alt_mat <- apply(gt_mat, 2, has_alt)
  
  # Boolean row index for LoF variants
  lof_rows <- anno_d$ExonicFunc %in% lof_classes
  
  data.frame(
    col_id     = present,
    cgnd_id    = col_to_id[present],
    n_damaging = colSums(alt_mat,                         na.rm=TRUE),
    n_lof      = colSums(alt_mat[lof_rows,  , drop=FALSE], na.rm=TRUE),
    n_nonsyn   = colSums(
      alt_mat[anno_d$ExonicFunc == "nonsynonymous_SNV", , drop=FALSE],
      na.rm=TRUE
    ),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# STEP 4: RUN ACROSS ALL CHROMOSOMES
# ============================================================

exon_files <- sort(list.files(".", 
                              pattern = "anno_report_exon\\.txt$",
                              full.names = TRUE))
cat("Chromosome files found:", length(exon_files), "\n")
print(basename(exon_files))

burden_list <- lapply(exon_files, parse_burden_chr,
                      use_cols        = use_cols,
                      col_to_id       = col_to_id,
                      damaging_classes = damaging_classes,
                      lof_classes     = lof_classes)

# Aggregate: sum across chromosomes per sample
burden_all <- bind_rows(burden_list) %>%
  group_by(col_id, cgnd_id) %>%
  summarise(
    total_damaging = sum(n_damaging, na.rm=TRUE),
    total_lof      = sum(n_lof,      na.rm=TRUE),
    total_nonsyn   = sum(n_nonsyn,   na.rm=TRUE),
    .groups = "drop"
  )

cat("\n=== Genome-wide burden summary (n =", nrow(burden_all), "samples) ===\n")
print(summary(burden_all[, c("total_damaging","total_lof","total_nonsyn")]))

write.csv(burden_all, "genomic_burden_per_sample_all_chrom.csv",
          row.names=FALSE)

# ============================================================
# STEP 5: MERGE WITH SUBTYPE AND CLINICAL DATA
# expr_sample_id (CASE-NEU) → cluster_final
# ============================================================

# ============================================================
# FIXED ID CHAIN
# ============================================================

# Step 1: Add Participant_ID to cluster_key by stripping Sample_ID suffix
# CASE-NEUAA599TMX-5310-T → CASE-NEUAA599TMX
cluster_key <- ax_clean_df[, c("sample_id", "cluster_clean")] %>%
  mutate(
    Participant_ID = sub("-[0-9]+-T$", "", sample_id)
  )

# Verify stripping worked
cat("Sample ID examples before/after stripping:\n")
print(head(cluster_key[, c("sample_id","Participant_ID")], 5))
# Should show: CASE-NEUAA599TMX-5310-T → CASE-NEUAA599TMX

# Step 2: Verify crosswalk uses same short form
cat("\nCrosswalk expr_sample_id examples:\n")
print(head(crosswalk$expr_sample_id, 5))
# Should show: CASE-NEUAA599TMX

# Step 3: Check overlap now
overlap_check <- intersect(cluster_key$Participant_ID, crosswalk$expr_sample_id)
cat("\nOverlap after stripping:", length(overlap_check), "samples\n")

# Step 4: Add Participant_ID to df_nfl_uniq the same way
df_nfl_uniq <- df_nfl_uniq %>%
  mutate(
    Participant_ID = sub("-[0-9]+-T$", "", sample_id)
  )

# Step 5: Clean join chain — everything joins on Participant_ID
burden_subtype <- burden_all %>%
  # cgnd_id → Participant_ID (short form)
  left_join(
    crosswalk %>% rename(Participant_ID = expr_sample_id),
    by = "cgnd_id"
  ) %>%
  # Participant_ID → cluster_final
  left_join(
    cluster_key %>% dplyr::select(Participant_ID, cluster_clean, sample_id),
    by = "Participant_ID"
  ) %>%
  # Participant_ID → clinical data
  left_join(
    df_nfl_uniq %>%
      dplyr::select(Participant_ID, logNfL,
                    ALSFRS_R_PROGRESSION_SLOPE,
                    AGE_AT_SYMPTOM_ONSET, SEX),
    by = "Participant_ID"
  ) %>%
  filter(!is.na(cluster_clean))

cat("\n=== Final merged table ===\n")
cat("Rows:", nrow(burden_subtype), "\n")
print(table(burden_subtype$cluster_clean))

# Sanity check — no duplicate samples
n_dup <- sum(duplicated(burden_subtype$cgnd_id))
cat("Duplicate CGND IDs:", n_dup, "\n")
# Should be 0

# Preview the joined table
cat("\nExample rows:\n")
print(head(burden_subtype[, c("cgnd_id","Participant_ID",
                              "cluster_clean","total_damaging",
                              "logNfL")], 5))


cat("\n=== Sporadic samples with burden + subtype ===\n")
print(table(burden_subtype$cluster_clean))
cat("Total:", nrow(burden_subtype), "\n")

# Sanity check — burden distributions look reasonable?
cat("\nMedian damaging variants per sample:",
    median(burden_subtype$total_damaging), "\n")
cat("Range:", range(burden_subtype$total_damaging), "\n")

# ============================================================
# BURDEN ANALYSIS — ALL CHECKS PASSED, PROCEED
# Chr2 LoF median = 19 → genome-wide ~237 (expected range 100-400) ✓
# Chr2 damaging median = 748 → genome-wide ~12,470 ✓
# ============================================================

subtype_colors <- c("S-ALS-1"="#E64B35",
                    "S-ALS-2"="#4DBBD5",
                    "S-ALS-3"="#00A087")

# ============================================================
# STEP 6: BURDEN BY SUBTYPE
# ============================================================

cat("=== Genome-wide LoF quantiles (all chromosomes) ===\n")
print(quantile(burden_subtype$total_lof, probs=c(0,.25,.5,.75,1)))

cat("\n=== Burden by subtype ===\n")
burden_subtype %>%
  group_by(cluster_clean) %>%
  summarise(
    n                  = n(),
    mean_damaging      = round(mean(total_damaging), 0),
    median_damaging    = median(total_damaging),
    mean_lof           = round(mean(total_lof), 1),
    median_lof         = median(total_lof),
    n_with_lof         = sum(total_lof > 0),
    .groups = "drop"
  ) %>%
  print()

# ANOVA: total damaging
aov_dam <- aov(total_damaging ~ cluster_clean, data=burden_subtype)
cat("\nANOVA — total damaging burden ~ subtype:\n")
print(summary(aov_dam))

# ANOVA: LoF
aov_lof <- aov(total_lof ~ cluster_clean, data=burden_subtype)
cat("\nANOVA — LoF burden ~ subtype:\n")
print(summary(aov_lof))

# Tukey if significant
if (summary(aov_dam)[[1]]$`Pr(>F)`[1] < 0.05) {
  cat("\nTukey HSD — total damaging:\n")
  print(TukeyHSD(aov_dam))
}
if (summary(aov_lof)[[1]]$`Pr(>F)`[1] < 0.05) {
  cat("\nTukey HSD — LoF burden:\n")
  print(TukeyHSD(aov_lof))
}

# Kruskal-Wallis as non-parametric check
# (burden counts are right-skewed, KW is more robust)
cat("\nKruskal-Wallis — total damaging:\n")
print(kruskal.test(total_damaging ~ cluster_clean, data=burden_subtype))

cat("\nKruskal-Wallis — LoF burden:\n")
print(kruskal.test(total_lof ~ cluster_clean, data=burden_subtype))

# ============================================================
# STEP 7: BURDEN VS CLINICAL OUTCOMES
# ============================================================

df_bn <- burden_subtype %>% filter(!is.na(logNfL))
df_bs <- burden_subtype %>% filter(!is.na(ALSFRS_R_PROGRESSION_SLOPE))

cat("\n=== n for clinical models ===\n")
cat("NfL:", nrow(df_bn), "\n")
cat("ALSFRS-R slope:", nrow(df_bs), "\n")

# NfL ~ damaging burden
m_nfl_dam <- lm(logNfL ~ log1p(total_damaging) + AGE_AT_SYMPTOM_ONSET + SEX,
                data = df_bn)
cat("\n=== NfL ~ total damaging burden ===\n")
print(summary(m_nfl_dam)$coefficients)

# NfL ~ LoF burden
m_nfl_lof <- lm(logNfL ~ log1p(total_lof) + AGE_AT_SYMPTOM_ONSET + SEX,
                data = df_bn)
cat("\n=== NfL ~ LoF burden ===\n")
print(summary(m_nfl_lof)$coefficients)

# ALSFRS-R ~ damaging burden
m_slope_dam <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ log1p(total_damaging) +
                    AGE_AT_SYMPTOM_ONSET + SEX,
                  data = df_bs)
cat("\n=== ALSFRS-R ~ total damaging burden ===\n")
print(summary(m_slope_dam)$coefficients)

# ALSFRS-R ~ LoF burden
m_slope_lof <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ log1p(total_lof) +
                    AGE_AT_SYMPTOM_ONSET + SEX,
                  data = df_bs)
cat("\n=== ALSFRS-R ~ LoF burden ===\n")
print(summary(m_slope_lof)$coefficients)

# FDR correction across all 4 clinical burden tests
burden_clinical_pvals <- c(
  summary(m_nfl_dam)$coefficients["log1p(total_damaging)", "Pr(>|t|)"],
  summary(m_nfl_lof)$coefficients["log1p(total_lof)",      "Pr(>|t|)"],
  summary(m_slope_dam)$coefficients["log1p(total_damaging)","Pr(>|t|)"],
  summary(m_slope_lof)$coefficients["log1p(total_lof)",    "Pr(>|t|)"]
)
burden_fdr <- p.adjust(burden_clinical_pvals, method="BH")

cat("\n=== Clinical burden tests — FDR summary ===\n")
print(data.frame(
  test     = c("NfL~damaging","NfL~LoF","ALSFRS~damaging","ALSFRS~LoF"),
  raw_p    = round(burden_clinical_pvals, 4),
  FDR      = round(burden_fdr, 4)
))

# ============================================================
# STEP 8: VISUALIZATION
# ============================================================

pdf("genomic_burden_by_subtype.pdf", width=10, height=4)
par(mfrow=c(1,3), mar=c(5,5,3,2))

for (v in c("total_damaging","total_lof","total_nonsyn")) {
  lab <- c(total_damaging="Total damaging",
           total_lof="Loss-of-function",
           total_nonsyn="Nonsynonymous SNV")[v]
  boxplot(
    log1p(burden_subtype[[v]]) ~ burden_subtype$cluster_clean,
    col     = subtype_colors[levels(burden_subtype$cluster_clean)],
    ylab    = paste0("log(1+n) ", lab),
    xlab    = "", las=2, main=lab, outline=FALSE
  )
  stripchart(
    log1p(burden_subtype[[v]]) ~ burden_subtype$cluster_clean,
    vertical=TRUE, method="jitter",
    pch=19, cex=0.4, col="black", add=TRUE
  )
}
dev.off()

# NfL vs LoF scatter coloured by subtype
pdf("nfl_vs_lof_burden_scatter.pdf", width=5, height=4)
par(mar=c(5,5,2,2))
plot(
  log1p(df_bn$total_lof),
  df_bn$logNfL,
  pch  = 19, cex = 0.7,
  col  = subtype_colors[as.character(df_bn$cluster_clean)],
  xlab = "log(1 + LoF variant count)",
  ylab = "log10 NfL (pg/mL)",
  main = "NfL vs LoF burden"
)
abline(lm(logNfL ~ log1p(total_lof), data=df_bn), lwd=2)
legend("topright", legend=names(subtype_colors),
       col=subtype_colors, pch=19, bty="n", cex=0.8)
dev.off()

# ============================================================
# SENSITIVITY: Is the negative LoF~NfL driven by outliers?
# ============================================================

# Check for high-burden outliers
cat("High LoF burden samples (>700):\n")
print(df_bn %>% filter(total_lof > 700) %>%
        dplyr::select(cgnd_id, cluster_clean, total_lof, logNfL))

# Rerun without top 5% LoF burden
lof_95 <- quantile(df_bn$total_lof, 0.95)
df_bn_trim <- df_bn %>% filter(total_lof <= lof_95)

m_nfl_lof_trim <- lm(
  logNfL ~ log1p(total_lof) + AGE_AT_SYMPTOM_ONSET + SEX,
  data = df_bn_trim
)
cat("\n=== NfL ~ LoF (top 5% removed, n=",
    nrow(df_bn_trim), ") ===\n")
print(summary(m_nfl_lof_trim)$coefficients["log1p(total_lof)", ])

# If direction and significance hold after trimming → result is robust
# If it flips or becomes NS → driven by outliers, report cautiously

# Does LoF ~ NfL association differ by subtype?
m_nfl_lof_int <- lm(
  logNfL ~ log1p(total_lof) * cluster_clean +
    AGE_AT_SYMPTOM_ONSET + SEX,
  data = df_bn
)

cat("=== Interaction: LoF burden × subtype on NfL ===\n")
print(summary(m_nfl_lof_int)$coefficients)

# ANOVA to test whether interaction term improves fit
m_nfl_lof_main <- lm(
  logNfL ~ log1p(total_lof) + cluster_clean +
    AGE_AT_SYMPTOM_ONSET + SEX,
  data = df_bn
)

cat("\n=== Model comparison: main effects vs interaction ===\n")
print(anova(m_nfl_lof_main, m_nfl_lof_int))

# Subtype-stratified slopes
cat("\n=== LoF ~ NfL slope within each subtype ===\n")
for (s in c("S-ALS-1","S-ALS-2","S-ALS-3")) {
  df_sub <- df_bn %>% filter(cluster_clean == s)
  m_sub  <- lm(logNfL ~ log1p(total_lof) + AGE_AT_SYMPTOM_ONSET + SEX,
               data = df_sub)
  coef_row <- summary(m_sub)$coefficients["log1p(total_lof)", ]
  cat(sprintf("  %s (n=%d): β=%.3f, p=%.4f\n",
              s, nrow(df_sub),
              coef_row["Estimate"],
              coef_row["Pr(>|t|)"]))
}





# Rare variant burden analysis
# ============================================================
# NEXT SESSION — START HERE
# Rare variant filter + genomic-only analyses
# ============================================================

# ============================================================
# STEP 1: APPLY MAF FILTER TO ANNOTATION FILE
# Columns: "Max Freq of 3" and "Freq:1KG/ESP/ExAC"
# ============================================================

# Check frequency column format on chr2
anno_chr2 <- fread("ALS.chr2.anno_report_exon.txt",
                   header=TRUE, sep="\t", fill=TRUE, quote="",
                   showProgress=FALSE)

freq_cols <- grep("Freq|freq|MAF|maf|1KG|ESP|ExAC",
                  colnames(anno_chr2), ignore.case=TRUE, value=TRUE)
cat("Frequency columns found:", paste(freq_cols, collapse=", "), "\n")

# Check example values
cat("\nExample Max Freq values:\n")
print(head(anno_chr2[["Max Freq of 3"]], 20))

# The "Max Freq of 3" column is the maximum allele frequency
# across 1000 Genomes, ESP, and ExAC
# Values like "0.0003", ".", "0" — need to parse as numeric

# Parse and check distribution
max_freq <- suppressWarnings(
  as.numeric(anno_chr2[["Max Freq of 3"]])
)
cat("\nFrequency distribution:\n")
print(quantile(max_freq, probs=c(0,.25,.5,.75,.9,.95,.99,1), na.rm=TRUE))

# How many variants pass MAF < 0.01?
cat("\nVariants with Max Freq < 0.01 (rare):",
    sum(max_freq < 0.01 | is.na(max_freq)), "\n")
cat("Variants with Max Freq >= 0.01 (common):",
    sum(max_freq >= 0.01, na.rm=TRUE), "\n")

# ============================================================
# STEP 2: RARE VARIANT BURDEN FUNCTION
# Same as parse_burden_chr but with MAF filter added
# ============================================================

# ============================================================
# STEP 2: RUN RARE BURDEN LOOP ACROSS ALL CHROMOSOMES
# Confirmed: "." → NA → treated as rare. Format is clean.
# ============================================================

parse_rare_burden_chr <- function(filepath, use_cols, col_to_id,
                                  damaging_classes, lof_classes,
                                  maf_threshold = 0.01) {
  cat("  Reading:", basename(filepath), "... ")
  
  anno <- fread(filepath, header=TRUE, sep="\t",
                fill=TRUE, quote="", showProgress=FALSE)
  
  # Parse frequency
  freq_col <- grep("Max Freq", colnames(anno), value=TRUE)[1]
  anno$max_freq_num <- suppressWarnings(as.numeric(anno[[freq_col]]))
  
  # Filter: damaging AND (rare OR unknown frequency)
  anno_rare <- anno[
    ExonicFunc %in% damaging_classes &
      (is.na(max_freq_num) | max_freq_num < maf_threshold)
  ]
  cat(nrow(anno_rare), "rare damaging variants\n")
  
  if (nrow(anno_rare) == 0) return(NULL)
  
  present <- intersect(use_cols, colnames(anno_rare))
  if (length(present) == 0) return(NULL)
  
  gt_mat  <- as.matrix(anno_rare[, present, with=FALSE])
  alt_mat <- apply(gt_mat, 2, has_alt)
  
  lof_rows <- anno_rare$ExonicFunc %in% lof_classes
  
  data.frame(
    col_id          = present,
    cgnd_id         = col_to_id[present],
    n_rare_damaging = colSums(alt_mat, na.rm=TRUE),
    n_rare_lof      = colSums(alt_mat[lof_rows, , drop=FALSE], na.rm=TRUE),
    n_rare_nonsyn   = colSums(
      alt_mat[anno_rare$ExonicFunc == "nonsynonymous_SNV", , drop=FALSE],
      na.rm=TRUE
    ),
    stringsAsFactors = FALSE
  )
}

# Run across all chromosomes — expect ~5–15 min
cat("Starting rare burden loop...\n")
rare_burden_list <- lapply(exon_files, parse_rare_burden_chr,
                           use_cols         = use_cols,
                           col_to_id        = col_to_id,
                           damaging_classes = damaging_classes,
                           lof_classes      = lof_classes,
                           maf_threshold    = 0.01)

# Aggregate across chromosomes per sample
rare_burden_all <- bind_rows(rare_burden_list) %>%
  group_by(col_id, cgnd_id) %>%
  summarise(
    rare_damaging = sum(n_rare_damaging, na.rm=TRUE),
    rare_lof      = sum(n_rare_lof,      na.rm=TRUE),
    rare_nonsyn   = sum(n_rare_nonsyn,   na.rm=TRUE),
    .groups = "drop"
  )

cat("\n=== Rare variant burden summary (all chromosomes) ===\n")
print(summary(rare_burden_all[, c("rare_damaging","rare_lof","rare_nonsyn")]))
cat("n samples:", nrow(rare_burden_all), "\n")


# ============================================================
# STEP 3: RARE VARIANT BURDEN BY SUBTYPE
# ============================================================

rare_burden_sub <- rare_burden_all %>%
  left_join(crosswalk %>% rename(Participant_ID = expr_sample_id),
            by = "cgnd_id") %>%
  left_join(ax_clean_df %>%
              mutate(Participant_ID = sub("-[0-9]+-T$", "", sample_id)) %>%
              dplyr::select(Participant_ID, cluster_clean),
            by = "Participant_ID") %>%
  left_join(df_nfl_uniq %>%
              mutate(Participant_ID = sub("-[0-9]+-T$", "", sample_id)) %>%
              dplyr::select(Participant_ID, logNfL,
                            ALSFRS_R_PROGRESSION_SLOPE,
                            AGE_AT_SYMPTOM_ONSET, SEX),
            by = "Participant_ID") %>%
  filter(!is.na(cluster_clean)) %>%
  distinct(cgnd_id, .keep_all = TRUE)

cat("n after merge:", nrow(rare_burden_sub), "\n")
cat("cluster_final distribution:\n")
print(table(rare_burden_sub$cluster_clean, useNA = "ifany"))

# --- Burden by subtype ---
cat("\n=== RARE variant burden by subtype ===\n")
rare_burden_sub %>%
  group_by(cluster_clean) %>%
  summarise(
    n               = n(),
    mean_rare_dam   = round(mean(rare_damaging), 1),
    median_rare_dam = median(rare_damaging),
    mean_rare_lof   = round(mean(rare_lof), 2),
    median_rare_lof = median(rare_lof)
  ) %>%
  print()

# --- ANOVA ---
aov_rare_dam <- aov(rare_damaging ~ cluster_clean, data = rare_burden_sub)
aov_rare_lof <- aov(rare_lof      ~ cluster_clean, data = rare_burden_sub)

cat("\nANOVA rare damaging p:", summary(aov_rare_dam)[[1]]$`Pr(>F)`[1], "\n")
cat("ANOVA rare LoF p:",       summary(aov_rare_lof)[[1]]$`Pr(>F)`[1], "\n")

# Kruskal-Wallis as non-parametric check
cat("\nKruskal-Wallis rare damaging p:",
    kruskal.test(rare_damaging ~ cluster_clean, data = rare_burden_sub)$p.value, "\n")
cat("Kruskal-Wallis rare LoF p:",
    kruskal.test(rare_lof ~ cluster_clean, data = rare_burden_sub)$p.value, "\n")

# Post-hoc only if significant
if (summary(aov_rare_dam)[[1]]$`Pr(>F)`[1] < 0.05) {
  cat("\nTukey HSD — rare damaging:\n")
  print(TukeyHSD(aov_rare_dam))
}
if (summary(aov_rare_lof)[[1]]$`Pr(>F)`[1] < 0.05) {
  cat("\nTukey HSD — rare LoF:\n")
  print(TukeyHSD(aov_rare_lof))
}

# --- Outlier check ---
cat("\nTop 5 rare_damaging outliers:\n")
rare_burden_sub %>%
  arrange(desc(rare_damaging)) %>%
  dplyr::select(cgnd_id, cluster_clean, rare_damaging, rare_lof) %>%
  head(5) %>% print()


# ============================================================
# STEP 4: Paul Q1 — Rare variant burden vs clinical
# No cluster labels; genome-wide association with NfL / ALSFRS-R
# ============================================================

library(dunn.test)   # install.packages("dunn.test") if needed

# --- Dunn post-hoc on marginal KW (rare damaging) ---
cat("=== Dunn post-hoc: rare damaging by subtype ===\n")
dunn.test(rare_burden_sub$rare_damaging,
          rare_burden_sub$cluster_clean,
          method = "bh")

# --- Paul Q1: Rare LoF ~ NfL (no cluster covariate) ---
df_rare_nfl   <- rare_burden_sub %>% filter(!is.na(logNfL))
df_rare_slope <- rare_burden_sub %>% filter(!is.na(ALSFRS_R_PROGRESSION_SLOPE))

cat("\n=== Rare LoF ~ NfL (pooled, no cluster) n =", nrow(df_rare_nfl), "===\n")
m_rare_lof_nfl <- lm(logNfL ~ log1p(rare_lof) + AGE_AT_SYMPTOM_ONSET + SEX,
                     data = df_rare_nfl)
print(summary(m_rare_lof_nfl)$coefficients)

cat("\n=== Rare damaging ~ NfL (pooled) ===\n")
m_rare_dam_nfl <- lm(logNfL ~ log1p(rare_damaging) + AGE_AT_SYMPTOM_ONSET + SEX,
                     data = df_rare_nfl)
print(summary(m_rare_dam_nfl)$coefficients)

# --- Outlier sensitivity: remove top 1% rare_damaging ---
q99 <- quantile(df_rare_nfl$rare_damaging, 0.99)
df_rare_nfl_clean <- df_rare_nfl %>% filter(rare_damaging <= q99)
cat("\n=== Rare LoF ~ NfL — outliers removed (n =", nrow(df_rare_nfl_clean), ") ===\n")
m_rare_lof_nfl_clean <- lm(logNfL ~ log1p(rare_lof) + AGE_AT_SYMPTOM_ONSET + SEX,
                           data = df_rare_nfl_clean)
print(summary(m_rare_lof_nfl_clean)$coefficients)

# --- ALSFRS-R ---
cat("\n=== Rare LoF ~ ALSFRS-R slope (pooled, n =", nrow(df_rare_slope), ") ===\n")
m_rare_lof_slope <- lm(ALSFRS_R_PROGRESSION_SLOPE ~ log1p(rare_lof) +
                         AGE_AT_SYMPTOM_ONSET + SEX,
                       data = df_rare_slope)
print(summary(m_rare_lof_slope)$coefficients)

# ============================================================
# STEP 5: Paul Q2 — Genomic-only subtyping
# Use ALL damaging variants (no MAF filter) in sporadic samples
# Gene x sample burden matrix → ConsensusClusterPlus
# ============================================================

# NOTE: Using all_damaging per your clarification — not rare-only

parse_gene_burden_chr <- function(filepath, use_cols, col_to_id,
                                  damaging_classes) {
  cat("  Gene burden:", basename(filepath), "... ")
  anno <- fread(filepath, header=TRUE, sep="\t",
                fill=TRUE, quote="", showProgress=FALSE)
  
  anno_dam <- anno[ExonicFunc %in% damaging_classes]
  cat(nrow(anno_dam), "damaging variants\n")
  
  if (nrow(anno_dam) == 0) return(NULL)
  
  present <- intersect(use_cols, colnames(anno_dam))
  if (length(present) == 0) return(NULL)
  
  gt_mat  <- as.matrix(anno_dam[, present, with=FALSE])
  alt_mat <- apply(gt_mat, 2, has_alt)   # your existing function
  
  # Gene column
  gene_col <- intersect(c("Gene","gene","Gene.refGene"), colnames(anno_dam))[1]
  genes <- anno_dam[[gene_col]]
  
  # Sum per gene per sample
  gene_burden <- rowsum(alt_mat, group = genes)
  
  data.frame(
    gene    = rownames(gene_burden),
    gene_burden,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

cat("\n=== Building gene-level burden matrix (all damaging, sporadic) ===\n")
gene_list <- lapply(exon_files, parse_gene_burden_chr,
                    use_cols         = use_cols,
                    col_to_id        = col_to_id,
                    damaging_classes = damaging_classes)

# Combine across chromosomes: sum per gene
gene_list_valid <- Filter(Negate(is.null), gene_list)
gene_burden_combined <- bind_rows(gene_list_valid) %>%
  group_by(gene) %>%
  summarise(across(everything(), sum, na.rm=TRUE)) %>%
  tibble::column_to_rownames("gene")

cat("Gene x sample matrix:", nrow(gene_burden_combined), "genes x",
    ncol(gene_burden_combined), "samples\n")

# --- Filter: keep genes with damaging variant in ≥5% of samples ---
min_samples <- round(0.05 * ncol(gene_burden_combined))
gene_prev   <- rowSums(gene_burden_combined > 0)
gene_burden_filt <- gene_burden_combined[gene_prev >= min_samples, ]
cat("After 5% prevalence filter:", nrow(gene_burden_filt), "genes\n")

# Transpose for clustering (samples as rows)
burden_mat_t <- t(as.matrix(gene_burden_filt))

# ============================================================
# FIX PCA + DIMENSIONALITY REDUCTION FOR GENOMIC CLUSTERING
# ============================================================

# Step 1: Remove zero-variance genes (same count in every sample)
gene_var <- apply(burden_mat_t, 2, var)
burden_mat_nozero <- burden_mat_t[, gene_var > 0]
cat("After zero-variance removal:", ncol(burden_mat_nozero), "genes\n")

# Step 2: PCA on filtered matrix
pca_geno <- prcomp(burden_mat_nozero, scale. = TRUE)

# Variance explained
var_exp <- summary(pca_geno)$importance[2, 1:20] * 100
cat("\nVariance explained PC1-20:\n")
print(round(var_exp, 2))

cum_var <- summary(pca_geno)$importance[3, 1:20] * 100
cat("\nCumulative variance PC1-20:\n")
print(round(cum_var, 1))

# How many PCs to reach 50% variance?
n_pcs_50 <- which(cum_var >= 50)[1]
cat("\nPCs needed for 50% variance:", n_pcs_50, "\n")

# ============================================================
# GENOMIC SUBTYPING RESULT — REPORT WITHOUT CCP
# The PCA structure already answers the question
# ============================================================

# 1. Confirm top PCs don't separate transcriptomic subtypes
cat("=== ANOVA: genomic PC1-10 vs transcriptomic subtype ===\n")
pc_df <- data.frame(
  col_id = rownames(pca_geno$x),
  pca_geno$x[, 1:10]
) %>%
  mutate(cgnd_id = col_to_id[col_id]) %>%
  left_join(
    rare_burden_sub %>% dplyr::select(cgnd_id, cluster_clean),
    by = "cgnd_id"
  ) %>%
  filter(!is.na(cluster_clean))

for (pc in paste0("PC", 1:10)) {
  pval <- summary(aov(as.formula(paste(pc, "~ cluster_clean")),
                      data = pc_df))[[1]]$`Pr(>F)`[1]
  cat(pc, ": p =", round(pval, 4),
      ifelse(pval < 0.05, " *", ""), "\n")
}

# 2. Silhouette at k=3 using top 30 PCs (direct test, no CCP needed)
# If transcriptomic subtypes impose no structure on genomic space,
# silhouette will be near zero
library(cluster)

pc30 <- pca_geno$x[, 1:30]

# Match sample order to transcriptomic labels
pc30_matched <- pc30[rownames(pc30) %in% pc_df$col_id, ]
labels_matched <- pc_df$cluster_clean[
  match(rownames(pc30_matched), pc_df$col_id)
]

dist_geno <- dist(pc30_matched, method = "euclidean")
sil_geno  <- silhouette(as.integer(labels_matched), dist_geno)
cat("\n=== Silhouette of TRANSCRIPTOMIC subtypes in GENOMIC space ===\n")
cat("Mean silhouette:", round(mean(sil_geno[, 3]), 3), "\n")
cat("(Near 0 = transcriptomic subtypes have no genomic structure)\n")
cat("(Compare: transcriptomic silhouette was 0.33)\n")

# Per-subtype silhouette
sil_df <- data.frame(
  cluster   = labels_matched,
  sil_width = sil_geno[, 3]
)
cat("\nPer-subtype mean silhouette in genomic space:\n")
print(sil_df %>%
        group_by(cluster) %>%
        summarise(mean_sil = round(mean(sil_width), 3), n = n()))

# 3. ARI: k-means at k=3 on genomic PCs vs transcriptomic subtypes
set.seed(42)
km_geno <- kmeans(pc30_matched, centers = 3, nstart = 50)

ari_geno <- adjustedRandIndex(km_geno$cluster, as.integer(labels_matched))
cat("\n=== ARI: k-means genomic k=3 vs transcriptomic subtypes ===\n")
cat("ARI:", round(ari_geno, 3), "\n")
cat("(ARI near 0 = genomic and transcriptomic subtypes are independent)\n")


# Blind clustering for mutation data
# ============================================================
# BLIND GENOMIC CLUSTERING — No reference to transcriptomic labels
# ============================================================

library(ConsensusClusterPlus)

# Use top 30 PCs (captures ~14% of variance — best available)
pc30 <- pca_geno$x[, 1:30]

cat("=== Blind ConsensusClusterPlus on genomic PC scores ===\n")
cat("Input:", nrow(pc30), "samples x", ncol(pc30), "PCs\n\n")

results_blind <- ConsensusClusterPlus(
  t(pc30),          # features x samples
  maxK       = 6,
  reps       = 500,
  pItem      = 0.80,
  pFeature   = 1.0,
  clusterAlg = "hc",
  distance   = "euclidean",
  seed       = 42,
  plot       = NULL,
  verbose    = FALSE
)

# PAC
calcPAC <- function(res, k) {
  cm <- res[[k]]$consensusMatrix
  length(cm[cm > 0.1 & cm < 0.9]) / length(cm)
}
pac_blind <- sapply(2:6, function(k) calcPAC(results_blind, k))

cat("=== PAC by k — blind genomic clustering ===\n")
print(data.frame(k = 2:6, PAC = round(pac_blind, 3)))

best_k <- which.min(pac_blind) + 1
cat("Minimum PAC at k =", best_k, "(PAC =", round(min(pac_blind), 3), ")\n")
cat("\nFor reference: transcriptomic clustering PAC at k=3 was ~0.12\n")
cat("PAC > 0.25 = no stable structure\n")


# Does the k=2 "stable" genomic split associate with anything real?
labels_k2 <- results_blind[[2]]$consensusClass

geno_k2_df <- data.frame(
  col_id      = names(labels_k2),
  geno_k2     = labels_k2
) %>%
  mutate(cgnd_id = col_to_id[col_id]) %>%
  left_join(rare_burden_sub %>%
              dplyr::select(cgnd_id, cluster_clean, logNfL,
                            ALSFRS_R_PROGRESSION_SLOPE,
                            rare_damaging, rare_lof),
            by = "cgnd_id") %>%
  filter(!is.na(cluster_clean))

cat("=== k=2 genomic split: sample sizes ===\n")
print(table(geno_k2_df$geno_k2))

cat("\n=== k=2 vs transcriptomic subtype ===\n")
print(table(geno_k2_df$geno_k2, geno_k2_df$cluster_clean,
            dnn = c("Genomic k=2", "Transcriptomic")))
ari_k2 <- adjustedRandIndex(geno_k2_df$geno_k2,
                            as.integer(geno_k2_df$cluster_clean))
cat("ARI:", round(ari_k2, 3), "\n")

cat("\n=== k=2 genomic split vs NfL ===\n")
t.test(logNfL ~ geno_k2, data = geno_k2_df %>% filter(!is.na(logNfL)))

cat("\n=== k=2 genomic split vs rare burden ===\n")
cat("Rare damaging — group means:\n")
print(tapply(geno_k2_df$rare_damaging, geno_k2_df$geno_k2, mean, na.rm=TRUE))
cat("t-test p:", t.test(rare_damaging ~ geno_k2, data=geno_k2_df)$p.value, "\n")

# ============================================================
# ALS-GENE-SPECIFIC BURDEN ANALYSIS
# Uses the "Known ALS Variant" or gene-of-interest column
# in the ANNOVAR annotation files
# ============================================================

# First check what ALS-gene column looks like in chr2
file_chr2 <- exon_files[grep("\\bchr2\\b", exon_files)]

anno_chr2 <- fread(file_chr2,
                   header=TRUE, sep="\t", fill=TRUE,
                   quote="", showProgress=FALSE)

# Find ALS-relevant columns
als_cols <- grep("ALS|als|Gene|gene|known|Known|disease|Disease",
                 colnames(anno_chr2), ignore.case=TRUE, value=TRUE)
cat("Potentially relevant columns:\n")
print(als_cols)

# Define canonical ALS genes (top ~30 by evidence)
als_genes <- c("SOD1","C9orf72","TARDBP","FUS","ANG","ATXN2",
               "CHCHD10","HNRNPA1","HNRNPA2B1","KIF5A","MATR3",
               "NEK1","OPTN","PFN1","SETX","SIGMAR1","SQSTM1",
               "TBK1","UBQLN2","VAPB","VCP","DCTN1","ERBB4",
               "GLE1","LGALSL","MOBP","SCFD1","TIA1","TUBA4A",
               "ANXA11","GLT8D1","NEFH","CAV1","DNAJC7")

cat("\nALS gene list: n =", length(als_genes), "genes\n")

# ============================================================
# Parse ALS-gene burden per chromosome
# ============================================================

parse_als_burden_chr <- function(filepath, use_cols, col_to_id,
                                 damaging_classes, als_genes,
                                 maf_threshold = 0.01) {
  anno <- fread(filepath, header=TRUE, sep="\t",
                fill=TRUE, quote="", showProgress=FALSE)
  
  # Gene column
  gene_col <- intersect(c("Gene.refGene","Gene","gene"), colnames(anno))[1]
  if (is.na(gene_col)) return(NULL)
  
  # MAF filter
  freq_col <- grep("Max Freq", colnames(anno), value=TRUE)[1]
  anno$max_freq_num <- suppressWarnings(as.numeric(anno[[freq_col]]))
  
  # Filter: ALS gene + damaging + rare
  anno_als <- anno[
    get(gene_col) %in% als_genes &
      ExonicFunc %in% damaging_classes &
      (is.na(max_freq_num) | max_freq_num < maf_threshold)
  ]
  
  if (nrow(anno_als) == 0) return(NULL)
  
  present <- intersect(use_cols, colnames(anno_als))
  if (length(present) == 0) return(NULL)
  
  gt_mat  <- as.matrix(anno_als[, present, with=FALSE])
  alt_mat <- apply(gt_mat, 2, has_alt)
  
  genes <- anno_als[[gene_col]]
  
  # Per-gene per-sample burden
  gene_burden <- rowsum(alt_mat, group = genes)
  
  data.frame(
    gene    = rownames(gene_burden),
    gene_burden,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

cat("\n=== Parsing ALS-gene rare variant burden ===\n")
als_gene_list <- lapply(exon_files, parse_als_burden_chr,
                        use_cols         = use_cols,
                        col_to_id        = col_to_id,
                        damaging_classes = damaging_classes,
                        als_genes        = als_genes,
                        maf_threshold    = 0.01)

# Combine across chromosomes
als_valid <- Filter(Negate(is.null), als_gene_list)
als_combined <- bind_rows(als_valid) %>%
  group_by(gene) %>%
  summarise(across(everything(), \(x) sum(x, na.rm=TRUE))) %>%
  tibble::column_to_rownames("gene")

cat("ALS genes with rare damaging variants found:",
    nrow(als_combined), "/", length(als_genes), "\n")

# Transpose: samples x genes
als_mat <- t(as.matrix(als_combined))

# ============================================================
# Summary: how many sporadic patients carry rare ALS-gene variants?
# ============================================================

# Binary carrier status per gene
carrier_mat <- (als_mat > 0) * 1

cat("\n=== Carrier rate per ALS gene (% of sporadic samples) ===\n")
carrier_rates <- sort(colMeans(carrier_mat) * 100, decreasing=TRUE)
print(round(carrier_rates, 1))

cat("\n=== Total ALS-gene rare variant carriers (any gene) ===\n")
any_carrier <- rowSums(carrier_mat) > 0
cat("Carriers:", sum(any_carrier), "/", length(any_carrier),
    sprintf("(%.1f%%)\n", mean(any_carrier)*100))

# ============================================================
# Merge with subtype and clinical
# ============================================================

als_burden_df <- data.frame(
  col_id     = rownames(als_mat),
  als_mat,
  any_carrier = any_carrier,
  check.names = FALSE
) %>%
  mutate(cgnd_id = col_to_id[col_id]) %>%
  left_join(rare_burden_sub %>%
              dplyr::select(cgnd_id, cluster_clean, logNfL,
                            ALSFRS_R_PROGRESSION_SLOPE,
                            AGE_AT_SYMPTOM_ONSET, SEX),
            by = "cgnd_id") %>%
  filter(!is.na(cluster_clean))

# ============================================================
# Test 1: ALS-gene carrier status by subtype
# ============================================================

cat("\n=== ALS-gene rare variant carrier rate by subtype ===\n")
als_burden_df %>%
  group_by(cluster_clean) %>%
  summarise(
    n          = n(),
    n_carrier  = sum(any_carrier, na.rm=TRUE),
    pct        = round(mean(any_carrier, na.rm=TRUE)*100, 1)
  ) %>% print()

chisq_carrier <- chisq.test(
  table(als_burden_df$any_carrier, als_burden_df$cluster_clean)
)
cat("Chi-square p:", round(chisq_carrier$p.value, 4), "\n")

# ============================================================
# Test 2: Total ALS-gene rare burden by subtype
# ============================================================

als_burden_df$total_als_rare <- rowSums(als_mat[
  rownames(als_mat) %in% als_burden_df$col_id, ], na.rm=TRUE)

cat("\n=== Total ALS-gene rare burden by subtype ===\n")
als_burden_df %>%
  group_by(cluster_clean) %>%
  summarise(mean = round(mean(total_als_rare), 3),
            median = median(total_als_rare)) %>% print()

kruskal.test(total_als_rare ~ cluster_clean, data = als_burden_df)

# ============================================================
# Test 3: ALS-gene rare burden → NfL (no cluster label)
# ============================================================

df_als_nfl <- als_burden_df %>% filter(!is.na(logNfL))

cat("\n=== ALS-gene rare LoF burden → NfL ===\n")
m_als_nfl <- lm(logNfL ~ log1p(total_als_rare) +
                  AGE_AT_SYMPTOM_ONSET + SEX,
                data = df_als_nfl)
print(summary(m_als_nfl)$coefficients)

# ============================================================
# RERUN CARRIER ANALYSIS — NEFH EXCLUDED
# NEFH has hypervariable repeat region causing inflated calls
# ============================================================

# Exclude NEFH from carrier matrix
carrier_mat_clean <- carrier_mat[, colnames(carrier_mat) != "NEFH"]

any_carrier_clean <- rowSums(carrier_mat_clean) > 0
cat("=== Carrier rate WITHOUT NEFH ===\n")
cat("Carriers:", sum(any_carrier_clean), "/", length(any_carrier_clean),
    sprintf("(%.1f%%)\n", mean(any_carrier_clean)*100))

# Update als_burden_df
als_burden_df$any_carrier_clean    <- any_carrier_clean[als_burden_df$col_id]
als_burden_df$total_als_rare_clean <- rowSums(
  carrier_mat_clean[rownames(carrier_mat_clean) %in% als_burden_df$col_id, ],
  na.rm=TRUE
)

cat("\n=== Carrier rate by subtype (NEFH excluded) ===\n")
als_burden_df %>%
  group_by(cluster_clean) %>%
  summarise(
    n         = n(),
    n_carrier = sum(any_carrier_clean, na.rm=TRUE),
    pct       = round(mean(any_carrier_clean, na.rm=TRUE)*100, 1)
  ) %>% print()

chisq_clean <- chisq.test(
  table(als_burden_df$any_carrier_clean, als_burden_df$cluster_clean)
)
cat("Chi-square p:", round(chisq_clean$p.value, 4), "\n")

# NfL model without NEFH
m_als_nfl_clean <- lm(logNfL ~ log1p(total_als_rare_clean) +
                        AGE_AT_SYMPTOM_ONSET + SEX,
                      data = als_burden_df %>% filter(!is.na(logNfL)))
cat("\n=== ALS-gene burden (NEFH excluded) → NfL ===\n")
print(summary(m_als_nfl_clean)$coefficients)

# Gene-level carrier rates without NEFH (clean summary for paper)
cat("\n=== Per-gene carrier rates (NEFH excluded, top genes only) ===\n")
top_genes <- carrier_rates[carrier_rates >= 1 & names(carrier_rates) != "NEFH"]
print(round(top_genes, 1))


# ============================================================
# DISEASE DURATION / SURVIVAL ANALYSIS
# ============================================================

# Compute disease duration for patients with both values
survival_df <- df_nfl_uniq %>%
  mutate(Participant_ID = sub("-[0-9]+-T$", "", sample_id)) %>%
  left_join(ax_clean_df %>%
              mutate(Participant_ID = sub("-[0-9]+-T$", "", sample_id)) %>%
              dplyr::select(Participant_ID, cluster_clean),
            by = "Participant_ID") %>%
  filter(!is.na(cluster_clean)) %>%
  distinct(Participant_ID, .keep_all = TRUE)

# Check availability
cat("=== AGE_AT_DEATH availability ===\n")
cat("Total sporadic samples:", nrow(survival_df), "\n")
cat("With AGE_AT_DEATH:", sum(!is.na(survival_df$AGE_AT_DEATH)), "\n")
cat("With both onset and death:",
    sum(!is.na(survival_df$AGE_AT_DEATH) &
          !is.na(survival_df$AGE_AT_SYMPTOM_ONSET)), "\n")

# Compute disease duration
survival_df <- survival_df %>%
  mutate(disease_duration = AGE_AT_DEATH - AGE_AT_SYMPTOM_ONSET)

# Missingness by subtype
cat("\n=== AGE_AT_DEATH missingness by subtype ===\n")
survival_df %>%
  group_by(cluster_clean) %>%
  summarise(
    n            = n(),
    n_with_death = sum(!is.na(AGE_AT_DEATH)),
    pct_complete = round(mean(!is.na(AGE_AT_DEATH))*100, 1)
  ) %>% print()

# Chi-square: is missingness subtype-dependent?
chisq_surv <- chisq.test(
  table(!is.na(survival_df$AGE_AT_DEATH), survival_df$cluster_clean)
)
cat("Missingness chi-sq p:", round(chisq_surv$p.value, 4), "\n")

# ============================================================
# Analysis 1: Disease duration by subtype (deceased only)
# ============================================================

surv_complete <- survival_df %>%
  filter(!is.na(disease_duration))

cat("\n=== Disease duration by subtype (deceased patients only) ===\n")
surv_complete %>%
  group_by(cluster_clean) %>%
  summarise(
    n      = n(),
    mean   = round(mean(disease_duration), 2),
    median = round(median(disease_duration), 2),
    sd     = round(sd(disease_duration), 2)
  ) %>% print()

aov_dur <- aov(disease_duration ~ cluster_clean, data = surv_complete)
cat("ANOVA p:", summary(aov_dur)[[1]]$`Pr(>F)`[1], "\n")
kw_dur <- kruskal.test(disease_duration ~ cluster_clean, data = surv_complete)
cat("Kruskal-Wallis p:", round(kw_dur$p.value, 4), "\n")

if (summary(aov_dur)[[1]]$`Pr(>F)`[1] < 0.1) {
  cat("\nTukey HSD:\n")
  print(TukeyHSD(aov_dur))
  
  library(dunn.test)
  dunn.test(surv_complete$disease_duration,
            surv_complete$cluster_final, method="bh")
}

# ============================================================
# Analysis 2: Age at death adjusted for age at onset
# ============================================================

cat("\n=== Age at death ~ subtype + age at onset + sex ===\n")
m_death <- lm(AGE_AT_DEATH ~ cluster_clean + AGE_AT_SYMPTOM_ONSET + SEX,
              data = surv_complete)
print(summary(m_death)$coefficients)

# ============================================================
# Analysis 3: Disease duration ~ NfL (does NfL predict survival?)
# ============================================================

surv_nfl <- surv_complete %>% filter(!is.na(logNfL))
cat("\n=== Disease duration ~ NfL (deceased patients) n =",
    nrow(surv_nfl), "===\n")
m_nfl_dur <- lm(disease_duration ~ logNfL + AGE_AT_SYMPTOM_ONSET + SEX,
                data = surv_nfl)
print(summary(m_nfl_dur)$coefficients)

