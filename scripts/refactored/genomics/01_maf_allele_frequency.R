# ============================================================
# Refactored entrypoint for original script: MAF_allele_frequency.R
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
if (isTRUE(ALS_ARGS$help)) { als_print_help("01_maf_allele_frequency"); quit(save = "no", status = 0) }
CFG <- als_load_config(ALS_ARGS$config)
ALSPROJ_BASE_DIR <- als_choose(ALS_ARGS$base_dir, CFG, "project.base_dir", getwd())
ALSPROJ_OUTDIR <- als_choose(ALS_ARGS$outdir, CFG, "outputs.root", file.path(getwd(), "results", "01_maf_allele_frequency"))
ALSPROJ_SCRIPT_ID <- "01_maf_allele_frequency"
als_init(ALSPROJ_SCRIPT_ID, ALSPROJ_BASE_DIR, ALSPROJ_OUTDIR, ALS_ARGS$redirect_outputs)
als_register_output_wrappers()
als_load_objects(ALS_ARGS$objects)
on.exit(als_write_session_info(ALSPROJ_SCRIPT_ID), add = TRUE)
# ---- End bootstrap ----

# ============================================================
# MAF DISTRIBUTION — from GATK joint VCF AF field
# ============================================================

library(data.table)
library(ggplot2)

cat("Reading AF data...\n")
af_df <- fread(
  "all_chr_AF.txt",
  header = FALSE,
  col.names = c("CHROM","POS","REF","ALT","AF","AN"),
  colClasses = c("character","integer","character",
                 "character","numeric","integer")
)

cat("Raw variants loaded:", nrow(af_df), "\n")

# Remove sites where AF failed to parse (multiallelic edge cases)
af_df <- af_df[!is.na(AF) & AF > 0 & AF <= 1]

# Compute MAF = min(AF, 1-AF)
# AF from GATK is ALT allele frequency; MAF folds it
af_df[, MAF := pmin(AF, 1 - AF)]

# Filter: remove invariant (shouldn't exist post-split but just in case)
af_df <- af_df[MAF > 0]

cat("Biallelic variants with MAF > 0:", nrow(af_df), "\n")

# AN distribution — shows genotyping completeness
cat("\nAN (alleles called) distribution:\n")
print(quantile(af_df$AN, probs=c(0,.1,.25,.5,.75,.9,1), na.rm=TRUE))
max_AN <- max(af_df$AN, na.rm=TRUE)
cat("Max possible AN (2 x n_samples):", max_AN, "\n")
cat("Estimated n_samples:", max_AN/2, "\n")

# MAF quantiles
cat("\nMAF distribution:\n")
print(quantile(af_df$MAF,
               probs=c(0,.01,.05,.10,.25,.50,.75,.90,.95,.99,1)))

# Variant counts by tier
tiers <- list(
  "Singleton/very rare (MAF<0.001)"   = af_df$MAF < 0.001,
  "Rare (MAF 0.001–0.01)"             = af_df$MAF >= 0.001 & af_df$MAF < 0.01,
  "Low frequency (MAF 0.01–0.05)"     = af_df$MAF >= 0.01  & af_df$MAF < 0.05,
  "Common (MAF ≥ 0.05)"               = af_df$MAF >= 0.05
)
cat("\nVariant counts by frequency tier:\n")
for (nm in names(tiers)) {
  n <- sum(tiers[[nm]])
  pct <- mean(tiers[[nm]]) * 100
  cat(sprintf("  %-40s %8s  (%.1f%%)\n",
              nm, format(n, big.mark=","), pct))
}

# ============================================================
# PLOT 1: Full SFS — log scale x-axis
# Classic site frequency spectrum shape
# ============================================================

p1 <- ggplot(af_df, aes(x = MAF)) +
  geom_histogram(bins = 120, fill = "#2C5F8A", color = NA, alpha = 0.85) +
  scale_x_log10(
    breaks = c(0.0001, 0.001, 0.01, 0.05, 0.10, 0.25, 0.50),
    labels = c("0.0001","0.001","0.01","0.05","0.10","0.25","0.50"),
    limits = c(0.00005, 0.5)
  ) +
  geom_vline(xintercept = 0.01, color = "#E05A3A",
             linewidth = 0.9, linetype = "dashed") +
  annotate("text", x = 0.007, y = Inf,
           label = "MAF < 0.01\nrare variant\nthreshold",
           hjust = 1, vjust = 1.4, color = "#E05A3A", size = 3.5) +
  labs(
    title = "Site Frequency Spectrum — AnswerALS Joint VCF",
    subtitle = paste0(
      format(nrow(af_df), big.mark=","),
      " biallelic variants | all autosomes + chrX | n ≈ ",
      round(max_AN/2), " samples"
    ),
    x = "Minor Allele Frequency (log scale)",
    y = "Number of variants"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

# ============================================================
# PLOT 2: Rare variant zoom (MAF < 0.05)
# Validates the MAF < 0.01 cutoff used in burden analysis
# ============================================================

p2 <- ggplot(af_df[MAF < 0.05], aes(x = MAF)) +
  geom_histogram(bins = 80, fill = "#2C5F8A", color = NA, alpha = 0.85) +
  geom_vline(xintercept = 0.01, color = "#E05A3A",
             linewidth = 0.9, linetype = "dashed") +
  annotate("text", x = 0.012, y = Inf,
           label = "MAF = 0.01\nburden cutoff",
           hjust = 0, vjust = 1.4, color = "#E05A3A", size = 3.5) +
  labs(
    title = "Rare Variant Frequency Spectrum (MAF < 0.05)",
    subtitle = paste0(
      format(sum(af_df$MAF < 0.05), big.mark=","),
      " variants shown"
    ),
    x = "Minor Allele Frequency",
    y = "Number of variants"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

# ============================================================
# PLOT 3: Tier composition — for slide summary panel
# ============================================================

tier_df <- data.frame(
  tier  = factor(
    c("Singleton/very rare\n(MAF<0.001)",
      "Rare\n(0.001–0.01)",
      "Low frequency\n(0.01–0.05)",
      "Common\n(MAF≥0.05)"),
    levels = c("Common\n(MAF≥0.05)",
               "Low frequency\n(0.01–0.05)",
               "Rare\n(0.001–0.01)",
               "Singleton/very rare\n(MAF<0.001)")
  ),
  count = c(sum(af_df$MAF <  0.001),
            sum(af_df$MAF >= 0.001 & af_df$MAF < 0.01),
            sum(af_df$MAF >= 0.01  & af_df$MAF < 0.05),
            sum(af_df$MAF >= 0.05)),
  fill  = c("#1A3A5C","#2C5F8A","#7BAFD4","#C8DFF0")
)
tier_df$pct <- tier_df$count / sum(tier_df$count) * 100

p3 <- ggplot(tier_df, aes(x = 1, y = pct, fill = tier)) +
  geom_bar(stat = "identity", width = 0.55) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_stack(vjust = 0.5),
            size = 4.2, color = "white", fontface = "bold") +
  scale_fill_manual(values = setNames(tier_df$fill, tier_df$tier)) +
  coord_flip() +
  labs(
    title = "Variant Composition by Frequency Tier",
    x = NULL, y = "Percentage of all variants", fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10)
  ) +
  guides(fill = guide_legend(reverse = TRUE))

# Save
ggsave("MAF_full_SFS.pdf",          p1, width = 8,  height = 5)
ggsave("MAF_rare_zoom.pdf",         p2, width = 8,  height = 5)
ggsave("MAF_tier_composition.pdf",  p3, width = 9,  height = 3.5)

cat("\nAll plots saved.\n")
