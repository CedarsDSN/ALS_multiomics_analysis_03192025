# =====================================================================
# make_tables.R
# Builds Table 1, Table 2 and Supplementary Tables S1-S4 from pinned
# objects. Writes CSVs you can paste into Word.
#
# Supplementary tables are numbered in order of first citation in the
# manuscript (Nature Communications requirement):
#   S1  One-vs-rest LASSO classifier panels      (Results, CCT2 section)
#   S2  Germline tests                            (Results, germline section)
#   S3  Clinical endpoint screen, 23 endpoints    (Results, clinical screen)
#   S4  Axis composition (+ S4b gene list)        (Methods, axis definition)
#
# Values are pulled from saved objects wherever they exist. Where a value
# has to be typed (because no saved object holds it), the console prints
# "TYPED:" so you can check it against Table 3 and the Results text.
# Run from ~/Documents/DSN/ALS in a fresh R session.
# =====================================================================

suppressPackageStartupMessages({ library(dplyr) })

find1 <- function(f, dirs = c(".", "~/Documents/DSN/ALS", "~/Documents",
                              "~/Documents/DSN/ALS/proteomics",
                              "~/Documents/DSN/ALS/targetALS",
                              "~/Documents/DSN/ALS/clinical_endpoints")) {
  for (d in path.expand(dirs)) {
    p <- file.path(d, f); if (file.exists(p)) return(p)
    if (dir.exists(d)) {
      g <- list.files(d, pattern = paste0("^", f, "$"), recursive = TRUE, full.names = TRUE)
      if (length(g)) return(g[1])
    }
  }
  stop("not found: ", f)
}

b   <- readRDS(find1("backbone_FINAL.rds"))
pl  <- readRDS(find1("patient_level_FINAL.rds"))
m   <- as.data.frame(b$meta)
pat <- as.data.frame(pl$pat)

fmt_msd <- function(x, dp = 1) if (all(is.na(x))) NA_character_ else
  sprintf(paste0("%.", dp, "f (%.", dp, "f)"),
          mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_n   <- function(x, tot) sprintf("%d (%.0f%%)", x, 100 * x / tot)

## ---------------------------------------------------------------
## TABLE 1  Cohort characteristics
## ---------------------------------------------------------------
ppl <- m[!duplicated(m$Participant_ID), ]
ppl$grp <- factor(ppl$Mutation_Subtype,
                  levels = c("Control", "Sporadic", "C9orf72", "SOD1"))
m$grp   <- factor(m$Mutation_Subtype,
                  levels = c("Control", "Sporadic", "C9orf72", "SOD1"))

t1 <- lapply(levels(ppl$grp), function(g) {
  p <- ppl[ppl$grp == g, ]; n <- nrow(p)
  data.frame(
    Characteristic = c("Participants", "RNA-seq samples", "Female", "Male",
                       "Age at symptom onset, years, mean (SD)",
                       "Age at first collection, years, mean (SD)"),
    Value = c(as.character(n),
              as.character(sum(m$grp == g, na.rm = TRUE)),
              fmt_n(sum(p$SEX == "Female", na.rm = TRUE), n),
              fmt_n(sum(p$SEX == "Male",   na.rm = TRUE), n),
              fmt_msd(suppressWarnings(as.numeric(p$AGE_AT_SYMPTOM_ONSET))),
              if ("Age_at_First_PBMC_Collection" %in% names(p))
                fmt_msd(suppressWarnings(as.numeric(p$Age_at_First_PBMC_Collection)))
              else NA_character_),
    Group = g, stringsAsFactors = FALSE)
})
t1 <- do.call(rbind, t1) |> tidyr::pivot_wider(names_from = Group, values_from = Value)
cat("\n===== TABLE 1: cohort characteristics =====\n"); print(as.data.frame(t1))
write.csv(t1, "Table1_cohort.csv", row.names = FALSE)

## ---------------------------------------------------------------
## TABLE 2  Subtype characteristics (sporadic patients)
## ---------------------------------------------------------------
## CHANGED: match differentiation quality by sample (the earliest-timepoint
## sample that carries each patient's subtype label), not by participant,
## because dq is sample-level and match() on Participant_ID takes whichever
## sample happens to come first.
dq <- try(readRDS(find1("diff_quality_sensitivity_FINAL.rds"))$meta, silent = TRUE)
if (!inherits(dq, "try-error") && "diff_quality" %in% names(dq)) {
  sid <- intersect(c("sample_id", "Sample_ID"), names(dq))
  if (length(sid) && "sample_id" %in% names(pat)) {
    pat$diff_quality <- dq$diff_quality[match(pat$sample_id, dq[[sid[1]]])]
    cat("diff_quality matched by sample:", sum(!is.na(pat$diff_quality)), "of", nrow(pat), "\n")
  } else {
    pat$diff_quality <- dq$diff_quality[match(pat$Participant_ID, dq$Participant_ID)]
    cat("NOTE: no sample ID column in dq; diff_quality matched by participant\n")
  }
}

## CHANGED: axis rows at 2 dp so they match the Results text.
t2 <- pat %>% group_by(subtype) %>% summarise(
  `Patients, n (%)`               = fmt_n(n(), nrow(pat)),
  `Metabolic axis, mean (SD)`     = fmt_msd(AXIS_METABOLIC, 2),
  `Inflammation axis, mean (SD)`  = fmt_msd(AXIS_INFLAMMATION, 2),
  `Cell cycle axis, mean (SD)`    = fmt_msd(AXIS_CELL_CYCLE, 2),
  `CCT2, VST, mean (SD)`          = fmt_msd(CCT2, 2),
  `Age at onset, years, mean (SD)` = fmt_msd(AGE_AT_SYMPTOM_ONSET),
  `Female, n (%)`                 = fmt_n(sum(SEX == "Female", na.rm = TRUE), n()),
  `log10 NfL, mean (SD)`          = fmt_msd(logNfL),
  `ALSFRS-R baseline, mean (SD)`  = fmt_msd(ALSFRS_R_Baseline_Value),
  `ALSFRS-R slope, points/month, mean (SD)` = fmt_msd(ALSFRS_R_PROGRESSION_SLOPE),
  `Differentiation quality, mean (SD)` =
    if ("diff_quality" %in% names(pat)) fmt_msd(diff_quality) else NA_character_,
  .groups = "drop")
t2 <- as.data.frame(t(t2)); colnames(t2) <- t2[1, ]; t2 <- t2[-1, , drop = FALSE]
t2 <- cbind(Characteristic = rownames(t2), t2)
cat("\n===== TABLE 2: subtype characteristics =====\n"); print(t2, row.names = FALSE)
write.csv(t2, "Table2_subtypes.csv", row.names = FALSE)

## ---------------------------------------------------------------
## SUPPLEMENTARY TABLE S1  One-vs-rest LASSO classifier panels
## (cited first, in the Results CCT2 section)
## ---------------------------------------------------------------
## Finds the object saved by onevsrestLASSO.R (a list with res, summary, genes).
lasso_file <- NULL
cand <- list.files(c(".", path.expand("~/Documents/DSN/ALS")),
                   pattern = "(lasso|ovr|onevsrest).*\\.rds$",
                   ignore.case = TRUE, recursive = TRUE, full.names = TRUE)
for (f in unique(cand)) {
  o <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.list(o) && all(c("summary", "genes") %in% names(o))) { lasso_file <- f; break }
}
if (is.null(lasso_file)) {
  cat("\nNOTE: could not find the onevsrestLASSO.R output. Check line 149 of",
      "onevsrestLASSO.R for its file name and set lasso_file by hand.\n")
} else {
  lz <- readRDS(lasso_file)
  cat("\n===== S1 source:", lasso_file, "=====\n")
  str(lz[c("summary", "genes")], max.level = 2)
  
  summ <- as.data.frame(lz$summary)
  ## Gene panels: genes may be a named list (one vector per subtype) or a data frame.
  if (is.list(lz$genes) && !is.data.frame(lz$genes)) {
    panels <- data.frame(Subtype = names(lz$genes),
                         `Panel size` = sapply(lz$genes, length),
                         `Selected genes` = sapply(lz$genes, paste, collapse = ", "),
                         check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    panels <- as.data.frame(lz$genes)
  }
  write.csv(summ,   "TableS1_lasso_summary.csv", row.names = FALSE)
  write.csv(panels, "TableS1_lasso_panels.csv",  row.names = FALSE)
  cat("Wrote TableS1_lasso_summary.csv and TableS1_lasso_panels.csv.",
      "Check: panel sizes 12 to 16, nested AUC 0.936 to 0.970,",
      "CCT2 bootstrap selection 74.7% for S-ALS-3.\n")
}

## ---------------------------------------------------------------
## SUPPLEMENTARY TABLE S2  Germline tests
## (cited second, in the Results germline section)
## ---------------------------------------------------------------
gb <- readRDS(find1("genomic_burden_FINAL.rds"))
L  <- gb$landmarks
cat("\nAvailable landmarks in genomic_burden_FINAL.rds:\n"); print(names(L))

## Use the saved value if the landmark exists; otherwise use the typed value
## and say so on the console.
lm_or <- function(key, fmt, typed) {
  if (!is.null(L[[key]])) sprintf(fmt, L[[key]])
  else { cat("TYPED:", typed, "(landmark '", key, "' not found)\n", sep = ""); typed }
}

s3 <- data.frame(
  Analysis = c("Genome-wide rare damaging burden",
               "Genome-wide rare damaging burden",
               "Genome-wide rare LoF burden",
               "CCT2 rare damaging carriers", "CCT2 cis-eQTL",
               "CCT2 cis-eQTL in ALS GWAS", "TRiC complex carrier rate",
               "TRiC complex burden", "MYC target pathway burden",
               "MYC target carrier rate", "MYC target per-gene scan",
               "Genomic k-means vs subtype", "Blind genomic consensus, k = 2"),
  Test = c("ANOVA", "Kruskal-Wallis", "ANOVA", "carrier count",
           "linear model, 3,782 variants", "lookup, top 50 variants",
           "chi-square", "ANOVA", "Kruskal-Wallis", "chi-square",
           "Fisher, 128 genes", "adjusted Rand index", "adjusted Rand index"),
  Result = c(sprintf("P = %.2f", L$anova_rare_dam_p),
             lm_or("kw_rare_dam_p", "P = %.2f", "P = 0.05"),        # CHANGED: KW row added
             sprintf("P = %.2f", L$anova_rare_lof_p),
             "0 of 360 patients",
             lm_or("eqtl_min_q",   "minimum FDR q = %.3f", "minimum FDR q = 0.712"),
             "42 of 50 matched, minimum P = 0.116",
             lm_or("tric_carrier_p", "P = %.2f", "P = 0.19"),
             lm_or("tric_burden_p",  "P = %.2f", "P = 0.35"),
             lm_or("myc_burden_p",   "P = %.2f", "P = 0.85"),
             lm_or("myc_carrier_p",  "P = %.2f", "P = 0.73"),
             "top hit XPOT, q = 0.362",
             sprintf("%.3f", L$ari_kmeans),
             sprintf("%.3f (359/1 split)", L$ari_blind_k2)),
  stringsAsFactors = FALSE)
cat("TYPED: GWAS lookup (42 of 50, P = 0.116), XPOT q = 0.362, 0 of 360 carriers.",
    "Check these against Table 3 and the Results.\n")
write.csv(s3, "TableS2_genomic_tests.csv", row.names = FALSE)

## ---------------------------------------------------------------
## SUPPLEMENTARY TABLE S3  Clinical endpoint screen
## (cited third, in the Results clinical screen section)
## ---------------------------------------------------------------
sc <- as.data.frame(readRDS(find1("clinical_screen_combined_FINAL.rds")))
cat("\nScreen rows before filtering:", nrow(sc), "\n")
print(sc[, intersect(c("endpoint", "n", "note"), names(sc))])

## CHANGED: keep only the 23 modeled endpoints (the text says 25 screened,
## 2 with fewer than 20 evaluable patients not modeled).
sc <- sc[!is.na(sc$p_S3) & sc$n >= 20, ]
if (nrow(sc) != 23) warning("Expected 23 modeled endpoints, found ", nrow(sc),
                            ". Check the filter before using this table.")
cat("Minimum subtype FDR q:", round(min(sc$p_S3_FDR_all, na.rm = TRUE), 3),
    "(the Fig. 8 legend says 0.12)\n")

s2 <- sc[order(sc$p_S3), c("endpoint", "n", "beta_S3", "p_S3", "p_S3_FDR_all",
                           "beta_CCT2_uni", "p_CCT2_uni", "p_CCT2_FDR_all")]
names(s2) <- c("Endpoint", "n", "S-ALS-3 vs S-ALS-2 beta", "Subtype P", "Subtype FDR q",
               "CCT2 beta (per SD)", "CCT2 P", "CCT2 FDR q")
write.csv(s2, "TableS3_clinical_screen.csv", row.names = FALSE)

## ---------------------------------------------------------------
## SUPPLEMENTARY TABLE S4  Axis composition
## (cited last, in the Methods axis definition)
## ---------------------------------------------------------------
s1 <- data.frame(
  Axis = rep(c("Metabolic", "Cell cycle", "Inflammation"), c(3, 4, 3)),
  `Hallmark gene set` = c(
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_FATTY_ACID_METABOLISM",
    "HALLMARK_GLYCOLYSIS", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1", "HALLMARK_DNA_REPAIR",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_COMPLEMENT"),
  `Genes in axis` = rep(sapply(b$axis_sets[c("AXIS_METABOLIC", "AXIS_CELL_CYCLE",
                                             "AXIS_INFLAMMATION")], length), c(3, 4, 3)),
  check.names = FALSE)
write.csv(s1, "TableS4_axis_composition.csv", row.names = FALSE)

## Optional companion: the actual genes in each axis, one row per gene.
sl <- stack(b$axis_sets[c("AXIS_METABOLIC", "AXIS_CELL_CYCLE", "AXIS_INFLAMMATION")])
names(sl) <- c("Gene", "Axis")
write.csv(sl, "TableS4b_axis_genes.csv", row.names = FALSE)
cat("\nAxis gene counts:", paste(names(table(sl$Axis)), table(sl$Axis), collapse = "; "),
    "(expect 512 / 610 / 538)\n")

cat("\nWrote Table1_cohort.csv, Table2_subtypes.csv,\n",
    "TableS1_lasso_summary.csv + TableS1_lasso_panels.csv (if the LASSO file was found),\n",
    "TableS2_genomic_tests.csv, TableS3_clinical_screen.csv,\n",
    "TableS4_axis_composition.csv + TableS4b_axis_genes.csv\n")