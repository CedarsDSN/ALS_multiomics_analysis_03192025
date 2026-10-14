################################################################################
## 06_clinical_screen.R   (AUDITED 2026-08-22)
##
## ONE runner for the full clinical endpoint screen. Previously a 7-step
## console workflow (preamble + hand-sourcing five wave scripts + a manual
## combine_screen() call); the combined whole-screen FDR that backs the
## manuscript's "uniformly null" claim only existed if the last call was
## remembered. This file does all of it top to bottom.
##
## Inputs : backbone_FINAL.rds, patient_level_FINAL.rds, portal clinical
##          CSVs under CLIN_DIR, and the wave monolith AllClinicalVariableTesting_newlabels.R.
## Output : clinical_screen_combined_FINAL.rds / .csv (verification-gated).
##
## Landmarks: anchor 306 patients (S-ALS-2 = 204, S-ALS-3 = 102) |
## ~25 linear endpoints | 0 significant at whole-screen FDR < 0.05 |
## combined table must match the saved rds (endpoint set, n, p_S3).
##
## Audit fixes: empty-table diagnostic (cluster_final -> subtype); explicit
## covariate-column guard; the wave monolith sourced by the runner (wd pinned
## to OUT_DIR so its mixed relative saveRDS paths land in one place);
## combine_screen() called automatically; verification gate added.
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(survival)
})

## ---------------------------------------------------------------------------
## CONFIG
## ---------------------------------------------------------------------------
CLIN_DIR <- "~/Documents/DSN/ALS/metadata/clinical"
OUT_DIR  <- "~/Documents/DSN/ALS/clinical_endpoints"
if (!dir.exists(path.expand(OUT_DIR))) dir.create(path.expand(OUT_DIR), recursive = TRUE)

DAYS_PER_YEAR <- 365.25
MIN_VISITS    <- 3
MIN_SPAN_DAYS <- 90
COVARS <- "AGE_AT_SYMPTOM_ONSET + SEX + SV1 + SV2 + SV3 + SV4 + SV5 + SV6"

## ---------------------------------------------------------------------------
## HELPERS
## ---------------------------------------------------------------------------
to_pid <- function(x) sub("-[0-9]+-T$", "", x)

## max of the valid (finite) values, NA if none. Used for VC trials, grip
## trials and HHD trials. This was missing from the scripts as circulated.
maxvalid <- function(...) {
  v <- suppressWarnings(as.numeric(c(...)))
  v <- v[is.finite(v)]
  if (length(v) == 0) NA_real_ else max(v)
}

slope_per_pid <- function(dt, id_col, day_col, val_col) {
  dt <- as.data.table(dt)
  dt <- dt[is.finite(dt[[day_col]]) & is.finite(dt[[val_col]])]
  dt[, .(
    n_visits = .N,
    span_d   = max(get(day_col)) - min(get(day_col)),
    slope_yr = if (uniqueN(get(day_col)) >= MIN_VISITS &&
                   (max(get(day_col)) - min(get(day_col))) >= MIN_SPAN_DAYS) {
      as.numeric(coef(lm(get(val_col) ~ get(day_col)))[2]) * DAYS_PER_YEAR
    } else NA_real_,
    baseline = get(val_col)[which.min(get(day_col))]
  ), by = id_col]
}

fit_linear_endpoint <- function(dat, endpoint, label) {
  dat <- as.data.table(dat)
  dat <- dat[is.finite(dat[[endpoint]])]
  n <- nrow(dat)
  if (n < 20) {
    return(data.table(endpoint = label, n = n,
                      beta_S3 = NA_real_, p_S3 = NA_real_,
                      beta_CCT2_uni = NA_real_, p_CCT2_uni = NA_real_,
                      beta_CCT2_adj = NA_real_, p_CCT2_adj = NA_real_,
                      note = "n<20, not modelled"))
  }
  f_sub <- as.formula(sprintf("%s ~ is_S3 + %s", endpoint, COVARS))
  f_uni <- as.formula(sprintf("%s ~ CCT2_z + %s", endpoint, COVARS))
  f_adj <- as.formula(sprintf("%s ~ CCT2_z + is_S3 + %s", endpoint, COVARS))
  m_sub <- lm(f_sub, dat); m_uni <- lm(f_uni, dat); m_adj <- lm(f_adj, dat)
  cs <- function(m, term) {
    s <- summary(m)$coefficients
    if (term %in% rownames(s)) s[term, c("Estimate", "Pr(>|t|)")] else c(NA, NA)
  }
  b_s3 <- cs(m_sub, "is_S3"); b_un <- cs(m_uni, "CCT2_z"); b_aj <- cs(m_adj, "CCT2_z")
  data.table(endpoint = label, n = n,
             beta_S3 = b_s3[1], p_S3 = b_s3[2],
             beta_CCT2_uni = b_un[1], p_CCT2_uni = b_un[2],
             beta_CCT2_adj = b_aj[1], p_CCT2_adj = b_aj[2],
             note = "+CCT2 collapse vs stable = state-marker vs graded")
}

## ---------------------------------------------------------------------------
## ANCHOR  from the rebuilt pipeline
## ---------------------------------------------------------------------------
bb <- readRDS("backbone_FINAL.rds")
pl <- readRDS("patient_level_FINAL.rds")

pat  <- as.data.table(pl$pat)          # earliest timepoint, 391 patients
meta <- as.data.table(bb$meta)
if (!"sample_id" %in% names(meta)) meta[, sample_id := Sample_ID]

cat("subtype counts, all sporadic patients:\n")
print(table(pat$subtype))

## covariates: SV1..SV6 and SEX from the metadata, one row per participant
sv_cols <- paste0("SV", 1:6)
stopifnot(all(c(sv_cols, "SEX", "AGE_AT_SYMPTOM_ONSET") %in% names(meta)))
cov_tab <- unique(meta[, c("Participant_ID", "SEX", "AGE_AT_SYMPTOM_ONSET",
                           sv_cols), with = FALSE], by = "Participant_ID")

lab_col <- base::intersect(c("cluster_final", "subtype", "cluster_r1"), names(pat))[1]
cat("using subtype column:", lab_col, "\n")
print(table(pat[[lab_col]]))

anchor <- merge(
  pat[, .(Participant_ID = as.character(Participant_ID),
          sample_id,
          subtype = as.character(get(lab_col)),
          CCT2)],
  cov_tab[, Participant_ID := as.character(Participant_ID)],
  by = "Participant_ID", all.x = TRUE)

anchor[, group := "Sporadic"]

## S-ALS-2 as reference, S-ALS-3 as the group of interest, matching the original
anchor <- anchor[subtype %in% c("S-ALS-2", "S-ALS-3")]
anchor[, subtype := factor(subtype, levels = c("S-ALS-2", "S-ALS-3"))]
anchor[, is_S3   := as.integer(subtype == "S-ALS-3")]
anchor[, CCT2_z  := as.numeric(scale(CCT2))]

cat(sprintf("\nAnchor: %d patients (S-ALS-2 = %d, S-ALS-3 = %d)\n",
            nrow(anchor), sum(anchor$subtype == "S-ALS-2"),
            sum(anchor$subtype == "S-ALS-3")))
cat("  previously: 226 and 74\n")

miss <- anchor[!is.finite(AGE_AT_SYMPTOM_ONSET) | is.na(SEX) | !is.finite(SV1), .N]
if (miss > 0) cat(sprintf("  WARNING: %d anchor rows missing a covariate\n", miss))

## sanity: the covariate model must be estimable
cat("\ncovariate completeness:\n")
for (v in c("AGE_AT_SYMPTOM_ONSET", "SEX", sv_cols))
  cat(sprintf("  %-22s %d of %d\n", v, sum(!is.na(anchor[[v]])), nrow(anchor)))

## ---------------------------------------------------------------------------
## After running v1 to v5, combine and apply FDR across the WHOLE screen.
## The original scripts correct within each wave of 5, which is more lenient
## than the manuscript's stated "correction across the linear screen".
## Run this block at the end.
## ---------------------------------------------------------------------------
combine_screen <- function() {
  waves <- list(get0("linear_results"), get0("linear_v2"), get0("linear_v3"),
                get0("linear_v4"), get0("linear_v5"))
  waves <- waves[!vapply(waves, is.null, logical(1))]
  if (length(waves) == 0) { cat("no wave results in session\n"); return(invisible(NULL)) }
  all <- rbindlist(waves, use.names = TRUE, fill = TRUE)
  all[, p_S3_FDR_wave := p_S3_FDR]
  all[, p_S3_FDR_all  := p.adjust(p_S3, method = "BH")]
  all[, p_CCT2_FDR_all := p.adjust(p_CCT2_uni, method = "BH")]
  setorder(all, p_S3)
  cat("\n=========== COMBINED LINEAR SCREEN, FDR across all tests ===========\n")
  print(all[, .(endpoint, n,
                beta_S3 = round(beta_S3, 3), p_S3 = signif(p_S3, 3),
                FDR_wave = signif(p_S3_FDR_wave, 3),
                FDR_all  = signif(p_S3_FDR_all, 3))])
  cat(sprintf("\n%d linear endpoints tested. Significant at FDR_all < 0.05: %d\n",
              nrow(all), sum(all$p_S3_FDR_all < 0.05, na.rm = TRUE)))
  ## saving is handled by the verification gate in the runner
  invisible(all)
}

## ---------------------------------------------------------------------------
## RUN THE WAVES, THEN COMBINE (no manual steps)
## The five waves live as sections of ONE monolith
## (AllClinicalVariableTesting_newlabels.R), which expects `anchor` and the
## helpers defined above and creates linear_results / linear_v2..v5 in the
## session. Its own intermediate saveRDS calls use mixed relative paths, so
## the working directory is set to OUT_DIR first to corral them.
## ---------------------------------------------------------------------------
WAVES_FILE <- "~/Documents/DSN/ALS/AllClinicalVariableTesting_newlabels.R"
if (!file.exists(path.expand(WAVES_FILE)))
  stop("Waves file not found: ", WAVES_FILE, " -- edit WAVES_FILE.")
owd <- getwd(); setwd(path.expand(OUT_DIR))
cat("\n----- sourcing", basename(WAVES_FILE), "(waves v1-v5) -----\n")
source(path.expand(WAVES_FILE))
setwd(owd)

all <- combine_screen()

## ---------------------------------------------------------------------------
## VERIFICATION GATE vs the saved combined screen
## ---------------------------------------------------------------------------
final_rds <- file.path(OUT_DIR, "clinical_screen_combined_FINAL.rds")
if (file.exists(path.expand(final_rds)) && !is.null(all)) {
  oldall <- as.data.table(readRDS(path.expand(final_rds)))
  ## order-and-compare (duplicate-label safe: "Creatinine baseline" appears
  ## in waves 1 and 4, so a merge on endpoint cross-joins)
  setorder(all, endpoint, n); setorder(oldall, endpoint, n)
  same_set <- base::setequal(all$endpoint, oldall$endpoint)
  ok <- isTRUE(all.equal(all[, .(endpoint, n, p_S3, beta_S3)],
                         oldall[, .(endpoint, n, p_S3, beta_S3)]))
  cat(sprintf("\nverification: endpoint set %s | rows %d vs %d | tables equal %s\n",
              same_set, nrow(all), nrow(oldall), ok))
  if (same_set && ok) {
    cat("EXACT MATCH: combined screen reproduces the saved results. Script VERIFIED.\n")
  } else {
    fwrite(all, file.path(OUT_DIR, "clinical_screen_AUDITCHECK.csv"))
    saveRDS(all, file.path(OUT_DIR, "clinical_screen_AUDITCHECK.rds"))
    cat("MISMATCH: saved as clinical_screen_AUDITCHECK; report back before",
        "the FINAL files are replaced.\n")
  }
}
sessionInfo()
################################################################################