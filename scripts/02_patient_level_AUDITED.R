# =====================================================================
# 02_patient_level.R   (AUDITED 2026-08-22)
# Backbone -> patient-level table (391 sporadic patients).
# Dedup: earliest timepoint per participant (slice_min on the numeric
# visit id in sample_id; with_ties = FALSE, deterministic given the
# backbone's row order). Consumes ax_round1 / cluster_r1, which are the
# CANONICAL labels (round 1; see audited script 01 header).
#
# AUDIT ADDITIONS:
#   - setwd pinned; NfL file is REQUIRED (previously an if-exists guard
#     could silently skip the merge and save a table without logNfL).
#   - Non-destructive verification against the saved
#     patient_level_FINAL.rds before overwriting.
#
# Landmarks: 391 patients | subtype 85/204/102 | AUC 0.895 (CI ~0.859-0.928,
# seed 42) | headline coefs: NfL 0.031/-0.001, slope -0.131/-0.098,
# baseline 0.240/0.261 | subtype agreement with saved table = 1.0
#
# ---- original header ----
# Switch to earliest-timepoint deduplication
#
# Your convention (PRSetTesting.R):
#   ax_clean_dedup <- ax_clean_df %>%
#     group_by(Participant_ID) %>%
#     slice_min(timepoint, n = 1, with_ties = FALSE)
#
# My figure scripts averaged across a patient's samples instead.
# This rebuilds the patient-level table your way, saves it once, and
# reports whether any figure number moves.
#
# Output: patient_level_FINAL.rds  (use this everywhere from now on)
# =====================================================================

suppressPackageStartupMessages({ library(dplyr) })
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::mutate,
                               dplyr::slice, base::setdiff, base::intersect,
                               base::union, base::unname, .quiet = TRUE)
}
base::set.seed(42)
base::setwd("~/Documents/DSN/ALS")

b    <- base::readRDS("backbone_FINAL.rds")
d    <- b$ax_round1
meta <- b$meta
E    <- b$expr_symbol
AX   <- c("AXIS_METABOLIC", "AXIS_INFLAMMATION", "AXIS_CELL_CYCLE")

if (base::is.null(d$Participant_ID))
  d$Participant_ID <- meta$Participant_ID[base::match(d$sample_id, meta$sample_id)]
d$subtype   <- d$cluster_r1
d$CCT2      <- E["CCT2", d$sample_id]
d$timepoint <- base::as.integer(base::sub(".*-(\\d+)-T$", "\\1", d$sample_id))

base::cat("samples: ", base::nrow(d),
          "   patients: ", base::length(base::unique(d$Participant_ID)),
          "   with repeats: ",
          base::sum(base::table(d$Participant_ID) > 1), "\n", sep = "")

# ---- your convention: earliest timepoint -----------------------------
pat <- d %>%
  dplyr::group_by(Participant_ID) %>%
  dplyr::slice_min(timepoint, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  base::as.data.frame()

base::cat("after earliest-timepoint dedup: ", base::nrow(pat), " (expect 391)\n", sep = "")
base::cat("duplicate IDs: ", base::sum(base::duplicated(pat$Participant_ID)), "\n", sep = "")
base::print(base::table(pat$subtype))

# ---- what my averaging approach gave, for comparison -----------------
avg <- stats::aggregate(base::cbind(CCT2) ~ Participant_ID, data = d, FUN = base::mean)
lb  <- base::tapply(base::as.character(d$subtype), d$Participant_ID, function(x) {
  t <- base::table(x); if (base::length(t) > 1) NA_character_ else base::names(t)[1] })
avg$subtype <- base::factor(lb[avg$Participant_ID], levels = base::paste0("S-ALS-", 1:3))

base::cat("\naveraging approach dropped ", base::sum(base::is.na(avg$subtype)),
          " ambiguous patients; earliest-timepoint keeps all\n", sep = "")

cmp <- base::merge(pat[, base::c("Participant_ID", "subtype", "CCT2")],
                   avg, by = "Participant_ID", suffixes = base::c("_first", "_mean"))
base::cat("subtype label disagreements: ",
          base::sum(base::as.character(cmp$subtype_first) !=
                      base::as.character(cmp$subtype_mean), na.rm = TRUE), "\n", sep = "")
base::cat("CCT2 correlation between methods: ",
          base::round(stats::cor(cmp$CCT2_first, cmp$CCT2_mean), 5), "\n", sep = "")

# ---- attach clinical -------------------------------------------------
cc <- base::c("Participant_ID", "SEX", "AGE_AT_SYMPTOM_ONSET", "NYGC_CGND_ID",
              "ALSFRS_R_PROGRESSION_SLOPE", "ALSFRS_R_Baseline_Value")
cc <- cc[cc %in% base::colnames(meta)]
pat <- base::merge(pat,
                   meta[!base::duplicated(meta$Participant_ID), cc, drop = FALSE],
                   by = "Participant_ID", all.x = TRUE)

nfl_path <- "answerals/NfL/NfL_Seimens_Atellica_Immunoassay_Proteomics_2022_07.csv"
if (!base::file.exists(nfl_path))
  base::stop("NfL file not found at ", nfl_path,
             " -- required; a patient table without logNfL breaks script 05.")
{
  nf <- utils::read.csv(nfl_path, stringsAsFactors = FALSE)
  tc <- base::grep("Timepoint", base::colnames(nf), value = TRUE)[1]
  if (!base::is.na(tc)) nf <- nf[nf[[tc]] == 1, ]
  nf <- nf[, base::c("SubjID", "Result.Value")]
  base::colnames(nf) <- base::c("SubjID", "NfL")
  nf$NfL <- base::as.numeric(base::gsub(",", "", nf$NfL))
  nf <- nf[!base::is.na(nf$NfL) & nf$NfL > 0 & !base::duplicated(nf$SubjID), ]
  nf$logNfL <- base::log10(nf$NfL)
  pat$SubjID <- base::sub("^CASE-", "", pat$Participant_ID)
  pat <- base::merge(pat, nf[, base::c("SubjID", "logNfL")], by = "SubjID", all.x = TRUE)
}
pat$CCT2z <- base::as.numeric(base::scale(pat$CCT2))
pat$SEX   <- base::factor(pat$SEX)

# ---- rerun the headline numbers --------------------------------------
base::cat("\n===== FIGURE 5: CCT2 classifier =====\n")
auc_fn <- function(s, p) {
  r <- base::rank(s); n1 <- base::sum(p); n0 <- base::sum(!p)
  (base::sum(r[p]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
sb  <- pat[pat$subtype %in% base::c("S-ALS-2", "S-ALS-3"), ]
sb$subtype <- base::droplevels(sb$subtype)
ps  <- sb$subtype == "S-ALS-3"
auc <- auc_fn(sb$CCT2, ps)
i3 <- base::which(ps); i2 <- base::which(!ps)
bt <- base::replicate(1000, {
  s <- base::c(base::sample(i3, base::length(i3), TRUE),
               base::sample(i2, base::length(i2), TRUE))
  auc_fn(sb$CCT2[s], ps[s]) })
ci <- stats::quantile(bt, base::c(0.025, 0.975))
base::cat("AUC = ", base::round(auc, 4), "  95% CI ", base::round(ci[1], 4),
          " to ", base::round(ci[2], 4), "   (was 0.895, 0.859-0.928)\n", sep = "")
base::cat("n = ", base::sum(ps), " S-ALS-3 vs ", base::sum(!ps), " S-ALS-2\n", sep = "")
base::cat("above S-ALS-2 max: ",
          base::round(100 * base::mean(sb$CCT2[ps] > base::max(sb$CCT2[!ps])), 1),
          "%   (was 28.7%)\n", sep = "")

base::cat("\n===== FIGURE 8: clinical =====\n")
fit2 <- function(y, lab) {
  df <- pat[stats::complete.cases(pat[, base::c(y, "CCT2z", "SEX",
                                                "AGE_AT_SYMPTOM_ONSET", "subtype")]), ]
  if (base::nrow(df) < 30) { base::cat(lab, ": too few\n"); return(NULL) }
  f1 <- stats::lm(stats::as.formula(base::paste(y, "~ CCT2z + SEX + AGE_AT_SYMPTOM_ONSET")), df)
  f2 <- stats::lm(stats::as.formula(base::paste(y, "~ CCT2z + subtype + SEX + AGE_AT_SYMPTOM_ONSET")), df)
  g <- function(f, m) base::data.frame(endpoint = lab, model = m, n = stats::nobs(f),
                                       beta = base::summary(f)$coefficients["CCT2z", 1],
                                       p = base::summary(f)$coefficients["CCT2z", 4],
                                       lo = stats::confint(f, "CCT2z")[1], hi = stats::confint(f, "CCT2z")[2])
  base::rbind(g(f1, "CCT2 alone"), g(f2, "CCT2 + subtype"))
}
res <- base::do.call(base::rbind, base::list(
  fit2("logNfL", "log10 NfL"),
  fit2("ALSFRS_R_PROGRESSION_SLOPE", "ALSFRS-R slope"),
  fit2("ALSFRS_R_Baseline_Value", "ALSFRS-R baseline")))
out <- base::data.frame(endpoint = res$endpoint, model = res$model, n = res$n,
                        beta = base::round(res$beta, 4),
                        p = base::signif(res$p, 3), stringsAsFactors = FALSE)
base::print(out, row.names = FALSE)
base::cat("previously: NfL 0.031/-0.001, slope -0.131/-0.098, baseline 0.240/0.261\n")

out_obj <- base::list(pat = pat, auc = auc, ci = ci, boot = bt,
                      coefs = res, dedup = "earliest timepoint")

base::cat("\nverifying against saved patient_level_FINAL.rds ...\n")
if (base::file.exists("patient_level_FINAL.rds")) {
  oldp <- base::readRDS("patient_level_FINAL.rds")
  oldpat <- if (base::is.data.frame(oldp)) oldp else oldp$pat
  v <- base::merge(pat[, base::c("Participant_ID", "subtype")],
                   oldpat[, base::c("Participant_ID", "subtype")],
                   by = "Participant_ID", suffixes = base::c("_new", "_old"))
  agree <- base::mean(base::as.character(v$subtype_new) ==
                      base::as.character(v$subtype_old))
  base::cat("patients matched:", base::nrow(v), "| subtype agreement:",
            base::round(agree, 4), "\n")
  nfl_ok <- "logNfL" %in% base::names(pat) &&
            base::sum(!base::is.na(pat$logNfL)) >=
            base::sum(!base::is.na(oldpat$logNfL))
  if (agree == 1 && base::nrow(v) == 391 && nfl_ok) {
    base::saveRDS(out_obj, "patient_level_FINAL.rds")
    base::cat("EXACT MATCH: patient_level_FINAL.rds refreshed. Script VERIFIED.\n")
  } else {
    base::saveRDS(out_obj, "patient_level_AUDITCHECK.rds")
    base::cat("MISMATCH (agreement ", agree, ", n ", base::nrow(v),
              ", nfl_ok ", nfl_ok, "): saved as patient_level_AUDITCHECK.rds;",
              " existing file untouched. Report back.\n", sep = "")
  }
} else {
  base::saveRDS(out_obj, "patient_level_FINAL.rds")
  base::cat("no prior patient table found; saved patient_level_FINAL.rds.\n")
}