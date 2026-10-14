# =====================================================================
# 05_downstream_ladders.R   (AUDITED 2026-08-22, SV-patched)
#
# AUDIT CHANGES from the circulated v2:
#   - SV1..SV6 (742-sample anchor SVA, via backbone meta) added to every
#     covariate-adjusted rung of L1-L4/L8 and to all six BIC models,
#     honouring the Methods text. Attached to the patient frame by
#     sample_id with a hard 391/391 match check.
#   - PRS ladders (L5-L7) intentionally keep their original spec
#     (SEX only), matching the PRSet analysis; the sex-adjustment
#     sensitivity for PRS is tracked separately.
#   - Within-subtype rung minimum raised 20 -> 40 (nine parameters).
#   - Verification gate: the L1/L3 covariate and +subtype rungs must
#     reproduce the verified 2026-08-19 console run (p-values to 1e-6:
#     NfL 0.045169 -> 0.751872; ALSFRS 0.0047696 -> 0.0393746;
#     n 297 / 244) before FINAL files are written.
#   - setwd pinned; conflicted preferences declared.
#
# ---- original header ----
# DOWNSTREAM ANALYSES, rerun on final labels   (complete, v2)
#
#   L1  CCT2 -> logNfL                    L5  each PRS -> CCT2
#   L2  AXIS_METABOLIC -> logNfL          L6  PRSet 3x3 matrix, clean subset
#   L3  CCT2 -> ALSFRS-R slope            L7  PRS -> axis under glial + DQ
#   L4  AXIS_METABOLIC -> ALSFRS-R slope  L8  CCT2 -> disease duration
#
# Each ladder: univariate -> +age+sex -> +subtype -> within-subtype,
# plus head-to-head and BIC comparison per outcome.
#
# Run diff_quality_sensitivity_final.R FIRST (supplies composites).
#
# Fixed in v2: CGND column arrives suffixed .x / .y after merging, so
# the PRS join now resolves whichever copy is populated.
#
# Output: downstream_analyses_FINAL.rds, downstream_ladders_FINAL.csv
# =====================================================================

suppressPackageStartupMessages({ library(dplyr); library(data.table) })
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::mutate,
                               dplyr::slice, base::setdiff, base::intersect,
                               base::union, base::unname, base::setequal,
                               .quiet = TRUE)
}
base::set.seed(42)
base::setwd("~/Documents/DSN/ALS")

b  <- base::readRDS("backbone_FINAL.rds")
pl <- base::readRDS("patient_level_FINAL.rds")
meta <- base::as.data.frame(b$meta)
AX <- c("AXIS_METABOLIC", "AXIS_INFLAMMATION", "AXIS_CELL_CYCLE")

## ---------------------------------------------------------------------
## master frame
## ---------------------------------------------------------------------
df <- base::as.data.frame(pl$pat)
df$cluster_final <- base::factor(df$subtype, levels = base::paste0("S-ALS-", 1:3))
df$CCT2_expr <- df$CCT2
df$SEX <- base::factor(df$SEX)

## SV1..SV6 from the backbone metadata (anchor SVA), matched by sample_id
SVS <- base::paste0("SV", 1:6)
stopifnot(base::all(SVS %in% base::colnames(meta)))
sv_idx <- base::match(df$sample_id, meta$sample_id)
stopifnot(!base::any(base::is.na(sv_idx)))
for (s in SVS) df[[s]] <- meta[[s]][sv_idx]
stopifnot(!base::any(base::is.na(df[, SVS])))
COV  <- base::paste("AGE_AT_SYMPTOM_ONSET + SEX +",
                    base::paste(SVS, collapse = " + "))
base::cat("SVs attached to", base::nrow(df), "patients; COV =", COV, "\n")

## composites
if (base::file.exists("diff_quality_sensitivity_FINAL.rds")) {
  s <- base::readRDS("diff_quality_sensitivity_FINAL.rds")$meta
  cc <- base::c("astro_score","micro_score","oligo_score","maturity_score",
                "progenitor_score","diff_quality","prolif_score")
  cc <- cc[cc %in% base::colnames(s)]
  df <- base::merge(df, s[, base::c("Participant_ID", cc)],
                    by = "Participant_ID", all.x = TRUE)
  base::cat("composites merged: ", base::paste(cc, collapse = ", "), "\n", sep = "")
} else base::cat("NOTE: diff_quality_sensitivity_FINAL.rds absent; L7 will skip\n")

## age at death for disease duration
add_cols <- base::intersect(base::c("AGE_AT_DEATH"), base::colnames(meta))
if (base::length(add_cols)) {
  mm <- meta[!base::duplicated(meta$Participant_ID),
             base::c("Participant_ID", add_cols), drop = FALSE]
  df <- base::merge(df, mm, by = "Participant_ID", all.x = TRUE)
  df$disease_duration <- base::as.numeric(df$AGE_AT_DEATH) -
    base::as.numeric(df$AGE_AT_SYMPTOM_ONSET)
}

## ---- CGND -> PRS, resolving suffixed duplicates ----------------------
cg_cols <- base::grep("CGND", base::colnames(df), value = TRUE, ignore.case = TRUE)
if (base::length(cg_cols) == 0) {
  ## not on df, pull from meta
  cgm <- base::grep("CGND", base::colnames(meta), value = TRUE, ignore.case = TRUE)[1]
  if (!base::is.na(cgm)) {
    xw <- meta[!base::duplicated(meta$Participant_ID),
               base::c("Participant_ID", cgm)]
    base::colnames(xw) <- base::c("Participant_ID", "cgnd_raw")
    df <- base::merge(df, xw, by = "Participant_ID", all.x = TRUE)
  } else df$cgnd_raw <- NA_character_
} else {
  vals <- base::lapply(cg_cols, function(k) base::as.character(df[[k]]))
  df$cgnd_raw <- base::Reduce(function(a, bb)
    base::ifelse(!base::is.na(a) & a != "", a, bb), vals)
  base::cat("CGND columns found: ", base::paste(cg_cols, collapse = ", "), "\n", sep = "")
}
df$cgnd <- base::sub("-b[0-9]+$", "", df$cgnd_raw)
base::cat("non-missing CGND: ",
          base::sum(!base::is.na(df$cgnd) & df$cgnd != ""), " of ",
          base::nrow(df), "\n", sep = "")

prs_f <- "prset/prset_axis_metabolic.best"
if (base::file.exists(prs_f)) {
  p <- utils::read.table(prs_f, header = TRUE, stringsAsFactors = FALSE)
  p <- p[p$In_Regression == "Yes", ]
  prs <- base::data.frame(cgnd = base::sub("-b[0-9]+$", "", p$IID),
                          PRS_METABOLIC = p$AXIS_METABOLIC,
                          PRS_CELL_CYCLE = p$AXIS_CELL_CYCLE,
                          PRS_INFLAMMATION = p$AXIS_INFLAMMATION,
                          stringsAsFactors = FALSE)
  df <- base::merge(df, prs, by = "cgnd", all.x = TRUE)
  base::cat("PRS linked: ", base::sum(!base::is.na(df$PRS_CELL_CYCLE)),
            " of ", base::nrow(df), "   (expect ~355)\n", sep = "")
} else base::cat("PRS .best file not found; L5-L7 will skip\n")

base::cat("\nmaster frame: ", base::nrow(df), " patients\n", sep = "")
for (v in base::c("logNfL","ALSFRS_R_PROGRESSION_SLOPE","disease_duration",
                  "PRS_CELL_CYCLE","diff_quality","CCT2_expr"))
  if (v %in% base::colnames(df))
    base::cat("  ", base::sprintf("%-30s %d\n", v, base::sum(!base::is.na(df[[v]]))))

dc <- df[df$cluster_final %in% base::c("S-ALS-2","S-ALS-3"), ]
dc$cluster_final <- base::droplevels(dc$cluster_final)
base::cat("\nclean subset: ", base::nrow(dc), "\n", sep = "")
base::print(base::table(dc$cluster_final))

## ---------------------------------------------------------------------
## helpers
## ---------------------------------------------------------------------
grab <- function(fit, term, model, outcome, predictor) {
  s <- base::summary(fit)$coefficients
  if (!term %in% base::rownames(s)) return(NULL)
  base::data.frame(outcome = outcome, predictor = predictor, model = model,
                   n = stats::nobs(fit), beta = s[term, 1], se = s[term, 2],
                   p = s[term, 4], adjR2 = base::summary(fit)$adj.r.squared,
                   stringsAsFactors = FALSE)
}

ladder <- function(dat, outcome, predictor) {
  need <- base::c(outcome, predictor, "AGE_AT_SYMPTOM_ONSET", "SEX",
                  "cluster_final", base::paste0("SV", 1:6))
  if (!base::all(need %in% base::colnames(dat))) return(NULL)
  d <- dat[stats::complete.cases(dat[, need]), ]
  if (base::nrow(d) < 25) {
    base::cat("  ", outcome, " ~ ", predictor, ": n = ", base::nrow(d),
              ", skipped\n", sep = ""); return(NULL) }
  f <- function(rhs) stats::lm(stats::as.formula(
    base::paste(outcome, "~", rhs)), data = d)
  out <- base::list(
    grab(f(predictor), predictor, "univariate", outcome, predictor),
    grab(f(base::paste(predictor, "+", COV)),
         predictor, "+ age + sex + SVs", outcome, predictor),
    grab(f(base::paste(predictor, "+", COV, "+ cluster_final")),
         predictor, "+ subtype", outcome, predictor))
  for (lv in base::levels(d$cluster_final)) {
    ds <- d[d$cluster_final == lv, ]
    if (base::nrow(ds) < 40) next
    fit <- stats::lm(stats::as.formula(base::paste(
      outcome, "~", predictor, "+", COV)), data = ds)
    out <- base::c(out, base::list(grab(fit, predictor,
                                        base::paste0("within ", lv), outcome, predictor)))
  }
  base::do.call(base::rbind, out)
}

bic_tab <- function(dat, outcome) {
  need <- base::c(outcome, "CCT2_expr", "AXIS_METABOLIC",
                  "AGE_AT_SYMPTOM_ONSET", "SEX", "cluster_final",
                  base::paste0("SV", 1:6))
  if (!base::all(need %in% base::colnames(dat))) return(NULL)
  d <- dat[stats::complete.cases(dat[, need]), ]
  if (base::nrow(d) < 25) return(NULL)
  mods <- base::list(
    "subtype only"   = base::paste("cluster_final +", COV),
    "CCT2 only"      = base::paste("CCT2_expr +", COV),
    "AXIS only"      = base::paste("AXIS_METABOLIC +", COV),
    "subtype + CCT2" = base::paste("CCT2_expr + cluster_final +", COV),
    "subtype + AXIS" = base::paste("AXIS_METABOLIC + cluster_final +", COV),
    "all three"      = base::paste("CCT2_expr + AXIS_METABOLIC + cluster_final +", COV))
  r <- base::do.call(base::rbind, base::lapply(base::names(mods), function(k) {
    fit <- stats::lm(stats::as.formula(base::paste(outcome, "~", mods[[k]])), d)
    base::data.frame(outcome = outcome, model = k, n = stats::nobs(fit),
                     BIC = stats::BIC(fit),
                     adjR2 = base::summary(fit)$adj.r.squared,
                     stringsAsFactors = FALSE) }))
  r$dBIC <- r$BIC - base::min(r$BIC)
  r[base::order(r$BIC), ]
}

## ---------------------------------------------------------------------
## L1 to L4 and L8
## ---------------------------------------------------------------------
res <- base::list()
## SCOPE DECISION (audit, 2026-08-22): outcome ladders, BIC, and
## head-to-head run on ALL 391 sporadic patients with three-level
## subtype, matching the verified console run and the mediation claim
## (omnibus subtype F). The old Aug-13 file used S2/S3 only without SVs,
## under which the adjusted ALSFRS rung was p=0.115; with the
## Methods-compliant SVs the result is robust in BOTH scopes
## (all-391 p=0.039, S2/S3 p=0.032). PRS sections keep the S2/S3
## clean subset per the PRSet spec.
for (o in base::c("logNfL", "ALSFRS_R_PROGRESSION_SLOPE", "disease_duration")) {
  if (!o %in% base::colnames(df)) next
  if (base::sum(!base::is.na(df[[o]])) < 25) {
    base::cat("\n", o, ": too few observations, skipped\n", sep = ""); next }
  base::cat("\n===================================================\n")
  base::cat(o, "\n")
  base::cat("===================================================\n")
  for (pr in base::c("CCT2_expr", "AXIS_METABOLIC")) {
    r <- ladder(df, o, pr)
    if (!base::is.null(r)) {
      base::cat("\n", pr, ":\n", sep = "")
      base::print(base::data.frame(model = r$model, n = r$n,
                                   beta = base::round(r$beta, 4), p = base::signif(r$p, 3)),
                  row.names = FALSE)
      res[[base::paste(o, pr)]] <- r
    }
  }
  bt <- bic_tab(df, o)
  if (!base::is.null(bt)) {
    base::cat("\nBIC comparison:\n")
    base::print(base::data.frame(model = bt$model, n = bt$n,
                                 BIC = base::round(bt$BIC, 1), dBIC = base::round(bt$dBIC, 1),
                                 adjR2 = base::round(bt$adjR2, 4)), row.names = FALSE)
    res[[base::paste0("BIC_", o)]] <- bt
  }
  need <- base::c(o, "CCT2_expr", "AXIS_METABOLIC", "AGE_AT_SYMPTOM_ONSET",
                  "SEX", "cluster_final", base::paste0("SV", 1:6))
  d <- df[stats::complete.cases(df[, need]), ]
  if (base::nrow(d) >= 25) {
    hh <- stats::lm(stats::as.formula(base::paste(o,
      "~ CCT2_expr + AXIS_METABOLIC + cluster_final +", COV)), d)
    base::cat("\nhead to head:\n")
    base::print(base::round(base::summary(hh)$coefficients, 4))
    if (base::requireNamespace("car", quietly = TRUE)) {
      base::cat("VIF:\n"); base::print(base::round(car::vif(hh), 2)) }
  }
}

## ---------------------------------------------------------------------
## L5 to L7
## ---------------------------------------------------------------------
if ("PRS_CELL_CYCLE" %in% base::colnames(dc) &&
    base::sum(!base::is.na(dc$PRS_CELL_CYCLE)) >= 50) {
  base::cat("\n===================================================\n")
  base::cat("PRS analyses\n")
  base::cat("===================================================\n")
  
  base::cat("\nL5  each PRS -> CCT2 (sex-adjusted, clean subset)\n")
  for (pr in base::c("PRS_METABOLIC","PRS_CELL_CYCLE","PRS_INFLAMMATION")) {
    fit <- stats::lm(stats::as.formula(base::paste("CCT2_expr ~", pr, "+ SEX")), dc)
    co <- base::summary(fit)$coefficients[pr, ]
    base::cat(base::sprintf("  %-18s beta %+.4f  p %.4f  n %d\n",
                            pr, co[1], co[4], stats::nobs(fit)))
  }
  base::cat("\n  PRS_CELL_CYCLE -> CCT2, + subtype:\n")
  base::print(base::round(base::summary(stats::lm(
    CCT2_expr ~ PRS_CELL_CYCLE + SEX + cluster_final, dc))$coefficients, 4))
  
  base::cat("\nL6  3x3 PRS-axis matrix, clean subset, sex-adjusted\n")
  pv <- base::c("PRS_METABOLIC","PRS_CELL_CYCLE","PRS_INFLAMMATION")
  mat <- base::do.call(base::rbind, base::lapply(pv, function(pp)
    base::do.call(base::rbind, base::lapply(AX, function(aa) {
      fit <- stats::lm(stats::as.formula(base::paste(aa, "~", pp, "+ SEX")), dc)
      s <- base::summary(fit)
      base::data.frame(PRS = pp, Axis = aa, R2 = s$r.squared,
                       beta = s$coefficients[pp, 1], p = s$coefficients[pp, 4],
                       stringsAsFactors = FALSE) }))))
  base::print(base::data.frame(PRS = mat$PRS, Axis = mat$Axis,
                               R2 = base::round(mat$R2, 4), p = base::signif(mat$p, 3)), row.names = FALSE)
  base::cat("  mean off-diagonal R2: ",
            base::round(base::mean(mat$R2[base::sub("PRS_", "", mat$PRS) !=
                                            base::sub("AXIS_", "", mat$Axis)]), 5), "\n", sep = "")
  res$prs_matrix <- mat
  
  if ("diff_quality" %in% base::colnames(dc)) {
    base::cat("\nL7  PRS_CELL_CYCLE -> AXIS_CELL_CYCLE, progressive adjustment\n")
    steps <- base::c(
      "baseline"         = "PRS_CELL_CYCLE + SEX",
      "+ glial"          = "PRS_CELL_CYCLE + SEX + astro_score + micro_score + oligo_score",
      "+ glial + DQ"     = "PRS_CELL_CYCLE + SEX + astro_score + micro_score + oligo_score + diff_quality",
      "+ all composites" = "PRS_CELL_CYCLE + SEX + astro_score + micro_score + oligo_score + diff_quality + maturity_score + progenitor_score")
    pr7 <- base::do.call(base::rbind, base::lapply(base::names(steps), function(k) {
      fit <- stats::lm(stats::as.formula(
        base::paste("AXIS_CELL_CYCLE ~", steps[[k]])), dc)
      s <- base::summary(fit)
      base::data.frame(step = k, beta = s$coefficients["PRS_CELL_CYCLE", 1],
                       p = s$coefficients["PRS_CELL_CYCLE", 4],
                       adjR2 = s$adj.r.squared, n = stats::nobs(fit),
                       stringsAsFactors = FALSE) }))
    base::print(base::data.frame(step = pr7$step, n = pr7$n,
                                 beta = base::round(pr7$beta, 3), p = base::signif(pr7$p, 3),
                                 adjR2 = base::round(pr7$adjR2, 4)), row.names = FALSE)
    res$prs_adjusted <- pr7
    
    base::cat("\n  same progression on CCT2:\n")
    c7 <- base::do.call(base::rbind, base::lapply(base::names(steps), function(k) {
      fit <- stats::lm(stats::as.formula(
        base::paste("CCT2_expr ~", steps[[k]])), dc)
      s <- base::summary(fit)
      base::data.frame(step = k, beta = s$coefficients["PRS_CELL_CYCLE", 1],
                       p = s$coefficients["PRS_CELL_CYCLE", 4],
                       stringsAsFactors = FALSE) }))
    base::print(base::data.frame(step = c7$step, beta = base::round(c7$beta, 3),
                                 p = base::signif(c7$p, 3)), row.names = FALSE)
    res$prs_cct2_adjusted <- c7
  }
} else base::cat("\nPRS not linked; L5 to L7 skipped\n")

## ---------------------------------------------------------------------
## save
## ---------------------------------------------------------------------
lad_names <- base::grep(" ", base::names(res), value = TRUE)
if (base::length(lad_names)) {
  lad <- base::do.call(base::rbind, res[lad_names])
  utils::write.csv(lad, "downstream_ladders_FINAL.csv", row.names = FALSE)
  base::cat("\nwrote downstream_ladders_FINAL.csv (", base::nrow(lad),
            " rows)\n", sep = "")
}
bic_names <- base::grep("^BIC_", base::names(res), value = TRUE)
if (base::length(bic_names)) {
  utils::write.csv(base::do.call(base::rbind, res[bic_names]),
                   "downstream_BIC_FINAL.csv", row.names = FALSE)
  base::cat("wrote downstream_BIC_FINAL.csv\n")
}
## ---- verification gate: L1/L3 must reproduce the verified console run ----
lad_all <- base::do.call(base::rbind, res[base::grep(" ", base::names(res), value = TRUE)])
pick <- function(outc, mod) {
  r <- lad_all[lad_all$outcome == outc & lad_all$predictor == "CCT2_expr" &
                 lad_all$model == mod, ]
  if (base::nrow(r) == 1) r else NULL
}
targets <- base::list(
  base::list(pick("logNfL", "+ age + sex + SVs"), p = 0.04516912, n = 297),
  base::list(pick("logNfL", "+ subtype"),         p = 0.75187214, n = 297),
  base::list(pick("ALSFRS_R_PROGRESSION_SLOPE", "+ age + sex + SVs"),
             p = 0.00476965, n = 244),
  base::list(pick("ALSFRS_R_PROGRESSION_SLOPE", "+ subtype"),
             p = 0.03937459, n = 244))
ok <- base::vapply(targets, function(t)
  !base::is.null(t[[1]]) && t[[1]]$n == t$n &&
    base::abs(t[[1]]$p - t$p) < 1e-6, logical(1))
base::cat("\nverification (L1/L3 vs verified console run):",
          base::sum(ok), "of 4 rungs match\n")
if (base::all(ok)) {
  if (base::file.exists("downstream_analyses_FINAL.rds") &&
      !base::file.exists("downstream_analyses_spotcheck_SV.rds"))
    base::file.rename("downstream_analyses_FINAL.rds",
                      "downstream_analyses_spotcheck_SV.rds")
  base::saveRDS(base::c(res, base::list(df = df, clean = dc, covariates = COV)),
                "downstream_analyses_FINAL.rds")
  base::cat("EXACT MATCH: saved downstream_analyses_FINAL.rds (full ladders,",
            "SV-adjusted). Script VERIFIED.\n")
} else {
  base::saveRDS(base::c(res, base::list(df = df, clean = dc, covariates = COV)),
                "downstream_analyses_AUDITCHECK.rds")
  base::cat("MISMATCH: saved as downstream_analyses_AUDITCHECK.rds; the",
            "existing FINAL file is untouched. Print the four rungs and",
            "report back.\n")
}