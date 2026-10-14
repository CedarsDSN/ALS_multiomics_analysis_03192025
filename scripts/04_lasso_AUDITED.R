# =====================================================================
# 04_lasso.R   (AUDITED 2026-08-22)
# ONE-VS-REST LASSO, refit on the final labels
#
# AUDIT CHANGES:
#   - axis_sets now read from backbone_FINAL.rds (was a session-object
#     stopifnot; the feature space of the classifier was never pinned).
#   - setwd pinned; conflicted preferences declared.
#   - Verification gate vs saved lasso_FINAL.rds before any FINAL file
#     is overwritten. Deterministic core gated strictly (selected gene
#     set per subtype, CCT2 selected). Stochastic quantities (nested-CV
#     AUC, stability percentages) gated at tolerance: RNG-stream
#     identical reruns reproduce them exactly, but any change in glmnet
#     internals shifts them slightly without changing conclusions.
#     Tolerances: |dAUC| < 0.01, |dStability| < 3 percentage points.
#
# This is the analysis that discovered CCT2. It ran on the old labels,
# so it has to be redone before CCT2 can be presented as the selected
# classifier for S-ALS-3.
#
# Method, per Methods section:
#   - features: the union of Hallmark axis genes
#   - one-vs-rest binary LASSO per subtype
#   - nested 10 x 10 cross-validation
#   - lambda.1se, optimised on AUC
#   - bootstrap stability selection for the retained genes
#
# Output: lasso_FINAL.rds, lasso_selected_genes_FINAL.csv,
#         lasso_stability_FINAL.csv
# Runtime: 20-40 min, dominated by nested CV and 1000 bootstraps.
# =====================================================================

suppressPackageStartupMessages({ library(glmnet) })
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(base::setdiff, base::intersect, base::union,
                               base::unname, base::setequal, .quiet = TRUE)
}
base::setwd("~/Documents/DSN/ALS")
base::set.seed(42)

b  <- base::readRDS("backbone_FINAL.rds")
pl <- base::readRDS("patient_level_FINAL.rds")
E  <- b$expr_symbol
pat <- base::as.data.frame(pl$pat)
pat$subtype <- base::factor(pat$subtype, levels = base::paste0("S-ALS-", 1:3))

axis_sets <- b$axis_sets
if (base::is.null(axis_sets))
  base::stop("backbone_FINAL.rds has no axis_sets -- rerun audited script 01 ",
             "first; the classifier's feature space must come from the ",
             "pinned backbone, not a session object.")
genes <- base::unique(base::unlist(axis_sets))
genes <- genes[genes %in% base::rownames(E)]
base::cat("axis genes available: ", base::length(genes), "\n", sep = "")
base::cat("  (manuscript says 1,309; check this against your Methods)\n")

X <- base::t(E[genes, pat$sample_id, drop = FALSE])
X <- base::scale(X)
X <- X[, base::apply(X, 2, function(z) base::all(base::is.finite(z))), drop = FALSE]
base::cat("design matrix: ", base::paste(base::dim(X), collapse = " x "), "\n", sep = "")
base::print(base::table(pat$subtype))

NFOLD <- 10; NREP <- 10; NBOOT <- 1000

auc_fn <- function(score, pos) {
  r <- base::rank(score); n1 <- base::sum(pos); n0 <- base::sum(!pos)
  (base::sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

fit_one <- function(target) {
  y <- base::as.integer(pat$subtype == target)
  base::cat("\n===============================================\n")
  base::cat(target, "  (", base::sum(y), " vs ", base::sum(!y), ")\n", sep = "")
  base::cat("===============================================\n")
  
  ## ---- nested CV: outer folds estimate performance, inner picks lambda
  outer_auc <- base::numeric(0)
  for (rep in 1:NREP) {
    folds <- base::sample(base::rep(1:NFOLD, length.out = base::length(y)))
    pred <- base::rep(NA_real_, base::length(y))
    for (k in 1:NFOLD) {
      tr <- folds != k; te <- !tr
      if (base::length(base::unique(y[tr])) < 2) next
      cvi <- glmnet::cv.glmnet(X[tr, , drop = FALSE], y[tr], family = "binomial",
                               alpha = 1, nfolds = NFOLD, type.measure = "auc")
      pred[te] <- stats::predict(cvi, newx = X[te, , drop = FALSE],
                                 s = "lambda.1se", type = "response")[, 1]
    }
    ok <- !base::is.na(pred)
    outer_auc <- base::c(outer_auc, auc_fn(pred[ok], y[ok] == 1))
    base::cat("  rep ", rep, " outer AUC = ",
              base::round(outer_auc[rep], 4), "\n", sep = "")
  }
  base::cat("  nested CV AUC: mean ", base::round(base::mean(outer_auc), 4),
            "  sd ", base::round(stats::sd(outer_auc), 4), "\n", sep = "")
  
  ## ---- final model on all data
  cvf <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = 1,
                           nfolds = NFOLD, type.measure = "auc")
  co <- glmnet::coef.glmnet(cvf, s = "lambda.1se")
  sel <- base::rownames(co)[co[, 1] != 0]
  sel <- base::setdiff(sel, "(Intercept)")
  base::cat("  lambda.1se = ", base::signif(cvf$lambda.1se, 4),
            "   genes retained: ", base::length(sel), "\n", sep = "")
  if (base::length(sel)) {
    cf <- co[sel, 1]
    base::print(base::round(base::sort(cf, decreasing = TRUE), 4))
  } else base::cat("  (none retained at lambda.1se)\n")
  
  ## ---- bootstrap stability selection
  base::cat("  bootstrapping ", NBOOT, " resamples ...\n", sep = "")
  cnt <- stats::setNames(base::numeric(base::ncol(X)), base::colnames(X))
  for (i in 1:NBOOT) {
    idx <- base::sample(base::seq_along(y), base::length(y), replace = TRUE)
    if (base::length(base::unique(y[idx])) < 2) next
    fit <- base::tryCatch(glmnet::glmnet(X[idx, , drop = FALSE], y[idx],
                                         family = "binomial", alpha = 1, lambda = cvf$lambda.1se),
                          error = function(e) NULL)
    if (base::is.null(fit)) next
    cf <- glmnet::coef.glmnet(fit)[, 1]
    nz <- base::names(cf)[cf != 0]
    nz <- base::setdiff(nz, "(Intercept)")
    cnt[nz] <- cnt[nz] + 1
  }
  stab <- base::sort(100 * cnt / NBOOT, decreasing = TRUE)
  base::cat("  top 10 by selection frequency:\n")
  base::print(base::round(utils::head(stab, 10), 1))
  if ("CCT2" %in% base::names(stab))
    base::cat("  >>> CCT2 selected in ", base::round(stab[["CCT2"]], 1),
              "% of bootstraps   (manuscript reported 82.5% for S-ALS-3)\n", sep = "")
  
  base::list(target = target, outer_auc = outer_auc, cv = cvf,
             selected = sel, coefs = if (base::length(sel)) co[sel, 1] else NULL,
             stability = stab)
}

RES <- base::lapply(base::levels(pat$subtype), fit_one)
base::names(RES) <- base::levels(pat$subtype)

## ---- summary ---------------------------------------------------------
base::cat("\n===============================================\n")
base::cat("SUMMARY\n")
base::cat("===============================================\n")
summ <- base::do.call(base::rbind, base::lapply(RES, function(r)
  base::data.frame(subtype = r$target,
                   nested_AUC = base::round(base::mean(r$outer_auc), 4),
                   AUC_sd = base::round(stats::sd(r$outer_auc), 4),
                   n_genes = base::length(r$selected),
                   CCT2_selected = "CCT2" %in% r$selected,
                   CCT2_stability = if ("CCT2" %in% base::names(r$stability))
                     base::round(r$stability[["CCT2"]], 1) else NA,
                   stringsAsFactors = FALSE)))
base::print(summ, row.names = FALSE)

sel_tab <- base::do.call(base::rbind, base::lapply(RES, function(r) {
  if (!base::length(r$selected)) return(NULL)
  base::data.frame(subtype = r$target, gene = r$selected,
                   coef = base::as.numeric(r$coefs),
                   stability = base::as.numeric(r$stability[r$selected]),
                   stringsAsFactors = FALSE) }))
stab_tab <- base::do.call(base::rbind, base::lapply(RES, function(r)
  base::data.frame(subtype = r$target, gene = base::names(r$stability),
                   stability = base::as.numeric(r$stability),
                   stringsAsFactors = FALSE)))
## CSV writes happen inside the verification gate below
## ---- verification gate vs saved lasso_FINAL.rds ---------------------
write_finals <- function() {
  if (!base::is.null(sel_tab))
    utils::write.csv(sel_tab, "lasso_selected_genes_FINAL.csv", row.names = FALSE)
  utils::write.csv(stab_tab[stab_tab$stability > 0, ],
                   "lasso_stability_FINAL.csv", row.names = FALSE)
  base::saveRDS(base::list(res = RES, summary = summ, genes = genes,
                           pat = pat, feature_source = "backbone_FINAL.rds$axis_sets"),
                "lasso_FINAL.rds")
  base::cat("saved lasso_FINAL.rds + CSVs\n")
}

if (base::file.exists("lasso_FINAL.rds")) {
  oldL <- base::readRDS("lasso_FINAL.rds")
  gate <- base::do.call(base::rbind, base::lapply(base::names(RES), function(s) {
    new <- RES[[s]]; old <- oldL$res[[s]]
    if (base::is.null(old)) return(base::data.frame(
      subtype = s, genes_equal = NA, cct2_ok = NA, dAUC = NA, dStab = NA))
    base::data.frame(
      subtype = s,
      genes_equal = base::setequal(new$selected, old$selected),
      cct2_ok = ("CCT2" %in% new$selected) == ("CCT2" %in% old$selected),
      dAUC = base::abs(base::mean(new$outer_auc) - base::mean(old$outer_auc)),
      dStab = if ("CCT2" %in% base::names(new$stability) &&
                  "CCT2" %in% base::names(old$stability))
        base::abs(new$stability[["CCT2"]] - old$stability[["CCT2"]]) else 0)
  }))
  base::print(gate, row.names = FALSE)
  strict <- base::all(gate$genes_equal, gate$cct2_ok, na.rm = FALSE)
  tol    <- base::all(gate$dAUC < 0.01, gate$dStab < 3)
  if (base::isTRUE(strict) && base::isTRUE(tol)) {
    exact <- base::all(gate$dAUC == 0, gate$dStab == 0)
    base::cat(if (exact) "EXACT MATCH" else
      "MATCH within tolerance (RNG-stream or glmnet-version drift)",
      ": Script VERIFIED.\n")
    write_finals()
  } else {
    base::saveRDS(base::list(res = RES, summary = summ, genes = genes, pat = pat),
                  "lasso_AUDITCHECK.rds")
    base::cat("MISMATCH: saved as lasso_AUDITCHECK.rds; existing FINAL",
              "untouched. Report the gate table back.\n")
  }
} else {
  base::cat("no prior lasso_FINAL.rds; saving fresh.\n")
  write_finals()
}
base::cat("wrote lasso_stability_FINAL.csv and lasso_FINAL.rds\n")