################################################################################
## clinical_endpoint_screen.R
##
## Paul's request: beyond NfL and ALSFRS-R, which AnswerALS clinical variables
## carry a subtype / CCT2 signal?
##
## Screens five new endpoint families against the sporadic subtypes and CCT2:
##   (A) Vital capacity (SVC/FVC % predicted) slope        -- continuous, longitudinal
##   (B) Composite survival: death OR trach OR PAV (Cox)   -- proper time-to-event
##   (C) Grip strength slope                               -- motor-specific, supporting
##   (D) Cognitive (CBS total) + behavioral (carbeh) slope -- FTD axis, supporting
##   (E) Site of onset (bulbar vs limb) x subtype          -- prognostic stratifier xtab
##
## Conventions (matching df_nfl_uniq / df_als_only workflow):
##   - join on sample_id via Participant_ID crosswalk: sub("-[0-9]+-T$","",sample_id)
##   - covariates: AGE_AT_SYMPTOM_ONSET + SEX + SV1..SV6
##   - primary contrast: S-ALS-2 (ref) vs S-ALS-3, sporadic only
##   - CCT2: univariate, then + subtype adjustment (parallel-tracks reading)
##   - honest framing: this is a SCREEN. Report power (n) for every test; do not
##     bury underpowered nulls.
##
## Each form lives in its own CSV named after the Form Name column, e.g.
##   Vital_Capacity.csv, Grip_Strength_Testing.csv, ALS_CBS.csv,
##   CNS_Lability_Scale.csv, Mortality.csv, Tracheostomy.csv,
##   Permanent_Assisted_Ventilation.csv, AALSHXFX.csv, Demographics.csv
## Set CLIN_DIR to wherever those live.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(survival)
})

## ---------------------------------------------------------------------------
## 0. CONFIG -- edit these two lines to your environment
## ---------------------------------------------------------------------------
CLIN_DIR <- "~/Documents/DSN/ALS/metadata/clinical"
OUT_DIR  <- "~/Documents/DSN/ALS/clinical_endpoints"   # where to write the summary table / RDS

DAYS_PER_YEAR <- 365.25

## Helper: strip sample suffix to Participant_ID (the established crosswalk)
to_pid <- function(x) sub("-[0-9]+-T$", "", x)

## Helper: robust longitudinal slope per participant via lm(value ~ day).
## Returns annualized slope; requires >= MIN_VISITS distinct timepoints
## spanning >= MIN_SPAN_DAYS, else NA (so single-visit patients are dropped,
## not forced to slope 0 -- this is exactly the kind of silent bias to avoid).
MIN_VISITS    <- 3
MIN_SPAN_DAYS <- 90
slope_per_pid <- function(dt, id_col, day_col, val_col) {
  dt <- dt[is.finite(dt[[day_col]]) & is.finite(dt[[val_col]])]
  setDT(dt)
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

## ---------------------------------------------------------------------------
## 1. ANCHOR TABLE -- subtype + CCT2 + covariates, one row per sample/patient
##
## We reuse the established objects. df_nfl_uniq already carries:
##   sample_id, group, logNfL, ALSFRS_R_PROGRESSION_SLOPE,
##   COGNITIVE_PROGRESSION_SLOPE, AGE_AT_SYMPTOM_ONSET, SEX, SV1..SV6
## and cluster_clean / ax_clean_df carry the FINAL subtype labels.
## cct2_df carries CCT2 VST (needs distinct(Participant_ID) due to longits).
##
## Build `anchor`: sample_id, Participant_ID, subtype, CCT2, covariates,
## restricted to sporadic S-ALS-2 / S-ALS-3.
## ---------------------------------------------------------------------------

## stopifnot(exists("df_nfl_uniq"), exists("cluster_clean"),
##           exists("ax_clean_df"), exists("cct2_df"))
## 
## ## subtype labels keyed by sample_id (cluster_clean aligns to ax_clean_df rows)
## subtype_tab <- data.table(
##   sample_id = rownames(ax_clean_df),
##   subtype   = paste0("S-ALS-", as.integer(cluster_clean))
## )
## 
## ## CCT2, de-duplicated to one row per participant (established trap)
## cct2_tab <- as.data.table(cct2_df) %>%
##   distinct(Participant_ID, .keep_all = TRUE) %>%
##   transmute(Participant_ID = as.character(Participant_ID),
##             CCT2 = CCT2)   # CCT2 column = VST value; rename here if yours differs
## 
## ## covariates + group from the established clinical anchor
## cov_tab <- as.data.table(df_nfl_uniq)[, .(
##   sample_id,
##   group,
##   AGE_AT_SYMPTOM_ONSET,
##   SEX,
##   SV1, SV2, SV3, SV4, SV5, SV6
## )]
## 
## anchor <- merge(cov_tab, subtype_tab, by = "sample_id", all.x = TRUE)
## anchor[, Participant_ID := to_pid(sample_id)]
## anchor <- merge(anchor, cct2_tab, by = "Participant_ID", all.x = TRUE)
## 
## ## sporadic only + clean half (S-ALS-2 vs S-ALS-3); S-ALS-2 as reference
## anchor <- anchor[group == "Sporadic" & subtype %in% c("S-ALS-2", "S-ALS-3")]
## anchor[, is_S3   := as.integer(subtype == "S-ALS-3")]
## anchor[, subtype := factor(subtype, levels = c("S-ALS-2", "S-ALS-3"))]
## anchor[, CCT2_z  := as.numeric(scale(CCT2))]   # z-score CCT2 for comparable betas

cat(sprintf("Anchor: %d patients (S-ALS-2=%d, S-ALS-3=%d)\n",
            nrow(anchor), sum(anchor$subtype=="S-ALS-2"),
            sum(anchor$subtype=="S-ALS-3")))

## ---------------------------------------------------------------------------
## 2. ENDPOINT MODELLING HELPERS
## ---------------------------------------------------------------------------
COVARS <- "AGE_AT_SYMPTOM_ONSET + SEX + SV1 + SV2 + SV3 + SV4 + SV5 + SV6"

## Linear endpoint: subtype model, CCT2-univariate, CCT2+subtype (mediation read)
fit_linear_endpoint <- function(dat, endpoint, label) {
  dat <- dat[is.finite(dat[[endpoint]])]
  n   <- nrow(dat)
  if (n < 20) {
    return(data.table(endpoint = label, n = n,
                      beta_S3 = NA, p_S3 = NA,
                      beta_CCT2_uni = NA, p_CCT2_uni = NA,
                      beta_CCT2_adj = NA, p_CCT2_adj = NA,
                      note = "n<20, not modelled"))
  }
  f_sub <- as.formula(sprintf("%s ~ is_S3 + %s", endpoint, COVARS))
  f_uni <- as.formula(sprintf("%s ~ CCT2_z + %s", endpoint, COVARS))
  f_adj <- as.formula(sprintf("%s ~ CCT2_z + is_S3 + %s", endpoint, COVARS))
  m_sub <- lm(f_sub, dat); m_uni <- lm(f_uni, dat); m_adj <- lm(f_adj, dat)
  cs <- function(m, term) {
    s <- summary(m)$coefficients
    if (term %in% rownames(s)) s[term, c("Estimate","Pr(>|t|)")] else c(NA,NA)
  }
  b_s3  <- cs(m_sub, "is_S3")
  b_un  <- cs(m_uni, "CCT2_z")
  b_aj  <- cs(m_adj, "CCT2_z")
  data.table(
    endpoint = label, n = n,
    beta_S3 = b_s3[1], p_S3 = b_s3[2],
    beta_CCT2_uni = b_un[1], p_CCT2_uni = b_un[2],
    beta_CCT2_adj = b_aj[1], p_CCT2_adj = b_aj[2],
    note = "+CCT2 collapse vs stable = state-marker vs graded"
  )
}

## ---------------------------------------------------------------------------
## 3. (A) VITAL CAPACITY SLOPE  -- SVC/FVC % predicted
##
## vctyp: 1=SVC, 2=FVC. svcp1/2/3 and fvcp1/2/3 = % predicted per trial.
## ALS convention: prefer SVC; take the best (max) valid trial per visit;
## upright preferred (vcpos==2) but keep all and let slope absorb it.
## Filter svcdn==1 (Done). Then slope of %predicted vs Visit_Date (days).
## ---------------------------------------------------------------------------
vc <- fread(file.path(CLIN_DIR, "Vital_Capacity.csv"))
vc[, Participant_ID := as.character(Participant_ID)]

## --- (A) VITAL CAPACITY: coerce all %pred trial columns to numeric first ---
pct_cols <- c("svcp1","svcp2","svcp3","fvcp1","fvcp2","fvcp3")
vc[, (pct_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(as.character(x)))),
   .SDcols = pct_cols]

vc[, svc_pct := mapply(maxvalid, svcp1, svcp2, svcp3)]
vc[, fvc_pct := mapply(maxvalid, fvcp1, fvcp2, fvcp3)]
vc[, vc_pct  := fifelse(is.finite(svc_pct), svc_pct, fvc_pct)]   # SVC preferred
vc_use <- vc[is.finite(vc_pct) & is.finite(Visit_Date)]

vc_slope <- slope_per_pid(vc_use, "Participant_ID", "Visit_Date", "vc_pct")
setnames(vc_slope, "slope_yr", "VC_slope")
vc_dat <- merge(anchor, vc_slope[, .(Participant_ID, VC_slope, n_visits)],
                by = "Participant_ID", all.x = TRUE)
res_VC <- fit_linear_endpoint(vc_dat, "VC_slope", "VC %pred slope (/yr)")

## ---------------------------------------------------------------------------
## 4. (C) GRIP STRENGTH SLOPE  -- max of L/R, both trials, pounds
##   grpstat==1 Done; *grpntck==1 means "not tested" -> exclude that side
## ---------------------------------------------------------------------------
gs <- fread(file.path(CLIN_DIR, "Grip_Strength_Testing.csv"))
gs[, Participant_ID := as.character(Participant_ID)]
gs[lgrpntck == 1, c("lgrpt1","lgrpt2") := NA]
gs[rgrpntck == 1, c("rgrpt1","rgrpt2") := NA]
gs[, grip := mapply(maxvalid, lgrpt1, lgrpt2, rgrpt1, rgrpt2)]
gs_use <- gs[is.finite(grip) & is.finite(Visit_Date)]
gs_slope <- slope_per_pid(gs_use, "Participant_ID", "Visit_Date", "grip")
setnames(gs_slope, "slope_yr", "GRIP_slope")
gs_dat <- merge(anchor, gs_slope[, .(Participant_ID, GRIP_slope, n_visits)],
                by = "Participant_ID", all.x = TRUE)
res_GRIP <- fit_linear_endpoint(gs_dat, "GRIP_slope", "Grip strength slope (lb/yr)")

## ---------------------------------------------------------------------------
## 5. (D) COGNITIVE (CBS) + BEHAVIORAL (carbeh) SLOPE
##   cbstot = CBS total; carbeh = behavioral total. cbsdn==1 Collected.
## ---------------------------------------------------------------------------
cbs <- fread(file.path(CLIN_DIR, "ALS_CBS.csv"))
cbs[, Participant_ID := as.character(Participant_ID)]
cbs_use <- cbs[is.finite(Visit_Date)]

cbs_cog <- slope_per_pid(cbs_use[is.finite(cbstot)],
                         "Participant_ID","Visit_Date","cbstot")
setnames(cbs_cog, "slope_yr", "CBS_slope")
cog_dat <- merge(anchor, cbs_cog[, .(Participant_ID, CBS_slope)],
                 by = "Participant_ID", all.x = TRUE)
res_CBS <- fit_linear_endpoint(cog_dat, "CBS_slope", "CBS cognitive slope (/yr)")

cbs_beh <- slope_per_pid(cbs_use[is.finite(carbeh)],
                         "Participant_ID","Visit_Date","carbeh")
setnames(cbs_beh, "slope_yr", "BEH_slope")
beh_dat <- merge(anchor, cbs_beh[, .(Participant_ID, BEH_slope)],
                 by = "Participant_ID", all.x = TRUE)
res_BEH <- fit_linear_endpoint(beh_dat, "BEH_slope", "CBS behavioral slope (/yr)")

## CNS lability (pseudobulbar affect) total -- baseline value, not slope
## (typically too sparse longitudinally; test baseline as a quick look)
cns <- fread(file.path(CLIN_DIR, "CNS_Lability_Scale.csv"))
cns[, Participant_ID := as.character(Participant_ID)]
cns_base <- cns[is.finite(cnslstot) & is.finite(Visit_Date)][
  order(Participant_ID, Visit_Date),
  .(CNS_baseline = cnslstot[1]), by = Participant_ID]
cns_dat <- merge(anchor, cns_base, by = "Participant_ID", all.x = TRUE)
res_CNS <- fit_linear_endpoint(cns_dat, "CNS_baseline", "CNS lability baseline")

## ---------------------------------------------------------------------------
## 6. (B) COMPOSITE SURVIVAL  -- time to (death OR tracheostomy OR PAV)
##   The proper time-to-event answer to the disease-duration / survivorship
##   problem. Time origin = symptom onset (AALSHXFX onsetdt is days from
##   baseline visit; Visit_Date fields are also days-from-baseline, so all
##   event dates are on the same day-from-baseline axis). Censor at last
##   observed clinical visit if no event.
## ---------------------------------------------------------------------------
mort  <- fread(file.path(CLIN_DIR, "Mortality.csv"))[, .(
  Participant_ID = as.character(Participant_ID), event_day = dieddt)]
trach <- fread(file.path(CLIN_DIR, "Tracheostomy.csv"))[, .(
  Participant_ID = as.character(Participant_ID), event_day = trchdt)]
pav   <- fread(file.path(CLIN_DIR, "Permanent_Assisted_Ventilation.csv"))[, .(
  Participant_ID = as.character(Participant_ID), event_day = pavdt)]

events <- rbindlist(list(mort, trach, pav))[is.finite(event_day)]
first_event <- events[, .(event_day = min(event_day)), by = Participant_ID]

## last-seen day = max Visit_Date across the longitudinal forms we loaded
last_seen <- rbindlist(list(
  vc_use[,  .(Participant_ID, d = Visit_Date)],
  gs_use[,  .(Participant_ID, d = Visit_Date)],
  cbs_use[, .(Participant_ID, d = Visit_Date)]
))[is.finite(d), .(last_day = max(d)), by = Participant_ID]

surv_tab <- merge(anchor[, .(Participant_ID, sample_id, subtype, is_S3, CCT2_z,
                             AGE_AT_SYMPTOM_ONSET, SEX,
                             SV1,SV2,SV3,SV4,SV5,SV6)],
                  last_seen, by = "Participant_ID", all.x = TRUE)
surv_tab <- merge(surv_tab, first_event, by = "Participant_ID", all.x = TRUE)
surv_tab[, status := as.integer(is.finite(event_day))]
surv_tab[, time_d := fifelse(status == 1, event_day, last_day)]
surv_tab <- surv_tab[is.finite(time_d) & time_d > 0]

## Cox: subtype, then CCT2; covariate-adjusted
cox_subtype <- coxph(Surv(time_d, status) ~ is_S3 +
                       AGE_AT_SYMPTOM_ONSET + SEX, data = surv_tab)
cox_cct2    <- coxph(Surv(time_d, status) ~ CCT2_z +
                       AGE_AT_SYMPTOM_ONSET + SEX, data = surv_tab)
cox_both    <- coxph(Surv(time_d, status) ~ CCT2_z + is_S3 +
                       AGE_AT_SYMPTOM_ONSET + SEX, data = surv_tab)

cox_row <- function(m, term, label) {
  s <- summary(m)$coefficients
  if (!term %in% rownames(s)) return(data.table(model=label, HR=NA, p=NA))
  data.table(model = label,
             HR = exp(s[term,"coef"]),
             p  = s[term,"Pr(>|z|)"])
}
res_SURV <- rbindlist(list(
  cox_row(cox_subtype, "is_S3",  "S-ALS-3 (subtype only)"),
  cox_row(cox_cct2,    "CCT2_z", "CCT2 (univariate)"),
  cox_row(cox_both,    "CCT2_z", "CCT2 (+ subtype)"),
  cox_row(cox_both,    "is_S3",  "S-ALS-3 (+ CCT2)")
))
res_SURV[, n_events := sum(surv_tab$status)]
res_SURV[, n_total  := nrow(surv_tab)]

## ---------------------------------------------------------------------------
## 7. (E) SITE OF ONSET x SUBTYPE  -- bulbar vs limb prognostic stratifier
## ---------------------------------------------------------------------------
hx <- fread(file.path(CLIN_DIR, "AALSHXFX.csv"))
hx[, Participant_ID := as.character(Participant_ID)]
hx_site <- hx[, .(bulbar = max(hxblb, na.rm=TRUE),
                  limb   = max(hxli,  na.rm=TRUE)), by = Participant_ID]
hx_site[, onset_site := fifelse(bulbar==1 & limb!=1, "bulbar",
                                fifelse(limb==1 & bulbar!=1, "limb", NA_character_))]
site_dat <- merge(anchor, hx_site[, .(Participant_ID, onset_site)],
                  by = "Participant_ID", all.x = TRUE)
site_tab <- table(site_dat$subtype, site_dat$onset_site)
site_test <- if (all(dim(site_tab) >= 2)) fisher.test(site_tab) else NULL

## ---------------------------------------------------------------------------
## 8. ASSEMBLE + WRITE
## ---------------------------------------------------------------------------
linear_results <- rbindlist(list(res_VC, res_GRIP, res_CBS, res_BEH, res_CNS),
                            use.names = TRUE)

## BH-FDR across the linear subtype tests (honest multiple-testing control)
linear_results[, p_S3_FDR := p.adjust(p_S3, method = "BH")]

cat("\n================ LINEAR ENDPOINT SCREEN ================\n")
print(linear_results)
cat("\n================ COMPOSITE SURVIVAL (Cox) =============\n")
print(res_SURV)
cat(sprintf("\nSurvival: %d patients, %d events (death/trach/PAV)\n",
            res_SURV$n_total[1], res_SURV$n_events[1]))
cat("\n================ ONSET SITE x SUBTYPE =================\n")
print(site_tab)
if (!is.null(site_test))
  cat(sprintf("Fisher p = %.3g\n", site_test$p.value))

saveRDS(list(linear = linear_results, survival = res_SURV,
             site_tab = site_tab, site_test = site_test,
             anchor = anchor),
        file.path(OUT_DIR, "clinical_endpoint_screen_results.rds"))

cat("\nDone. Results saved to clinical_endpoint_screen_results.rds\n")
################################################################################


# Let's do one more round with more clinical variables
################################################################################
## clinical_endpoint_screen_v2.R
##
## Second-wave clinical endpoints for the sporadic subtypes (S-ALS-2 vs S-ALS-3)
## and CCT2, following the v1 screen. Adds:
##   (F) HHD megascore slope  -- z-scored 18-muscle composite (objective motor,
##                               better-powered upgrade of the grip result)
##   (G) Body weight slope + BMI slope (Vital_Signs) -- catabolic/metabolic axis,
##                               mechanistically aligned with S-ALS-3 metabolic state
##   (H) Serum creatinine + CK (Auxiliary_Chemistry) -- muscle-mass / breakdown
##   (I) Survival composite v2: death OR trach OR PAV OR gastrostomy (more events)
##
## Reuses the v1 anchor + fit_linear_endpoint + slope_per_pid. RUN v1 FIRST so
## `anchor`, `slope_per_pid`, `fit_linear_endpoint`, `COVARS`, `to_pid`,
## `maxvalid`, `last_seen`, `first_event` already exist. If running standalone,
## copy those helpers over.
##
## Conventions identical to v1: sporadic only, S-ALS-2 ref vs S-ALS-3,
## covariates AGE_AT_SYMPTOM_ONSET + SEX + SV1..6, CCT2 z-scored, BH-FDR.
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(survival)
})

stopifnot(exists("anchor"), exists("fit_linear_endpoint"),
          exists("slope_per_pid"), exists("maxvalid"), exists("to_pid"))

CLIN_DIR <- get0("CLIN_DIR",
                 ifnotfound = "/common/compbiomed-dsn/projects/ALS_multiomics_analysis_03192025/clinical")
DAYS_PER_YEAR <- 365.25

## ---------------------------------------------------------------------------
## (F) HAND-HELD DYNAMOMETRY MEGASCORE
##
## 18 muscle stems, each with <stem>t1/t2/t3 (trials) and <stem>nt (not-tested).
## Per visit, per muscle: take max valid trial (force); if nt==1, that muscle
## is NA for the visit. Then z-score EACH muscle across all visit-rows (so the
## different force ranges across muscles are comparable), and the megascore for
## a visit = mean of available z-scored muscle values. Slope of megascore vs day.
## ---------------------------------------------------------------------------
hhd <- fread(file.path(CLIN_DIR, "Hand_Held_Dynamometry.csv"))
hhd[, Participant_ID := as.character(Participant_ID)]

stems <- c("lad","lee","lef","lfdic","lhf","lke","lkf","lsf","lwe",
           "rad","ree","ref","rfdic","rhf","rke","rkf","rsf","rwe")

## per-muscle max valid trial, respecting the not-tested checkbox
for (s in stems) {
  t1 <- paste0(s, "t1"); t2 <- paste0(s, "t2"); t3 <- paste0(s, "t3")
  nt <- paste0(s, "nt")
  # coerce trials to numeric defensively
  hhd[, (c(t1,t2,t3)) := lapply(.SD, function(x) suppressWarnings(as.numeric(as.character(x)))),
      .SDcols = c(t1,t2,t3)]
  hhd[, (s) := mapply(maxvalid, get(t1), get(t2), get(t3))]
  if (nt %in% names(hhd)) hhd[get(nt) == 1, (s) := NA_real_]
}

## z-score each muscle across all visit-rows (centre/scale per muscle)
for (s in stems) {
  mu <- mean(hhd[[s]], na.rm = TRUE); sdv <- sd(hhd[[s]], na.rm = TRUE)
  hhd[, (paste0(s, "_z")) := if (is.finite(sdv) && sdv > 0) (get(s) - mu)/sdv else NA_real_]
}

## megascore = row mean of available z-scored muscles; require >= 6 muscles
z_cols <- paste0(stems, "_z")
hhd[, n_musc := rowSums(!is.na(.SD)), .SDcols = z_cols]
hhd[, hhd_mega := rowMeans(.SD, na.rm = TRUE), .SDcols = z_cols]
hhd[n_musc < 6, hhd_mega := NA_real_]

hhd_use <- hhd[is.finite(hhd_mega) & is.finite(Visit_Date)]
hhd_slope <- slope_per_pid(hhd_use, "Participant_ID", "Visit_Date", "hhd_mega")
setnames(hhd_slope, "slope_yr", "HHD_slope")
hhd_dat <- merge(anchor, hhd_slope[, .(Participant_ID, HHD_slope, n_visits)],
                 by = "Participant_ID", all.x = TRUE)
res_HHD <- fit_linear_endpoint(hhd_dat, "HHD_slope", "HHD megascore slope (z/yr)")

## ---------------------------------------------------------------------------
## (G) BODY WEIGHT + BMI SLOPE  (Vital_Signs)
##   weight units vary (weightu); convert lb->kg so the slope is interpretable.
##   1 = (commonly) lb, 2 = kg in AnswerALS coding -- confirm; we convert any
##   value flagged lb. BMI is unitless and safe as-is.
## ---------------------------------------------------------------------------
vs <- fread(file.path(CLIN_DIR, "Vital_Signs.csv"))
vs[, Participant_ID := as.character(Participant_ID)]
vs[, weight := suppressWarnings(as.numeric(as.character(weight)))]
vs[, bmi    := suppressWarnings(as.numeric(as.character(bmi)))]

## normalise weight to kg where flagged as lb (weightu==1 assumed lb); guard if absent
if ("weightu" %in% names(vs)) {
  vs[, weight_kg := fifelse(weightu == 1 & is.finite(weight), weight * 0.453592, weight)]
} else {
  vs[, weight_kg := weight]   # assume already kg
}
## plausibility filter
vs[!is.finite(weight_kg) | weight_kg < 25 | weight_kg > 250, weight_kg := NA_real_]
vs[!is.finite(bmi)       | bmi < 10       | bmi > 70,        bmi := NA_real_]

wt_use <- vs[is.finite(weight_kg) & is.finite(Visit_Date)]
wt_slope <- slope_per_pid(wt_use, "Participant_ID", "Visit_Date", "weight_kg")
setnames(wt_slope, "slope_yr", "WT_slope")
wt_dat <- merge(anchor, wt_slope[, .(Participant_ID, WT_slope)],
                by = "Participant_ID", all.x = TRUE)
res_WT <- fit_linear_endpoint(wt_dat, "WT_slope", "Body weight slope (kg/yr)")

bmi_use <- vs[is.finite(bmi) & is.finite(Visit_Date)]
bmi_slope <- slope_per_pid(bmi_use, "Participant_ID", "Visit_Date", "bmi")
setnames(bmi_slope, "slope_yr", "BMI_slope")
bmi_dat <- merge(anchor, bmi_slope[, .(Participant_ID, BMI_slope)],
                 by = "Participant_ID", all.x = TRUE)
res_BMI <- fit_linear_endpoint(bmi_dat, "BMI_slope", "BMI slope (/yr)")

## ---------------------------------------------------------------------------
## (H) SERUM CREATININE + CK  (Auxiliary_Chemistry)
##   Values live in the TEXT result fields (accrrslt, acckrslt). Parse to numeric.
##   Test as BASELINE (earliest visit) -- these are sparse longitudinally.
##   Creatinine is a muscle-mass prognostic marker; CK reflects muscle breakdown.
## ---------------------------------------------------------------------------
ac <- fread(file.path(CLIN_DIR, "Auxiliary_Chemistry.csv"))
ac[, Participant_ID := as.character(Participant_ID)]
ac[, creat := suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(accrrslt))))]
ac[, ck    := suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(acckrslt))))]
## plausibility (serum creatinine ~0.2-3 mg/dL; CK can be high but cap absurd)
ac[!is.finite(creat) | creat <= 0 | creat > 20,   creat := NA_real_]
ac[!is.finite(ck)    | ck    <= 0 | ck    > 50000, ck   := NA_real_]

creat_base <- ac[is.finite(creat) & is.finite(Visit_Date)][
  order(Participant_ID, Visit_Date), .(creat_base = creat[1]), by = Participant_ID]
creat_dat <- merge(anchor, creat_base, by = "Participant_ID", all.x = TRUE)
res_CREAT <- fit_linear_endpoint(creat_dat, "creat_base", "Creatinine baseline (mg/dL)")

ck_base <- ac[is.finite(ck) & is.finite(Visit_Date)][
  order(Participant_ID, Visit_Date), .(ck_base = log10(ck[1])), by = Participant_ID]  # log CK
ck_dat <- merge(anchor, ck_base, by = "Participant_ID", all.x = TRUE)
res_CK <- fit_linear_endpoint(ck_dat, "ck_base", "log10 CK baseline")

## ---------------------------------------------------------------------------
## (I) SURVIVAL COMPOSITE v2: death / trach / PAV / gastrostomy
##   Augments v1 first_event with feeding-tube placement to raise event count.
## ---------------------------------------------------------------------------
gtube <- fread(file.path(CLIN_DIR, "Feeding_Tube_Placement.csv"))[, .(
  Participant_ID = as.character(Participant_ID), event_day = ftpaccdt)]
gtube <- gtube[is.finite(event_day)]

## first_event from v1 already pools death/trach/PAV; add gastrostomy
stopifnot(exists("first_event"), exists("last_seen"))
events_v2 <- rbindlist(list(first_event[, .(Participant_ID, event_day)],
                            gtube), use.names = TRUE)
first_event_v2 <- events_v2[, .(event_day = min(event_day)), by = Participant_ID]

surv_v2 <- merge(anchor[, .(Participant_ID, is_S3, CCT2_z,
                            AGE_AT_SYMPTOM_ONSET, SEX)],
                 last_seen, by = "Participant_ID", all.x = TRUE)
surv_v2 <- merge(surv_v2, first_event_v2, by = "Participant_ID", all.x = TRUE)
surv_v2[, status := as.integer(is.finite(event_day))]
surv_v2[, time_d := fifelse(status == 1, event_day, last_day)]
surv_v2 <- surv_v2[is.finite(time_d) & time_d > 0]

cox2_sub  <- coxph(Surv(time_d, status) ~ is_S3  + AGE_AT_SYMPTOM_ONSET + SEX, surv_v2)
cox2_cct2 <- coxph(Surv(time_d, status) ~ CCT2_z + AGE_AT_SYMPTOM_ONSET + SEX, surv_v2)
cox2_both <- coxph(Surv(time_d, status) ~ CCT2_z + is_S3 + AGE_AT_SYMPTOM_ONSET + SEX, surv_v2)
cox_row <- function(m, term, label) {
  s <- summary(m)$coefficients
  if (!term %in% rownames(s)) return(data.table(model=label, HR=NA, p=NA))
  data.table(model = label, HR = exp(s[term,"coef"]), p = s[term,"Pr(>|z|)"])
}
res_SURV2 <- rbindlist(list(
  cox_row(cox2_sub,  "is_S3",  "S-ALS-3 (subtype only)"),
  cox_row(cox2_cct2, "CCT2_z", "CCT2 (univariate)"),
  cox_row(cox2_both, "CCT2_z", "CCT2 (+ subtype)"),
  cox_row(cox2_both, "is_S3",  "S-ALS-3 (+ CCT2)")
))
res_SURV2[, n_events := sum(surv_v2$status)]
res_SURV2[, n_total  := nrow(surv_v2)]

## ---------------------------------------------------------------------------
## ASSEMBLE
## ---------------------------------------------------------------------------
linear_v2 <- rbindlist(list(res_HHD, res_WT, res_BMI, res_CREAT, res_CK),
                       use.names = TRUE)
linear_v2[, p_S3_FDR := p.adjust(p_S3, method = "BH")]

cat("\n============== LINEAR ENDPOINT SCREEN v2 ==============\n")
print(linear_v2)
cat("\n============== SURVIVAL COMPOSITE v2 (death/trach/PAV/gtube) =======\n")
print(res_SURV2)
cat(sprintf("\nSurvival v2: %d patients, %d events\n",
            res_SURV2$n_total[1], res_SURV2$n_events[1]))

## quick coverage report (so power is explicit, per the screen philosophy)
cat("\n--- endpoint coverage (n with valid value in anchor) ---\n")
cat(sprintf("  HHD megascore slope : %d\n", nrow(hhd_dat[is.finite(HHD_slope)])))
cat(sprintf("  Weight slope        : %d\n", nrow(wt_dat[is.finite(WT_slope)])))
cat(sprintf("  BMI slope           : %d\n", nrow(bmi_dat[is.finite(BMI_slope)])))
cat(sprintf("  Creatinine baseline : %d\n", nrow(creat_dat[is.finite(creat_base)])))
cat(sprintf("  log CK baseline     : %d\n", nrow(ck_dat[is.finite(ck_base)])))

saveRDS(list(linear = linear_v2, survival_v2 = res_SURV2),
        "clinical_endpoint_screen_v2_results.rds")
cat("\nDone. Saved clinical_endpoint_screen_v2_results.rds\n")
################################################################################


# Hhd fdi test

################################################################################
## hhd_fdi_test.R
##
## Adjudicate the grip-vs-HHD non-replication.
##   Grip strength (distal hand) separates subtypes (p=7e-4, n=54).
##   HHD whole-body megascore does NOT (n=81, null).
## Question: is grip a regional distal-hand signal that the proximal-heavy
##   megascore dilutes, or a fragile small-n effect?
##
## Test the HHD muscle closest to grip: first dorsal interosseous (FDI),
##   left+right (lfdic / rfdic). If FDI tracks S-ALS-3 like grip -> regional
##   reading holds. If FDI is also null -> grip stands alone, soften it.
##
## Run AFTER v1 + v2 (needs anchor, slope_per_pid, slope_per_pid_2 if used,
## fit_linear_endpoint, maxvalid, hhd object already loaded & parsed).
## Reuses the hhd data.table built in v2 (per-muscle max-trial cols present).
################################################################################

suppressPackageStartupMessages({ library(data.table); library(dplyr) })
stopifnot(exists("anchor"), exists("fit_linear_endpoint"), exists("hhd"))

DAYS_PER_YEAR <- 365.25

## relaxed slope helper (>=2 visits, >=90 days) so FDI is compared to HHD on the
## SAME footing as the n=81 megascore run (apples-to-apples threshold)
slope2 <- function(dt, id_col, day_col, val_col) {
  dt <- as.data.table(dt)[is.finite(get(day_col)) & is.finite(get(val_col))]
  dt[, .(n_visits = .N,
         slope_yr = if (uniqueN(get(day_col)) >= 2 &&
                        (max(get(day_col)) - min(get(day_col))) >= 90)
           as.numeric(coef(lm(get(val_col) ~ get(day_col)))[2]) * DAYS_PER_YEAR
         else NA_real_),
     by = id_col]
}

## also a >=3 version to mirror grip's threshold exactly
slope3 <- function(dt, id_col, day_col, val_col) {
  dt <- as.data.table(dt)[is.finite(get(day_col)) & is.finite(get(val_col))]
  dt[, .(n_visits = .N,
         slope_yr = if (uniqueN(get(day_col)) >= 3 &&
                        (max(get(day_col)) - min(get(day_col))) >= 90)
           as.numeric(coef(lm(get(val_col) ~ get(day_col)))[2]) * DAYS_PER_YEAR
         else NA_real_),
     by = id_col]
}

## ---------------------------------------------------------------------------
## Build FDI-only score per visit.
## hhd already has per-muscle max-trial cols lfdic / rfdic (built in v2),
## already NA'd where not-tested. FDI score = mean of available L/R, z-scored.
## ---------------------------------------------------------------------------
stopifnot(all(c("lfdic","rfdic") %in% names(hhd)))

## z-score each side across visit-rows (so L/R on comparable scale), then average
for (s in c("lfdic","rfdic")) {
  mu <- mean(hhd[[s]], na.rm = TRUE); sdv <- sd(hhd[[s]], na.rm = TRUE)
  hhd[, (paste0(s, "_z")) := if (is.finite(sdv) && sdv > 0) (get(s) - mu)/sdv else NA_real_]
}
hhd[, fdi_n := rowSums(!is.na(.SD)), .SDcols = c("lfdic_z","rfdic_z")]
hhd[, fdi_score := rowMeans(.SD, na.rm = TRUE), .SDcols = c("lfdic_z","rfdic_z")]
hhd[fdi_n < 1, fdi_score := NA_real_]

fdi_use <- hhd[is.finite(fdi_score) & is.finite(Visit_Date)]

## ---- FDI slope at both thresholds ----
fdi2 <- slope2(fdi_use, "Participant_ID", "Visit_Date", "fdi_score")
setnames(fdi2, "slope_yr", "FDI_slope")
fdi_dat2 <- merge(anchor, unique(fdi2, by="Participant_ID")[, .(Participant_ID, FDI_slope)],
                  by = "Participant_ID", all.x = TRUE)
res_FDI2 <- fit_linear_endpoint(fdi_dat2, "FDI_slope", "FDI slope (z/yr, >=2 visits)")

fdi3 <- slope3(fdi_use, "Participant_ID", "Visit_Date", "fdi_score")
setnames(fdi3, "slope_yr", "FDI_slope")
fdi_dat3 <- merge(anchor, unique(fdi3, by="Participant_ID")[, .(Participant_ID, FDI_slope)],
                  by = "Participant_ID", all.x = TRUE)
res_FDI3 <- fit_linear_endpoint(fdi_dat3, "FDI_slope", "FDI slope (z/yr, >=3 visits)")

## ---------------------------------------------------------------------------
## Direct grip-vs-FDI agreement on shared patients (do the two distal-hand
## measures even agree with each other?)  -- sanity on the regional reading.
## ---------------------------------------------------------------------------
agree_txt <- "gs_slope not in session; skipped grip-FDI correlation"
if (exists("gs_slope")) {
  gv <- unique(as.data.table(gs_slope)[, .(Participant_ID, GRIP_slope)], by="Participant_ID")
  fv <- unique(fdi2[, .(Participant_ID, FDI_slope)], by="Participant_ID")
  chk <- merge(gv, fv, by="Participant_ID")
  if (nrow(chk) >= 10) {
    rr <- cor(chk$GRIP_slope, chk$FDI_slope, use="complete.obs")
    agree_txt <- sprintf("grip vs FDI slope: r = %.3f (n = %d shared)", rr, nrow(chk))
  }
}

## ---------------------------------------------------------------------------
## REPORT
## ---------------------------------------------------------------------------
cat("\n================ FDI-SPECIFIC HHD TEST ================\n")
cat("Reference results:\n")
cat("  Grip strength slope : beta_S3 = -21.78, p = 7.1e-4 (n=54)  [SIGNIFICANT]\n")
cat("  HHD megascore slope : beta_S3 = +0.03,  p = 0.87   (n=81)  [NULL]\n\n")
cat(">>> FDI (first dorsal interosseous, distal hand) -- closest HHD muscle to grip:\n\n")
print(res_FDI2)
cat("\n")
print(res_FDI3)
cat("\n", agree_txt, "\n", sep="")

cat("\n--- INTERPRETATION GUIDE ---\n")
cat("If FDI beta_S3 is NEGATIVE and significant (like grip):\n")
cat("   -> regional reading holds: distal-hand decline tracks S-ALS-3; the\n")
cat("      whole-body megascore dilutes it. Grip + FDI = coherent distal signal.\n")
cat("If FDI is NULL too:\n")
cat("   -> grip stands alone; objective-motor signal is grip-specific / fragile.\n")
cat("      Soften the grip claim, disclose non-replication across HHD + FDI.\n")
cat("Also check grip-FDI r: high positive r supports that they measure the\n")
cat("   same construct; near-zero r undercuts the regional rescue.\n")
################################################################################

# clinical endpoint v3
################################################################################
## clinical_endpoint_screen_v3.R
##
## Third-wave endpoints:
##   (J) Auxiliary_Chemistry_Labs -- SECOND chemistry form (separate from the
##       sparse Auxiliary_Chemistry). Test urate (uric acid) + creatinine baseline.
##       Urate is an established ALS prognostic antioxidant biomarker.
##   (K) UMN / LMN burden (AALSDXFX) -- count of body regions with UMN / LMN
##       signs at diagnostic exam. Mechanistic: do the motor-neuron-intrinsic
##       subtypes map onto the UMN-vs-LMN-predominant clinical axis?
##   (L) El Escorial diagnostic certainty (ordinal) by subtype -- quick check.
##
## Run AFTER v1 (needs anchor, fit_linear_endpoint, to_pid, COVARS).
## UMN/LMN coding (AALSDXFX): 1=Yes, 2=No, 90=Not Done  per region.
## Regions: bulbar, L/R upper ext, L/R lower ext, trunk.
################################################################################

suppressPackageStartupMessages({ library(data.table); library(dplyr) })
stopifnot(exists("anchor"), exists("fit_linear_endpoint"), exists("to_pid"))

## ---------------------------------------------------------------------------
## (J) AUXILIARY_CHEMISTRY_LABS  -- urate + creatinine baseline
##   First: COVERAGE CHECK. If this form is as sparse as the other one, stop.
## ---------------------------------------------------------------------------
acl <- fread(file.path(CLIN_DIR, "Auxiliary_Chemistry_Labs.csv"))
acl[, Participant_ID := as.character(Participant_ID)]

cat("=== Auxiliary_Chemistry_Labs coverage ===\n")
cat("  rows:", nrow(acl), "\n")
cat("  creatinine text populated:", sum(!is.na(acl$accrrslt)), "\n")
cat("  uric acid text populated :", sum(!is.na(acl$acuarslt)), "\n")
cat("  CK text populated        :", sum(!is.na(acl$acckrslt)), "\n")

parse_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(x))))
acl[, creat := parse_num(accrrslt)]
acl[, urate := parse_num(acuarslt)]
acl[, ck    := parse_num(acckrslt)]
## plausibility
acl[!is.finite(creat) | creat <= 0 | creat > 20,    creat := NA_real_]
acl[!is.finite(urate) | urate <= 0 | urate > 20,    urate := NA_real_]   # mg/dL
acl[!is.finite(ck)    | ck    <= 0 | ck    > 50000,  ck   := NA_real_]

date_col <- if ("Visit_Date" %in% names(acl)) "Visit_Date" else "labdt"
baseline_of <- function(dt, val) {
  d <- dt[is.finite(get(val)) & is.finite(get(date_col))]
  d[order(Participant_ID, get(date_col)), .(v = get(val)[1]), by = Participant_ID]
}
creat_b <- baseline_of(acl, "creat"); setnames(creat_b, "v", "creat_base")
urate_b <- baseline_of(acl, "urate"); setnames(urate_b, "v", "urate_base")

creat_dat <- merge(anchor, creat_b, by = "Participant_ID", all.x = TRUE)
urate_dat <- merge(anchor, urate_b, by = "Participant_ID", all.x = TRUE)
res_CREAT <- fit_linear_endpoint(creat_dat, "creat_base", "Creatinine baseline (mg/dL)")
res_URATE <- fit_linear_endpoint(urate_dat, "urate_base", "Uric acid baseline (mg/dL)")

## ---------------------------------------------------------------------------
## (K) UMN / LMN BURDEN  (AALSDXFX)
##   Per region UMN col (*cumn) and LMN col (*clmn): 1=Yes, 2=No, 90=Not Done.
##   Burden = count of regions scored Yes(==1). Take earliest exam per patient.
## ---------------------------------------------------------------------------
dx <- fread(file.path(CLIN_DIR, "AALSDXFX.csv"))
dx[, Participant_ID := as.character(Participant_ID)]

umn_cols <- c("blbcumn","luecumn","ruecumn","llecumn","rlecumn","trnkcumn")
lmn_cols <- c("blbclmn","lueclmn","rueclmn","lleclmn","rleclmn","trnkclmn")
umn_cols <- umn_cols[umn_cols %in% names(dx)]
lmn_cols <- lmn_cols[lmn_cols %in% names(dx)]

## count "Yes" (==1) across regions; treat 90/NA as not-positive
yes_count <- function(DT, cols) rowSums(sapply(cols, function(c) DT[[c]] == 1), na.rm = TRUE)
dx[, umn_burden := yes_count(dx, umn_cols)]
dx[, lmn_burden := yes_count(dx, lmn_cols)]
dx[, n_umn_scored := rowSums(sapply(umn_cols, function(c) dx[[c]] %in% c(1,2)))]
dx[, n_lmn_scored := rowSums(sapply(lmn_cols, function(c) dx[[c]] %in% c(1,2)))]

date_dx <- if ("Visit_Date" %in% names(dx)) "Visit_Date" else "alsdxdt"
dx_base <- dx[is.finite(get(date_dx))][order(Participant_ID, get(date_dx)),
                                       .(umn_burden = umn_burden[1], lmn_burden = lmn_burden[1],
                                         n_umn_scored = n_umn_scored[1], n_lmn_scored = n_lmn_scored[1],
                                         elescr = if ("elescrlr" %in% names(dx)) elescrlr[1] else NA_real_),
                                       by = Participant_ID]
## only keep patients where at least some regions were actually scored
dx_base[n_umn_scored == 0, umn_burden := NA_real_]
dx_base[n_lmn_scored == 0, lmn_burden := NA_real_]

umn_dat <- merge(anchor, dx_base[, .(Participant_ID, umn_burden)], by="Participant_ID", all.x=TRUE)
lmn_dat <- merge(anchor, dx_base[, .(Participant_ID, lmn_burden)], by="Participant_ID", all.x=TRUE)
res_UMN <- fit_linear_endpoint(umn_dat, "umn_burden", "UMN burden (n regions)")
res_LMN <- fit_linear_endpoint(lmn_dat, "lmn_burden", "LMN burden (n regions)")

## ---------------------------------------------------------------------------
## (L) El Escorial certainty by subtype (ordinal 1-5) -- quick look
## ---------------------------------------------------------------------------
el_dat <- merge(anchor, dx_base[, .(Participant_ID, elescr)], by="Participant_ID", all.x=TRUE)
res_EL  <- fit_linear_endpoint(el_dat, "elescr", "El Escorial certainty (1-5)")

## ---------------------------------------------------------------------------
## ASSEMBLE
## ---------------------------------------------------------------------------
linear_v3 <- rbindlist(list(res_CREAT, res_URATE, res_UMN, res_LMN, res_EL),
                       use.names = TRUE)
linear_v3[, p_S3_FDR := p.adjust(p_S3, method = "BH")]

cat("\n============== LINEAR ENDPOINT SCREEN v3 ==============\n")
print(linear_v3)

cat("\n--- coverage (n in anchor) ---\n")
cat(sprintf("  Creatinine baseline : %d\n", nrow(creat_dat[is.finite(creat_base)])))
cat(sprintf("  Uric acid baseline  : %d\n", nrow(urate_dat[is.finite(urate_base)])))
cat(sprintf("  UMN burden          : %d\n", nrow(umn_dat[is.finite(umn_burden)])))
cat(sprintf("  LMN burden          : %d\n", nrow(lmn_dat[is.finite(lmn_burden)])))
cat(sprintf("  El Escorial         : %d\n", nrow(el_dat[is.finite(elescr)])))

## descriptive: mean UMN/LMN burden per subtype (interpretation aid)
cat("\n--- UMN / LMN burden by subtype (means) ---\n")
print(umn_dat[is.finite(umn_burden), .(mean_UMN = round(mean(umn_burden),2), n=.N), by=subtype])
print(lmn_dat[is.finite(lmn_burden), .(mean_LMN = round(mean(lmn_burden),2), n=.N), by=subtype])

saveRDS(linear_v3, "clinical_endpoint_screen_v3_results.rds")
cat("\nDone. Saved clinical_endpoint_screen_v3_results.rds\n")
################################################################################

# clinical endpoint v4
################################################################################
## clinical_endpoint_screen_v4.R
##
## Final two clinical measures, for completeness (both index UMN/spasticity,
## which the AALSDXFX UMN-burden test already covered as a null at n=234):
##   (M) Ashworth spasticity -- mean tone across 4 limbs; baseline + slope
##   (N) Reflexes UMN-sign burden -- count of pathological signs present
##       (Babinski, Hoffman, clonus, jaw jerk, palmomental, crossed adduction);
##       baseline + slope. Plus a brisk-DTR count (deep-tendon reflexes >=3).
##
## Run AFTER v1 (needs anchor, fit_linear_endpoint, slope_per_pid, to_pid).
##
## Codings:
##   Ashworth *arml/*armr/*legl/*legr: 1-5 severity, 6/7 = not tested -> NA
##   Reflex pathological signs (*bsscr Babinski, *hsscr Hoffman, *cclnscr cervical
##     clonus, *lclnscr lumbosacral clonus, cjjscr jaw jerk, cpsscr palmomental,
##     *cadscr crossed adduction): 1=Absent, 2=Present
##   Deep-tendon reflexes (biceps/triceps/brachiorad/patellar/ankle/finger-flex):
##     0-4 ordinal; >=3 = brisk (UMN-ish)
################################################################################

suppressPackageStartupMessages({ library(data.table); library(dplyr) })
stopifnot(exists("anchor"), exists("fit_linear_endpoint"), exists("slope_per_pid"))

## ---------------------------------------------------------------------------
## (M) ASHWORTH SPASTICITY
## ---------------------------------------------------------------------------
ash <- fread(file.path(CLIN_DIR, "Ashworth_Spasticity_Scale.csv"))
ash[, Participant_ID := as.character(Participant_ID)]
ash_cols <- c("asharml","asharmr","ashlegl","ashlegr")
ash_cols <- ash_cols[ash_cols %in% names(ash)]
## recode: 6/7 (not tested) -> NA; keep 1-5
for (c in ash_cols) {
  ash[, (c) := suppressWarnings(as.numeric(as.character(get(c))))]
  ash[get(c) %in% c(6,7), (c) := NA_real_]
}
ash[, ash_n   := rowSums(!is.na(.SD)), .SDcols = ash_cols]
ash[, ash_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = ash_cols]
ash[ash_n < 1, ash_mean := NA_real_]

date_ash <- if ("Visit_Date" %in% names(ash)) "Visit_Date" else "ashdt"
ash_use <- ash[is.finite(ash_mean) & is.finite(get(date_ash))]

## baseline
ash_base <- ash_use[order(Participant_ID, get(date_ash)),
                    .(ash_base = ash_mean[1]), by = Participant_ID]
ashb_dat <- merge(anchor, ash_base, by = "Participant_ID", all.x = TRUE)
res_ASH_b <- fit_linear_endpoint(ashb_dat, "ash_base", "Ashworth mean baseline (1-5)")

## slope
ash_slope <- slope_per_pid(ash_use, "Participant_ID", date_ash, "ash_mean")
setnames(ash_slope, "slope_yr", "ASH_slope")
ashs_dat <- merge(anchor, unique(ash_slope, by="Participant_ID")[, .(Participant_ID, ASH_slope)],
                  by = "Participant_ID", all.x = TRUE)
res_ASH_s <- fit_linear_endpoint(ashs_dat, "ASH_slope", "Ashworth slope (/yr)")

## ---------------------------------------------------------------------------
## (N) REFLEXES UMN-SIGN BURDEN
## ---------------------------------------------------------------------------
rfx <- fread(file.path(CLIN_DIR, "Reflexes.csv"))
rfx[, Participant_ID := as.character(Participant_ID)]

## pathological signs: coded 1=Absent, 2=Present -> present = (==2)
path_cols <- c("lbsscr","rbsscr",      # Babinski
               "lhsscr","rhsscr",      # Hoffman
               "lcclnscr","rcclnscr",  # cervical clonus
               "llclnscr","rlclnscr",  # lumbosacral clonus
               "lcadscr","rcadscr",    # crossed adduction
               "cjjscr","cpsscr","crcfrscr")  # jaw jerk, palmomental, facial
path_cols <- path_cols[path_cols %in% names(rfx)]
for (c in path_cols) rfx[, (c) := suppressWarnings(as.numeric(as.character(get(c))))]

rfx[, umn_present := rowSums(sapply(path_cols, function(c) rfx[[c]] == 2), na.rm = TRUE)]
rfx[, umn_scored  := rowSums(sapply(path_cols, function(c) rfx[[c]] %in% c(1,2)))]
rfx[umn_scored == 0, umn_present := NA_real_]

## brisk deep-tendon reflexes (0-4 ordinal): count >= 3
dtr_cols <- c("lbcpscr","rbcpscr","ltrcpscr","rtrcpscr","lbrchscr","rbrchscr",
              "lptlrscr","rptlrscr","lanklscr","ranklscr","lffscr","rffscr")
dtr_cols <- dtr_cols[dtr_cols %in% names(rfx)]
for (c in dtr_cols) rfx[, (c) := suppressWarnings(as.numeric(as.character(get(c))))]
rfx[, brisk_dtr  := rowSums(sapply(dtr_cols, function(c) rfx[[c]] >= 3), na.rm = TRUE)]
rfx[, dtr_scored := rowSums(sapply(dtr_cols, function(c) is.finite(rfx[[c]])))]
rfx[dtr_scored == 0, brisk_dtr := NA_real_]

date_rfx <- if ("Visit_Date" %in% names(rfx)) "Visit_Date" else "rflxdt"

## baseline UMN-sign burden
rfx_use_u <- rfx[is.finite(umn_present) & is.finite(get(date_rfx))]
umn_base <- rfx_use_u[order(Participant_ID, get(date_rfx)),
                      .(umn_sign_base = umn_present[1]), by = Participant_ID]
umnb_dat <- merge(anchor, umn_base, by = "Participant_ID", all.x = TRUE)
res_UMNSIGN_b <- fit_linear_endpoint(umnb_dat, "umn_sign_base", "UMN-sign count baseline")

## slope of UMN-sign burden
umn_slope <- slope_per_pid(rfx_use_u, "Participant_ID", date_rfx, "umn_present")
setnames(umn_slope, "slope_yr", "UMNSIGN_slope")
umns_dat <- merge(anchor, unique(umn_slope, by="Participant_ID")[, .(Participant_ID, UMNSIGN_slope)],
                  by = "Participant_ID", all.x = TRUE)
res_UMNSIGN_s <- fit_linear_endpoint(umns_dat, "UMNSIGN_slope", "UMN-sign slope (/yr)")

## brisk-DTR baseline
rfx_use_d <- rfx[is.finite(brisk_dtr) & is.finite(get(date_rfx))]
dtr_base <- rfx_use_d[order(Participant_ID, get(date_rfx)),
                      .(brisk_base = brisk_dtr[1]), by = Participant_ID]
dtrb_dat <- merge(anchor, dtr_base, by = "Participant_ID", all.x = TRUE)
res_DTR_b <- fit_linear_endpoint(dtrb_dat, "brisk_base", "Brisk-DTR count baseline")

## ---------------------------------------------------------------------------
## ASSEMBLE
## ---------------------------------------------------------------------------
linear_v4 <- rbindlist(list(res_ASH_b, res_ASH_s, res_UMNSIGN_b,
                            res_UMNSIGN_s, res_DTR_b), use.names = TRUE)
linear_v4[, p_S3_FDR := p.adjust(p_S3, method = "BH")]

cat("\n============== LINEAR ENDPOINT SCREEN v4 (Ashworth + Reflexes) ==============\n")
print(linear_v4)

cat("\n--- coverage (n in anchor) ---\n")
cat(sprintf("  Ashworth baseline    : %d\n", nrow(ashb_dat[is.finite(ash_base)])))
cat(sprintf("  Ashworth slope       : %d\n", nrow(ashs_dat[is.finite(ASH_slope)])))
cat(sprintf("  UMN-sign baseline    : %d\n", nrow(umnb_dat[is.finite(umn_sign_base)])))
cat(sprintf("  UMN-sign slope       : %d\n", nrow(umns_dat[is.finite(UMNSIGN_slope)])))
cat(sprintf("  Brisk-DTR baseline   : %d\n", nrow(dtrb_dat[is.finite(brisk_base)])))

cat("\n--- means by subtype (interpretation aid) ---\n")
print(umnb_dat[is.finite(umn_sign_base), .(mean_UMNsign = round(mean(umn_sign_base),2), n=.N), by=subtype])
print(ashb_dat[is.finite(ash_base),     .(mean_Ashworth = round(mean(ash_base),2),    n=.N), by=subtype])

saveRDS(linear_v4, "clinical_endpoint_screen_v4_results.rds")
cat("\nDone. Saved clinical_endpoint_screen_v4_results.rds\n")
################################################################################



################################################################################
## clinical_endpoint_screen_v5.R
##
## Uses the AnswerALS data-portal MASTER table (one row per participant,
## consortium-curated derived values) -- aals_dataportal_datatable.csv.
## Join key Participant_ID is already CASE-NEU... format (matches anchor).
##
##   (O) ALSFRS-R baseline severity (n~890)  -- NEW (we only tested the slope)
##   (P) CBS baseline cognition (n~659)      -- higher-powered re-test of the
##                                              weak n=88 cognitive slope null
##   (Q) Portal ALSFRS-R progression slope (n~542) -- RECONCILE against the
##       ~186 slope in df_nfl_uniq; if concordant, a power-upgraded re-test of
##       the anchor endpoint
##   (R) Revised El Escorial certainty (n~891, full coverage) -- re-test
##   (S) Number_of_Visits by subtype -- dropout/retention robustness check
##
## Run AFTER v1 (needs anchor, fit_linear_endpoint). df_nfl_uniq for reconcile.
################################################################################

suppressPackageStartupMessages({ library(data.table); library(dplyr) })
stopifnot(exists("anchor"), exists("fit_linear_endpoint"))

## PORTAL_PATH <- "/path/to/aals_dataportal_datatable.csv"

port <- fread("~/Documents/DSN/ALS/metadata/aals_dataportal_datatable.csv")
port[, Participant_ID := as.character(Participant_ID)]

## ---- merge selected portal columns onto the anchor ----
keep <- c("Participant_ID","ALSFRS_R_Baseline_Value","ALSFRS_R_Latest_Value",
          "ALSFRS_R_PROGRESSION_SLOPE","CBS_Baseline_Value","CBS_Latest_Value",
          "COGNITIVE_PROGRESSION_SLOPE","REVISED_EL_ESCORIAL_CRITERIA",
          "Number_of_Visits","AGE_AT_DEATH")
keep <- keep[keep %in% names(port)]
portsel <- port[, ..keep]

## numeric coercions
numcols <- setdiff(keep, c("Participant_ID","REVISED_EL_ESCORIAL_CRITERIA"))
portsel[, (numcols) := lapply(.SD, function(x) suppressWarnings(as.numeric(as.character(x)))),
        .SDcols = numcols]

## El Escorial text -> ordinal 1-5
el_map <- c("Suspected"=1, "Possible"=2, "Probable Laboratory Supported"=3,
            "Probable"=4, "Definite"=5)
portsel[, elescr_portal := el_map[as.character(REVISED_EL_ESCORIAL_CRITERIA)]]

dat <- merge(anchor, portsel, by = "Participant_ID", all.x = TRUE)

## ---------------------------------------------------------------------------
## (O) ALSFRS-R baseline severity
## ---------------------------------------------------------------------------
res_ALSFRSbase <- fit_linear_endpoint(dat, "ALSFRS_R_Baseline_Value",
                                      "ALSFRS-R baseline (0-48)")

## ---------------------------------------------------------------------------
## (P) CBS baseline cognition (higher power than the longitudinal slope)
## ---------------------------------------------------------------------------
res_CBSbase <- fit_linear_endpoint(dat, "CBS_Baseline_Value",
                                   "CBS baseline (0-20)")

## ---------------------------------------------------------------------------
## (Q) Portal ALSFRS-R slope -- RECONCILE first, then test
## ---------------------------------------------------------------------------
cat("=== ALSFRS-R slope reconciliation: portal vs df_nfl_uniq ===\n")
if (exists("df_nfl_uniq") && "ALSFRS_R_PROGRESSION_SLOPE" %in% names(df_nfl_uniq)) {
  own <- as.data.table(df_nfl_uniq)[, .(sample_id,
                                        own_slope = ALSFRS_R_PROGRESSION_SLOPE)]
  own[, Participant_ID := sub("-[0-9]+-T$", "", sample_id)]
  rec <- merge(own[, .(Participant_ID, own_slope)],
               portsel[, .(Participant_ID, portal_slope = ALSFRS_R_PROGRESSION_SLOPE)],
               by = "Participant_ID")
  rec <- rec[is.finite(own_slope) & is.finite(portal_slope)]
  cat(sprintf("  shared patients with both slopes: %d\n", nrow(rec)))
  if (nrow(rec) >= 10) {
    cat(sprintf("  correlation: r = %.3f\n", cor(rec$own_slope, rec$portal_slope)))
    cat(sprintf("  mean own=%.3f  mean portal=%.3f\n",
                mean(rec$own_slope), mean(rec$portal_slope)))
    cat("  -> r>0.9 = same definition, safe to use portal (more coverage)\n")
    cat("     r 0.6-0.9 = related but different windowing; report cautiously\n")
    cat("     r<0.6 = different metric; do NOT swap, investigate\n")
  }
} else {
  cat("  df_nfl_uniq not available; skipping reconciliation\n")
}

res_ALSFRSslope_portal <- fit_linear_endpoint(dat, "ALSFRS_R_PROGRESSION_SLOPE",
                                              "ALSFRS-R slope PORTAL (/mo)")

## ---------------------------------------------------------------------------
## (R) Revised El Escorial certainty (full coverage)
## ---------------------------------------------------------------------------
res_EL <- fit_linear_endpoint(dat, "elescr_portal", "El Escorial certainty PORTAL (1-5)")

## ---------------------------------------------------------------------------
## (S) Number of visits -- dropout / retention check (not a disease endpoint)
## ---------------------------------------------------------------------------
res_VISITS <- fit_linear_endpoint(dat, "Number_of_Visits", "Number of visits (retention)")

## ---------------------------------------------------------------------------
## ASSEMBLE
## ---------------------------------------------------------------------------
linear_v5 <- rbindlist(list(res_ALSFRSbase, res_CBSbase,
                            res_ALSFRSslope_portal, res_EL, res_VISITS),
                       use.names = TRUE)
linear_v5[, p_S3_FDR := p.adjust(p_S3, method = "BH")]

cat("\n============== LINEAR ENDPOINT SCREEN v5 (portal master table) ==============\n")
print(linear_v5)

cat("\n--- coverage (n in anchor) ---\n")
for (v in c("ALSFRS_R_Baseline_Value","CBS_Baseline_Value",
            "ALSFRS_R_PROGRESSION_SLOPE","elescr_portal","Number_of_Visits"))
  cat(sprintf("  %-28s : %d\n", v, nrow(dat[is.finite(get(v))])))

cat("\n--- means by subtype (interpretation aid) ---\n")
print(dat[is.finite(ALSFRS_R_Baseline_Value),
          .(mean_ALSFRS_base = round(mean(ALSFRS_R_Baseline_Value),2), n=.N), by=subtype])
print(dat[is.finite(CBS_Baseline_Value),
          .(mean_CBS_base = round(mean(CBS_Baseline_Value),2), n=.N), by=subtype])
print(dat[is.finite(Number_of_Visits),
          .(mean_visits = round(mean(Number_of_Visits),2), n=.N), by=subtype])

saveRDS(linear_v5, "clinical_endpoint_screen_v5_results.rds")
cat("\nDone. Saved clinical_endpoint_screen_v5_results.rds\n")
################################################################################




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