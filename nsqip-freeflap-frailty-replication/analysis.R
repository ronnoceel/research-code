# analysis.R -- Age vs. frailty after microsurgical free-flap reconstruction
# ===========================================================================
# Targeted replication of:
#   Jubbal KT, Zavlin D, Suliman A. "The effect of age on microsurgical free flap
#   outcomes: An analysis of 5,951 cases." Microsurgery 2017;37(8):858-864.
#   PMID 28573680. DOI 10.1002/micr.30189.
#
# Data: ACS-NSQIP Participant Use Files (PUF) 2022-2024. Free tissue transfers by
#   primary CPT. Case-level PUF data are NOT included (ACS-NSQIP Data Use Agreement
#   prohibits redistribution); place the tab-delimited PUFs where nsqip_load.R's
#   nsqip_path() expects them, or point CACHE at a prebuilt column cache.
#
# Pipeline (all results are OUTCOMES, reported as-is; never tuned to a target):
#   1. Cohort, descriptives, complication incidence, univariate age gradient.
#   2. Multivariable models: is age independent after adjustment? (per year + strata)
#   3. ASA>=3 prevalence audit (overall / by year / by flap type).
#   4. mFI-5 frailty: nested Model A (no frailty) vs B (+mFI-5); age log-odds
#      attenuation; discrimination (AUC) and AIC incl. a frailty-only (age-dropped)
#      model; categorical mFI-5; head-to-head (mFI-5 vs age/BMI/ASA/smoking).
#   5. Preoperative hypoalbuminemia: multiple imputation (mice, m=50; White's rule
#      for ~47% missingness) as the PRIMARY combined model, a missing-indicator
#      model and a complete-case model as robustness/sensitivity, and the age
#      log-odds attenuation after mFI-5 (reported as attenuation only -- a raw
#      log-odds difference is not a mediation metric because of odds-ratio
#      non-collapsibility; Schuster 2021, BMC Med Res Methodol 21:136).
#   6. 30-day mortality (~27 events): Firth's penalized logistic regression
#      (logistf), since events-per-variable is too low for conventional MLE.
#
# Run (from this directory):  Rscript analysis.R
# Requires: R (>=4.2) with tidyverse, broom, mice, logistf.
# ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(broom); library(mice); library(logistf)
})
source("nsqip_load.R")
dir.create("outputs", showWarnings = FALSE)

# AUC (c-statistic) via the Mann-Whitney rank statistic (no external dependency)
auc_manual <- function(y, p) {
  ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# --- Columns needed across the three PUF years -----------------------------
NEED <- c("Age","SEX","CPT","HEIGHT","WEIGHT","ASACLAS","SMOKE","DIABETES","HYPERMED",
          "STEROID","BLEEDDIS","FNSTATUS2","OPTIME","TOTHLOS","DOpertoD",
          "SUPINFEC","WNDINFD","ORGSPCSSI","DEHIS","RETURNOR","REOPERATION1",
          "OUPNEUMO","REINTUB","FAILWEAN","PULEMBOL","OTHDVT","CDMI","CDARREST",
          "RENAFAIL","RENAINSF","URNINFEC","OTHSYSEP",
          "OTHBLEED","TRANSFUS",
          "HXCHF","HXCOPD","PRALBUM")
CACHE <- file.path("data","raw_cache.rds")   # optional prebuilt cache (git-ignored)
if (file.exists(CACHE)) {
  raw <- readRDS(CACHE)
} else {
  raw <- load_nsqip(NEED); dir.create("data", showWarnings = FALSE); saveRDS(raw, CACHE)
}

# occurrence fields null-out as "No Complication"; a few (e.g. RENAFAIL) use "No"
comp    <- function(x) as.integer(!is.na(x) & !(x %in% c("No Complication","No")))
FLAP_CPT <- c("15756","15757","15758","15842","15845",
              "19364","20955","20962","20969","43496")

# ===========================================================================
# Cohort construction
# ===========================================================================
dat <- raw %>%
  filter(CPT %in% FLAP_CPT) %>%
  mutate(
    age    = suppressWarnings(as.integer(if_else(Age == "90+", "90", Age))),
    female = SEX == "female",
    ht     = suppressWarnings(as.numeric(HEIGHT)),
    wt     = suppressWarnings(as.numeric(WEIGHT)),
    ht     = if_else(ht > 0, ht, NA_real_),
    wt     = if_else(wt > 0, wt, NA_real_),
    bmi    = 703 * wt / (ht^2),                          # inches/lb -> kg/m^2
    bmi    = if_else(bmi >= 12 & bmi <= 80, bmi, NA_real_),
    flap_site = case_when(
      CPT == "19364"                      ~ "Breast",
      CPT %in% c("15842","15845")         ~ "Facial",
      CPT %in% c("20955","20962","20969") ~ "Bone/osteocutaneous",
      CPT == "43496"                      ~ "Jejunum",
      TRUE                                ~ "Soft-tissue (muscle/skin/fascial)"),
    asa3     = as.integer(ASACLAS %in% c("3-Severe Disturb","4-Life Threat","5-Moribund")),
    bleeddis = as.integer(BLEEDDIS == "Yes"),
    diab     = as.integer(DIABETES %in% c("INSULIN","NON-INSULIN")),
    htn      = as.integer(HYPERMED == "Yes"),
    smoke    = as.integer(SMOKE == "Yes"),
    dep      = as.integer(FNSTATUS2 %in% c("Partially Dependent","Totally Dependent")),
    optime   = suppressWarnings(as.numeric(OPTIME)),
    optime   = if_else(optime > 0, optime, NA_real_),
    los      = suppressWarnings(as.numeric(TOTHLOS)),
    dopertod = suppressWarnings(as.numeric(DOpertoD)),
    # occurrence composites (bleeding excluded from the primary surgical composite)
    reop     = as.integer((!is.na(RETURNOR) & RETURNOR == "Yes") |
                          (!is.na(REOPERATION1) & REOPERATION1 == "Yes")),
    surg     = as.integer(comp(SUPINFEC) | comp(WNDINFD) | comp(ORGSPCSSI) |
                          comp(DEHIS) | reop),
    medical  = as.integer(comp(OUPNEUMO) | comp(REINTUB) | comp(FAILWEAN) | comp(PULEMBOL) |
                          comp(OTHDVT) | comp(CDMI) | comp(CDARREST) | comp(RENAFAIL) |
                          comp(RENAINSF) | comp(URNINFEC) | comp(OTHSYSEP)),
    bleed    = as.integer(comp(OTHBLEED) | (!is.na(TRANSFUS) & TRANSFUS == "Yes")),
    mort30   = as.integer(!is.na(dopertod) & dopertod >= 0 & dopertod <= 30),
    # flap failure has no field in the modern PUF -> proxied by return-to-OR
    flapfail_proxy = reop,
    # mFI-5 components
    mfi_chf  = as.integer(HXCHF == "Yes"),
    mfi_dm   = as.integer(DIABETES %in% c("INSULIN","NON-INSULIN")),
    mfi_copd = as.integer(HXCOPD == "Yes"),
    mfi_htn  = as.integer(HYPERMED == "Yes"),
    mfi_dep  = as.integer(FNSTATUS2 %in% c("Partially Dependent","Totally Dependent")),
    m_chf_ok  = !is.na(HXCHF)    & HXCHF    %in% c("Yes","No"),
    m_copd_ok = !is.na(HXCOPD)   & HXCOPD   %in% c("Yes","No"),
    m_htn_ok  = !is.na(HYPERMED) & HYPERMED %in% c("Yes","No"),
    m_dm_ok   = !is.na(DIABETES) & DIABETES %in% c("INSULIN","NON-INSULIN","NO"),
    m_dep_ok  = !is.na(FNSTATUS2)& FNSTATUS2 %in% c("Independent","Partially Dependent","Totally Dependent"),
    # preoperative albumin (-99 sentinel / implausible -> missing)
    albumin  = suppressWarnings(as.numeric(PRALBUM)),
    albumin  = if_else(albumin > 0 & albumin < 8, albumin, NA_real_),
    hypoalb  = if_else(is.na(albumin), NA_integer_, as.integer(albumin < 3.5))
  ) %>%
  filter(!is.na(age), age >= 18)

# prolonged operative time = top decile within cohort
op_thresh   <- quantile(dat$optime, 0.90, na.rm = TRUE)
dat$prolong <- as.integer(dat$optime >= op_thresh)
dat$surg_incl_bld <- as.integer(dat$surg | dat$bleed)
dat$bmi5    <- dat$bmi / 5
dat$age_grp <- relevel(factor(cut(dat$age, c(-Inf,50,60,70,80,Inf),
                     labels = c("18-49","50-59","60-69","70-79","80+"), right = FALSE)), ref = "18-49")

# mFI-5 = sum of 5 components (require all 5 codeable)
dat$mfi_ok <- dat$m_chf_ok & dat$m_copd_ok & dat$m_htn_ok & dat$m_dm_ok & dat$m_dep_ok
dat$mfi5   <- with(dat, ifelse(mfi_ok, mfi_chf + mfi_dm + mfi_copd + mfi_htn + mfi_dep, NA_integer_))

cat(sprintf("Cohort n = %d free tissue transfers (NSQIP 2022-2024)\n", nrow(dat)))
cat(sprintf("  mean age %.1f; %.1f%% female; median BMI %.1f; %.1f%% ASA>=3; albumin missing %.1f%%\n",
            mean(dat$age), 100*mean(dat$female), median(dat$bmi, na.rm=TRUE),
            100*mean(dat$asa3), 100*mean(is.na(dat$albumin))))

# ===========================================================================
# 1. Descriptives, incidence, univariate age gradient
# ===========================================================================
site_tab <- dat %>% count(flap_site) %>% mutate(pct = 100*n/sum(n)) %>% arrange(desc(n))
write_csv(site_tab, "outputs/flap_site_composition.csv")

inc <- tibble(
  outcome = c("Surgical (excl bleeding)","Surgical (incl bleeding)","Medical",
              "Reoperation/return-to-OR","30-day mortality"),
  pct = 100*c(mean(dat$surg), mean(dat$surg_incl_bld), mean(dat$medical),
              mean(dat$reop), mean(dat$mort30)))
write_csv(inc, "outputs/incidence_overall.csv")

by_age <- dat %>% group_by(age_grp) %>%
  summarise(n=n(), surgical=100*mean(surg), medical=100*mean(medical),
            reoperation=100*mean(reop), mortality=100*mean(mort30),
            asa3=100*mean(asa3), bmi=median(bmi,na.rm=TRUE),
            optime=median(optime,na.rm=TRUE), los=median(los,na.rm=TRUE), .groups="drop")
write_csv(by_age, "outputs/complications_by_age_stratum.csv")

trend_p <- function(y) {
  tb <- dat %>% mutate(y = .data[[y]]) %>% group_by(age_grp) %>%
    summarise(x = sum(y), n = n(), .groups = "drop")
  if (sum(tb$x) < 2 || all(tb$x == 0)) return(NA_real_)
  tryCatch(suppressWarnings(prop.trend.test(tb$x, tb$n)$p.value), error = function(e) NA_real_)
}
trend <- tibble(outcome = c("surg","medical","reop","mort30"),
                p_trend = map_dbl(c("surg","medical","reop","mort30"), trend_p))
write_csv(trend, "outputs/univariate_age_trend_pvalues.csv")

# ===========================================================================
# 2. Multivariable models -- is AGE independent after adjustment?
#    Full comorbidity adjustment (the base reproduction of Jubbal).
# ===========================================================================
RHS_yr  <- "age + female + bmi5 + asa3 + diab + smoke + htn + dep + prolong"
RHS_grp <- "age_grp + female + bmi5 + asa3 + diab + smoke + htn + dep + prolong"

fit_or <- function(rhs, outcome, label, data = dat) {   # Wald CIs (large-n)
  m <- glm(reformulate(rhs, outcome), binomial(), data = data)
  broom::tidy(m) %>% filter(term != "(Intercept)") %>%
    transmute(model = label, term, OR = exp(estimate),
              lo = exp(estimate - 1.96*std.error), hi = exp(estimate + 1.96*std.error),
              beta = estimate, p = p.value)
}

models_yr <- bind_rows(
  fit_or(RHS_yr, "surg",    "Surgical"),
  fit_or(RHS_yr, "medical", "Medical"),
  fit_or(RHS_yr, "reop",    "Reoperation (flap-fail proxy)"),
  fit_or(RHS_yr, "mort30",  "30-day mortality"))
write_csv(models_yr, "outputs/adjusted_or_age_per_year.csv")

models_grp <- bind_rows(
  fit_or(RHS_grp, "surg",    "Surgical"),
  fit_or(RHS_grp, "medical", "Medical"),
  fit_or(RHS_grp, "reop",    "Reoperation (flap-fail proxy)"),
  fit_or(RHS_grp, "mort30",  "30-day mortality"))
write_csv(models_grp, "outputs/adjusted_or_age_strata.csv")

verdict <- models_yr %>% filter(term %in% c("age","asa3","bmi5","prolong")) %>%
  mutate(result = ifelse(lo>1 & p<0.05, "sig +", ifelse(hi<1 & p<0.05, "sig -", "ns")))
write_csv(verdict, "outputs/reproduction_verdict.csv")

# ===========================================================================
# 3. ASA>=3 audit (overall / by year / by flap type) + ASACLAS distribution
# ===========================================================================
asa_dist   <- dat %>% count(ASACLAS) %>% mutate(pct = 100*n/sum(n)) %>% arrange(desc(n))
asa_overall<- tibble(stratum="Overall", n=nrow(dat), asa3_pct=100*mean(dat$asa3))
asa_byyear <- dat %>% group_by(PUFYEAR) %>% summarise(n=n(), asa3_pct=100*mean(asa3), .groups="drop") %>%
  transmute(stratum=paste0("PUF ",PUFYEAR), n, asa3_pct)
asa_bytype <- dat %>% mutate(type=if_else(CPT=="19364","Breast (19364)","Non-breast")) %>%
  group_by(type) %>% summarise(n=n(), asa3_pct=100*mean(asa3), .groups="drop") %>%
  transmute(stratum=type, n, asa3_pct)
asa_audit_out <- bind_rows(asa_overall, asa_byyear, asa_bytype,
  asa_dist %>% transmute(stratum=paste0("ASACLAS: ",ASACLAS), n, asa3_pct=pct))
write_csv(asa_audit_out, "outputs/asa_audit.csv")

# ===========================================================================
# 4. mFI-5 frailty sensitivity
# ===========================================================================
n_excl_mfi <- sum(is.na(dat$mfi5))
mfi_dist   <- dat %>% filter(!is.na(mfi5)) %>% count(mfi5) %>% mutate(pct=100*n/sum(n))
frail_prev <- dat %>% filter(!is.na(mfi5)) %>% summarise(frail_ge2_pct=100*mean(mfi5>=2)) %>% pull()
mfi_dist_out <- bind_rows(
  mfi_dist %>% transmute(metric=paste0("mFI-5 = ",mfi5), n, pct),
  tibble(metric="Frail (mFI-5 >= 2), %", n=NA_integer_, pct=frail_prev),
  tibble(metric="Excluded (incomplete components), n", n=n_excl_mfi, pct=100*n_excl_mfi/nrow(dat)),
  tibble(metric="mean mFI-5 (complete cases)", n=NA_integer_, pct=mean(dat$mfi5,na.rm=TRUE)))
write_csv(mfi_dist_out, "outputs/mfi5_distribution.csv")

# analysis set: complete mFI-5 + non-frailty covariates for Models A/B
mset <- dat %>% filter(!is.na(mfi5), !is.na(bmi), !is.na(prolong), !is.na(asa3))
mset$mfi5_cat <- relevel(factor(cut(mset$mfi5, c(-Inf,0.5,1.5,Inf), labels=c("0","1",">=2"))), ref="0")

# Model A: age + non-frailty covariates (diabetes/HTN/dependence enter only via mFI-5)
# Model B: A + mFI-5
RHS_A_yr  <- "age + female + bmi5 + asa3 + smoke + prolong"
RHS_B_yr  <- "age + female + bmi5 + asa3 + smoke + prolong + mfi5"
RHS_A_grp <- "age_grp + female + bmi5 + asa3 + smoke + prolong"
RHS_B_grp <- "age_grp + female + bmi5 + asa3 + smoke + prolong + mfi5"
OUTS   <- c("surg","medical","reop","mort30")
OUTLAB <- c(surg="Surgical", medical="Medical",
            reop="Reoperation (flap-fail proxy)", mort30="30-day mortality")

nested <- list(); fitrows <- list()
for (o in OUTS) {
  a_yr <- fit_or(RHS_A_yr, o, paste0(OUTLAB[o],": A (age/yr, no frailty)"), mset)
  b_yr <- fit_or(RHS_B_yr, o, paste0(OUTLAB[o],": B (+mFI-5)"), mset)
  a_gp <- fit_or(RHS_A_grp, o, paste0(OUTLAB[o],": A strata (no frailty)"), mset)
  b_gp <- fit_or(RHS_B_grp, o, paste0(OUTLAB[o],": B strata (+mFI-5)"), mset)
  nested[[o]] <- bind_rows(a_yr, b_yr, a_gp, b_gp) %>% mutate(outcome=OUTLAB[o])

  bA <- a_yr$beta[a_yr$term=="age"]; bB <- b_yr$beta[b_yr$term=="age"]
  bA80 <- a_gp$beta[a_gp$term=="age_grp80+"]; bB80 <- b_gp$beta[b_gp$term=="age_grp80+"]
  mA <- glm(reformulate(RHS_A_yr, o), binomial(), data=mset)
  mB <- glm(reformulate(RHS_B_yr, o), binomial(), data=mset)
  mF <- glm(reformulate("female + bmi5 + asa3 + smoke + prolong + mfi5", o), binomial(), data=mset)
  aucf <- function(m) auc_manual(mset[[o]], predict(m, type="response"))
  fitrows[[o]] <- tibble(
    outcome=OUTLAB[o],
    age_OR_A=exp(bA), age_OR_B=exp(bB), age_pct_atten=100*(bA-bB)/bA,
    age80_OR_A=exp(bA80), age80_OR_B=exp(bB80), age80_pct_atten=100*(bA80-bB80)/bA80,
    age_p_A=a_yr$p[a_yr$term=="age"], age_p_B=b_yr$p[b_yr$term=="age"],
    AUC_A=aucf(mA), AUC_B=aucf(mB), AUC_frailtyOnly=aucf(mF),
    AIC_A=AIC(mA), AIC_B=AIC(mB), AIC_frailtyOnly=AIC(mF))
}
write_csv(bind_rows(nested),  "outputs/mfi5_nested_models.csv")
write_csv(bind_rows(fitrows), "outputs/mfi5_model_fit.csv")

nested_cat <- list()
for (o in OUTS)
  nested_cat[[o]] <- fit_or("age + female + bmi5 + asa3 + smoke + prolong + mfi5_cat", o,
                            paste0(OUTLAB[o],": B (+mFI-5 categorical)"), mset) %>% mutate(outcome=OUTLAB[o])
write_csv(bind_rows(nested_cat), "outputs/mfi5_nested_categorical.csv")

h2h <- list()
for (o in OUTS)
  h2h[[o]] <- fit_or("mfi5 + age + bmi5 + asa3 + smoke", o, OUTLAB[o], mset) %>% mutate(outcome=OUTLAB[o])
write_csv(bind_rows(h2h), "outputs/mfi5_headtohead.csv")

# ===========================================================================
# 5. Frailty + preoperative hypoalbuminemia
#    Base set: complete mFI-5 + covariates; only albumin is missing -> imputed.
# ===========================================================================
base <- dat %>%
  filter(!is.na(mfi5), !is.na(bmi5), !is.na(prolong), !is.na(asa3)) %>%
  transmute(albumin, age, female, bmi5, asa3, smoke, prolong, mfi5, surg, medical)
pct_miss <- 100*mean(is.na(base$albumin))
cat(sprintf("\nAlbumin-combined base set n = %d; albumin missing %.1f%% (White's rule -> m ~= %d)\n",
            nrow(base), pct_miss, ceiling(pct_miss)))

# (5a) PRIMARY: multiple imputation of albumin (mice, m=50, PMM, outcomes in model)
M <- 50
meth <- make.method(base); meth[] <- ""; meth["albumin"] <- "pmm"
pred <- make.predictorMatrix(base); pred[,] <- 0
pred["albumin", c("age","female","bmi5","asa3","smoke","prolong","mfi5","surg","medical")] <- 1
imp  <- mice(base, m=M, method=meth, predictorMatrix=pred, maxit=10, seed=20260912, printFlag=FALSE)
long <- complete(imp, "long", include=TRUE)
long$hypoalb <- as.integer(long$albumin < 3.5)
long$age_grp <- relevel(factor(cut(long$age, c(-Inf,50,60,70,80,Inf),
                     labels=c("18-49","50-59","60-69","70-79","80+"), right=FALSE)), ref="18-49")
imp2 <- as.mids(long)

pool_or <- function(mids_obj, expr_glm, outcome_label, analysis_label) {
  fit <- with(mids_obj, eval(expr_glm))
  ps  <- summary(pool(fit), conf.int=TRUE)
  tibble(analysis=analysis_label, outcome=outcome_label, term=as.character(ps$term),
         OR=exp(ps$estimate), lo=exp(ps$`2.5 %`), hi=exp(ps$`97.5 %`),
         beta=ps$estimate, p=ps$p.value)
}
mi_rows <- bind_rows(
  pool_or(imp2, quote(glm(surg    ~ age + mfi5 + hypoalb + female + bmi5 + asa3 + smoke + prolong, family=binomial)), "Surgical", "MI_m50_continuous_age"),
  pool_or(imp2, quote(glm(medical ~ age + mfi5 + hypoalb + female + bmi5 + asa3 + smoke + prolong, family=binomial)), "Medical",  "MI_m50_continuous_age"),
  pool_or(imp2, quote(glm(surg    ~ age_grp + mfi5 + hypoalb + female + bmi5 + asa3 + smoke + prolong, family=binomial)), "Surgical", "MI_m50_strata_age"),
  pool_or(imp2, quote(glm(medical ~ age_grp + mfi5 + hypoalb + female + bmi5 + asa3 + smoke + prolong, family=binomial)), "Medical",  "MI_m50_strata_age"))
write_csv(mi_rows, "outputs/albumin_mi_m50.csv")

# (5b) ROBUSTNESS: missing-albumin indicator (normal / hypoalbuminemic / missing)
ind <- base %>%
  mutate(alb3 = relevel(factor(case_when(is.na(albumin) ~ "missing",
                                         albumin < 3.5  ~ "hypoalbuminemic",
                                         TRUE           ~ "normal")), ref="normal"),
         age_grp = relevel(factor(cut(age, c(-Inf,50,60,70,80,Inf),
                     labels=c("18-49","50-59","60-69","70-79","80+"), right=FALSE)), ref="18-49"))
ind_or <- function(f, outcome_label, analysis_label) {
  m <- glm(as.formula(f), binomial(), data=ind)
  broom::tidy(m) %>% filter(term != "(Intercept)") %>%
    transmute(analysis=analysis_label, outcome=outcome_label, term,
              OR=exp(estimate), lo=exp(estimate-1.96*std.error), hi=exp(estimate+1.96*std.error),
              beta=estimate, p=p.value)
}
ind_rows <- bind_rows(
  ind_or("surg    ~ age + mfi5 + alb3 + female + bmi5 + asa3 + smoke + prolong", "Surgical", "MissingIndicator_continuous_age"),
  ind_or("medical ~ age + mfi5 + alb3 + female + bmi5 + asa3 + smoke + prolong", "Medical",  "MissingIndicator_continuous_age"),
  ind_or("surg    ~ age_grp + mfi5 + alb3 + female + bmi5 + asa3 + smoke + prolong", "Surgical", "MissingIndicator_strata_age"),
  ind_or("medical ~ age_grp + mfi5 + alb3 + female + bmi5 + asa3 + smoke + prolong", "Medical",  "MissingIndicator_strata_age"))
write_csv(ind_rows, "outputs/albumin_missing_indicator.csv")

# (5c) SENSITIVITY: complete-case combined model
cc <- base %>% filter(!is.na(albumin)) %>%
  mutate(hypoalb = as.integer(albumin < 3.5),
         age_grp = relevel(factor(cut(age, c(-Inf,50,60,70,80,Inf),
                     labels=c("18-49","50-59","60-69","70-79","80+"), right=FALSE)), ref="18-49"))
cc_or <- function(f, outcome_label, analysis_label) {
  m <- glm(as.formula(f), binomial(), data=cc)
  broom::tidy(m) %>% filter(term != "(Intercept)") %>%
    transmute(analysis=analysis_label, outcome=outcome_label, term,
              OR=exp(estimate), lo=exp(estimate-1.96*std.error), hi=exp(estimate+1.96*std.error),
              beta=estimate, p=p.value)
}
cc_rows <- bind_rows(
  cc_or("surg    ~ age + mfi5 + hypoalb + female + bmi5 + asa3 + smoke + prolong", "Surgical", "CompleteCase_continuous_age"),
  cc_or("medical ~ age + mfi5 + hypoalb + female + bmi5 + asa3 + smoke + prolong", "Medical",  "CompleteCase_continuous_age"))
write_csv(cc_rows, "outputs/albumin_complete_case.csv")

# (5d) Age log-odds ATTENUATION after adding mFI-5 (attenuation only; not mediation)
att_base <- dat %>% filter(!is.na(mfi5), !is.na(bmi5), !is.na(prolong), !is.na(asa3))
attenuation <- function(outcome) {
  mA <- glm(reformulate(c("age","female","bmi5","asa3","smoke","prolong"), outcome), binomial(), att_base)
  mB <- glm(reformulate(c("age","female","bmi5","asa3","smoke","prolong","mfi5"), outcome), binomial(), att_base)
  bA <- coef(mA)["age"]; bB <- coef(mB)["age"]
  tibble(outcome=outcome, OR_age_A=exp(bA), OR_age_B=exp(bB), logodds_attenuation_pct=100*(bA-bB)/bA)
}
write_csv(bind_rows(attenuation("medical"), attenuation("surg")), "outputs/albumin_age_attenuation.csv")

# ===========================================================================
# 6. 30-day mortality -- Firth-penalized logistic regression (few events)
# ===========================================================================
n_ev <- sum(mset$mort30)
firth_tidy <- function(form, label) {
  m <- logistf(form, data = mset)
  tibble(model=label, term=names(coef(m)), OR=exp(coef(m)),
         lo=exp(m$ci.lower), hi=exp(m$ci.upper), beta=as.numeric(coef(m)), p=m$prob) %>%
    filter(term != "(Intercept)")
}
firth_all <- bind_rows(
  firth_tidy(mort30 ~ age + female + bmi5 + asa3 + smoke + prolong + mfi5,
             "Firth 30-day mortality: age per year (+mFI-5)"),
  firth_tidy(mort30 ~ age_grp + female + bmi5 + asa3 + smoke + prolong + mfi5,
             "Firth 30-day mortality: age strata (+mFI-5)")) %>%
  mutate(events = n_ev, n = nrow(mset))
write_csv(firth_all, "outputs/mortality_firth.csv")
cat(sprintf("\n30-day deaths = %d of %d (%.3f%%); mortality modeled with Firth penalization.\n",
            n_ev, nrow(mset), 100*mean(mset$mort30)))

cat("\nDone. Outputs written to outputs/.\n")
