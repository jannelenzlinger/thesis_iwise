# =============================================================================
# 07_models.R
# First modelling attempt: multilevel ordinal logistic regression
#   Outcomes : FLI_3item (d.fl) first, then INCOME_5 (d.inc)
#   M0       : unadjusted     = iwisescore
#   M1       : fully adjusted = rhs_ctx (individual + country-level covariates)
#   Weights  : w_fit
#   Section 3, step by step:
#     3a  Mundlak check (country-year mean needed?)    random intercept only
#     3b  Random-intercept models (without the mean)
#     3c  Does the centring work?                      random intercept only
#     3d  Random IWISE slope: fit and compare with 3b
#     3e  Final models (switch: use_slope) -> used by all checks
# Then checks 1-4: proportional odds, linearity, influence, goodness of fit.
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
library(tidyverse)                   # data handling (dplyr, purrr, stringr) and ggplot2
library(ordinal)                     # clm() and clmm(): cumulative link (mixed) models
library(MASS, exclude = "select")    # polr(), needed by brant(); keeps dplyr's select()
library(brant)                       # Brant test of proportional odds
library(splines)                     # ns(): natural splines for the linearity check
library(generalhoslem)               # lipsitz.test(), pulkrob.chisq(), pulkrob.deviance()


# ---- 1. DATA ----------------------------------------------------------------
list2env(readRDS("data/iwise_data_prepared.rds"), envir = globalenv())  # restores d.fl, d.inc, rhs_ctx, v.* objects

prep <- function(data, y) {                                # builds the analytic sample for one outcome
  data %>%
    filter(cc) %>%                                         # complete cases = exactly the rows the model uses
    mutate(.id          = row_number(),                    # row id, to trace influential observations later
           country_year = factor(.data[[v.cluster]]),      # grouping factor for the random intercept
           iwise_cat    = factor(iwise_cat, ordered = FALSE),  # unordered: avoids polynomial contrasts
           "{y}" := factor(.data[[y]], ordered = TRUE))    # clm()/clmm() need an ordered factor outcome
}

d.fl.cc  <- prep(d.fl,  "FLI_3item")                       # FLI analytic sample
d.inc.cc <- prep(d.inc, "INCOME_5")                        # income analytic sample

table(d.fl.cc$FLI_3item)                                   # check: 5 levels, lowest first
table(d.inc.cc$INCOME_5)                                   # check: 5 levels, poorest quintile first

v.watch <- c("iwisescore", "log_gdp_c",                    # terms to monitor in checks 3 and LOO:
             paste0(v.wgi.main, "_c"), v.jmp,              # IWISE + country-level covariates
             paste0(v.gbd.model, "_c"))                    # (centred names, see 05 Section 5.2)


# ---- 2. FORMULAS AND HELPERS ------------------------------------------------
# with-mean versions are built here explicitly for the check in Section 3a.
rhs_unadj      <- "iwisescore"                             # M0: IWISE only
rhs_unadj_mean <- "iwisescore + iwise_mean"                # M0 + country-year mean (Mundlak)
rhs_ctx_mean   <- paste(rhs_ctx, "+ iwise_mean")           # M1 + country-year mean (Mundlak)

re_term   <- "(1 | country_year)"                          # random intercept: one shift per country-year
re_slope  <- "(1 + iwisescore | country_year)"             # + random IWISE slope (correlated):
                                                           # each country-year its own IWISE effect

# the random-effects structure is an argument. Section 3 always
# states it explicitly (re = re_term / re = re_slope). Without it, the
# FINAL structure re_final is used (set in 3e), so all checks follow
# the decision made in 3d automatically.

make_f <- function(y, rhs, random = TRUE, re = re_final) { # turns text into a model formula
  as.formula(paste(y, "~", rhs,                            # outcome ~ predictors
                   if (random) paste("+", re) else ""))    # + random effects (or none)
}
#function helps for loops and makes sure the settings are the same everywhere
fit_clmm <- function(data, y, rhs, re = re_final) {        # same model everywhere, used in loops
  clmm(make_f(y, rhs, re = re), data = data, weights = w_fit,  # weighted multilevel ordinal logit
       link = "logit", Hess = TRUE, nAGQ = 1)              # logit link, keep Hessian, Laplace approx.
}                                                          # (random slopes need nAGQ = 1)

tidy_or <- function(m, label) {                            # odds-ratio table for the fixed effects
  co <- summary(m)$coefficients[names(m$beta), , drop = FALSE]  # betas only (thresholds dropped)
  tibble(model = label,                                    # model label
         term  = rownames(co),                             # predictor
         b     = co[, 1],                                  # log-odds coefficient
         se    = co[, 2],                                  # standard error
         OR    = exp(b),                                   # odds ratio
         lo    = exp(b - 1.96 * se),                       # lower 95% CI (Wald)
         hi    = exp(b + 1.96 * se),                       # upper 95% CI (Wald)
         p     = co[, 4])                                  # p-value
}

load_models <- function(file) {                            # restores every model in an .rds file
  obj <- readRDS(file)                                     # named list of models
  list2env(obj, envir = globalenv())                       # each name becomes an object again
  message("Loaded from ", file, ": ", paste(names(obj), collapse = ", "))
  invisible(names(obj))
}

# with a random slope, VarCorr() is a 2x2 matrix; [1, 1] is the
# intercept variance. The ICC then applies at iwisescore = 0 (no water
# insecurity), because between-country variance changes with IWISE.
icc <- function(m) {                                       # latent-scale ICC from the fitted model
  s2 <- VarCorr(m)[[1]][1, 1]                              # random-intercept variance
  s2 / (s2 + pi^2 / 3)                                     # pi^2/3 = level-1 variance of the logistic
}


# ---- RELOAD: use instead of refitting Section 3 ------------------------------
use_slope <- TRUE                                          # must match your decision in 3e
re_final  <- if (use_slope) re_slope else re_term          # default structure for fit_clmm()

load_models("models_mundlak.rds")                          # 3a: m.*.nm, m.*.wm
m.fl.0.ri  <- m.fl.0.nm                                    # 3b: instant, just renaming
m.fl.1.ri  <- m.fl.1.nm
m.inc.0.ri <- m.inc.0.nm
m.inc.1.ri <- m.inc.1.nm

load_models("models_random_slope.rds")                     # 3d: m.*.rs
load_models("models_final.rds")                            # 3e: m.fl.0, m.fl.1, m.inc.0, m.inc.1
load_models("models_random_slope.rds")                     # 3d: m.*.rs
load_models("models_slope_check.rds")                      # 3d-ii: m.fl.1.rs.chk
# load_models("models_unnormalised.rds")                   # Check 5, if needed
# load_models("loo_attempt1.rds")                          # Check 3c: loo.fl, loo.inc


# ---- 3a. IS THE COUNTRY-YEAR MEAN NEEDED? (MUNDLAK) -------------------------
# Random intercept only (re = re_term): the decision below is based on these.
# iwise_mean = contextual effect (between minus within). If it is ~0 and the
# iwisescore OR does not change, the mean can be dropped.
# Without the mean, iwisescore is a blend of within and between effects.
m.fl.0.nm  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_unadj,      re = re_term)  # FLI, unadjusted, without mean
m.fl.0.wm  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_unadj_mean, re = re_term)  # FLI, unadjusted, with mean
m.fl.1.nm  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx,        re = re_term)  # FLI, adjusted, without mean
m.fl.1.wm  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx_mean,   re = re_term)  # FLI, adjusted, with mean
m.inc.0.nm <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_unadj,      re = re_term)  # income, unadjusted, without mean
m.inc.0.wm <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_unadj_mean, re = re_term)  # income, unadjusted, with mean
m.inc.1.nm <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx,        re = re_term)  # income, adjusted, without mean
m.inc.1.wm <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx_mean,   re = re_term)  # income, adjusted, with mean

# 3a-i. IWISE odds ratios side by side
bind_rows(tidy_or(m.fl.0.wm,  "FLI M0 with mean"), tidy_or(m.fl.0.nm,  "FLI M0 without mean"),
          tidy_or(m.fl.1.wm,  "FLI M1 with mean"), tidy_or(m.fl.1.nm,  "FLI M1 without mean"),
          tidy_or(m.inc.0.wm, "INC M0 with mean"), tidy_or(m.inc.0.nm, "INC M0 without mean"),
          tidy_or(m.inc.1.wm, "INC M1 with mean"), tidy_or(m.inc.1.nm, "INC M1 without mean")) %>%
  filter(term %in% c("iwisescore", "iwise_mean")) %>%      # IWISE rows only
  select(model, term, OR, lo, hi, p) %>%                   # OR, 95% CI, p-value
  print()                                                  # similar iwisescore ORs = the mean barely matters

# 3a-ii. Likelihood-ratio test: does adding the mean improve the fit?
#        (Hedeker 2015, p. 9). Small p = within and between effects differ.
anova(m.fl.0.nm,  m.fl.0.wm)                               # FLI, unadjusted   (p = 0.001)
anova(m.fl.1.nm,  m.fl.1.wm)                               # FLI, adjusted     (p = 0.9346)
anova(m.inc.0.nm, m.inc.0.wm)                              # income, unadjusted (p < 0.001)
anova(m.inc.1.nm, m.inc.1.wm)                              # income, adjusted  (p = 0.013, dAIC = 4)

saveRDS(mget(c("m.fl.0.nm", "m.fl.0.wm", "m.fl.1.nm", "m.fl.1.wm",       # mget() returns a named list
               "m.inc.0.nm", "m.inc.0.wm", "m.inc.1.nm", "m.inc.1.wm")),
        "models_mundlak.rds")
# =============================================================================
# DECISION: EXCLUDE MEAN
# The iwisescore OR is unchanged in the adjusted models (FLI 0.967, income
# 0.982, with and without), fit improves not at all for FLI and only
# marginally for income. The with-mean models are reported as sensitivity.
# Made under the random-intercept model; the random-slope model with the
# mean did not converge (NaN SEs, implausible slope SD of 0.45).
# =============================================================================


# ---- 3b. RANDOM-INTERCEPT MODELS (WITHOUT THE MEAN) -------------------------
# These are the "without mean" models from 3a, reused instead of refitted.
m.fl.0.ri  <- m.fl.0.nm                                    # FLI, unadjusted
m.fl.1.ri  <- m.fl.1.nm                                    # FLI, adjusted
m.inc.0.ri <- m.inc.0.nm                                   # income, unadjusted
m.inc.1.ri <- m.inc.1.nm                                   # income, adjusted

# Estimation quality
sapply(list(fl0 = m.fl.0.ri, fl1 = m.fl.1.ri, inc0 = m.inc.0.ri, inc1 = m.inc.1.ri),
       function(m) c(condHess = summary(m)$condHess,       # < ~1e5 fine, > ~1e6 problematic
                     max_grad = max(abs(m$gradient))))     # close to 0 = converged
bind_rows(tidy_or(m.fl.1.ri, "FLI M1"), tidy_or(m.inc.1.ri, "INC M1")) %>%
  filter(is.na(se))                                        # 0 rows = all standard errors defined

sapply(list(fl0 = m.fl.0.ri, fl1 = m.fl.1.ri,              # ICC per model: share of variance
            inc0 = m.inc.0.ri, inc1 = m.inc.1.ri), icc)    # between country-years


# ---- 3c. DOES THE CENTRING WORK? (random intercept only) --------------------
# The models use centred covariates (05, Section 5.2). Here the adjusted
# models are refitted on the RAW covariates and compared.
# First run (uncentred FLI M1): condHess NaN, SEs undefined, country-level
# betas off (e.g. bas_rate_pp +0.068 vs -0.046 centred) = not converged.
# Centred: condHess 67,948 (FLI) and 52,164 (income), all SEs defined.
rhs_ctx_raw <- reduce(v.centre,                            # swap each centred name back to raw:
                      ~ str_replace_all(.x, paste0("\\b", .y, "_c\\b"), .y),
                      .init = rhs_ctx)                     # e.g. log_gdp_c -> log_gdp
rhs_ctx_raw                                                # check: no "_c" names left

m.fl.1.raw  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx_raw, re = re_term)  # FLI, uncentred
m.inc.1.raw <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx_raw, re = re_term)  # income, uncentred
# (warnings about the variance-covariance matrix are expected for the raw fits)

# 3c-i. Condition number: centred should be clearly lower (and not NaN)
sapply(list(fl_raw  = m.fl.1.raw,  fl_centred  = m.fl.1.ri,
            inc_raw = m.inc.1.raw, inc_centred = m.inc.1.ri),
       function(m) summary(m)$condHess)

# 3c-ii. Log-likelihood: centred equal or higher = raw fit had not fully converged
c(fl_raw  = logLik(m.fl.1.raw),  fl_centred  = logLik(m.fl.1.ri),
  inc_raw = logLik(m.inc.1.raw), inc_centred = logLik(m.inc.1.ri))

# 3c-iii. Betas side by side (same term order; row names from the raw model)
cbind(raw = m.fl.1.raw$beta,  centred = m.fl.1.ri$beta)    # FLI
cbind(raw = m.inc.1.raw$beta, centred = m.inc.1.ri$beta)   # income

saveRDS(mget(c("m.fl.1.raw", "m.inc.1.raw")), "models_uncentred.rds")

# ---- 3d. RANDOM IWISE SLOPE --------------------------------------------------
# Same models as 3b, plus a random IWISE slope (re = re_slope): each
# country-year has its own IWISE effect, spread around the average (fixed)
# effect with SD sigma_slope. Adds 2 parameters: slope variance and the
# intercept-slope correlation. Slower; may warn (see 3d-i).

# ---- 3d. RANDOM IWISE SLOPE --------------------------------------------------
# Question: does the IWISE effect differ between country-years?
# A random slope gives each country-year its own IWISE effect, spread around
# the average (fixed) effect. Compared with 3b it adds 2 parameters: the slope
# SD and the intercept-slope correlation.
#
#   3d-1  Fit the four slope models, or load them if already saved
#   3d-2  Did the fits work?      estimation quality + plausible values
#   3d-3  Are the fits stable?    same result with another optimiser?
#   3d-4  Is the slope needed?    boundary-corrected LR test + AIC   <- THE test
#   3d-5  What does it change?    average OR, spread, plot           (description)
#   3d-6  Decision                -> used in 3e
#
# 3d-2 and 3d-3 only check that the slope models can be trusted. Only 3d-4
# tests whether the slope is needed. 3d-5 describes the result.
# Needs the random-intercept models from 3b (m.*.ri) in memory.


# m.fl.0.rs  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_unadj, re = re_slope)  # FLI, unadjusted
# m.fl.1.rs  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx,   re = re_slope)  # FLI, adjusted
# m.inc.0.rs <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_unadj, re = re_slope)  # income, unadjusted
# m.inc.1.rs <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx,   re = re_slope)  # income, adjusted

# saveRDS(list(m.fl.0.rs = m.fl.0.rs, m.fl.1.rs = m.fl.1.rs,  # save: slow to refit
#             m.inc.0.rs = m.inc.0.rs, m.inc.1.rs = m.inc.1.rs),
#        "models_random_slope.rds")



stopifnot(exists("m.fl.0.ri"), exists("m.fl.1.ri"),        # 3b models must be loaded
          exists("m.inc.0.ri"), exists("m.inc.1.ri"))


# ---- 3d-1. FIT OR LOAD -------------------------------------------------------
# Loads the saved models if the file exists, then fits only what is missing.
# To force a refit: delete the .rds file and rm() the model objects.
f.rs <- "models_random_slope.rds"
if (file.exists(f.rs)) load_models(f.rs)

if (!exists("m.fl.0.rs"))  m.fl.0.rs  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_unadj, re = re_slope)
if (!exists("m.fl.1.rs"))  m.fl.1.rs  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx,   re = re_slope)
if (!exists("m.inc.0.rs")) m.inc.0.rs <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_unadj, re = re_slope)
if (!exists("m.inc.1.rs")) m.inc.1.rs <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx,   re = re_slope)

saveRDS(mget(c("m.fl.0.rs", "m.fl.1.rs", "m.inc.0.rs", "m.inc.1.rs")), f.rs)

# Each slope model next to its intercept-only version (used in 3d-2 to 3d-5)
pairs <- list(fl0  = list(ri = m.fl.0.ri,  rs = m.fl.0.rs),
              fl1  = list(ri = m.fl.1.ri,  rs = m.fl.1.rs),
              inc0 = list(ri = m.inc.0.ri, rs = m.inc.0.rs),
              inc1 = list(ri = m.inc.1.ri, rs = m.inc.1.rs))


# ---- 3d-2. DID THE FITS WORK? ------------------------------------------------
# One column per model. A fit can only be trusted if ALL rows look fine.
fit_quality <- function(p) {
  vc  <- VarCorr(p$rs)[[1]]                                # 2x2 random-effects matrix
  sds <- attr(vc, "stddev")
  c(condHess    = suppressWarnings(summary(p$rs)$condHess),  # NaN or > 1e6 = not identified
    max_grad    = max(abs(p$rs$gradient)),                 # large = may not have converged
    logLik_gain = as.numeric(logLik(p$rs) - logLik(p$ri)), # must be >= 0 (ri is nested in rs)
    sd_int_ri   = attr(VarCorr(p$ri)[[1]], "stddev")[[1]], # intercept SD without slope ...
    sd_int_rs   = sds[[1]],                                # ... and with: should be similar
    sd_slope    = sds[["iwisescore"]],                     # plausible: ~0.01-0.05 per point
    corr        = attr(vc, "correlation")[1, 2])           # intercept-slope correlation
}
round(sapply(pairs, fit_quality), 4)

# RESULT: fl0, fl1, inc1 fine (slope SD 0.019 FLI, 0.016 income).
#   inc0 FAILED: condHess NaN, gradient 158, slope SD 0.53, intercept SD
#   0.10 -> 0.56. Implausible; unadjusted income has almost no between-
#   country variance, so the slope cannot be identified. -> inc0 dropped below.


# ---- 3d-3. ARE THE FITS STABLE? ----------------------------------------------
# Refit the adjusted models with another optimiser (ucminf). If both land on
# the same logLik, betas and SDs, the original fit is a real maximum.
fit_ucminf <- function(data, y, rhs, re = re_slope) {      # same as fit_clmm(), other optimiser
  clmm(make_f(y, rhs, re = re), data = data, weights = w_fit, link = "logit",
       Hess = TRUE, nAGQ = 1, control = clmm.control(method = "ucminf"))
}

f.chk <- "models_slope_check.rds"                          # fit or load, as in 3d-1
if (file.exists(f.chk)) load_models(f.chk)
if (!exists("m.fl.1.rs.chk"))  m.fl.1.rs.chk  <- fit_ucminf(d.fl.cc,  "FLI_3item", rhs_ctx)
if (!exists("m.inc.1.rs.chk")) m.inc.1.rs.chk <- fit_ucminf(d.inc.cc, "INCOME_5",  rhs_ctx)
saveRDS(mget(c("m.fl.1.rs.chk", "m.inc.1.rs.chk")), f.chk)

cmp_fits <- function(m1, m2) {                             # how far apart are two fits?
  c(logLik_diff   = as.numeric(logLik(m2) - logLik(m1)),   # ~0 = same maximum
    max_beta_diff = max(abs(m1$beta - m2$beta)),           # < ~0.001 = same estimates
    max_sd_diff   = max(abs(attr(VarCorr(m1)[[1]], "stddev") -
                            attr(VarCorr(m2)[[1]], "stddev"))))  # < ~0.001 = same SDs
}
round(cbind(fl1  = cmp_fits(m.fl.1.rs,  m.fl.1.rs.chk),
            inc1 = cmp_fits(m.inc.1.rs, m.inc.1.rs.chk)), 5)

# RESULT: 0 difference, results stable


# ---- 3d-4. IS THE SLOPE NEEDED? -----------------------------------
# LR test: slope model vs intercept-only model (2 extra parameters).
# Under H0 the slope variance is 0, the edge of its possible values, so the
# usual chi2(2) p-value is too large. Corrected p = 50:50 mix of chi2(1) and
# chi2(2) (Snijders & Bosker). dAIC < 0 = slope model is better.
lr_slope <- function(p) {
  lr <- as.numeric(2 * (logLik(p$rs) - logLik(p$ri)))      # LR statistic
  c(LR      = lr,
    p_naive = pchisq(lr, df = 2, lower.tail = FALSE),      # standard test (too conservative)
    p_mix   = 0.5 * pchisq(lr, 1, lower.tail = FALSE) +    # boundary-corrected p-value
              0.5 * pchisq(lr, 2, lower.tail = FALSE),
    dAIC    = AIC(p$rs) - AIC(p$ri))                       # < 0 = slope better
}
round(sapply(pairs[c("fl0", "fl1", "inc1")], lr_slope), 4)  # inc0 left out: fit failed (3d-2)

# RESULT: random slope clearly improves fit in fl0, fl1, inc1
#   (LR 189-259, boundary-corrected p < 0.001, dAIC -183 to -255).
#   DECISION: keep the slope (use_slope <- TRUE); inc0 intercept only.


# ---- 3d-5. WHAT DOES THE SLOPE CHANGE? (DESCRIPTION, NO TEST) -----------------
# a. Average IWISE OR with and without the slope. CI usually wider with it.
imap_dfr(pairs[c("fl1", "inc1")], function(p, nm)
  bind_rows(tidy_or(p$ri, paste(nm, "intercept only")),
            tidy_or(p$rs, paste(nm, "+ slope")))) %>%
  filter(term == "iwisescore") %>%
  select(model, OR, lo, hi, p)
# RESULT: average IWISE OR stable with the slope (FLI 0.967 -> 0.964,
#   income 0.982 -> 0.981); CIs ~3x wider, as they now include between-
#   country variation in the effect.


# b. Spread: range in which 95% of the country-year IWISE ORs lie
slope_range <- function(m) {
  b  <- m$beta[["iwisescore"]]                             # average (fixed) IWISE effect
  sd <- attr(VarCorr(m)[[1]], "stddev")[["iwisescore"]]    # SD of the country-year slopes
  exp(c(average = b, lower95 = b - 1.96 * sd, upper95 = b + 1.96 * sd))  # ORs per IWISE point
}
round(sapply(list(fl1 = m.fl.1.rs, inc1 = m.inc.1.rs), slope_range), 3)
# RESULT: 95% of country-year IWISE ORs lie within 0.929-0.999 (FLI) and
#   0.950-1.012 (income). FLI: negative everywhere, strength varies.
#   Income: mostly negative, around zero in some country-years.


# c. Caterpillar plot: each country-year's own IWISE OR (average + deviation).
#    ranef() values are shrunken toward 0, so the plot slightly understates
#    the spread; b. is the better summary.
plot_cy_slopes <- function(m, data, title) {
  re <- ranef(m)$country_year                              # deviations per country-year
  cy <- tibble(country_year = rownames(re),
               OR = exp(m$beta[["iwisescore"]] + re[["iwisescore"]])) %>%
    left_join(data %>% distinct(country_year, countrynew) %>%   # add country names
                mutate(country_year = as.character(country_year)),
              by = "country_year")
  ggplot(cy, aes(OR, reorder(countrynew, OR))) +
    geom_vline(xintercept = exp(m$beta[["iwisescore"]]), colour = "red") +  # average OR
    geom_vline(xintercept = 1, colour = "grey60", linetype = 2) +          # no effect
    geom_point(size = 1) +
    labs(title = title, x = "OR per IWISE point", y = NULL) +
    theme_minimal() + theme(axis.text.y = element_text(size = 5))
}
plot_cy_slopes(m.fl.1.rs,  d.fl.cc,  "FLI: IWISE OR per country-year")
plot_cy_slopes(m.inc.1.rs, d.inc.cc, "Income: IWISE OR per country-year")


# ---- 3e. FINAL MODELS  --------
# Keep the random slope for fl0, fl1 and inc1, inc0 always stays intercept only (slope not identifiable).
use_slope <- TRUE                                          # FALSE = random intercept only
re_final  <- if (use_slope) re_slope else re_term          # used by fit_clmm() from here on

m.fl.0  <- if (use_slope) m.fl.0.rs  else m.fl.0.ri        # FLI, unadjusted
m.fl.1  <- if (use_slope) m.fl.1.rs  else m.fl.1.ri        # FLI, adjusted
m.inc.0 <- m.inc.0.ri                                      # income, unadjusted: always intercept only (3d-2)
m.inc.1 <- if (use_slope) m.inc.1.rs else m.inc.1.ri       # income, adjusted

saveRDS(list(m.fl.0 = m.fl.0, m.fl.1 = m.fl.1,             # save final fits
             m.inc.0 = m.inc.0, m.inc.1 = m.inc.1,
             m.fl.1.wm = m.fl.1.wm, m.inc.1.wm = m.inc.1.wm),  # + sensitivity: with the mean
        "models_final.rds")

summary(m.fl.1)                                            # full output
summary(m.inc.1)                                           # same for income

res <- bind_rows(tidy_or(m.fl.0,  "FLI M0"), tidy_or(m.fl.1,  "FLI M1"),   # all ORs in one table
                 tidy_or(m.inc.0, "INC M0"), tidy_or(m.inc.1, "INC M1"))
res %>% filter(term == "iwisescore") %>% print()           # IWISE rows

bind_rows(tidy_or(m.fl.1, "FLI M1"), tidy_or(m.inc.1, "INC M1")) %>%
  filter(term %in% v.watch[-1]) %>%                        # country-level ORs
  select(model, term, OR, lo, hi, p)

sapply(list(fl0 = m.fl.0, fl1 = m.fl.1,                    # ICC (with random slope: at IWISE = 0)
            inc0 = m.inc.0, inc1 = m.inc.1), icc)



# =============================================================================
# CHECKS
# brant(), nominal_test(), Cook's distance and the Lipsitz / Pulkstenis-
# Robinson tests do not accept clmm() objects. They are run on single-level
# stand-ins with the SAME fixed part (rhs_ctx), sample and (where possible)
# weights. Refits (checks 2 and 3) use the real multilevel model with the
# FINAL random-effects structure (re_final, 3e), except the iwise_cat
# comparison, which is intercept-only on both sides.
# =============================================================================

#each test gets a single-level stand-in: the same predictors, sample (and where possible) weights, just without the ransdom intercept
f.fl.fix  <- make_f("FLI_3item", rhs_ctx, random = FALSE)  # FLI formula without random intercept
f.inc.fix <- make_f("INCOME_5",  rhs_ctx, random = FALSE)  # income formula without random intercept

c.fl  <- clm(f.fl.fix,  data = d.fl.cc,  weights = w_fit, link = "logit")  # weighted clm stand-in, FLI
c.inc <- clm(f.inc.fix, data = d.inc.cc, weights = w_fit, link = "logit")  # weighted clm stand-in, income

p.fl  <- polr(f.fl.fix,  data = d.fl.cc,  method = "logistic", Hess = TRUE)  # polr stand-in (unweighted:
p.inc <- polr(f.inc.fix, data = d.inc.cc, method = "logistic", Hess = TRUE)  # brant/generalhoslem ignore weights)

saveRDS(mget(c("c.fl", "c.inc", "p.fl", "p.inc", "cuts.fl")), "models_standins.rds")

# load_models("models_standins.rds")                         # checks: c.*, p.*, cuts.*

# ---- CHECK 1. WEIGHT NORMALISATION: Sensivitity check -------------------------------
# w_fit = v.weight divided by its mean in the complete cases (05, Section 9),
# so the weights sum to the model N. Here the final adjusted models are
# refitted with the weight as it is (normalised over the full dataset, not
# the model sample) and compared with the final models from 3e.
# Expectation: clmm() reads weights as counts. Multiplying all weights by a
#   constant k acts like a sample of k * N, so SEs change by 1/sqrt(k). In a
#   multilevel model it also shifts the balance between the data and the
#   random effects, so ORs and random-effect SDs can move slightly.
# Uses the final structure (re_final). For speed, add re = re_term.
# logLik/AIC are NOT comparable here: the weights differ, so the scale differs.

v.weight                                                   # check: the weight w_fit was built from

# 3f-i. How far is the weight from normalised in each model sample?
w_scale <- function(data) {                                # one row per model sample
  data %>% summarise(n       = n(),                        # model N
                     sum_raw = sum(.data[[v.weight]]),     # sum of the weight as it is
                     sum_fit = sum(w_fit),                 # sum of w_fit (should equal n)
                     k       = sum_raw / n)                # scale factor: 1 = no difference
}
bind_rows(FLI = w_scale(d.fl.cc), income = w_scale(d.inc.cc), .id = "outcome")

# 3f-ii. Refit with the unnormalised weight
d.fl.un  <- d.fl.cc  %>% mutate(w_fit = .data[[v.weight]]) # copy with w_fit replaced by the raw weight
d.inc.un <- d.inc.cc %>% mutate(w_fit = .data[[v.weight]]) # (fit_clmm() always uses the column w_fit)

m.fl.1.un  <- fit_clmm(d.fl.un,  "FLI_3item", rhs_ctx)     # FLI, adjusted, unnormalised weight
m.inc.1.un <- fit_clmm(d.inc.un, "INCOME_5",  rhs_ctx)     # income, adjusted, unnormalised weight

saveRDS(list(m.fl.1.un = m.fl.1.un, m.inc.1.un = m.inc.1.un),  # save: slow to refit
        "models_unnormalised.rds")

# 3f-iii. ORs and SEs side by side
cmp_weights <- function(m.norm, m.un, data) {
  k <- sum(data[[v.weight]]) / nrow(data)                  # scale factor from 3f-i
  inner_join(tidy_or(m.norm, "norm") %>% select(term, OR_norm = OR, se_norm = se),  # normalised
             tidy_or(m.un,   "un")   %>% select(term, OR_un   = OR, se_un   = se),  # unnormalised
             by = "term") %>%
    mutate(OR_change_pct = 100 * (OR_un / OR_norm - 1),    # % change in the OR
           se_ratio      = se_un / se_norm,                # observed SE ratio
           se_expected   = 1 / sqrt(k))                    # expected ratio if only N changed
}
cmp.fl  <- cmp_weights(m.fl.1,  m.fl.1.un,  d.fl.cc)       # FLI
cmp.inc <- cmp_weights(m.inc.1, m.inc.1.un, d.inc.cc)      # income

cmp.fl  %>% filter(term %in% v.watch)                      # key terms, FLI
cmp.inc %>% filter(term %in% v.watch)                      # key terms, income

bind_rows(FLI = cmp.fl, income = cmp.inc, .id = "outcome") %>%   # overview across all terms
  group_by(outcome) %>%
  summarise(max_OR_change_pct = max(abs(OR_change_pct)),   # largest OR shift
            se_ratio_min      = min(se_ratio),             # range of SE ratios ...
            se_ratio_max      = max(se_ratio),
            se_expected       = first(se_expected))        # ... vs the expected ratio

# 3f-iv. Random-effect SDs and ICC
sapply(list(fl_norm  = m.fl.1,  fl_un  = m.fl.1.un,        # SD of intercept (and slope)
            inc_norm = m.inc.1, inc_un = m.inc.1.un),
       function(m) attr(VarCorr(m)[[1]], "stddev"))
sapply(list(fl_norm  = m.fl.1,  fl_un  = m.fl.1.un,        # ICC
            inc_norm = m.inc.1, inc_un = m.inc.1.un), icc)

# DECISION (fill in after running):
# k close to 1, ORs change < ~1%, SE ratios close to se_expected
#   -> normalisation makes no practical difference; keep w_fit (SEs match the real N).
# k far from 1 -> normalisation matters for the SEs; keep w_fit, report this check.



# ---- CHECK 2: PROPORTIONAL ODDS ---------------------------------------------
# 1a. Brant test: omnibus + one test per predictor. Significant = PO violated.
brant(p.fl)                                                # FLI
brant(p.inc)                                               # income
# If brant() errors (possible with many factor levels), rely on 1b and 1c.

# 1b. Likelihood-ratio test per predictor: lets each one vary across cuts.
nominal_test(c.fl)                                         # FLI; small p = that predictor violates PO
nominal_test(c.inc)                                        # income

# 1c. Coefficients across cuts: one binary logit per split of the outcome.
fit_cuts <- function(data, y) {                            # fits K binary models (K = levels - 1)
  K <- nlevels(data[[y]]) - 1                              # 5 levels -> 4 cuts
  map(set_names(seq_len(K), paste0("above_", seq_len(K))), function(k) {
    d.k <- mutate(data, y_bin = as.integer(as.integer(.data[[y]]) > k))  # 1 = above category k
    glm(as.formula(paste("y_bin ~", rhs_ctx)),             # same fixed part as M1
        data = d.k, weights = w_fit,                       # same sample and weights
        family = quasibinomial)                            # quasi: no warnings for non-integer weights
  })
}

plot_cuts <- function(cuts, m.clm, title) {                # plots each coefficient per cut
  est <- imap_dfr(cuts, function(m, cut) {                 # collect coefficients from each cut model
    co <- summary(m)$coefficients[-1, , drop = FALSE]      # drop the intercept
    tibble(cut = cut, term = rownames(co), b = co[, 1], se = co[, 2])
  })
  ref <- tibble(term = names(m.clm$beta), b = m.clm$beta)  # single PO estimate from clm()
  ggplot(est, aes(cut, b)) +
    geom_hline(yintercept = 0, colour = "grey60") +        # no-effect line
    geom_hline(data = ref, aes(yintercept = b),            # red dashed = PO estimate
               colour = "red", linetype = 2) +
    geom_pointrange(aes(ymin = b - 1.96 * se,              # cut-specific estimate with 95% CI
                        ymax = b + 1.96 * se)) +
    facet_wrap(~ term, scales = "free_y") +                # one panel per predictor
    labs(title = title, x = "Binary split", y = "Log-odds (95% CI)") +
    theme_minimal()
}

cuts.fl  <- fit_cuts(d.fl.cc,  "FLI_3item")                # FLI cut models (reused in check 3)
cuts.inc <- fit_cuts(d.inc.cc, "INCOME_5")                 # income cut models (reused in check 3)
plot_cuts(cuts.fl,  c.fl,  "FLI: coefficients across cuts")     # flat points near red line = PO fine
plot_cuts(cuts.inc, c.inc, "Income: coefficients across cuts")  # look for trends, not p-values

saveRDS(mget(c("cuts.fl", "cuts.inc")), "models_cuts.rds") # save: reload with load_models()
#load_models("models_cuts.rds")                             # cuts.fl, cuts.inc (Check 2)

# ---- CHECK 3: LINEARITY OF THE LOGIT ----------------------------------------
spline_check <- function(data, y, m.lin, var, df) {        # compares linear vs spline term in M1
  rhs_s <- str_replace(rhs_ctx, paste0("\\b", var, "\\b"), # swap the linear term ...
                       sprintf("ns(%s, df = %d)", var, df))  # ... for a natural spline
  m.s   <- fit_clmm(data, y, rhs_s)                        # refit the multilevel model with the spline
  print(anova(m.lin, m.s))                                 # LR test: small p = non-linear
  b   <- ns(data[[var]], df = df)                          # same basis, to redraw the fitted curve
  x   <- seq(min(data[[var]]), max(data[[var]]), length.out = 100)  # grid over the observed range
  bet <- m.s$beta[str_detect(names(m.s$beta), fixed(paste0("ns(", var)))]  # spline coefficients
  tibble(x      = x,
         spline = drop(predict(b, x) %*% bet),             # spline effect, relative to the minimum
         linear = m.lin$beta[[var]] * (x - min(x))) %>%    # linear effect, same reference point
    pivot_longer(-x, names_to = "form") %>%                # long format for ggplot
    ggplot(aes(x, value, colour = form)) +
    geom_line() +                                          # curves overlap = linear is fine
    labs(title = paste(y, "-", var), x = var, y = "Log-odds vs minimum") +
    theme_minimal()
}

spline_check(d.fl.cc,  "FLI_3item", m.fl.1,  "iwisescore",  df = 4)  # IWISE (zero-inflated)
spline_check(d.fl.cc,  "FLI_3item", m.fl.1,  "bas_rate_pp", df = 3)  # JMP rate (can be negative)
spline_check(d.inc.cc, "INCOME_5",  m.inc.1, "iwisescore",  df = 4)  # same for income
spline_check(d.inc.cc, "INCOME_5",  m.inc.1, "bas_rate_pp", df = 3)

# Alternative: IWISE in categories (non-nested, so compare AIC; lower = better)
# >>> CHANGED: a random slope for continuous iwisescore does not fit a model
#     with categorical IWISE, so both sides of this comparison are
#     intercept-only (m.*.1.ri from Section 3d).
m.fl.cat  <- fit_clmm(d.fl.cc,  "FLI_3item", str_replace(rhs_ctx, "\\biwisescore\\b", "iwise_cat"),
                      re = re_term)
m.inc.cat <- fit_clmm(d.inc.cc, "INCOME_5",  str_replace(rhs_ctx, "\\biwisescore\\b", "iwise_cat"),
                      re = re_term)
AIC(m.fl.1.ri,  m.fl.cat)                                  # FLI: continuous vs categorical IWISE
AIC(m.inc.1.ri, m.inc.cat)                                 # income
# <<< END CHANGED

saveRDS(mget(c("m.fl.cat", "m.inc.cat")), "models_iwise_cat.rds")

# ---- CHECK 4: INFLUENTIAL OBSERVATIONS --------------------------------------
# 3a. Cook's D and standardised residuals from the cut models (check 1c).
#     Each respondent gets their worst value across cuts.
infl_table <- function(cuts) {
  cd <- sapply(cuts, cooks.distance)                       # Cook's D: rows = respondents, cols = cuts
  rs <- sapply(cuts, rstandard)                            # standardised deviance residuals
  tibble(.id      = seq_len(nrow(cd)),                     # matches data$.id (no rows dropped)
         cook_max = apply(cd, 1, max),                     # worst Cook's D across cuts
         rstd_max = apply(abs(rs), 1, max))                # largest absolute residual across cuts
}

inf.fl  <- infl_table(cuts.fl)                             # FLI
inf.inc <- infl_table(cuts.inc)                            # income

summ_infl <- function(inf) {                               # quick overview
  inf %>% summarise(n          = n(),                      # respondents
                    n_cook_4n  = sum(cook_max > 4 / n()),  # above the common 4/n cut-off
                    n_resid_3  = sum(rstd_max > 3),        # |standardised residual| > 3
                    cook_p99   = quantile(cook_max, .99),  # threshold for the top 1%
                    cook_top   = max(cook_max))            # single most influential value
}
summ_infl(inf.fl)
summ_infl(inf.inc)

ggplot(inf.fl, aes(.id, cook_max)) + geom_point(size = .3) +   # index plot: look for isolated spikes
  labs(title = "FLI: Cook's distance (max over cuts)", x = "Respondent", y = "Cook's D") +
  theme_minimal()

# 3b. Refit M1 without the most influential 1%.
refit_drop <- function(data, y, inf, m.full) {
  top <- inf %>% slice_max(cook_max, prop = 0.01) %>% pull(.id)  # ids of the top 1%
  d.k <- data %>% filter(!.id %in% top) %>%                # drop them
    mutate(w_fit = w_fit / mean(w_fit))                    # re-rescale weights to the new N
  bind_rows(tidy_or(m.full, "all"),                        # original estimates
            tidy_or(fit_clmm(d.k, y, rhs_ctx), "minus top 1%")) %>%  # estimates without top 1%
    filter(term %in% v.watch) %>%                          # key terms only
    select(term, model, OR, lo, hi) %>%
    arrange(term, model)                                   # side by side per term
}
refit_drop(d.fl.cc,  "FLI_3item", inf.fl,  m.fl.1)         # similar ORs = results not driven by a few cases
refit_drop(d.inc.cc, "INCOME_5",  inf.inc, m.inc.1)

# 3c. Leave one country-year out. SLOW: one refit per country-year
#     (slower still with the random slope).
loo_country <- function(data, y) {
  ids <- levels(droplevels(data$country_year))             # all country-years in the sample
  map_dfr(ids, function(cy) {                              # one refit per dropped country-year
    message("Leaving out ", cy)                            # progress in the console
    d.k <- data %>% filter(country_year != cy) %>%         # drop this country-year
      mutate(country_year = droplevels(country_year),      # remove its empty factor level
             w_fit = w_fit / mean(w_fit))                  # re-rescale weights
    m <- tryCatch(fit_clmm(d.k, y, rhs_ctx), error = function(e) NULL)  # skip fits that fail
    if (is.null(m)) return(NULL)
    tidy_or(m, cy) %>% filter(term %in% v.watch)           # keep the key terms
  })
}

loo.fl  <- loo_country(d.fl.cc,  "FLI_3item")              # FLI
loo.inc <- loo_country(d.inc.cc, "INCOME_5")               # income
saveRDS(list(loo.fl = loo.fl, loo.inc = loo.inc), "loo_attempt1.rds")  # save: too slow to rerun

plot_loo <- function(loo, m.full, title) {                 # each point = estimate without one country-year
  ref <- tidy_or(m.full, "all") %>% filter(term %in% v.watch)  # full-sample OR and CI
  ggplot(loo, aes(model, OR)) +
    geom_hline(data = ref, aes(yintercept = OR), colour = "red") +               # full-sample OR
    geom_hline(data = ref, aes(yintercept = lo), colour = "red", linetype = 2) + # full-sample CI
    geom_hline(data = ref, aes(yintercept = hi), colour = "red", linetype = 2) +
    geom_point(size = .8) +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = title, x = "Country-year left out", y = "Odds ratio") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, size = 5))
}
plot_loo(loo.fl,  m.fl.1,  "FLI: leave one country-year out")
plot_loo(loo.inc, m.inc.1, "Income: leave one country-year out")

outside_ci <- function(loo, m.full, data) {                # which country-years move an OR outside the CI?
  ref <- tidy_or(m.full, "all") %>% select(term, lo_full = lo, hi_full = hi)
  names <- data %>% distinct(country_year, countrynew) %>%  # country names for the ids
    mutate(country_year = as.character(country_year))
  loo %>% left_join(ref, by = "term") %>%
    filter(OR < lo_full | OR > hi_full) %>%                # outside the full-sample 95% CI
    left_join(names, by = c("model" = "country_year"))
}
outside_ci(loo.fl,  m.fl.1,  d.fl.cc)
outside_ci(loo.inc, m.inc.1, d.inc.cc)


# ---- CHECK 5: GOODNESS OF FIT -----------------------------------------------
# Small p = poor fit. With N ~ 80,000 both tests reject on tiny misfit,
# so read them together with checks 1-3.
lipsitz.test(p.fl,  g = 10)                                # Lipsitz: 10 groups of predicted score, FLI
lipsitz.test(p.inc, g = 10)                                # income

v.pr <- c("female", "urban_imp")                           # covariate patterns for Pulkstenis-Robinson;
                                                           # more variables = sparse cells
pulkrob.chisq(p.fl,     catvars = v.pr)                    # Pulkstenis-Robinson chi-square, FLI
pulkrob.deviance(p.fl,  catvars = v.pr)                    # deviance version, FLI
pulkrob.chisq(p.inc,    catvars = v.pr)                    # income
pulkrob.deviance(p.inc, catvars = v.pr)