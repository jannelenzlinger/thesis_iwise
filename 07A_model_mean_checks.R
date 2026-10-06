# =============================================================================
# 07_models.R
# First modelling attempt: multilevel ordinal logistic regression
#   Outcomes : FLI_3item (d.fl) first, then INCOME_5 (d.inc)
#   M0       : unadjusted     = iwisescore
#   M1       : fully adjusted = rhs_ctx (individual + country-level covariates)
#   Both     : random intercept per country_year, weights = w_fit
#   Section 3 first checks whether the country-year mean (Mundlak) is needed,
#   then sets the final models (decision: exclude the mean).
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

make_f <- function(y, rhs, random = TRUE) {                # turns text into a model formula
  as.formula(paste(y, "~", rhs,                            # outcome ~ predictors
                   if (random) paste("+", re_term) else "")) # + random intercept (or none)
}
#function helps for loops and makes sure the settings are the same everywhere
fit_clmm <- function(data, y, rhs) {                       # same model everywhere, used in loops
  clmm(make_f(y, rhs), data = data, weights = w_fit,       # weighted multilevel ordinal logit
       link = "logit", Hess = TRUE, nAGQ = 1)              # logit link, keep Hessian, Laplace approx.
}

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

icc <- function(m) {                                       # latent-scale ICC from the fitted model
  s2 <- as.numeric(VarCorr(m)[[1]])                        # random-intercept variance
  s2 / (s2 + pi^2 / 3)                                     # pi^2/3 = level-1 variance of the logistic
}


# ---- 3a. IS THE COUNTRY-YEAR MEAN NEEDED? (MUNDLAK) -------------------------  # >>> CHANGED: now first
# iwise_mean = contextual effect (between minus within). If it is ~0 and the
# iwisescore OR does not change, the mean can be dropped.
# Without the mean, iwisescore is a blend of within and between effects.
m.fl.0.nm  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_unadj)        # FLI, unadjusted, without mean
m.fl.0.wm  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_unadj_mean)   # FLI, unadjusted, with mean
m.fl.1.nm  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx)          # FLI, adjusted, without mean
m.fl.1.wm  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx_mean)     # FLI, adjusted, with mean
m.inc.0.nm <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_unadj)        # income, unadjusted, without mean
m.inc.0.wm <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_unadj_mean)   # income, unadjusted, with mean
m.inc.1.nm <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx)          # income, adjusted, without mean
m.inc.1.wm <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx_mean)     # income, adjusted, with mean

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
anova(m.fl.1.nm,  m.fl.1.wm)                               # FLI, adjusted     (p = 0.48)
anova(m.inc.0.nm, m.inc.0.wm)                              # income, unadjusted (p < 0.001)
anova(m.inc.1.nm, m.inc.1.wm)                              # income, adjusted  (p = 0.018, dAIC = 4)

# =============================================================================
# DECISION: EXCLUDE MEAN
# The iwisescore OR is unchanged in the adjusted models (FLI 0.967, income
# 0.982, with and without), fit improves not at all for FLI and only
# marginally for income. The with-mean models are reported as sensitivity.
# =============================================================================


# ---- 3b. FINAL MODELS (WITHOUT THE MEAN) ------------------------------------
# The final models are the "without mean" models from 3a. They are identical
# to a refit, so they are reused instead of fitted again.
m.fl.0  <- m.fl.0.nm                                       # FLI, unadjusted
m.fl.1  <- m.fl.1.nm                                       # FLI, adjusted
m.inc.0 <- m.inc.0.nm                                      # income, unadjusted
m.inc.1 <- m.inc.1.nm                                      # income, adjusted

saveRDS(list(m.fl.0 = m.fl.0, m.fl.1 = m.fl.1,             # save final fits
             m.inc.0 = m.inc.0, m.inc.1 = m.inc.1),  # + sensitivity models with the mean
        "models_attempt1.rds")

# Estimation quality of the final models
sapply(list(fl0 = m.fl.0, fl1 = m.fl.1,                    # condition number of the Hessian:
            inc0 = m.inc.0, inc1 = m.inc.1),               # < ~1e4 fine, > ~1e6 problematic
       function(m) summary(m)$condHess)
bind_rows(tidy_or(m.fl.1, "FLI M1"), tidy_or(m.inc.1, "INC M1")) %>%
  filter(is.na(se))                                        # 0 rows = all standard errors defined


# ---- 3c. DOES THE CENTRING WORK? --------------------------------------------
# The final models use centred covariates (05, Section 5.2). Here the same
# models are refitted on the RAW covariates and compared.
# First run (uncentred FLI M1): condHess NaN, SEs undefined, country-level
# betas off (e.g. bas_rate_pp +0.068 vs -0.046 centred) = not converged.
# Centred test model: condHess 71,488, individual-level betas unchanged.
rhs_ctx_raw <- reduce(v.centre,                            # swap each centred name back to raw:
                      ~ str_replace_all(.x, paste0("\\b", .y, "_c\\b"), .y),
                      .init = rhs_ctx)                     # e.g. log_gdp_c -> log_gdp
rhs_ctx_raw                                                # check: no "_c" names left

m.fl.1.raw  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx_raw)  # FLI, adjusted, uncentred
m.inc.1.raw <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx_raw)  # income, adjusted, uncentred
# (warnings about the variance-covariance matrix are expected for the raw fits)

# 3c-i. Condition number: centred should be clearly lower (and not NaN)
sapply(list(fl_raw  = m.fl.1.raw,  fl_centred  = m.fl.1,
            inc_raw = m.inc.1.raw, inc_centred = m.inc.1),
       function(m) summary(m)$condHess)

# 3c-ii. Log-likelihood: centred equal or higher = raw fit had not fully converged
c(fl_raw  = logLik(m.fl.1.raw),  fl_centred  = logLik(m.fl.1),
  inc_raw = logLik(m.inc.1.raw), inc_centred = logLik(m.inc.1))

# 3c-iii. Betas side by side (same term order; row names from the raw model).
#         Individual-level betas should match; country-level ones may differ
#         if the raw fit did not converge.
cbind(raw = m.fl.1.raw$beta,  centred = m.fl.1$beta)       # FLI
cbind(raw = m.inc.1.raw$beta, centred = m.inc.1$beta)      # income

# 3c-iv. Country-level ORs in the final (centred) models
bind_rows(tidy_or(m.fl.1, "FLI M1"), tidy_or(m.inc.1, "INC M1")) %>%
  filter(term %in% v.watch[-1]) %>%                        # country-level terms only
  select(model, term, OR, lo, hi, p)

summary(m.fl.1)                                            # full output
summary(m.inc.1)                                           # same for income

res <- bind_rows(tidy_or(m.fl.0,  "FLI M0"), tidy_or(m.fl.1,  "FLI M1"),   # all ORs in one table
                 tidy_or(m.inc.0, "INC M0"), tidy_or(m.inc.1, "INC M1"))
res %>% filter(term == "iwisescore") %>% print()           # IWISE rows

sapply(list(fl0 = m.fl.0, fl1 = m.fl.1,                    # ICC per model: share of variance
            inc0 = m.inc.0, inc1 = m.inc.1), icc)          # between country-years


# =============================================================================
# CHECKS
# brant(), nominal_test(), Cook's distance and the Lipsitz / Pulkstenis-
# Robinson tests do not accept clmm() objects. They are run on single-level
# stand-ins with the SAME fixed part (rhs_ctx), sample and (where possible)
# weights. Refits (checks 2 and 3) use the real multilevel model.
# =============================================================================

#each test gets a single-level stand-in: the same predictors, sample (and where possible) weights, just without the ransdom intercept
f.fl.fix  <- make_f("FLI_3item", rhs_ctx, random = FALSE)  # FLI formula without random intercept
f.inc.fix <- make_f("INCOME_5",  rhs_ctx, random = FALSE)  # income formula without random intercept

c.fl  <- clm(f.fl.fix,  data = d.fl.cc,  weights = w_fit, link = "logit")  # weighted clm stand-in, FLI
c.inc <- clm(f.inc.fix, data = d.inc.cc, weights = w_fit, link = "logit")  # weighted clm stand-in, income

p.fl  <- polr(f.fl.fix,  data = d.fl.cc,  method = "logistic", Hess = TRUE)  # polr stand-in (unweighted:
p.inc <- polr(f.inc.fix, data = d.inc.cc, method = "logistic", Hess = TRUE)  # brant/generalhoslem ignore weights)


# ---- CHECK 1: PROPORTIONAL ODDS ---------------------------------------------
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


# ---- CHECK 2: LINEARITY OF THE LOGIT ----------------------------------------
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
m.fl.cat  <- fit_clmm(d.fl.cc,  "FLI_3item", str_replace(rhs_ctx, "\\biwisescore\\b", "iwise_cat"))
m.inc.cat <- fit_clmm(d.inc.cc, "INCOME_5",  str_replace(rhs_ctx, "\\biwisescore\\b", "iwise_cat"))
AIC(m.fl.1,  m.fl.cat)                                     # FLI: continuous vs categorical IWISE
AIC(m.inc.1, m.inc.cat)                                    # income


# ---- CHECK 3: INFLUENTIAL OBSERVATIONS --------------------------------------
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

# 3c. Leave one country-year out. SLOW: one refit per country-year.
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


# ---- CHECK 4: GOODNESS OF FIT -----------------------------------------------
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