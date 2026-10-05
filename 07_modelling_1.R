# =============================================================================
# 07_models.R
# First modelling attempt: multilevel ordinal logistic regression
#   Outcomes : FLI_3item (d.fl) first, then INCOME_5 (d.inc)
#   M0       : unadjusted     = iwisescore + iwise_mean (Mundlak correction)
#   M1       : fully adjusted = rhs_ctx (individual + country-level covariates)
#   Both     : random intercept per country_year, weights = w_fit
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
dat <- readRDS("data/iwise_data_prepared.rds")                  # read saved list from 05
list2env(dat, envir = .GlobalEnv)                               # unpack list items into workspace
rm(dat)                                                         # remove the now-redundant list

prep <- function(data, y) {                                     # builds the analytic sample for one outcome
  data %>%
    filter(cc) %>%                                              # complete cases = exactly the rows the model uses
    mutate(.id          = row_number(),                         # row id, to trace influential observations later
           country_year = factor(.data[[v.cluster]]),           # grouping factor for the random intercept
           iwise_cat    = factor(iwise_cat, ordered = FALSE))   # unordered: avoids polynomial contrasts
}

d.fl.cc  <- prep(d.fl,  "FLI_3item")                       # FLI analytic sample
d.inc.cc <- prep(d.inc, "INCOME_5")                        # income analytic sample

table(d.fl.cc$FLI_3item)                                   # check: 5 levels, lowest first
table(d.inc.cc$INCOME_5)                                   # check: 5 levels, poorest quintile first

v.watch <- c("iwisescore", "iwise_mean", "log_gdp",        # terms to monitor in checks 3 and LOO:
             v.wgi.main, v.jmp, v.gbd)                     # IWISE terms + country-level covariates


# ---- 2. FORMULAS AND HELPERS ------------------------------------------------
rhs_unadj <- "iwisescore + iwise_mean"                     # M0: IWISE + its country-year mean (Mundlak)
re_term   <- "(1 | country_year)"                          # random intercept: one shift per country-year

make_f <- function(y, rhs, random = TRUE) {                # turns text into a model formula
  as.formula(paste(y, "~", rhs,                            # outcome ~ predictors
                   if (random) paste("+", re_term) else "")) # + random intercept (or none)
}

fit_clmm <- function(data, y, rhs) {                       # same model as Section 3, used in loops
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




sapply(d.fl.cc[c("iwisescore", "iwise_mean", "hhsize", "log_gdp", v.wgi.main, v.jmp, v.gbd)],
       function(x) round(c(mean = mean(x), sd = sd(x)), 2))       # 

# ---- 3. FIT THE MODELS ------------------------------------------------------
# Formula objects are created at top level so update() (used inside
# nominal_test) can find them later.
f.fl.0  <- make_f("FLI_3item", rhs_unadj)                  # FLI_3item ~ iwisescore + iwise_mean + (1 | country_year)
f.fl.1  <- make_f("FLI_3item", rhs_ctx)                    # FLI_3item ~ all covariates + (1 | country_year)
f.inc.0 <- make_f("INCOME_5",  rhs_unadj)                  # same for income
f.inc.1 <- make_f("INCOME_5",  rhs_ctx)                    # same for income

m.fl.0 <- clmm(f.fl.0,                                     # model formula
               data    = d.fl.cc,                          # complete-case FLI sample
               weights = w_fit,                            # regression weights, rescaled to sum to N
               link    = "logit",                          # cumulative logit = ordinal logistic
               Hess    = TRUE,                             # Hessian, needed for standard errors
               nAGQ    = 1)                                # Laplace approximation of the random effect
m.fl.1  <- clmm(f.fl.1,  data = d.fl.cc,  weights = w_fit, link = "logit", Hess = TRUE, nAGQ = 1)  # FLI, adjusted
m.inc.0 <- clmm(f.inc.0, data = d.inc.cc, weights = w_fit, link = "logit", Hess = TRUE, nAGQ = 1)  # income, unadjusted
m.inc.1 <- clmm(f.inc.1, data = d.inc.cc, weights = w_fit, link = "logit", Hess = TRUE, nAGQ = 1)  # income, adjusted

saveRDS(list(m.fl.0 = m.fl.0, m.fl.1 = m.fl.1,             # save fits so you don't refit every session
             m.inc.0 = m.inc.0, m.inc.1 = m.inc.1),
        "models_attempt1.rds")

summary(m.fl.1)                                            # full output; cond.H < 1e4 = no estimation problem
summary(m.inc.1)                                           # same for income

res <- bind_rows(tidy_or(m.fl.0,  "FLI M0"), tidy_or(m.fl.1,  "FLI M1"),   # all ORs in one table
                 tidy_or(m.inc.0, "INC M0"), tidy_or(m.inc.1, "INC M1"))
res %>% filter(str_detect(term, "iwise")) %>% print()      # IWISE rows: within (iwisescore) and contextual (iwise_mean)

sapply(list(fl0 = m.fl.0, fl1 = m.fl.1,                    # ICC per model: share of variance
            inc0 = m.inc.0, inc1 = m.inc.1), icc)          # between country-years

# ---- 3b. WITH VS WITHOUT THE COUNTRY-YEAR MEAN (MUNDLAK) --------------------  # >>> ADDED
# iwise_mean = contextual effect (between minus within). If it is ~0, the
# within and between effects are equal and the mean could be dropped.
# Without the mean, iwisescore is a blend of within and between effects.
rhs_unadj_nomean <- "iwisescore"                           # M0 without the mean: IWISE only
rhs_ctx_nomean   <- str_remove(rhs_ctx, "\\+ iwise_mean ") # M1 without the mean
rhs_ctx_nomean                                             # check: iwise_mean is gone, rest unchanged
 
m.fl.0.nm  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_unadj_nomean)  # FLI, unadjusted, no mean
m.fl.1.nm  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx_nomean)    # FLI, adjusted, no mean
m.inc.0.nm <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_unadj_nomean)  # income, unadjusted, no mean
m.inc.1.nm <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx_nomean)    # income, adjusted, no mean
 
# 3b-i. IWISE odds ratios side by side
bind_rows(tidy_or(m.fl.0,  "FLI M0 with mean"),  tidy_or(m.fl.0.nm,  "FLI M0 without mean"),
          tidy_or(m.fl.1,  "FLI M1 with mean"),  tidy_or(m.fl.1.nm,  "FLI M1 without mean"),
          tidy_or(m.inc.0, "INC M0 with mean"),  tidy_or(m.inc.0.nm, "INC M0 without mean"),
          tidy_or(m.inc.1, "INC M1 with mean"),  tidy_or(m.inc.1.nm, "INC M1 without mean")) %>%
  filter(term %in% c("iwisescore", "iwise_mean")) %>%      # IWISE rows only
  select(model, term, OR, lo, hi, p) %>%                   # OR, 95% CI, p-value
  print()                                                  # similar iwisescore ORs = the mean barely matters
 
# 3b-ii. Likelihood-ratio test: does adding the mean improve the fit?
#        Small p = within and between effects differ, so keep the mean
#        (Hedeker 2015, p. 9; equivalent to the iwise_mean Wald test).
anova(m.fl.0.nm,  m.fl.0)                                  # FLI, unadjusted
anova(m.fl.1.nm,  m.fl.1)                                  # FLI, adjusted
anova(m.inc.0.nm, m.inc.0)                                 # income, unadjusted
anova(m.inc.1.nm, m.inc.1)                                 # income, adjusted
 
# =============================================================================
# CHECKS
# brant(), nominal_test(), Cook's distance and the Lipsitz / Pulkstenis-
# Robinson tests do not accept clmm() objects. They are run on single-level
# stand-ins with the SAME fixed part (rhs_ctx), sample and (where possible)
# weights. Refits (checks 2 and 3) use the real multilevel model.
# =============================================================================

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