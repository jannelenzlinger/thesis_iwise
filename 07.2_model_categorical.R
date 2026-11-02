# =============================================================================
# 07.2_models_iwise_cat.R
# Final models from 07.1 (Section 3e) refitted with CATEGORICAL IWISE, and a
# check whether CONTINUOUS IWISE is adequate (Section 8).
#   iwise_cat replaces iwisescore. Everything else as in 07.1.
#   M0 : unadjusted     = iwise_cat
#   M1 : fully adjusted = rhs_ctx with iwisescore -> iwise_cat
#   Weights: w_fit. Outcomes: FLI_3item (d.fl), INCOME_5 (d.inc).
#
# Why two parts?
#   Part A (Sections 3-7) fits the categorical models. With categories, a
#   random IWISE slope becomes large (one slope per category), which caused
#   fitting problems.
#   Part B (Section 8) asks the question behind the choice: is the effect of
#   continuous IWISE close enough to a straight line on the log-odds scale
#   ("linearity of the logit")? If yes, continuous IWISE (07.1) can stay the
#   main model and the categorical models are a sensitivity analysis.
#
# Random effects in Part A: RANDOM IWISE SLOPE (needed in 07.1, 3d-4).
#   full   (1 + iwise_cat | country_year)  each category its own random
#          slope: 9 extra parameters vs random intercept
#   linear (1 + iwise_num | country_year)  one random slope across the
#          ordered categories, iwise_num = 0,1,2,3: 2 extra parameters,
#          same as in 07.1. Fixed part stays fully categorical.
#   Random-intercept fits (Section 3) are only the reference for the LR test.
#
#   PART A: CATEGORICAL MODELS          PART B: IS CONTINUOUS IWISE ADEQUATE?
#   1  Data                             8a  AIC: continuous vs categorical
#   2  Helpers                          8b  Spline + fine bins (direct check)
#   3  Random-intercept reference       8c  Piecewise linear (knots at cut-offs)
#   4  Random-slope models              8d  Decision
#   5  Checks, stability, LR tests
#   6  Final models
#   7  Results
# Saves: models_cat_ri.rds, models_cat_slope.rds, models_cat_final.rds
#        (-> used by 09), iwise_shape.rds, models_piecewise.rds
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
library(tidyverse)                   # dplyr, purrr, stringr, ggplot2
library(ordinal)                     # clmm()
library(splines)                     # ns(): natural splines (Section 8b)


# ---- 1. DATA  -------------------------
list2env(readRDS("data/iwise_data_prepared.rds"), envir = globalenv())  # d.fl, d.inc, rhs_ctx, v.*

prep <- function(data, y) {                                # analytic sample for one outcome
  data %>%
    filter(cc) %>%                                         # complete cases = rows the model uses
    mutate(.id          = row_number(),
           country_year = factor(.data[[v.cluster]]),      # grouping factor for random effects
           iwise_cat    = factor(iwise_cat, ordered = FALSE),  # unordered: one dummy per category
           iwise_num    = as.integer(iwise_cat) - 1,       # 0,1,2,3: used only for the linear random slope
           "{y}" := factor(.data[[y]], ordered = TRUE))    # clmm() needs an ordered outcome
}

d.fl.cc  <- prep(d.fl,  "FLI_3item")
d.inc.cc <- prep(d.inc, "INCOME_5")

# Levels must run from lowest to highest WI: the first level is the
# reference category AND iwise_num = 0. Check both:
levels(d.fl.cc$iwise_cat)
table(d.fl.cc$iwise_cat, d.fl.cc$iwise_num)                # check: diagonal, no sparse categories
table(d.inc.cc$iwise_cat)


# ---- 2. FORMULAS AND HELPERS ------------------------------------------------
rhs_unadj_cat <- "iwise_cat"                               # M0: categorical IWISE only
rhs_ctx_cat   <- str_replace(rhs_ctx, "\\biwisescore\\b", "iwise_cat")  # M1: swap the IWISE term
rhs_ctx_cat                                                # check: iwise_cat in, iwisescore out

re_term      <- "(1 | country_year)"                       # random intercept only
re_slope_cat <- "(1 + iwise_cat | country_year)"           # full: random slope per category
re_slope_num <- "(1 + iwise_num | country_year)"           # linear: one slope across categories
re_slope     <- "(1 + iwisescore | country_year)"          # continuous slope as in 07.1 (Section 8c)

make_f <- function(y, rhs, re) as.formula(paste(y, "~", rhs, "+", re))

fit_clmm <- function(data, y, rhs, re) {                   # same settings as 07.1
  clmm(make_f(y, rhs, re), data = data, weights = w_fit,
       link = "logit", Hess = TRUE, nAGQ = 1)              # random slopes need nAGQ = 1
}

tidy_or <- function(m, label) {                            # OR table for the fixed effects
  co <- summary(m)$coefficients[names(m$beta), , drop = FALSE]
  tibble(model = label, term = rownames(co),
         b = co[, 1], se = co[, 2], OR = exp(b),
         lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = co[, 4])
}

load_models <- function(file) {                            # restores every model in an .rds file
  obj <- readRDS(file)
  list2env(obj, envir = globalenv())
  message("Loaded from ", file, ": ", paste(names(obj), collapse = ", "))
  invisible(names(obj))
}

icc <- function(m) {                                       # latent-scale ICC; with a random slope
  s2 <- VarCorr(m)[[1]][1, 1]                              # it applies at the reference category
  s2 / (s2 + pi^2 / 3)
}


# =============================================================================
# PART A: CATEGORICAL MODELS
# =============================================================================

# ---- 3. RANDOM-INTERCEPT REFERENCE MODELS (NOT FINAL) -----------------------
# Purpose: the "without slope" version of each model. Needed for
#   - the LR test in Section 5d (is the random slope needed?)
#   - inc0, which keeps the random intercept only (see Section 4)
#   - the AIC comparison in Section 8a
m.fl.0.cat.ri  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_unadj_cat, re_term)  # FLI, unadjusted
m.fl.1.cat.ri  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx_cat,   re_term)  # FLI, adjusted
m.inc.0.cat.ri <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_unadj_cat, re_term)  # income, unadjusted
m.inc.1.cat.ri <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx_cat,   re_term)  # income, adjusted

saveRDS(mget(c("m.fl.0.cat.ri", "m.fl.1.cat.ri", "m.inc.0.cat.ri", "m.inc.1.cat.ri")),
        "models_cat_ri.rds")
# Reload instead of refitting: load_models("models_cat_ri.rds")


# ---- 4. RANDOM-SLOPE MODELS -------------------------------------------------
# Purpose: same models + random IWISE slope, so each country-year can have
# its own IWISE effect. SLOW.
# RESULT of the first run (Section 5a, all four with the planned slope):
#   fl0  full slope   FAILED  condHess NaN, 3 SEs NA, gradient 81
#   fl1  full slope   ok      condHess 9,422, SEs defined, gradient 8.2
#   inc0 linear slope FAILED  condHess NaN, 3 SEs NA, logLik BELOW the
#                             random-intercept model (= not converged)
#   inc1 full slope   ok      condHess 5,661, SEs defined, gradient 9.6
# DECISION: fl0 -> linear slope (fewer parameters); inc0 -> random intercept
#   only, as in 07.1 (3d-2): unadjusted income has too little between-country
#   variance to identify a slope. fl1, inc1 -> full slope, if 5c passes.
m.fl.0.cat.rs  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_unadj_cat, re_slope_num)  # FLI, unadjusted (linear)
m.fl.1.cat.rs  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx_cat,   re_slope_cat)  # FLI, adjusted (full)
m.inc.1.cat.rs <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx_cat,   re_slope_cat)  # income, adjusted (full)
# inc0: no slope model (random intercept only, Section 6)

slope_type <- c(fl0 = "linear", fl1 = "full", inc0 = "none", inc1 = "full")  # keep in sync with the fits

saveRDS(mget(c("m.fl.0.cat.rs", "m.fl.1.cat.rs", "m.inc.1.cat.rs", "slope_type")),
        "models_cat_slope.rds")
# Reload instead of refitting: load_models("models_cat_slope.rds")

# 4b. FALLBACK: if fl1 or inc1 fail 5a or 5c, refit with the linear slope
#     and update slope_type, then save again. Example for FLI adjusted:
# m.fl.1.cat.rs     <- fit_clmm(d.fl.cc, "FLI_3item", rhs_ctx_cat, re_slope_num)
# slope_type["fl1"] <- "linear"


# ---- 5. CHECKS, STABILITY, LR TESTS (as in 07.1, 3d-2 to 3d-4) --------------
# Each slope model next to its random-intercept version (inc0 has no slope model)
pairs <- list(fl0  = list(ri = m.fl.0.cat.ri,  rs = m.fl.0.cat.rs),
              fl1  = list(ri = m.fl.1.cat.ri,  rs = m.fl.1.cat.rs),
              inc1 = list(ri = m.inc.1.cat.ri, rs = m.inc.1.cat.rs))

# 5a. Did the fits work? A fit can only be trusted if ALL rows look fine.
fit_quality <- function(p) {
  c(condHess    = suppressWarnings(summary(p$rs)$condHess),   # NaN or > 1e6 = not identified
    max_grad    = max(abs(p$rs$gradient)),                    # large = may not have converged
    n_na_se     = sum(is.na(summary(p$rs)$coefficients[names(p$rs$beta), 2])),  # 0 = all SEs defined
    logLik_gain = as.numeric(logLik(p$rs) - logLik(p$ri)))    # must be >= 0 (ri is nested in rs)
}
round(sapply(pairs, fit_quality), 4)

VarCorr(m.fl.0.cat.rs)        # linear refit: "(Intercept)" and "iwise_num"
                              # old fit:      "(Intercept)" and three "iwise_cat..." rows

# 5b. Random-effect SDs and correlations: plausible? Correlations of +-1
#     mean the slope is not identified -> use the fallback (4b).
lapply(pairs, function(p) VarCorr(p$rs)[[1]])

# 5c. Are fl1 and inc1 stable? Their gradients (8-10) are not close to 0.
#     Refit with another optimiser: if logLik, betas and SDs are the same,
#     the original fit is a real maximum and the gradient does not matter.
fit_ucminf <- function(data, y, rhs, re) {                 # as fit_clmm(), other optimiser
  clmm(make_f(y, rhs, re), data = data, weights = w_fit, link = "logit",
       Hess = TRUE, nAGQ = 1, control = clmm.control(method = "ucminf"))
}
m.fl.1.cat.chk  <- fit_ucminf(d.fl.cc,  "FLI_3item", rhs_ctx_cat, re_slope_cat)
m.inc.1.cat.chk <- fit_ucminf(d.inc.cc, "INCOME_5",  rhs_ctx_cat, re_slope_cat)
saveRDS(mget(c("m.fl.1.cat.chk", "m.inc.1.cat.chk")), "models_cat_slope_check.rds")

cmp_fits <- function(m1, m2) {                             # how far apart are two fits?
  c(logLik_diff   = as.numeric(logLik(m2) - logLik(m1)),   # ~0 = same maximum
    max_beta_diff = max(abs(m1$beta - m2$beta)),           # < ~0.001 = same estimates
    max_sd_diff   = max(abs(attr(VarCorr(m1)[[1]], "stddev") -
                            attr(VarCorr(m2)[[1]], "stddev"))))  # < ~0.001 = same SDs
}
round(cbind(fl1  = cmp_fits(m.fl.1.cat.rs,  m.fl.1.cat.chk),
            inc1 = cmp_fits(m.inc.1.cat.rs, m.inc.1.cat.chk)), 5)
# If logLik_diff is clearly > 0, the ucminf fit is better: use it as m.*.1.cat.rs.
# If the fits differ a lot, the full slope is not stable: use the fallback (4b).

# 5d. Is the slope needed? LR test vs the random-intercept model.
#     linear: 50:50 mix of chi2(1) and chi2(2), as in 07.1 (3d-4).
#     full:   naive chi2 with k(k+1)/2 - 1 df; conservative (variances are
#             tested at their boundary, Snijders & Bosker), so read with dAIC.
k <- nlevels(d.fl.cc$iwise_cat)
lr_slope <- function(p, type) {
  lr <- as.numeric(2 * (logLik(p$rs) - logLik(p$ri)))
  df <- if (type == "full") k * (k + 1) / 2 - 1 else 2
  p_val <- if (type == "full") pchisq(lr, df, lower.tail = FALSE) else
             0.5 * pchisq(lr, 1, lower.tail = FALSE) + 0.5 * pchisq(lr, 2, lower.tail = FALSE)
  c(LR = lr, df = df, p = p_val, dAIC = AIC(p$rs) - AIC(p$ri))   # dAIC < 0 = slope better
}
round(sapply(names(pairs), function(nm) lr_slope(pairs[[nm]], slope_type[[nm]])), 4)


# ---- 6. FINAL CATEGORICAL MODELS --------------------------------------------
m.fl.0.cat  <- m.fl.0.cat.rs                               # FLI, unadjusted (linear slope)
m.fl.1.cat  <- m.fl.1.cat.rs                               # FLI, adjusted
m.inc.0.cat <- m.inc.0.cat.ri                              # income, unadjusted: random intercept only
m.inc.1.cat <- m.inc.1.cat.rs                              # income, adjusted

# Random-effects structure of the adjusted models, per outcome (used in 09,
# so the interaction models are nested in these)
re_of  <- c(full = re_slope_cat, linear = re_slope_num)
re_cat <- c(fl = re_of[[slope_type[["fl1"]]]], inc = re_of[[slope_type[["inc1"]]]])
re_cat

saveRDS(list(m.fl.0.cat = m.fl.0.cat, m.fl.1.cat = m.fl.1.cat,
             m.inc.0.cat = m.inc.0.cat, m.inc.1.cat = m.inc.1.cat,
             slope_type = slope_type, re_cat = re_cat),
        "models_cat_final.rds")


# ---- 7. RESULTS (CATEGORICAL) -----------------------------------------------
# 7a. IWISE ORs per category vs the reference category (average effects)
res.cat <- bind_rows(tidy_or(m.fl.0.cat,  "FLI M0"), tidy_or(m.fl.1.cat,  "FLI M1"),
                     tidy_or(m.inc.0.cat, "INC M0"), tidy_or(m.inc.1.cat, "INC M1"))
res.cat %>% filter(str_starts(term, "iwise_cat")) %>%
  select(model, term, OR, lo, hi, p) %>% print()

# 7b. Average IWISE ORs with vs without random slope (cf. 07.1, 3d-5a).
#     Expected: similar ORs, wider CIs with the slope.
bind_rows(tidy_or(m.fl.1.cat.ri,  "FLI M1 intercept only"), tidy_or(m.fl.1.cat,  "FLI M1 + slope"),
          tidy_or(m.inc.1.cat.ri, "INC M1 intercept only"), tidy_or(m.inc.1.cat, "INC M1 + slope")) %>%
  filter(str_starts(term, "iwise_cat")) %>%
  select(model, term, OR, lo, hi)

# 7c. Spread between country-years: range of 95% of the country-year ORs
#     per category (cf. 07.1, 3d-5b). full: SD per category; linear: the
#     random slope SD times iwise_num (= 1, 2, 3 for the categories).
cat_range <- function(m, type) {
  sds  <- attr(VarCorr(m)[[1]], "stddev")
  main <- grep("^iwise_cat", names(m$beta), value = TRUE)
  map_dfr(seq_along(main), function(j) {
    b  <- m$beta[[main[j]]]
    sd <- if (type == "full") sds[[main[j]]] else j * sds[["iwise_num"]]
    tibble(term = main[j], OR = exp(b),
           lower95 = exp(b - 1.96 * sd), upper95 = exp(b + 1.96 * sd))
  })
}
bind_rows(FLI = cat_range(m.fl.1.cat,  slope_type[["fl1"]]),
          INC = cat_range(m.inc.1.cat, slope_type[["inc1"]]), .id = "outcome")

# 7d. Joint Wald test: all iwise_cat terms = 0 (k - 1 df). No refit needed.
#     With this N it will almost certainly reject; the ORs and CIs matter more.
wald_iwise <- function(m) {
  bn <- grep("^iwise_cat", names(m$beta), value = TRUE)
  b  <- m$beta[bn]
  V  <- vcov(m)[bn, bn]
  W  <- drop(t(b) %*% solve(V) %*% b)
  c(W = W, df = length(bn), p = pchisq(W, length(bn), lower.tail = FALSE))
}
round(sapply(list(fl0 = m.fl.0.cat, fl1 = m.fl.1.cat,
                  inc0 = m.inc.0.cat, inc1 = m.inc.1.cat), wald_iwise), 4)

# 7e. ICC per model (at the reference category)
sapply(list(fl0 = m.fl.0.cat, fl1 = m.fl.1.cat,
            inc0 = m.inc.0.cat, inc1 = m.inc.1.cat), icc)

# 7f. Plot: adjusted ORs per category, both outcomes
res.cat %>%
  filter(str_starts(term, "iwise_cat"), str_detect(model, "M1")) %>%
  mutate(term = str_remove(term, "^iwise_cat")) %>%
  ggplot(aes(term, OR, ymin = lo, ymax = hi)) +
  geom_hline(yintercept = 1, colour = "grey60", linetype = 2) +
  geom_pointrange() +
  facet_wrap(~ model) +
  scale_y_log10() +
  labs(x = "WI category (vs reference)", y = "Adjusted OR (95% CI, log scale)") +
  theme_minimal()


# =============================================================================
# PART B: IS CONTINUOUS IWISE ADEQUATE? (LINEARITY OF THE LOGIT FOR IWISE)
# =============================================================================
# The continuous model (07.1) assumes that each extra IWISE point changes the
# log-odds of a higher outcome level by the same amount, wherever on the
# 0-36 scale it happens. 8a-8c test this in three ways:
#   8a  AIC: continuous vs categorical     indirect, whole models compared
#   8b  spline + fine bins                 direct: does the shape bend?
#   8c  piecewise, knots at the cut-offs   direct: do the slopes differ
#                                          between the IWISE categories?
# With N ~ 80,000, p-values will be small even for tiny deviations. Decide
# on the size of the deviation and the plots, not on p alone (8d).

load_models("models_final.rds")                            # 07.1 (3e): m.fl.1, m.inc.1 (random slope)
load_models("models_mundlak.rds")                          # 07.1 (3a): random-intercept versions
m.fl.1.ri  <- m.fl.1.nm                                    # FLI, adjusted, continuous, random intercept
m.inc.1.ri <- m.inc.1.nm                                   # income, adjusted, continuous, random intercept


# ---- 8a. AIC: CONTINUOUS VS CATEGORICAL -------------------------------------
# Non-nested models: compare AIC (lower = better).
# Random intercept on both sides = fair comparison of the IWISE term alone.
# Random slope on both sides = comparison of the final models, but the
# categorical slope has more parameters, which AIC also counts.
AIC(m.fl.1.ri,  m.fl.1.cat.ri)                             # FLI, random intercept
AIC(m.inc.1.ri, m.inc.1.cat.ri)                            # income, random intercept
AIC(m.fl.1,     m.fl.1.cat)                                # FLI, final models
AIC(m.inc.1,    m.inc.1.cat)                               # income, final models
# Categorical clearly lower = continuous may miss part of the shape -> 8b, 8c.
# Similar or continuous lower = the straight line is good enough.


# ---- 8b. SPLINE + FINE BINS (DIRECT CHECK) -----------------------------------
# Lets IWISE take any shape and compares it with the straight line:
#   spline: iwisescore -> ns(iwisescore, 3): a smooth curve. Linear is
#           nested in it -> LR test (2 df) and dAIC.
#   bins:   IWISE as many small categories (0, 1-3, 4-6, ...): no smoothing,
#           just the data. Do the bin estimates lie on the straight line?
# Random intercept for speed (as in 07.1, check 3); compared with m.*.1.ri.
# Everything is plotted relative to IWISE = 0.
iwise_shape <- function(data, y, m.lin, df = 3) {
  data <- data %>%
    mutate(iw_bin = cut(iwisescore, c(-Inf, 0, 3, 6, 9, 12, 15, 18, 21, 24, 30, Inf)))
  print(table(data$iw_bin))                                # check: no near-empty bins

  m.s <- fit_clmm(data, y, str_replace(rhs_ctx, "\\biwisescore\\b",          # spline model
                                       sprintf("ns(iwisescore, df = %d)", df)), re_term)
  m.b <- fit_clmm(data, y, str_replace(rhs_ctx, "\\biwisescore\\b", "iw_bin"), re_term)  # bin model

  # Straight line from the linear model
  x <- seq(0, max(data$iwisescore), length.out = 100)
  b <- m.lin$beta[["iwisescore"]]

  # Spline curve relative to IWISE = 0, with 95% CI (same basis as in the model)
  bn  <- grep("ns(iwisescore", names(m.s$beta), fixed = TRUE, value = TRUE)
  B   <- ns(data$iwisescore, df = df)
  X   <- sweep(predict(B, x), 2, predict(B, 0))
  fit <- drop(X %*% m.s$beta[bn])
  se  <- sqrt(rowSums((X %*% vcov(m.s)[bn, bn]) * X))

  # Bin estimates vs bin "0", placed at the mean score within each bin
  bins <- data %>% group_by(iw_bin) %>% summarise(x = mean(iwisescore), n = n()) %>%
    mutate(term = paste0("iw_bin", iw_bin),
           est  = c(0, m.b$beta[term[-1]]),
           se   = c(0, sqrt(diag(vcov(m.b))[term[-1]])))

  lr  <- 2 * as.numeric(logLik(m.s) - logLik(m.lin))
  tab <- tibble(outcome     = y,
                LR_spline   = lr,                                  # spline vs linear
                p_spline    = pchisq(lr, df - 1, lower.tail = FALSE),
                dAIC_spline = AIC(m.s) - AIC(m.lin),               # < 0 = spline better
                dAIC_bins   = AIC(m.b) - AIC(m.lin),               # < 0 = bins better
                max_dev_OR  = exp(max(abs(fit - b * x))),          # largest gap curve vs line, as OR
                OR_full_range = exp(b * max(x)))                   # total linear effect, 0 -> max, for scale

  p <- ggplot(tibble(x, fit, se), aes(x)) +
    geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), alpha = 0.15) +
    geom_line(aes(y = fit, colour = "spline")) +
    geom_line(aes(y = b * x, colour = "linear")) +
    geom_pointrange(data = bins, aes(x = x, y = est, ymin = est - 1.96 * se,
                                     ymax = est + 1.96 * se, colour = "bins"),
                    inherit.aes = FALSE, size = 0.3) +
    geom_vline(xintercept = c(2, 11, 23), linetype = 3, colour = "grey60") +  # category cut-offs
    geom_rug(data = data %>% distinct(iwisescore), aes(x = iwisescore), inherit.aes = FALSE) +
    labs(title = paste(y, "- shape of the IWISE effect"), x = "IWISE score",
         y = "Log-odds vs IWISE = 0", colour = NULL) +
    theme_minimal()

  list(table = tab, plot = p, m.spline = m.s, m.bins = m.b)
}

shape.fl  <- iwise_shape(d.fl.cc,  "FLI_3item", m.fl.1.ri)
shape.inc <- iwise_shape(d.inc.cc, "INCOME_5",  m.inc.1.ri)
saveRDS(list(shape.fl = shape.fl, shape.inc = shape.inc), "iwise_shape.rds")

bind_rows(shape.fl$table, shape.inc$table)
shape.fl$plot
shape.inc$plot
# Read: do the bins and the spline follow the straight line where most
# respondents are (rug)? Compare max_dev_OR with OR_full_range: a gap of
# 1.03 against a total effect of 0.30 is negligible; a gap of 1.3 is not.


# ---- 8c. PIECEWISE LINEAR (KNOTS AT THE CATEGORY CUT-OFFS) -------------------
# A straight line that may change slope at 2, 11 and 23, the upper limits of
# the IWISE categories (Frongillo et al. 2024: 0-2, 3-11, 12-23, 24-36).
# Answers: does the IWISE effect differ between the categories, without
# jumps at the cut-offs? Knots are fixed in advance, not chosen from the data.
# Hinge terms: iw_k1 = 0 up to 2, then rises 1 per point (etc.). The
# coefficient of each hinge = change in slope after that knot.
add_hinges <- function(d) mutate(d,
  iw_k1 = pmax(iwisescore - 2,  0),                        # extra slope above 2
  iw_k2 = pmax(iwisescore - 11, 0),                        # extra slope above 11
  iw_k3 = pmax(iwisescore - 23, 0))                        # extra slope above 23
d.fl.cc  <- add_hinges(d.fl.cc)
d.inc.cc <- add_hinges(d.inc.cc)

rhs_pw <- str_replace(rhs_ctx, "\\biwisescore\\b", "iwisescore + iw_k1 + iw_k2 + iw_k3")

# Same random part as the final continuous models (re_slope), so the linear
# model m.*.1 is nested in the piecewise one (= LR test valid).
m.fl.1.pw  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_pw, re_slope)
m.inc.1.pw <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_pw, re_slope)
saveRDS(mget(c("m.fl.1.pw", "m.inc.1.pw")), "models_piecewise.rds")
# Reload instead of refitting: load_models("models_piecewise.rds")

# Is piecewise better than linear? (nested, 3 extra fixed parameters)
lr3 <- function(m0, m1) {
  s <- 2 * as.numeric(logLik(m1) - logLik(m0))
  c(LR = s, p = pchisq(s, 3, lower.tail = FALSE), dAIC = AIC(m1) - AIC(m0))
}
round(rbind(fl = lr3(m.fl.1, m.fl.1.pw), inc = lr3(m.inc.1, m.inc.1.pw)), 4)

# OR between any two IWISE values (adds up the slopes of the segments crossed)
hinge <- function(x) c(iwisescore = x, iw_k1 = max(x - 2, 0),
                       iw_k2 = max(x - 11, 0), iw_k3 = max(x - 23, 0))
contrast_or <- function(m, x1, x2) {
  bn  <- names(hinge(0))
  L   <- hinge(x2) - hinge(x1)                             # weights of the 4 coefficients
  est <- sum(L * m$beta[bn])
  se  <- sqrt(drop(t(L) %*% vcov(m)[bn, bn] %*% L))
  c(x1 = x1, x2 = x2, OR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
}

# 8c-i. OR per IWISE point within each category range. Similar values = the
#       slope does not change between categories = linear is enough.
seg <- list(no_marginal = c(0, 1), low = c(5, 6), moderate = c(15, 16), high = c(28, 29))
t(sapply(seg, function(v) contrast_or(m.fl.1.pw,  v[1], v[2])))   # FLI
t(sapply(seg, function(v) contrast_or(m.inc.1.pw, v[1], v[2])))   # income

# 8c-ii. Each category vs no-to-marginal, at the median score per category
iw_ref <- function(data) {
  data %>%
    mutate(cat = cut(iwisescore, c(-Inf, 2, 11, 23, Inf),
                     labels = c("No-to-marginal", "Low", "Moderate", "High"))) %>%
    group_by(cat) %>%
    summarise(n = n(), iwise = median(iwisescore))
}
ref.fl  <- iw_ref(d.fl.cc)                                 # check n: a sparse "High" category
ref.inc <- iw_ref(d.inc.cc)                                # makes the last segment unstable
ref.fl
t(sapply(ref.fl$iwise,  function(x) contrast_or(m.fl.1.pw,  ref.fl$iwise[1],  x)))  # FLI
t(sapply(ref.inc$iwise, function(x) contrast_or(m.inc.1.pw, ref.inc$iwise[1], x)))  # income


# ---- 8d. DECISION -------------------------------------------------------------
# Linear is adequate if:
#   8a  categorical not clearly better than continuous (AIC)
#   8b  spline and bins close to the straight line; max_dev_OR small
#       relative to OR_full_range
#   8c  segment slopes similar; dAIC of piecewise small
#   -> continuous IWISE (07.1) stays the main model; categorical models
#      (Part A) are a sensitivity analysis.
# Clear bend in 8b and 8c:
#   -> piecewise (8c) or categorical (Part A) as the main model.
# DECISION: ...