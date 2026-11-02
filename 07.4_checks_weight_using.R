# =============================================================================
# 09_checks_wgt_using.R
# Output + assumption checks for the 8 models from 08_models_wgt_using.R.
# Run block by block: each model has its own section, so you can look at
#   (1) the model output, (2) proportional odds (PO), (3) linearity of the logit
# for one model before moving on to the next.
#
#   A  m.inc.0.wu       income, continuous IWISE, unadjusted
#   B  m.inc.1.wu       income, continuous IWISE, adjusted
#   C  m.fl.0.wu        FLI,    continuous IWISE, unadjusted
#   D  m.fl.1.wu        FLI,    continuous IWISE, adjusted
#   E  m.inc.0.cat.wu   income, categorical IWISE, unadjusted
#   F  m.inc.1.cat.wu   income, categorical IWISE, adjusted
#   G  m.fl.0.cat.wu    FLI,    categorical IWISE, unadjusted
#   H  m.fl.1.cat.wu    FLI,    categorical IWISE, adjusted
#
# Why the PO checks use models WITHOUT random effects:
#   brant() only works on polr() models and nominal_test() only on clm()
#   models; neither works on clmm(). So for each mixed model we fit the same
#   fixed part (same sample, same weights) without random effects and check
#   PO there. This is the standard workaround: if PO is violated in the
#   fixed-effects version, it is very likely violated in the mixed model too.
#
# Note on weights: with survey weights the log-likelihood is a
# pseudo-likelihood, so LR tests are approximate. Read them together with the
# plots and the change in the IWISE OR, not as strict yes/no tests.
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
library(tidyverse)
library(ordinal)                                           # clm(), clmm(), nominal_test()
library(brant)                                             # brant() test for polr models
library(splines)                                           # ns() natural splines
# MASS is NOT attached (it masks dplyr::select); we call MASS::polr() directly.


# ---- 1. DATA AND MODELS -----------------------------------------------------
# Same data preparation as in 08, so the samples are identical.
list2env(readRDS("data/iwise_data_prepared.rds"), envir = globalenv())  # d.fl, d.inc, rhs_ctx, v.*

prep <- function(data, y) {                                # analytic sample for one outcome
  data %>%
    filter(cc, !is.na(wgt_using)) %>%                      # complete cases; drops Palestine 2025
    mutate(country_year = factor(.data[[v.cluster]]),      # grouping factor
           iwise_cat    = factor(iwise_cat, ordered = FALSE),  # one dummy per category
           iwise_num    = as.integer(iwise_cat) - 1,       # 0,1,2,3: linear version of the categories
           "{y}" := factor(.data[[y]], ordered = TRUE),    # ordinal models need an ordered outcome
           w_fit        = wgt_using)                       # not rescaled
}

d.inc.cc <- prep(d.inc, "INCOME_5")
d.fl.cc  <- prep(d.fl,  "FLI_3item")

# Fixed parts and random parts, exactly as in 08
rhs_ctx_cat  <- str_replace(rhs_ctx, "\\biwisescore\\b", "iwise_cat")
re_term      <- "(1 | country_year)"
re_slope     <- "(1 + iwisescore | country_year)"
re_slope_cat <- "(1 + iwise_cat | country_year)"
re_slope_num <- "(1 + iwise_num | country_year)"

make_f   <- function(y, rhs, re) as.formula(paste(y, "~", rhs, "+", re))
fit_clmm <- function(data, y, rhs, re) {                   # same settings as 08
  clmm(make_f(y, rhs, re), data = data, weights = w_fit,
       link = "logit", Hess = TRUE, nAGQ = 1)
}

# Load the fitted models (no refitting needed)
list2env(readRDS("models_wgt_using_cont.rds"), envir = globalenv())
list2env(readRDS("models_wgt_using_cat.rds"),  envir = globalenv())


# ---- 2. HELPER FUNCTIONS ----------------------------------------------------

# 2a. Odds ratios with 95% Wald CI for the fixed effects of a clmm model.
or_table <- function(m) {
  co <- coef(summary(m))[names(m$beta), , drop = FALSE]    # keep predictors, drop thresholds
  tibble(term = rownames(co),
         OR   = exp(co[, 1]),                              # odds ratio
         lo   = exp(co[, 1] - 1.96 * co[, 2]),             # lower 95% limit
         hi   = exp(co[, 1] + 1.96 * co[, 2]),             # upper 95% limit
         p    = co[, 4])                                   # Wald p-value
}

# 2b. One binary logit per cut of the outcome ("above category k").
#     Under PO, every predictor should have roughly the same coefficient at
#     every cut. rhs = fixed part of the model being checked.
fit_cuts <- function(data, y, rhs) {
  K <- nlevels(data[[y]]) - 1                              # number of cuts = levels - 1
  map(set_names(seq_len(K), paste0("above_", seq_len(K))), function(k) {
    d.k <- mutate(data, y_bin = as.integer(as.integer(.data[[y]]) > k))  # 1 = above category k
    glm(as.formula(paste("y_bin ~", rhs)),                 # same fixed part as the model
        data = d.k, weights = w_fit,                       # same sample and weights
        family = quasibinomial)                            # quasi: no warnings for non-integer weights
  })
}

# 2c. Plot of the cut-specific coefficients against the single PO estimate.
#     Read: points roughly flat and around the red dashed line = PO fine.
#     A clear trend across cuts = PO violated for that predictor.
#     Look at size and direction of the trend, not at p-values (n is large).
plot_cuts <- function(cuts, m.clm, title) {
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
    scale_x_discrete(labels = function(x) sub("above_", "", x)) +  # works for any number of cuts
    labs(title = title, x = "Cut (above category k)", y = "Log-odds (95% CI)") +
    theme_minimal() +
    theme(strip.text = element_text(size = 7))
}

# 2d. Linearity check for ONE continuous variable.
#     Refits the mixed model with `var` replaced by a natural spline (df = 3),
#     keeping everything else (fixed part, random part, weights) the same.
#     Returns:
#       table: LR test spline vs linear (small p = non-linear) and the IWISE
#              OR(s) in both models (big % change = non-linearity matters
#              for the IWISE result; small change = can keep linear)
#       plot : spline shape vs straight line, both relative to the median
spline_check <- function(data, y, m.lin, rhs, re, var, df = 3) {
  rhs.s <- str_replace(rhs, paste0("\\b", var, "\\b"),    # replace var by ns(var, df)
                       sprintf("ns(%s, df = %d)", var, df))
  m.s   <- fit_clmm(data, y, rhs.s, re)                    # refit (SLOW)

  lr <- 2 * (as.numeric(logLik(m.s)) - as.numeric(logLik(m.lin)))  # LR by hand (anova() fails here)
  iw <- grep("^iwise", names(m.lin$beta), value = TRUE)    # IWISE term(s): iwisescore or iwise_cat dummies
  tab <- tibble(outcome = y, var = var, LR = lr,
                p       = pchisq(lr, df - 1, lower.tail = FALSE),  # df - 1 extra parameters
                term    = iw,                              # one row per IWISE term
                OR_lin  = exp(m.lin$beta[iw]),             # IWISE OR with var linear
                OR_spl  = exp(m.s$beta[iw]),               # IWISE OR with var as spline (NA if var is IWISE itself)
                OR_change_pct = 100 * (OR_spl / OR_lin - 1))  # % change in IWISE OR: the key number

  bn  <- grep(paste0("ns(", var), names(m.s$beta), fixed = TRUE, value = TRUE)  # spline coefficient names
  x   <- seq(min(data[[var]]), max(data[[var]]), length.out = 100)  # grid over observed range
  b   <- ns(data[[var]], df = df)                          # same spline basis as in the model
  X   <- sweep(predict(b, x), 2, predict(b, median(data[[var]])))  # basis at x minus basis at median
  fit <- drop(X %*% m.s$beta[bn])                          # spline effect relative to the median
  se  <- sqrt(rowSums((X %*% vcov(m.s)[bn, bn]) * X))      # its standard error at each x

  p <- ggplot(tibble(x, fit, se), aes(x)) +
    geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), alpha = 0.2) +  # 95% CI of spline
    geom_line(aes(y = fit, colour = "spline")) +
    geom_line(aes(y = m.lin$beta[[var]] * (x - median(data[[var]])), colour = "linear")) +
    geom_rug(data = tibble(x = unique(data[[var]])), aes(x), inherit.aes = FALSE) +  # where the data are
    labs(title = paste(y, "-", var), x = var, y = "Log-odds vs median", colour = NULL) +
    theme_minimal()

  list(table = tab, plot = p, model = m.s)                 # keep the spline model too: no refitting later
}

# 2e. Linearity check for CATEGORICAL IWISE.
#     Refits the model with iwise_cat replaced by iwise_num (0,1,2,3 = one
#     linear trend), same random part. The two models are nested, so an LR
#     test tells whether the categories deviate from a straight line.
#     The plot shows the category estimates next to the linear trend.
cat_linearity_check <- function(data, y, m.cat, rhs, re) {
  m.n <- fit_clmm(data, y, str_replace(rhs, "\\biwise_cat\\b", "iwise_num"), re)  # refit (SLOW)

  df  <- nlevels(data$iwise_cat) - 2                       # (categories - 1) dummies vs 1 slope
  lr  <- 2 * (as.numeric(logLik(m.cat)) - as.numeric(logLik(m.n)))
  tab <- tibble(outcome = y, LR = lr, df = df,
                p = pchisq(lr, df, lower.tail = FALSE))    # small p = categories not linear

  bn  <- grep("^iwise_cat", names(m.cat$beta), value = TRUE)  # category dummies
  co  <- coef(summary(m.cat))[bn, , drop = FALSE]
  est <- tibble(k  = 0:(length(bn)),                       # 0 = reference category
                b  = c(0, co[, 1]),                        # reference has log-odds 0
                se = c(0, co[, 2]))
  p <- ggplot(est, aes(k, b)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_abline(intercept = 0, slope = m.n$beta[["iwise_num"]],  # linear-trend model
                colour = "red", linetype = 2) +
    geom_pointrange(aes(ymin = b - 1.96 * se, ymax = b + 1.96 * se)) +  # category estimates
    scale_x_continuous(breaks = est$k, labels = levels(data$iwise_cat)) +
    labs(title = paste(y, "- IWISE categories vs linear trend (red)"),
         x = "IWISE category", y = "Log-odds vs reference") +
    theme_minimal()

  list(table = tab, plot = p, model = m.n)
}


# ---- 3. FORMULAS FOR THE PO CHECKS (fixed part only, no random effects) -----
# Defined here at top level on purpose: nominal_test() and brant() re-evaluate
# the model call, which fails if the formula only exists inside a function.
f.inc.0     <- INCOME_5  ~ iwisescore
f.inc.1     <- as.formula(paste("INCOME_5 ~",  rhs_ctx))
f.fl.0      <- FLI_3item ~ iwisescore
f.fl.1      <- as.formula(paste("FLI_3item ~", rhs_ctx))
f.inc.0.cat <- INCOME_5  ~ iwise_cat
f.inc.1.cat <- as.formula(paste("INCOME_5 ~",  rhs_ctx_cat))
f.fl.0.cat  <- FLI_3item ~ iwise_cat
f.fl.1.cat  <- as.formula(paste("FLI_3item ~", rhs_ctx_cat))


# =============================================================================
# A. m.inc.0.wu: INCOME, continuous IWISE, unadjusted, random intercept
# =============================================================================

# ---- A1. Output ----
summary(m.inc.0.wu)                                        # coefficients, thresholds, random-effect SD
or_table(m.inc.0.wu)                                       # IWISE OR with 95% CI

# ---- A2. Proportional odds ----
p.inc.0 <- MASS::polr(f.inc.0, data = d.inc.cc, weights = w_fit, Hess = TRUE)  # needed for brant()
brant(p.inc.0)                                             # significant = PO violated (very sensitive with large n)

c.inc.0 <- clm(f.inc.0, data = d.inc.cc, weights = w_fit, link = "logit")      # needed for nominal_test()
nominal_test(c.inc.0)                                      # LR test: lets IWISE vary across cuts; small p = PO violated

cuts.inc.0 <- fit_cuts(d.inc.cc, "INCOME_5", "iwisescore") # one binary logit per cut
plot_cuts(cuts.inc.0, c.inc.0, "A. Income M0 (cont.): coefficients across cuts")

# ---- A3. Linearity of the logit ----
# Only continuous predictor: iwisescore itself.
lin.inc.0.iwise <- spline_check(d.inc.cc, "INCOME_5", m.inc.0.wu, "iwisescore", re_term, "iwisescore")
lin.inc.0.iwise$table                                      # small p = IWISE effect not linear
lin.inc.0.iwise$plot                                       # spline close to straight line = linear OK


# =============================================================================
# B. m.inc.1.wu: INCOME, continuous IWISE, adjusted, random IWISE slope
# =============================================================================

# ---- B1. Output ----
summary(m.inc.1.wu)
or_table(m.inc.1.wu)                                       # all fixed effects as ORs
VarCorr(m.inc.1.wu)                                        # random intercept and IWISE slope SDs

# ---- B2. Proportional odds ----
p.inc.1 <- MASS::polr(f.inc.1, data = d.inc.cc, weights = w_fit, Hess = TRUE)
brant(p.inc.1)                                             # omnibus + one test per predictor
# If brant() errors (possible with many factor levels), rely on nominal_test and the plot.

c.inc.1 <- clm(f.inc.1, data = d.inc.cc, weights = w_fit, link = "logit")
nominal_test(c.inc.1)                                      # one LR test per predictor

cuts.inc.1 <- fit_cuts(d.inc.cc, "INCOME_5", rhs_ctx)
plot_cuts(cuts.inc.1, c.inc.1, "B. Income M1 (cont.): coefficients across cuts")

# ---- B3. Linearity of the logit (SLOW: one refit per variable) ----
# IWISE itself
lin.inc.1.iwise <- spline_check(d.inc.cc, "INCOME_5", m.inc.1.wu, rhs_ctx, re_slope, "iwisescore")
lin.inc.1.iwise$table
lin.inc.1.iwise$plot

# JMP basic sanitation rate
lin.inc.1.jmp <- spline_check(d.inc.cc, "INCOME_5", m.inc.1.wu, rhs_ctx, re_slope, "bas_rate_pp")
lin.inc.1.jmp$table                                        # look at p AND OR_change_pct
lin.inc.1.jmp$plot

# Household size
lin.inc.1.hhsize <- spline_check(d.inc.cc, "INCOME_5", m.inc.1.wu, rhs_ctx, re_slope, "hhsize_c")
lin.inc.1.hhsize$table
lin.inc.1.hhsize$plot

# Log GDP
lin.inc.1.gdp <- spline_check(d.inc.cc, "INCOME_5", m.inc.1.wu, rhs_ctx, re_slope, "log_gdp_c")
lin.inc.1.gdp$table
lin.inc.1.gdp$plot

# Government effectiveness
lin.inc.1.wgi <- spline_check(d.inc.cc, "INCOME_5", m.inc.1.wu, rhs_ctx, re_slope, "wgi_ge_sc_c")
lin.inc.1.wgi$table
lin.inc.1.wgi$plot

# Enteric DALYs
lin.inc.1.daly <- spline_check(d.inc.cc, "INCOME_5", m.inc.1.wu, rhs_ctx, re_slope, "daly_enteric_model_c")
lin.inc.1.daly$table
lin.inc.1.daly$plot

# Overview of all linearity checks for this model
bind_rows(lin.inc.1.iwise$table, lin.inc.1.jmp$table, lin.inc.1.hhsize$table,
          lin.inc.1.gdp$table,   lin.inc.1.wgi$table, lin.inc.1.daly$table)


# =============================================================================
# C. m.fl.0.wu: FLI, continuous IWISE, unadjusted, random IWISE slope
# =============================================================================

# ---- C1. Output ----
summary(m.fl.0.wu)
or_table(m.fl.0.wu)
VarCorr(m.fl.0.wu)

# ---- C2. Proportional odds ----
p.fl.0 <- MASS::polr(f.fl.0, data = d.fl.cc, weights = w_fit, Hess = TRUE)
brant(p.fl.0)

c.fl.0 <- clm(f.fl.0, data = d.fl.cc, weights = w_fit, link = "logit")
nominal_test(c.fl.0)

cuts.fl.0 <- fit_cuts(d.fl.cc, "FLI_3item", "iwisescore")
plot_cuts(cuts.fl.0, c.fl.0, "C. FLI M0 (cont.): coefficients across cuts")

# ---- C3. Linearity of the logit ----
lin.fl.0.iwise <- spline_check(d.fl.cc, "FLI_3item", m.fl.0.wu, "iwisescore", re_slope, "iwisescore")
lin.fl.0.iwise$table
lin.fl.0.iwise$plot


# =============================================================================
# D. m.fl.1.wu: FLI, continuous IWISE, adjusted, random IWISE slope
# =============================================================================

# ---- D1. Output ----
summary(m.fl.1.wu)
or_table(m.fl.1.wu)
VarCorr(m.fl.1.wu)

# ---- D2. Proportional odds ----
p.fl.1 <- MASS::polr(f.fl.1, data = d.fl.cc, weights = w_fit, Hess = TRUE)
brant(p.fl.1)

c.fl.1 <- clm(f.fl.1, data = d.fl.cc, weights = w_fit, link = "logit")
nominal_test(c.fl.1)

cuts.fl.1 <- fit_cuts(d.fl.cc, "FLI_3item", rhs_ctx)
plot_cuts(cuts.fl.1, c.fl.1, "D. FLI M1 (cont.): coefficients across cuts")

# ---- D3. Linearity of the logit (SLOW) ----
lin.fl.1.iwise <- spline_check(d.fl.cc, "FLI_3item", m.fl.1.wu, rhs_ctx, re_slope, "iwisescore")
lin.fl.1.iwise$table
lin.fl.1.iwise$plot

lin.fl.1.jmp <- spline_check(d.fl.cc, "FLI_3item", m.fl.1.wu, rhs_ctx, re_slope, "bas_rate_pp")
lin.fl.1.jmp$table
lin.fl.1.jmp$plot

lin.fl.1.hhsize <- spline_check(d.fl.cc, "FLI_3item", m.fl.1.wu, rhs_ctx, re_slope, "hhsize_c")
lin.fl.1.hhsize$table
lin.fl.1.hhsize$plot

lin.fl.1.gdp <- spline_check(d.fl.cc, "FLI_3item", m.fl.1.wu, rhs_ctx, re_slope, "log_gdp_c")
lin.fl.1.gdp$table
lin.fl.1.gdp$plot

lin.fl.1.wgi <- spline_check(d.fl.cc, "FLI_3item", m.fl.1.wu, rhs_ctx, re_slope, "wgi_ge_sc_c")
lin.fl.1.wgi$table
lin.fl.1.wgi$plot

lin.fl.1.daly <- spline_check(d.fl.cc, "FLI_3item", m.fl.1.wu, rhs_ctx, re_slope, "daly_enteric_model_c")
lin.fl.1.daly$table
lin.fl.1.daly$plot

bind_rows(lin.fl.1.iwise$table, lin.fl.1.jmp$table, lin.fl.1.hhsize$table,
          lin.fl.1.gdp$table,   lin.fl.1.wgi$table, lin.fl.1.daly$table)


# =============================================================================
# E. m.inc.0.cat.wu: INCOME, categorical IWISE, unadjusted, random intercept
# =============================================================================

# ---- E1. Output ----
summary(m.inc.0.cat.wu)
or_table(m.inc.0.cat.wu)                                   # OR of each category vs lowest WI

# ---- E2. Proportional odds ----
p.inc.0.cat <- MASS::polr(f.inc.0.cat, data = d.inc.cc, weights = w_fit, Hess = TRUE)
brant(p.inc.0.cat)

c.inc.0.cat <- clm(f.inc.0.cat, data = d.inc.cc, weights = w_fit, link = "logit")
nominal_test(c.inc.0.cat)                                  # one test for iwise_cat as a whole

cuts.inc.0.cat <- fit_cuts(d.inc.cc, "INCOME_5", "iwise_cat")
plot_cuts(cuts.inc.0.cat, c.inc.0.cat, "E. Income M0 (cat.): coefficients across cuts")

# ---- E3. Linearity ----
# No continuous predictor in this model. The question here is whether the
# IWISE categories follow a straight line (then continuous IWISE is fine).
lin.inc.0.cat <- cat_linearity_check(d.inc.cc, "INCOME_5", m.inc.0.cat.wu, "iwise_cat", re_term)
lin.inc.0.cat$table                                        # small p = categories deviate from a line
lin.inc.0.cat$plot                                         # points on the red line = linear trend OK


# =============================================================================
# F. m.inc.1.cat.wu: INCOME, categorical IWISE, adjusted, random slope per category
# =============================================================================

# ---- F1. Output ----
summary(m.inc.1.cat.wu)
or_table(m.inc.1.cat.wu)
VarCorr(m.inc.1.cat.wu)

# ---- F2. Proportional odds ----
p.inc.1.cat <- MASS::polr(f.inc.1.cat, data = d.inc.cc, weights = w_fit, Hess = TRUE)
brant(p.inc.1.cat)

c.inc.1.cat <- clm(f.inc.1.cat, data = d.inc.cc, weights = w_fit, link = "logit")
nominal_test(c.inc.1.cat)

cuts.inc.1.cat <- fit_cuts(d.inc.cc, "INCOME_5", rhs_ctx_cat)
plot_cuts(cuts.inc.1.cat, c.inc.1.cat, "F. Income M1 (cat.): coefficients across cuts")

# ---- F3. Linearity (SLOW) ----
# IWISE categories vs linear trend
lin.inc.1.cat <- cat_linearity_check(d.inc.cc, "INCOME_5", m.inc.1.cat.wu, rhs_ctx_cat, re_slope_cat)
lin.inc.1.cat$table
lin.inc.1.cat$plot

# Covariates: same checks as in B. Results should be very close to B; the
# table now shows the change in each IWISE category OR. Skip if time is short.
lin.inc.1.cat.jmp <- spline_check(d.inc.cc, "INCOME_5", m.inc.1.cat.wu, rhs_ctx_cat, re_slope_cat, "bas_rate_pp")
lin.inc.1.cat.jmp$table
lin.inc.1.cat.jmp$plot

lin.inc.1.cat.hhsize <- spline_check(d.inc.cc, "INCOME_5", m.inc.1.cat.wu, rhs_ctx_cat, re_slope_cat, "hhsize_c")
lin.inc.1.cat.hhsize$table
lin.inc.1.cat.hhsize$plot

lin.inc.1.cat.gdp <- spline_check(d.inc.cc, "INCOME_5", m.inc.1.cat.wu, rhs_ctx_cat, re_slope_cat, "log_gdp_c")
lin.inc.1.cat.gdp$table
lin.inc.1.cat.gdp$plot

lin.inc.1.cat.wgi <- spline_check(d.inc.cc, "INCOME_5", m.inc.1.cat.wu, rhs_ctx_cat, re_slope_cat, "wgi_ge_sc_c")
lin.inc.1.cat.wgi$table
lin.inc.1.cat.wgi$plot

lin.inc.1.cat.daly <- spline_check(d.inc.cc, "INCOME_5", m.inc.1.cat.wu, rhs_ctx_cat, re_slope_cat, "daly_enteric_model_c")
lin.inc.1.cat.daly$table
lin.inc.1.cat.daly$plot

bind_rows(lin.inc.1.cat.jmp$table, lin.inc.1.cat.hhsize$table, lin.inc.1.cat.gdp$table,
          lin.inc.1.cat.wgi$table, lin.inc.1.cat.daly$table)


# =============================================================================
# G. m.fl.0.cat.wu: FLI, categorical IWISE, unadjusted, linear random slope
# =============================================================================

# ---- G1. Output ----
summary(m.fl.0.cat.wu)
or_table(m.fl.0.cat.wu)
VarCorr(m.fl.0.cat.wu)

# ---- G2. Proportional odds ----
p.fl.0.cat <- MASS::polr(f.fl.0.cat, data = d.fl.cc, weights = w_fit, Hess = TRUE)
brant(p.fl.0.cat)

c.fl.0.cat <- clm(f.fl.0.cat, data = d.fl.cc, weights = w_fit, link = "logit")
nominal_test(c.fl.0.cat)

cuts.fl.0.cat <- fit_cuts(d.fl.cc, "FLI_3item", "iwise_cat")
plot_cuts(cuts.fl.0.cat, c.fl.0.cat, "G. FLI M0 (cat.): coefficients across cuts")

# ---- G3. Linearity ----
lin.fl.0.cat <- cat_linearity_check(d.fl.cc, "FLI_3item", m.fl.0.cat.wu, "iwise_cat", re_slope_num)
lin.fl.0.cat$table
lin.fl.0.cat$plot


# =============================================================================
# H. m.fl.1.cat.wu: FLI, categorical IWISE, adjusted, random slope per category
# =============================================================================

# ---- H1. Output ----
summary(m.fl.1.cat.wu)
or_table(m.fl.1.cat.wu)
VarCorr(m.fl.1.cat.wu)

# ---- H2. Proportional odds ----
p.fl.1.cat <- MASS::polr(f.fl.1.cat, data = d.fl.cc, weights = w_fit, Hess = TRUE)
brant(p.fl.1.cat)

c.fl.1.cat <- clm(f.fl.1.cat, data = d.fl.cc, weights = w_fit, link = "logit")
nominal_test(c.fl.1.cat)

cuts.fl.1.cat <- fit_cuts(d.fl.cc, "FLI_3item", rhs_ctx_cat)
plot_cuts(cuts.fl.1.cat, c.fl.1.cat, "H. FLI M1 (cat.): coefficients across cuts")

# ---- H3. Linearity (SLOW) ----
lin.fl.1.cat <- cat_linearity_check(d.fl.cc, "FLI_3item", m.fl.1.cat.wu, rhs_ctx_cat, re_slope_cat)
lin.fl.1.cat$table
lin.fl.1.cat$plot

# Covariates: should be very close to D. Skip if time is short.
lin.fl.1.cat.jmp <- spline_check(d.fl.cc, "FLI_3item", m.fl.1.cat.wu, rhs_ctx_cat, re_slope_cat, "bas_rate_pp")
lin.fl.1.cat.jmp$table
lin.fl.1.cat.jmp$plot

lin.fl.1.cat.hhsize <- spline_check(d.fl.cc, "FLI_3item", m.fl.1.cat.wu, rhs_ctx_cat, re_slope_cat, "hhsize_c")
lin.fl.1.cat.hhsize$table
lin.fl.1.cat.hhsize$plot

lin.fl.1.cat.gdp <- spline_check(d.fl.cc, "FLI_3item", m.fl.1.cat.wu, rhs_ctx_cat, re_slope_cat, "log_gdp_c")
lin.fl.1.cat.gdp$table
lin.fl.1.cat.gdp$plot

lin.fl.1.cat.wgi <- spline_check(d.fl.cc, "FLI_3item", m.fl.1.cat.wu, rhs_ctx_cat, re_slope_cat, "wgi_ge_sc_c")
lin.fl.1.cat.wgi$table
lin.fl.1.cat.wgi$plot

lin.fl.1.cat.daly <- spline_check(d.fl.cc, "FLI_3item", m.fl.1.cat.wu, rhs_ctx_cat, re_slope_cat, "daly_enteric_model_c")
lin.fl.1.cat.daly$table
lin.fl.1.cat.daly$plot

bind_rows(lin.fl.1.cat.jmp$table, lin.fl.1.cat.hhsize$table, lin.fl.1.cat.gdp$table,
          lin.fl.1.cat.wgi$table, lin.fl.1.cat.daly$table)


# =============================================================================
# 4. SAVE
# =============================================================================
# PO objects (fast to refit, but saved for convenience)
saveRDS(mget(ls(pattern = "^(p|c|cuts)\\.(inc|fl)\\.")), "checks_po_wgt_using.rds")
# Linearity results incl. the refitted spline / linear-trend models (slow to refit)
saveRDS(mget(ls(pattern = "^lin\\.")), "checks_linearity_wgt_using.rds")
# Reload later with: list2env(readRDS("checks_linearity_wgt_using.rds"), envir = globalenv())