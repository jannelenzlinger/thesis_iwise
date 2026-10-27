# =============================================================================
# 07.2_models_iwise_cat.R
# Final models from 07 (Section 3e) refitted with CATEGORICAL IWISE
#   iwise_cat replaces iwisescore. Everything else as in 07.
#   M0 : unadjusted     = iwise_cat
#   M1 : fully adjusted = rhs_ctx with iwisescore -> iwise_cat
#   Weights: w_fit. Outcomes: FLI_3item (d.fl), INCOME_5 (d.inc).
#
# Random effects: RANDOM IWISE SLOPE (as in 07.1, where it was clearly needed).
#   full   (1 + iwise_cat | country_year)  each category its own random
#          slope (9 extra parameters vs random intercept)
#   linear (1 + iwise_num | country_year)  one random slope across the
#          ordered categories, iwise_num = 0,1,2,3 (2 extra parameters,
#          same as in 07). Fixed part stays fully categorical.
#   fl0, fl1, inc1: full slope. inc0: linear slope (even one slope was
#   not identifiable in 07.1). If a fit fails the checks in Section 5,
#   refit it with the linear slope (Section 4b).
#   Random-intercept fits (Section 3) are only the reference for the LR test.
#
#   1  Data                         5  Checks and LR tests
#   2  Helpers                      6  Final models
#   3  Random-intercept reference   7  Results
#   4  Random-slope models
# Saves: models_cat_ri.rds, models_cat_slope.rds,
#        models_cat_final.rds (m.fl.0.cat, m.fl.1.cat, m.inc.0.cat,
#                              m.inc.1.cat, slope_type, re_cat) -> used by 09
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
library(tidyverse)                   # dplyr, purrr, stringr, ggplot2
library(ordinal)                     # clmm()


# ---- 1. DATA (same as 07.1, Section 1, + iwise_num) ---------------------------
list2env(readRDS("data/iwise_data_prepared.rds"), envir = globalenv())

prep <- function(data, y) {                                # analytic sample for one outcome
  data %>%
    filter(cc) %>%                                         # complete cases
    mutate(.id          = row_number(),
           country_year = factor(.data[[v.cluster]]),
           iwise_cat    = factor(iwise_cat, ordered = FALSE),  # unordered: treatment contrasts
           iwise_num    = as.integer(iwise_cat) - 1,       # 0,1,2,3 for the linear random slope
           "{y}" := factor(.data[[y]], ordered = TRUE))
}

d.fl.cc  <- prep(d.fl,  "FLI_3item")
d.inc.cc <- prep(d.inc, "INCOME_5")

# Levels must run from lowest to highest WI: the first level is the
# reference category AND iwise_num = 0. Check both:
levels(d.fl.cc$iwise_cat)
table(d.fl.cc$iwise_cat, d.fl.cc$iwise_num)                # check: diagonal, no sparse categories
table(d.inc.cc$iwise_cat)


# ---- 2. FORMULAS AND HELPERS ------------------------------------------------
rhs_unadj_cat <- "iwise_cat"
rhs_ctx_cat   <- str_replace(rhs_ctx, "\\biwisescore\\b", "iwise_cat")
rhs_ctx_cat                                                # check: iwise_cat in, iwisescore out

re_term      <- "(1 | country_year)"                       # random intercept (reference only)
re_slope_cat <- "(1 + iwise_cat | country_year)"           # full: random slope per category
re_slope_num <- "(1 + iwise_num | country_year)"           # linear: one random slope

make_f <- function(y, rhs, re) as.formula(paste(y, "~", rhs, "+", re))

fit_clmm <- function(data, y, rhs, re) {                   # same settings as 07
  clmm(make_f(y, rhs, re), data = data, weights = w_fit,
       link = "logit", Hess = TRUE, nAGQ = 1)              # random slopes need nAGQ = 1
}

tidy_or <- function(m, label) {                            # as in 07
  co <- summary(m)$coefficients[names(m$beta), , drop = FALSE]
  tibble(model = label, term = rownames(co),
         b = co[, 1], se = co[, 2], OR = exp(b),
         lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = co[, 4])
}

load_models <- function(file) {                            # as in 07
  obj <- readRDS(file)
  list2env(obj, envir = globalenv())
  message("Loaded from ", file, ": ", paste(names(obj), collapse = ", "))
  invisible(names(obj))
}

icc <- function(m) {                                       # latent-scale ICC; with a random slope
  s2 <- VarCorr(m)[[1]][1, 1]                              # it applies at the reference category
  s2 / (s2 + pi^2 / 3)
}


# ---- 3. RANDOM-INTERCEPT REFERENCE MODELS (NOT FINAL) -----------------------
# Only needed for the LR test in Section 5c (is the slope needed?).
# The adjusted ones already exist from 07, check 3: reused, not refitted.
load_models("models_iwise_cat.rds")                        # m.fl.cat, m.inc.cat (07, check 3)
m.fl.0.cat.ri  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_unadj_cat, re_term)  # FLI, unadjusted
m.fl.1.cat.ri  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx_cat,   re_term)  # FLI, adjusted
m.inc.0.cat.ri <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_unadj_cat, re_term)  # income, unadjusted
m.inc.1.cat.ri <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx_cat,   re_term)  # income, adjusted

saveRDS(mget(c("m.fl.0.cat.ri", "m.fl.1.cat.ri", "m.inc.0.cat.ri", "m.inc.1.cat.ri")),
        "models_cat_ri.rds")
# Reload instead of refitting: load_models("models_cat_ri.rds")


saveRDS(mget(c("m.fl.0.cat.ri", "m.fl.1.cat.ri", "m.inc.0.cat.ri", "m.inc.1.cat.ri")),
        "models_cat_ri.rds")
# Reload instead of refitting: load_models("models_cat_ri.rds")


# ---- 4. RANDOM-SLOPE MODELS -------------------------------------------------
# SLOW (one fit can take long).
m.fl.0.cat.rs  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_unadj_cat, re_slope_cat)  # FLI, unadjusted
m.fl.1.cat.rs  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx_cat,   re_slope_cat)  # FLI, adjusted
m.inc.0.cat.rs <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_unadj_cat, re_slope_num)  # income, unadjusted (linear)
m.inc.1.cat.rs <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx_cat,   re_slope_cat)  # income, adjusted

slope_type <- c(fl0 = "full", fl1 = "full", inc0 = "linear", inc1 = "full")  # keep in sync with the fits

saveRDS(mget(c("m.fl.0.cat.rs", "m.fl.1.cat.rs", "m.inc.0.cat.rs", "m.inc.1.cat.rs", "slope_type")),
        "models_cat_slope.rds")
# Reload instead of refitting: load_models("models_cat_slope.rds")

# 4b. FALLBACK: only for a model that fails the checks in Section 5.
#     Refit it with the linear slope and update slope_type, then save again.
#     Example for FLI adjusted:
# m.fl.1.cat.rs     <- fit_clmm(d.fl.cc, "FLI_3item", rhs_ctx_cat, re_slope_num)
# slope_type["fl1"] <- "linear"


# ---- 5. CHECKS AND LR TESTS (as in 07, 3d-2 and 3d-4) -----------------------
pairs <- list(fl0  = list(ri = m.fl.0.cat.ri,  rs = m.fl.0.cat.rs),
              fl1  = list(ri = m.fl.1.cat.ri,  rs = m.fl.1.cat.rs),
              inc0 = list(ri = m.inc.0.cat.ri, rs = m.inc.0.cat.rs),
              inc1 = list(ri = m.inc.1.cat.ri, rs = m.inc.1.cat.rs))

# 5a. Did the fits work? A fit can only be trusted if ALL rows look fine.
fit_quality <- function(p) {
  c(condHess    = suppressWarnings(summary(p$rs)$condHess),   # NaN or > 1e6 = not identified
    max_grad    = max(abs(p$rs$gradient)),                    # large = not converged
    n_na_se     = sum(is.na(summary(p$rs)$coefficients[names(p$rs$beta), 2])),  # 0 = all SEs defined
    logLik_gain = as.numeric(logLik(p$rs) - logLik(p$ri)))    # must be >= 0 (ri is nested)
}
round(sapply(pairs, fit_quality), 4)

# 5b. Random-effect SDs and correlations: plausible? Correlations of +-1
#     mean the slope is not identified -> use the fallback (4b).
lapply(pairs, function(p) VarCorr(p$rs)[[1]])

# 5c. Is the slope needed? LR test vs the random-intercept model.
#     linear: 50:50 mix of chi2(1) and chi2(2), as in 07 (3d-4).
#     full:   naive chi2, k(k+1)/2 - 1 df; conservative (variances tested at
#             their boundary, Snijders & Bosker), so read it with dAIC.
k <- nlevels(d.fl.cc$iwise_cat)
lr_slope <- function(p, type) {
  lr <- as.numeric(2 * (logLik(p$rs) - logLik(p$ri)))
  df <- if (type == "full") k * (k + 1) / 2 - 1 else 2
  p_val <- if (type == "full") pchisq(lr, df, lower.tail = FALSE) else
             0.5 * pchisq(lr, 1, lower.tail = FALSE) + 0.5 * pchisq(lr, 2, lower.tail = FALSE)
  c(LR = lr, df = df, p = p_val, dAIC = AIC(p$rs) - AIC(p$ri))   # dAIC < 0 = slope better
}
round(sapply(names(pairs), function(nm) lr_slope(pairs[[nm]], slope_type[[nm]])), 4)


# ---- 6. FINAL MODELS --------------------------------------------------------
m.fl.0.cat  <- m.fl.0.cat.rs                               # FLI, unadjusted
m.fl.1.cat  <- m.fl.1.cat.rs                               # FLI, adjusted
m.inc.0.cat <- m.inc.0.cat.rs                              # income, unadjusted
m.inc.1.cat <- m.inc.1.cat.rs                              # income, adjusted

# Random-effects structure of the adjusted models, per outcome (used in 09)
re_of  <- c(full = re_slope_cat, linear = re_slope_num)
re_cat <- c(fl = re_of[[slope_type[["fl1"]]]], inc = re_of[[slope_type[["inc1"]]]])
re_cat

saveRDS(list(m.fl.0.cat = m.fl.0.cat, m.fl.1.cat = m.fl.1.cat,
             m.inc.0.cat = m.inc.0.cat, m.inc.1.cat = m.inc.1.cat,
             slope_type = slope_type, re_cat = re_cat),
        "models_cat_final.rds")


# ---- 7. RESULTS -------------------------------------------------------------
# 7a. IWISE ORs per category (vs reference category): average effects
res.cat <- bind_rows(tidy_or(m.fl.0.cat,  "FLI M0"), tidy_or(m.fl.1.cat,  "FLI M1"),
                     tidy_or(m.inc.0.cat, "INC M0"), tidy_or(m.inc.1.cat, "INC M1"))
res.cat %>% filter(str_starts(term, "iwise_cat")) %>%
  select(model, term, OR, lo, hi, p) %>% print()

# 7b. Average IWISE ORs: random slope vs random intercept (cf. 07, 3d-5a).
#     CIs are expected to be wider with the slope.
bind_rows(tidy_or(m.fl.1.cat.ri,  "FLI M1 intercept only"), tidy_or(m.fl.1.cat,  "FLI M1 + slope"),
          tidy_or(m.inc.1.cat.ri, "INC M1 intercept only"), tidy_or(m.inc.1.cat, "INC M1 + slope")) %>%
  filter(str_starts(term, "iwise_cat")) %>%
  select(model, term, OR, lo, hi)

# 7c. Spread between country-years: range of 95% of the country-year ORs
#     per category (cf. 07, 3d-5b). full: SD per category; linear: the
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

# 7f. Continuous vs categorical IWISE (non-nested: AIC, lower = better).
#     Compared with the final continuous models from 07 (3e), which also
#     have a random slope. The random parts differ (the full slope has more
#     parameters), so this compares whole models.
load_models("models_final.rds")                            # m.fl.1, m.inc.1 (07, 3e)
AIC(m.fl.1,  m.fl.1.cat)
AIC(m.inc.1, m.inc.1.cat)

# 7g. Plot: adjusted ORs per category, both outcomes
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