# =============================================================================
# 07.8_final_interactions.R
#   Effect modification of the CONTINUOUS IWISE association (adjusted M1)
#   Predictor: iwisescore (continuous)   Outcomes: FLI_3item, INCOME_5
#   Modifiers: gender (female), employment, urbanicity (urban_imp)
#
#   Step A  ONE AT A TIME (primary): M1 + iwisescore:X, one model per X.
#           Each model answers one pre-specified question: "does the IWISE
#           association differ by X?" (effect modification, Knol &
#           VanderWeele 2012). Holm correction over the 3 tests, per outcome.
#   Step B  ALL TOGETHER (sensitivity): M1 + all three interactions.
#           First a global ("chunk") LR test of all interaction terms
#           (Harrell, RMS). Then each block adjusted for the others: an
#           included interaction can absorb an omitted one when modifiers
#           are correlated, e.g. gender and employment (Beiser-McGrath &
#           Beiser-McGrath 2020).
#   Report  Ratio of ORs (interaction terms) AND the IWISE OR within each
#           modifier group, with 95% CIs (Knol & VanderWeele 2012).
#           All ORs are per 1-point increase in the IWISE score (0-36).
#
#   Model:  clmm, logit link, weights = w_fit,
#           random intercept + random IWISE slope per country(-wave):
#           (1 + iwisescore | country_year).
#           M1 = rhs_ctx WITHOUT wgi_ge_sc_c (removed in Section 3).
#           The base model M1 is fitted here with the same data, weights and
#           random part, so every interaction model is nested in it.
#           SLOW: 5 random-slope fits per outcome (10 in total).
#   Notes:  - ORs are conditional on the random effects (a "typical"
#             country-wave), not population-averaged.
#           - With weights, logLik is a weighted pseudo-likelihood, so LR and
#             Wald p-values are approximate (model-based, not design-based;
#             Lumley & Scott 2014). Report this in the methods.
#   Saves:  models_cont_interactions.rds (m.fl.int, m.inc.int, spec.int)
#           interaction_results_cont.rds
#           figures/interactions_iwise_cont.png
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
library(tidyverse)   # dplyr, tidyr, purrr, stringr, ggplot2: data handling and plots
library(ordinal)     # clmm(): cumulative link mixed models for ordinal outcomes


# ---- 1. DATA ----------------------------------------------------------------
# Loads every object saved in the prepared file (d.fl, d.inc, v.cluster, rhs_ctx, ...)
# directly into the workspace.
list2env(readRDS("data/iwise_data_prepared.rds"), envir = globalenv())

prep <- function(data, y) {                                  # builds one analysis dataset for outcome y
  data %>%
    filter(cc) %>%                                           # keep complete cases only (cc = TRUE/FALSE flag)
    mutate(.id          = row_number(),                      # row id (not used later in this script)
           country_year = factor(.data[[v.cluster]]),        # cluster variable as factor, for the random part
           "{y}" := factor(.data[[y]], ordered = TRUE))      # outcome as ordered factor (required by clmm)
}

d.fl.cc  <- prep(d.fl,  "FLI_3item")                         # FLI analysis sample
d.inc.cc <- prep(d.inc, "INCOME_5")                          # income-quintile analysis sample
d.list   <- list(FLI = d.fl.cc, INC = d.inc.cc)              # both samples, for checks below

# Random part: random intercept + random IWISE slope per country-wave
# (correlated). Modifiers are individual-level, so they get no random slope.
re <- "(1 + iwisescore | country_year)"                      # random-effects term as text, pasted into formulas


# ---- 2. HELPERS -------------------------------------------------------------
make_f <- function(y, rhs, re) as.formula(paste(y, "~", rhs, "+", re))   # text -> formula "y ~ fixed + random"

fit_clmm <- function(data, y, rhs, re) {                     # fits one weighted ordinal mixed model
  clmm(make_f(y, rhs, re), data = data, weights = w_fit,     # w_fit = weight column in the data
       link = "logit", Hess = TRUE, nAGQ = 1)                # logit = proportional odds; Hess = needed for SEs;
}                                                            # nAGQ = 1 = Laplace (only option with random slopes)

tidy_or <- function(m, label) {                              # fixed effects of one model as a tidy OR table
  co <- summary(m)$coefficients[names(m$beta), , drop = FALSE]   # coefficient rows for fixed effects (no thresholds)
  tibble(model = label, term = rownames(co),                 # model label and term name
         b = co[, 1], se = co[, 2], OR = exp(b),             # log-OR, its SE, OR (per 1 point for iwisescore terms)
         lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = co[, 4])  # Wald 95% CI and p-value
}

load_models <- function(file) {                              # loads saved models into the workspace
  obj <- readRDS(file)                                       # named list of saved objects
  list2env(obj, envir = globalenv())                         # each element becomes an object (m.fl.int, m.inc.int, spec.int)
  message("Loaded from ", file, ": ", paste(names(obj), collapse = ", "))   # report what was loaded
  invisible(names(obj))                                      # return the names without printing
}

fit_quality <- function(m) {                                 # convergence diagnostics of one model
  c(condHess = suppressWarnings(summary(m)$condHess),        # condition number of Hessian (rule of thumb: < 1e4)
    max_grad = max(abs(m$gradient)),                         # largest absolute gradient (should be ~0, e.g. < 1e-3)
    n_na_se  = sum(is.na(summary(m)$coefficients[names(m$beta), 2])))   # fixed effects without SE (should be 0)
}

# Interaction coefficients of iwisescore with one modifier term. R may name
# them "iwisescore:female" or "female:iwisescore", so both parts are
# identified by name, not by position.
int_terms <- function(m, modterm) {
  nm <- names(m$beta)[str_detect(names(m$beta), fixed(":"))]   # all interaction coefficient names
  parts <- str_split(nm, fixed(":"))                         # split each name into its two parts
  tibble(name = nm,                                          # full coefficient name
         iw   = map_lgl(parts, ~ any(.x == "iwisescore")),   # TRUE if one part is iwisescore
         mod  = map_chr(parts, ~ .x[.x != "iwisescore"][1])) %>%   # the other part = modifier level
    filter(iw, !is.na(mod), startsWith(mod, modterm)) %>%    # keep IWISE x this modifier only
    select(-iw)                                              # return name + modifier level
}

# Readable group labels, e.g. "employment = Employed (ref)", "female = 1".
# `mods` = modifier parts of the coefficient names (from int_terms).
group_labels <- function(data, var, modterm, mods) {
  x <- data[[var]]                                           # the modifier variable
  categorical <- !is.numeric(x) || str_detect(modterm, fixed("factor("))   # factor/character/logical, or factor() in formula
  if (categorical) {
    ref  <- levels(factor(x))[1]                             # reference level = first factor level
    lev  <- str_remove(mods, fixed(modterm))                 # level name from coefficient, e.g. "factor(employment)2" -> "2"
  } else {
    ref  <- as.character(min(x, na.rm = TRUE))               # numeric 0/1: reference = lowest value
    lev  <- as.character(max(x, na.rm = TRUE))               # the single other value
  }
  c(paste0(var, " = ", ref, " (ref)"), paste0(var, " = ", lev))   # reference first, then the other groups
}

# IWISE OR (per 1 point) WITHIN each group of the modifier:
#   reference group : exp(b_iwise)
#   other group g   : exp(b_iwise + b_iwise:g), SE from vcov.
# In the all-together model, the other modifiers are at their reference level.
iwise_by_group <- function(m, mod_name, data, label) {
  modterm <- t.mod[[mod_name]]                               # modifier term as written in M1
  b  <- m$beta                                               # fixed-effect estimates
  V  <- vcov(m)[names(b), names(b)]                          # their covariance matrix
  it <- int_terms(m, modterm)                                # interaction terms of this modifier
  labs_g <- group_labels(data, v.mod[[mod_name]], modterm, it$mod)   # readable labels (reference first)
  map2_dfr(c("reference", it$mod), labs_g, function(g, lab) {   # loop: reference group, then each other level
    L <- setNames(numeric(length(b)), names(b))              # contrast vector, all zeros
    L["iwisescore"] <- 1                                     # IWISE main effect = slope in the reference group
    if (g != "reference") L[it$name[it$mod == g]] <- 1       # + interaction term for group g
    est <- sum(L * b)                                        # group-specific log-OR per IWISE point
    se  <- sqrt(drop(t(L) %*% V %*% L))                      # its SE: sqrt(L' V L)
    tibble(model = label, modifier = mod_name, group = lab,  # labels
           OR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))   # OR with 95% CI
  })
}

lr_test <- function(m0, m1) {                              # nested LR test by hand
  lr <- 2 * as.numeric(logLik(m1) - logLik(m0))            # (anova() fails for fit_clmm models)
  df <- length(m1$beta) - length(m0$beta)                  # extra fixed effects only
  c(LR = lr, df = df, p = pchisq(lr, df, lower.tail = FALSE))   # chi-square p-value (approximate with weights)
}

wald_block <- function(m, modterm) {                       # joint Wald test of one interaction
  bn <- int_terms(m, modterm)$name                         # block, inside a larger model
  b  <- m$beta[bn]                                         # estimates of the block
  V  <- vcov(m)[bn, bn, drop = FALSE]                      # their covariance matrix
  W  <- drop(t(b) %*% solve(V) %*% b)                      # Wald statistic b' V^-1 b
  c(W = W, df = length(bn), p = pchisq(W, length(bn), lower.tail = FALSE))   # chi-square test, df = no. of terms
}


# ---- 3. FORMULAS AND MODIFIERS ----------------------------------------------
# rhs_ctx = adjusted M1 right-hand side (from the prepared data), must contain
# the continuous iwisescore.
stopifnot(exists("rhs_ctx"), str_detect(rhs_ctx, "\\biwisescore\\b"))   # stop if M1 missing or has no iwisescore

# Drop wgi_ge_sc_c (government effectiveness) from M1 in this script.
v.drop  <- "wgi_ge_sc_c"                                   # covariate(s) to exclude
rhs_ctx <- str_trim(str_split(rhs_ctx, fixed("+"))[[1]]) %>%   # split M1 into single terms
  discard(~ str_detect(.x, paste0("\\b", v.drop, "\\b", collapse = "|"))) %>%   # remove every term containing it
  paste(collapse = " + ")                                  # rebuild the right-hand side
stopifnot(!str_detect(rhs_ctx, paste0("\\b", v.drop, "\\b", collapse = "|")))   # check: it is gone
rhs_ctx                                                    # print M1 to check

# Each modifier is used exactly as it appears in rhs_ctx (e.g. "female" or
# "factor(employment)"), so main effect and interaction match.
rhs_terms <- str_trim(str_split(rhs_ctx, fixed("+"))[[1]])   # split M1 into single terms
term_of <- function(var) {                                 # finds the M1 term that contains `var`
  hit <- rhs_terms[str_detect(rhs_terms, paste0("\\b", var, "\\b"))]   # terms with var as a whole word
  if (length(hit) != 1) stop(var, " must appear in exactly one term of rhs_ctx (found ", length(hit), ")")   # must be unique
  hit                                                      # return the term exactly as written
}
v.mod <- c(gender = "female", employment = "employment", urbanicity = "urban_imp")   # modifier variables (names = model labels)
t.mod <- map_chr(v.mod, term_of)                           # matching terms in M1
t.mod                                                      # check the terms

# Employment must be categorical, otherwise its interaction is linear (1 df).
# Checked in both samples; factor, character or factor() in the formula are OK.
sapply(d.list, function(d) sapply(v.mod, function(v) class(d[[v]])[1]))   # class of each modifier, per sample
emp_ok <- map_lgl(d.list, ~ !is.numeric(.x$employment) ||
                            str_detect(t.mod[["employment"]], fixed("factor(")))   # TRUE if employment is categorical
if (!all(emp_ok))
  stop("employment is numeric in ", paste(names(emp_ok)[!emp_ok], collapse = ", "),
       " and not wrapped in factor() in rhs_ctx: make it a factor first.")   # stop: a linear interaction would be wrong

# Formulas: base + one per modifier + all together
int_of  <- function(mods) paste0("iwisescore:", t.mod[mods], collapse = " + ")   # "iwisescore:term1 + iwisescore:term2 ..."
rhs_int <- c(base = rhs_ctx,                                                     # M1, no interaction
             map_chr(set_names(names(t.mod)), ~ paste(rhs_ctx, "+", int_of(.x))),   # M1 + one interaction each
             all  = paste(rhs_ctx, "+", int_of(names(t.mod))))                   # M1 + all three interactions
rhs_int                                                    # check: 5 formulas


# ---- 4. FIT OR LOAD ---------------------------------------------------------
# 5 models per outcome. Saved after every fit, so an interrupted run resumes.
# The file also stores the specification the models were fitted with
# (spec.int). If formulas, random part, sample or weights changed, the script
# stops instead of reusing old models. To force a refit: delete the .rds file.
f.int <- "models_cont_interactions.rds"                    # cache file for fitted models
spec.now <- list(rhs_int = rhs_int, re = re,                                  # current formulas and random part
                 n = map_int(d.list, nrow),                                   # current sample sizes
                 w = map_dbl(d.list, ~ sum(.x$w_fit)))                        # current sum of weights
m.fl.int <- list(); m.inc.int <- list(); spec.int <- NULL  # start empty (ignores old objects in the workspace)

if (file.exists(f.int)) {                                  # cache exists:
  load_models(f.int)                                       #   load models (and spec.int, if saved)
  if (is.null(spec.int)) {                                 #   older cache without a stored specification
    warning("Cache has no stored specification, so it cannot be verified. ",
            "Delete ", f.int, " if formulas or data changed since fitting.")
  } else if (!isTRUE(all.equal(spec.int, spec.now))) {     #   specification changed since fitting
    stop("Saved models do not match the current formulas/data. Delete ", f.int, " to refit.")
  }
}

for (nm in names(rhs_int)) {                               # loop over the 5 specifications
  if (is.null(m.fl.int[[nm]])) {                           # fit FLI model only if missing
    message("FLI: ", nm);    m.fl.int[[nm]]  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_int[[nm]], re)
  }
  if (is.null(m.inc.int[[nm]])) {                          # fit income model only if missing
    message("Income: ", nm); m.inc.int[[nm]] <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_int[[nm]], re)
  }
  saveRDS(list(m.fl.int = m.fl.int, m.inc.int = m.inc.int, spec.int = spec.now), f.int)   # save models + spec (resumable)
}

# Fit quality: all columns should look fine
round(sapply(m.fl.int,  fit_quality), 4)                   # diagnostics, FLI models (one column per model)
round(sapply(m.inc.int, fit_quality), 4)                   # diagnostics, income models

# Random intercept and slope variances (and their correlation) per outcome
map(list(FLI = m.fl.int$base, INC = m.inc.int$base), ~ VarCorr(.x))   # random-effect SDs + correlation, base models


# ---- 5. STEP A: ONE AT A TIME (PRIMARY) -------------------------------------
# LR test of each interaction block vs the base model; Holm over the 3 tests.
# With N ~ 80,000 small p-values are expected: judge the size of the ORs
# within groups (Section 7), not the p-value alone.
lr_tab <- function(ints, outcome) {                        # LR table for one outcome
  imap_dfr(ints[names(ints) != "base"], function(m, nm) {  # each interaction model vs base
    r <- lr_test(ints$base, m)                             # LR test
    tibble(outcome = outcome, model = nm, LR = r[["LR"]], df = r[["df"]], p = r[["p"]])   # one row per model
  }) %>%
    mutate(p_holm = p.adjust(replace(p, model == "all", NA), "holm"))  # "all" = global test (Step B); Holm per outcome
}
lr.int <- bind_rows(lr_tab(m.fl.int, "FLI"), lr_tab(m.inc.int, "INC"))   # both outcomes in one table
lr.int                                                     # print


# ---- 6. STEP B: ALL TOGETHER (SENSITIVITY) ----------------------------------
# 6a. Global chunk test = the "all" row of lr.int (all interaction terms at once).
lr.int %>% filter(model == "all")                          # global test, both outcomes

# 6b. Each block within the all-together model (adjusted for the others)
round(sapply(t.mod, function(tm) wald_block(m.fl.int$all,  tm)), 4)   # FLI
round(sapply(t.mod, function(tm) wald_block(m.inc.int$all, tm)), 4)   # income

# 6c. Do the ratio of ORs change between one-at-a-time and all-together?
#     Large changes = the modifications overlap (e.g. gender vs employment).
int_or <- function(ints, outcome) {                        # ratio-of-OR table for one outcome
  imap_dfr(ints[names(ints) != "base"], ~ tidy_or(.x, .y)) %>%   # fixed effects of all interaction models
    filter(str_detect(term, fixed(":")), str_detect(term, "iwisescore")) %>%   # keep IWISE interactions only
    mutate(outcome = outcome,                              # outcome label
           approach = if_else(model == "all", "all together", "one at a time"))   # which approach
}
ror <- bind_rows(int_or(m.fl.int, "FLI"), int_or(m.inc.int, "INC"))   # both outcomes
ror %>%                                                     # ratio of ORs, side by side
  select(outcome, term, approach, OR) %>%                  # columns needed
  pivot_wider(names_from = approach, values_from = OR) %>% # one column per approach
  mutate(change_pct = 100 * (`all together` / `one at a time` - 1)) %>%   # % change when adjusting for other interactions
  print(n = Inf)                                           # print all rows


# ---- 7. IWISE ORs WITHIN MODIFIER GROUPS ------------------------------------
# Main output for the paper (Knol & VanderWeele 2012): the IWISE OR within
# each group, plus the ratio of ORs (ror above).
mods <- names(t.mod)                                       # "gender", "employment", "urbanicity"
strat <- bind_rows(
  map_dfr(mods, ~ iwise_by_group(m.fl.int[[.x]],  .x, d.fl.cc,  "FLI, one at a time")),   # FLI, each modifier from its own model
  map_dfr(mods, ~ iwise_by_group(m.inc.int[[.x]], .x, d.inc.cc, "INC, one at a time")),   # income, same
  map_dfr(mods, ~ iwise_by_group(m.fl.int$all,    .x, d.fl.cc,  "FLI, all together")),    # FLI, all from the joint model
  map_dfr(mods, ~ iwise_by_group(m.inc.int$all,   .x, d.inc.cc, "INC, all together")))    # income, same
strat %>% print(n = Inf)                                   # print all rows

p.strat <- strat %>%
  filter(str_detect(model, "one at a time")) %>%           # plot primary (one-at-a-time) results only
  ggplot(aes(x = OR, xmin = lo, xmax = hi, y = group)) +   # OR and CI per group
  geom_vline(xintercept = 1, colour = "grey60", linetype = 2) +   # OR = 1: no association
  geom_pointrange(size = 0.3) +                            # point estimate + CI
  facet_grid(modifier ~ model, scales = "free_y", space = "free_y") +   # rows = modifier, columns = outcome
  scale_x_log10() +                                        # log scale: ORs symmetric around 1
  labs(x = "Adjusted OR per 1-point IWISE increase (95% CI, log scale)",
       y = "Modifier group") +                             # axis labels
  theme_minimal()                                          # clean theme
p.strat                                                    # show plot

dir.create("figures", showWarnings = FALSE)                # create figures folder if missing
ggsave("figures/interactions_iwise_cont.png", p.strat, width = 9, height = 7, dpi = 300, bg = "white")   # save figure
saveRDS(list(lr.int = lr.int, ror = ror, strat = strat), "interaction_results_cont.rds")   # save result tables