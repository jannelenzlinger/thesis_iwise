# =============================================================================
# PRE-MODEL CHECKS  |  Ordinal logistic regression
# Dataset: iwise_analysis, used here as d.iwise
# Outcomes: INDEX_FL_ord (9 levels), INCOME_5 (5 levels)
# Main predictor: iwisescore (0-36)
#
# THIS SCRIPT STOPS BEFORE THE MODEL IS FITTED.
# Everything here is a property of the DATA, not of a fitted model.
# See premodel_checks_explained.md for the reasoning behind each check.
# =============================================================================
 
 
# ---- 0. PACKAGES ------------------------------------------------------------
 
pkgs <- c("tidyverse",   # data wrangling and plotting
          "naniar",      # missing data summaries and plots
          "car",         # vif(): multicollinearity
          "lme4",        # lmer(): empty model for the ICC
          "performance", # icc()
          "detectseparation")  # separation check
 
missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing_pkgs)) install.packages(missing_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))
# invisible() stops R printing the list that lapply returns - we only want the
# side effect (packages attached), not the return value.
 
select <- dplyr::select   # MASS and others mask dplyr::select; be explicit
 
 
# ---- 1. LOAD DATA AND DEFINE VARIABLE GROUPS --------------------------------
# Work on a copy so the raw object stays clean and you can always restart.
 
d.iwise <- iwise_analysis
 
# Naming these groups once means you edit your variable list in ONE place
# instead of rewriting every command below.
 
v.outcomes    <- c("INDEX_FL_ord", "INCOME_5")
v.main        <- "iwisescore"      # main predictor, 0-36
v.ctrl        <- "iwise4"          # alternative 4-item version, 0-12
 
# The WGI come in two suffix families, one row per country-year:
#   _est = governance estimate, standard normal units, approx. -2.5 to 2.5
#   _sc  = the same estimate mapped onto a 0-100 scale
# They are transformations of each other, so NEVER use both for the same
# dimension in one model - that is perfect collinearity.
v.wgi.est <- c("wgi_va_est", "wgi_pv_est", "wgi_ge_est",
               "wgi_rq_est", "wgi_rl_est", "wgi_cc_est")
v.wgi.sc  <- c("wgi_va_sc",  "wgi_pv_sc",  "wgi_ge_sc",
               "wgi_rq_sc",  "wgi_rl_sc",  "wgi_cc_sc")
 
v.wgi <- v.wgi.est   # pick ONE family; _est reads as "per 1 SD of governance"
 
v.continuous  <- c("iwisescore", "hhsize", "gdp_pc_ppp", v.wgi)
v.categorical <- c("female", "urban_imp", "age_gp_profile",
                   "maritalstatus", "employment", "education")
v.indices     <- c("INDEX_PH", "INDEX_FS", "INDEX_CA", "INDEX_YD")
v.cluster     <- "country"
 
# IMPORTANT: build the variable list PER MODEL, not one global list.
# A respondent missing INCOME_5 should still count towards the INDEX_FL_ord
# model. Using one combined list is what drops them from both.
v.model.fl  <- c("INDEX_FL_ord", v.main, v.continuous, v.categorical, v.cluster)
v.model.inc <- c("INCOME_5",     v.main, v.continuous, v.categorical, v.cluster)
 
 
# =============================================================================
# CHECK 1: VARIABLE TYPES AND LEVEL ORDERING
# -----------------------------------------------------------------------------
# WHY: the model must know the outcome is ORDERED. A plain factor either errors
# or discards the ordering. Levels sorted alphabetically ("0","100","16.67")
# produce a model that runs fine and is completely wrong, with no warning.
# =============================================================================
 
# 1a. What does R currently think these variables are?
class(d.iwise$INDEX_FL_ord)      # want: "ordered" "factor"
class(d.iwise$INCOME_5)
levels(d.iwise$INDEX_FL_ord)     # want: ascending numeric order
levels(d.iwise$INCOME_5)
 
is.ordered(d.iwise$INDEX_FL_ord) # the check that actually matters
is.ordered(d.iwise$INCOME_5)
 
# NOTE: Positron's Variables pane may show "fct" even for ordered factors.
# is.ordered() is authoritative; the pane is a display layer.
 
# 1b. Recode. Safe to run even if already correct - it is idempotent.
#     as.character() -> the labels; as.numeric() -> real numbers;
#     factor() on a NUMERIC vector sorts levels numerically, which removes the
#     alphabetical hazard entirely.
d.iwise <- d.iwise %>%
  mutate(
    INDEX_FL_ord = factor(as.numeric(as.character(INDEX_FL_ord)), ordered = TRUE),
    INCOME_5     = factor(as.numeric(as.character(INCOME_5)),     ordered = TRUE),
 
    # Categorical predictors: labelled factors. The FIRST level becomes the
    # reference category the coefficients are compared against - choose it
    # deliberately rather than accepting R's alphabetical default.
    female         = factor(female, levels = c(0, 1),
                            labels = c("Male", "Female")),
    urban_imp      = factor(urban_imp, levels = c(0, 1),
                            labels = c("Rural", "Peri-urban/urban")),
    age_gp_profile = factor(age_gp_profile, levels = 0:3,
                            labels = c("15-24", "25-34", "35-49", "50+")),
    maritalstatus  = factor(maritalstatus, levels = 0:2,
                            labels = c("Single", "Married/partnered",
                                       "Sep/div/widowed")),
    employment     = factor(employment, levels = 0:3,
                            labels = c("Employed", "Underemployed",
                                       "Unemployed", "Out of workforce")),
    education      = factor(education, levels = 0:2,
                            labels = c("Elementary", "Secondary", "College")),
 
    across(all_of(v.continuous), as.numeric)
  )
 
# 1c. Confirm the recode did what you expect.
str(d.iwise %>% select(any_of(c(v.outcomes, v.continuous, v.categorical))))
 
# 1d. Direction check. R orders ascending, so level 1 = 0 (worst financial
#     life) and level 9 = 100 (best). Coefficients will then read as
#     "higher water insecurity -> odds of a HIGHER financial-life category".
#     If you want the reverse, use forcats::fct_rev() rather than flipping
#     signs by hand when interpreting.
levels(d.iwise$INDEX_FL_ord)
 
# 1e. Sanity-check the WGI scales before committing to a suffix family.
#     _est should sit roughly in [-2.5, 2.5]; _sc in [0, 100]. If reversed,
#     swap v.wgi.est and v.wgi.sc above.
d.iwise %>%
  select(any_of(c(v.wgi.est, v.wgi.sc))) %>%
  summarise(across(everything(),
                   list(min = ~min(.x, na.rm = TRUE),
                        max = ~max(.x, na.rm = TRUE)))) %>%
  pivot_longer(everything(), names_to = c("var", "stat"),
               names_pattern = "^(.*)_(min|max)$") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  print(n = Inf)
 
 
# =============================================================================
# CHECK 2: MISSING DATA
# -----------------------------------------------------------------------------
# WHY: the model silently deletes any row missing ANY variable in the formula.
# That is only harmless if the deleted rows are a random subset. If they differ
# systematically from the rows you keep, every coefficient is biased.
# =============================================================================
 
# 2a. How much is missing on each variable?
#     n_miss = count, pct_miss = percentage.
d.iwise %>%
  select(any_of(unique(c(v.model.fl, v.model.inc)))) %>%
  naniar::miss_var_summary() %>%
  print(n = 30)
 
# 2b. Which variables go missing TOGETHER?
#     The UpSet plot shows combinations. All bars are COUNTS, not percentages:
#     top bars = rows in that specific combination, left bars = total rows
#     missing on each variable, dots below = which variables are involved.
naniar::gg_miss_upset(d.iwise %>% select(any_of(unique(c(v.model.fl,
                                                         v.model.inc)))))
 
# 2c. How many rows does each model actually lose?
#     complete.cases() returns TRUE only if a row has NO missing value on any
#     listed variable. Doing this per model is the point: a row missing
#     INCOME_5 should not be dropped from the financial-life model.
d.iwise <- d.iwise %>%
  mutate(cc.fl  = complete.cases(select(., any_of(v.model.fl))),
         cc.inc = complete.cases(select(., any_of(v.model.inc))))
 
cat("\nRows retained - INDEX_FL_ord model:", sum(d.iwise$cc.fl),
    "of", nrow(d.iwise),
    sprintf("(%.1f%% lost)\n", 100 * mean(!d.iwise$cc.fl)))
cat("Rows retained - INCOME_5 model:    ", sum(d.iwise$cc.inc),
    "of", nrow(d.iwise),
    sprintf("(%.1f%% lost)\n", 100 * mean(!d.iwise$cc.inc)))
 
# 2d. Do the dropped rows DIFFER from the kept rows?
#     Under MCAR the two groups look alike. A large gap means listwise
#     deletion is not neutral and you should say so in the paper.
d.iwise %>%
  group_by(cc.fl) %>%
  summarise(n = n(),
            across(all_of(c("iwisescore", "hhsize", "gdp_pc_ppp")),
                   ~mean(.x, na.rm = TRUE)),
            .groups = "drop")
 
# 2e. Is missingness clustered by country?
#     This is the check that matters most for you: a country losing ~100% of
#     its respondents disappears from the analysis entirely, cutting your
#     cluster count and narrowing what your results generalise to.
d.iwise %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(n = n(),
            pct_dropped_fl  = round(100 * mean(!cc.fl),  1),
            pct_dropped_inc = round(100 * mean(!cc.inc), 1),
            .groups = "drop") %>%
  arrange(desc(pct_dropped_fl)) %>%
  print(n = Inf)
 
# 2f. For any country losing most of its cases, find WHICH variable is
#     responsible. The answer determines what to do:
#       missing OUTCOME            -> accept the loss, nothing to model
#       missing COUNTRY covariate  -> fix the merge, or drop that covariate
#       missing INDIVIDUAL covar.  -> consider multiple imputation (mice)
bad_countries <- d.iwise %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(pct = mean(!cc.fl), .groups = "drop") %>%
  filter(pct > 0.5) %>%
  pull(!!sym(v.cluster))
 
d.iwise %>%
  filter(.data[[v.cluster]] %in% bad_countries) %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(across(any_of(v.model.fl), ~round(100 * mean(is.na(.x)), 1)),
            .groups = "drop") %>%
  print(width = Inf)
 
 
# =============================================================================
# CHECK 3: OUTCOME DISTRIBUTION AND SPARSE CATEGORIES
# -----------------------------------------------------------------------------
# WHY: the model estimates one cut-point per threshold (8 for a 9-level
# outcome). Each is estimated from respondents near that threshold. A level
# holding very few people gives an unstable cut-point, and drains power from
# the proportional odds test later.
# =============================================================================
 
# 3a. Counts and proportions.
table(d.iwise$INDEX_FL_ord, useNA = "ifany")
round(100 * prop.table(table(d.iwise$INDEX_FL_ord)), 1)
 
table(d.iwise$INCOME_5, useNA = "ifany")
round(100 * prop.table(table(d.iwise$INCOME_5)), 1)
 
# LOOKING FOR: no level below about 2% of the sample.
 
# 3b. Visual version.
d.iwise %>%
  filter(!is.na(INDEX_FL_ord)) %>%
  count(INDEX_FL_ord) %>%
  ggplot(aes(INDEX_FL_ord, n)) +
  geom_col() +
  labs(title = "INDEX_FL_ord: respondents per level",
       x = "Financial Life Index", y = "n") +
  theme_minimal()
 
# 3c. IF sparse: collapse to a coarser but better-populated scale. Collapse on
#     substantive grounds, not on whatever gives the nicest p-value, and
#     report that you did it.
d.iwise <- d.iwise %>%
  mutate(INDEX_FL_c = cut(as.numeric(as.character(INDEX_FL_ord)),
                          breaks = c(-Inf, 25, 50, 75, Inf),
                          labels = c("Low", "Lower-mid", "Upper-mid", "High"),
                          ordered_result = TRUE))
table(d.iwise$INDEX_FL_c)
 
 
# =============================================================================
# CHECK 4: PREDICTOR DISTRIBUTIONS AND RARE CATEGORIES
# -----------------------------------------------------------------------------
# WHY: a category with a handful of respondents produces a coefficient too
# imprecise to interpret. Implausible values (a household of 97) are usually
# data errors. Heavy skew in gdp_pc_ppp distorts estimates.
# =============================================================================
 
# 4a. Categorical predictors: counts per level.
for (v in v.categorical) {
  cat("\n---", v, "---\n")
  print(table(d.iwise[[v]], useNA = "ifany"))
}
# LOOKING FOR: levels with very few cases -> collapse them.
 
# 4b. Continuous predictors: range, centre, spread, skew.
d.iwise %>%
  select(all_of(v.continuous)) %>%
  summarise(across(everything(),
                   list(min  = ~min(.x, na.rm = TRUE),
                        q25  = ~quantile(.x, .25, na.rm = TRUE),
                        med  = ~median(.x, na.rm = TRUE),
                        mean = ~mean(.x, na.rm = TRUE),
                        q75  = ~quantile(.x, .75, na.rm = TRUE),
                        max  = ~max(.x, na.rm = TRUE)))) %>%
  pivot_longer(everything(), names_to = c("var", "stat"),
               names_pattern = "^(.*)_(min|q25|med|mean|q75|max)$") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  print(n = Inf)
# LOOKING FOR: impossible values; mean far from median (= skew).
 
# 4c. Histograms.
d.iwise %>%
  select(all_of(v.continuous)) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 40) +
  facet_wrap(~name, scales = "free") +
  labs(title = "Distribution of continuous predictors") +
  theme_minimal()
 
# 4d. gdp_pc_ppp is income-like and usually right-skewed. Log transformation
#     is standard, and makes the coefficient read as "per 1% change in GDP".
d.iwise <- d.iwise %>% mutate(log_gdp = log(gdp_pc_ppp))
hist(d.iwise$log_gdp, breaks = 40, main = "log(GDP per capita)")
 
# 4e. iwisescore is often bunched at zero in water-insecurity data. Check the
#     shape - it affects whether you model it as continuous or categorical.
table(cut(d.iwise$iwisescore, breaks = c(-Inf, 0, 2, 11, 23, Inf),
          labels = c("0", "1-2", "3-11", "12-23", "24-36")))
# The standard IWISE cut-points are 0-2 no-to-marginal, 3-11 low,
# 12-23 moderate, 24-36 high (Young et al. 2019).
 
 
# =============================================================================
# CHECK 5: OUTCOME x CATEGORICAL PREDICTOR CROSS-TABS
# -----------------------------------------------------------------------------
# WHY: empty cells (zero respondents in a predictor-level x outcome-level
# combination) are the main cause of non-convergence and of separation.
# =============================================================================
 
for (v in v.categorical) {
  cat("\n=== INDEX_FL_ord x", v, "===\n")
  tb <- table(d.iwise[[v]], d.iwise$INDEX_FL_ord)
  print(tb)
  if (any(tb == 0)) cat(">>> WARNING: empty cell(s) present\n")
  if (any(tb < 5 & tb > 0)) cat(">>> NOTE: cell(s) with fewer than 5 cases\n")
}
 
for (v in v.categorical) {
  cat("\n=== INCOME_5 x", v, "===\n")
  tb <- table(d.iwise[[v]], d.iwise$INCOME_5)
  print(tb)
  if (any(tb == 0)) cat(">>> WARNING: empty cell(s) present\n")
}
# IF EMPTY CELLS: collapse the offending predictor category.
 
 
# =============================================================================
# CHECK 6: SAMPLE SIZE / EVENTS PER VARIABLE
# -----------------------------------------------------------------------------
# WHY: you need enough cases in the SMALLEST outcome category to estimate all
# your coefficients. Too few -> biased estimates, unreliable intervals.
# Rule of thumb: at least 10 events per estimated parameter.
# NOTE: a k-level factor costs k-1 parameters, not 1.
# =============================================================================
 
# Define the right-hand side ONCE here; reuse it for the model later.
rhs <- paste("iwisescore + female + urban_imp + age_gp_profile +",
             "maritalstatus + hhsize + employment + education")
rhs_ctx <- paste(rhs, "+ log_gdp + wgi_cc_est")   # ONE governance term only
 
# model.matrix() expands factors into dummies, so this counts real parameters.
n.params <- ncol(model.matrix(as.formula(paste("~", rhs_ctx)), data = d.iwise))
 
epv <- function(outcome, n.params) {
  tb <- table(d.iwise[[outcome]])
  cat(sprintf("\n%s: %d parameters | smallest category n = %d | EPV = %.1f %s\n",
              outcome, n.params, min(tb), min(tb) / n.params,
              ifelse(min(tb) / n.params < 10, "<- BELOW 10", "")))
}
epv("INDEX_FL_ord", n.params)
epv("INCOME_5",     n.params)
 
# IF BELOW 10: collapse outcome levels, drop covariates on theoretical
# grounds, or plan to use penalised estimation (brglm2).
 
 
# =============================================================================
# CHECK 7: MULTICOLLINEARITY
# -----------------------------------------------------------------------------
# WHY: if two predictors carry nearly the same information, the model cannot
# separate their effects. Coefficients become unstable - wide standard errors,
# signs that flip when you add or remove a variable.
#
# This is a property of the PREDICTORS ONLY. The outcome is irrelevant to it,
# which is why it can be checked before fitting the ordinal model.
# =============================================================================
 
# 7a. Correlation matrix. Catches PAIRWISE overlap only.
d.iwise %>%
  select(all_of(v.continuous), any_of(v.indices)) %>%
  cor(use = "pairwise.complete.obs") %>%
  round(2) %>%
  print()
# LOOKING FOR: |r| above ~0.8. Expect this among the six WGI indicators, and
# possibly between INDEX_YD and education (flagged in your codebook).
 
# 7b. Are the two WGI suffix families redundant with each other? They should
#     correlate at or near 1 - which is the evidence that you can use only one.
for (i in seq_along(v.wgi.est)) {
  if (all(c(v.wgi.est[i], v.wgi.sc[i]) %in% names(d.iwise))) {
    r <- cor(d.iwise[[v.wgi.est[i]]], d.iwise[[v.wgi.sc[i]]],
             use = "pairwise.complete.obs")
    cat(sprintf("%-12s vs %-12s r = %.3f\n", v.wgi.est[i], v.wgi.sc[i], r))
  }
}
 
# 7c. VIF. Catches overlap involving THREE OR MORE variables, which pairwise
#     correlation misses. VIF = 4 means this coefficient's standard error is
#     twice as wide as it would be with uncorrelated predictors.
#
#     VIF concerns only the predictor matrix, so a probe model with a random
#     outcome gives the identical answer - no need to fit the ordinal model.
set.seed(2024)
m.vifprobe <- lm(as.formula(paste("rnorm(nrow(d.iwise)) ~", rhs_ctx)),
                 data = d.iwise)
car::vif(m.vifprobe)
# For FACTORS this returns GVIF. Read the GVIF^(1/(2*Df)) column and SQUARE it
# before comparing to the usual thresholds of 5 or 10.
 
# 7d. All six WGI indicators together, to see how bad it is.
m.vifwgi <- lm(as.formula(paste("rnorm(nrow(d.iwise)) ~",
                                paste(v.wgi, collapse = " + "))),
               data = d.iwise)
car::vif(m.vifwgi)
 
# 7e. IF COLLINEAR: either pick ONE indicator on theoretical grounds, or
#     reduce all six to a single component. Run the PCA at COUNTRY level, not
#     respondent level, or the largest country dominates the components.
d.wgi.country <- d.iwise %>%
  select(all_of(c(v.cluster, v.wgi))) %>%
  distinct() %>%
  drop_na()
 
pca.wgi <- prcomp(d.wgi.country %>% select(all_of(v.wgi)), scale. = TRUE)
summary(pca.wgi)        # PC1 usually explains 80-90% of WGI variance
pca.wgi$rotation[, 1]   # loadings: all six should point the same direction
 
# Merge PC1 back on, then use wgi_pc1 in place of the six indicators.
d.iwise <- d.wgi.country %>%
  mutate(wgi_pc1 = pca.wgi$x[, 1]) %>%
  select(all_of(v.cluster), wgi_pc1) %>%
  right_join(d.iwise, by = v.cluster)
 
 
# =============================================================================
# CHECK 8: CLUSTERING / NON-INDEPENDENCE
# -----------------------------------------------------------------------------
# WHY: all regression assumes independent observations. Yours are nested in 76
# countries, and gdp_pc_ppp / wgi_* take ONE value per country. Ignoring this
# makes standard errors far too small, so you find effects that are not there.
#
# The Moulton problem: ~74,000 respondents, but only 76 independent
# observations of the country-level variables.
#
# The ICC needs an EMPTY model (no predictors), so it belongs here rather than
# after - it describes the data, not your model.
# =============================================================================
 
# 8a. How many clusters, and how big?
d.iwise %>%
  count(across(all_of(v.cluster))) %>%
  summarise(n_countries = n(),
            min_n = min(n), median_n = median(n), max_n = max(n))
# Fewer than ~40 clusters means you need a small-sample correction (CR2) later.
 
# 8b. ICC: what share of outcome variance sits BETWEEN countries?
m.iccprobe <- lme4::lmer(
  as.numeric(as.character(INDEX_FL_ord)) ~ 1 + (1 | country),
  data = d.iwise
)
performance::icc(m.iccprobe)
# LOOKING FOR: ICC above ~0.05 means clustering is NOT ignorable.
# In cross-national survey data it is frequently much higher.
 
# 8c. Visual: how much do country means differ?
d.iwise %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(mean_fl = mean(as.numeric(as.character(INDEX_FL_ord)),
                           na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  ggplot(aes(x = reorder(factor(country), mean_fl), y = mean_fl)) +
  geom_point() +
  coord_flip() +
  labs(title = "Country mean of INDEX_FL_ord",
       x = "Country", y = "Mean") +
  theme_minimal(base_size = 7)
 
# YOUR FIX (applied at the modelling stage, not here):
#   Option A - cluster-robust standard errors with a CR2 small-sample
#              correction (clubSandwich). Keeps coefficients, fixes inference.
#   Option B - multilevel ordinal model, ordinal::clmm() with (1 | country).
#              Preferred when you have country-level predictors and want to
#              say something about between-country variation.
 
 
# =============================================================================
# CHECK 9: SEPARATION
# -----------------------------------------------------------------------------
# WHY: separation means a predictor (or combination) perfectly predicts the
# outcome. Symptoms: coefficients above 15, standard errors in the thousands,
# convergence warnings. It is a DATA problem - more iterations never fix it.
#
# The ordinal model works by cutting the outcome at each threshold, so you can
# test each of those binary cuts in advance with a fast linear-programming
# check. No ordinal model needed.
# =============================================================================
 
check_separation <- function(data, outcome, rhs) {
  data <- data %>% drop_na(all_of(outcome))
  lv <- levels(droplevels(data[[outcome]]))
  cat("\n--- Separation check:", outcome, "---\n")
  for (k in seq_len(length(lv) - 1)) {
    d.k <- data %>%
      mutate(y_bin = as.integer(as.numeric(.data[[outcome]]) > k))
    fit <- glm(as.formula(paste("y_bin ~", rhs)), data = d.k,
               family = binomial,
               method = detectseparation::detect_separation)
    cat(sprintf("  cut above %-8s -> %s\n", lv[k],
                ifelse(fit$outcome, "SEPARATION DETECTED", "ok")))
  }
}
 
check_separation(d.iwise, "INDEX_FL_ord", rhs)
check_separation(d.iwise, "INCOME_5",     rhs)
 
# IF DETECTED: collapse the offending category, drop the variable, or plan to
# use Firth penalised likelihood (brglm2::brglmFit) at the modelling stage.
 
 
# =============================================================================
# SUMMARY
# =============================================================================
 
tribble(
  ~Check, ~What_it_tests,                              ~Section,
  1, "Variable types and level ordering",              "Check 1",
  2, "Missing data pattern and its consequences",      "Check 2",
  3, "Outcome levels adequately populated",            "Check 3",
  4, "Predictor distributions, rare categories",       "Check 4",
  5, "No empty cells in outcome x predictor tables",   "Check 5",
  6, "At least ~10 events per parameter",              "Check 6",
  7, "No severe multicollinearity",                    "Check 7",
  8, "Clustering quantified (ICC), fix chosen",        "Check 8",
  9, "No separation",                                  "Check 9"
) %>% print(n = Inf)
 
cat("
STILL TO DO, AFTER THE MODEL IS FITTED
--------------------------------------
These are properties of a FITTED MODEL - coefficients, predictions, residuals.
None of them exist yet, which is why they cannot be done here.
 
  10. PROPORTIONAL ODDS - the assumption specific to this model. Tests whether
      each predictor's effect is the same across every cut of the outcome. It
      is frequently violated. Three tests: Brant, per-variable LR tests
      (ordinal::nominal_test), and a plot of coefficients across cuts. With
      your sample size the formal tests reject on trivial departures, so the
      plot is the most honest of the three.
 
  11. LINEARITY OF THE LOGIT - whether continuous predictors relate linearly
      to the log-odds. The log-odds scale only exists once a model defines it.
      Also a PARTIAL relationship, so it depends on the other covariates.
 
  12. INFLUENTIAL OBSERVATIONS - whether a few extreme cases drive the result.
      Influence means 'how much do estimates change if I remove this case', so
      it needs estimates.
 
  13. GOODNESS OF FIT - whether predicted and observed frequencies agree.
 
Fit the model you intend to REPORT, then run these on it. Diagnostics are
conditional on the model: a variable can look fine alone and fail once the
others are in. Do NOT use these tests to select variables - variable selection
follows your research question, and the tests only tell you whether the model
you chose is trustworthy.
")
 
# =============================================================================
# END OF PRE-MODEL CHECKS
# =============================================================================