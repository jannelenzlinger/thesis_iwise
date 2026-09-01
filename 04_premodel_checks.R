# =============================================================================
# PRE-MODEL CHECKS  |  Ordinal logistic regression
# Dataset: iwise_analysis
# Outcomes: INDEX_FL_ord (Financial Life Index), INCOME_5 (Income Quintiles)
# Main predictor: iwisescore (Water Insecurity Experiences, 0-36)
# Written for Positron / R >= 4.2
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Nine checks, in order, on the DATA. It STOPS BEFORE THE MODEL IS FITTED.
# The four checks that require a fitted model are listed at the end.
# See premodel_checks_explained.md for the reasoning behind each one.
#
# TWO THINGS TO KNOW UP FRONT
# ---------------------------
# 1. Neither outcome is binary, so the model is ORDINAL logistic regression
#    (proportional odds / cumulative logit), not binary logistic. Every binary
#    assumption still applies, PLUS proportional odds - checked after fitting.
#
# 2. The Financial Life Index keeps ALL respondents with a computed index -
#    both 3-item and 4-item responders, per Gallup's own methodology. The
#    outcome therefore retains all 9 levels. Levels 25 and 75 are reachable
#    ONLY from 3-item responses, so they exist solely in that group; a
#    sensitivity check comparing the two groups is included in Check 3.
#
# STRUCTURE
# ---------
#   Section 0   Packages
#   Section 1   Load data, define variable groups
#   Check 1     Variable types and level ordering
#   [split]     Derived variables, then two analytic samples: d.fl and d.inc
#   Check 2     Missing data
#   Check 3     Outcome distribution and sparse categories
#   Check 4     Predictor distributions and rare categories
#   Check 5     Outcome x categorical predictor cross-tabs
#   Check 6     Sample size / events per variable
#   Check 7     Multicollinearity
#   Check 8     Clustering / non-independence
#   Check 9     Separation
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

# Load from FILE, not from an object already in memory, so the script runs
# standalone in a fresh session tomorrow. Use whichever matches your file.
d.iwise <- readRDS("data/iwise_analysis.rds")

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

v.wgi <- v.wgi.sc   # pick ONE family; _est reads as "per 1 SD of governance"

v.continuous  <- c("iwisescore", "hhsize", "gdp_pc_ppp", v.wgi)
v.categorical <- c("female", "urban_imp", "age_gp_profile",
                   "maritalstatus", "employment", "education")
# The four other Gallup indices are NOT used as predictors. They are missing
# for 30-73% of respondents (INDEX_YD 70-73%, INDEX_CA 40-47%, INDEX_PH 29-31%;
# INDEX_FS is complete), so including them would cut the analytic sample to
# roughly a quarter. Kept here only for the correlation matrix in Check 7.
v.indices     <- c("INDEX_PH", "INDEX_FS", "INDEX_CA", "INDEX_YD")
v.cluster     <- "iso3c"

# NOTE: the per-model variable lists (v.model.fl, v.model.inc) are defined
# AFTER the split below, once the derived variables exist.


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


# ---- SETUP: DERIVED VARIABLES FIRST -----------------------------------------
# IMPORTANT: create every derived variable on d.iwise BEFORE splitting.
# Anything you mutate afterwards will not propagate into d.fl / d.inc.

v.fli.items <- c("WP2319", "WP30", "WP31", "WP88")

d.iwise <- d.iwise %>%
  mutate(
    # How many FLI items did each respondent answer? Gallup computes the index
    # from 3 or 4 items; nobody in this file answered fewer than 2, and the
    # 2-item group has no index at all. Kept as a flag so you can run the
    # 3-item vs 4-item sensitivity check in Check 3 without refiltering.
    n_fli_answered = rowSums(!is.na(across(all_of(v.fli.items)))),
    fli_4item      = n_fli_answered == 4,

    # gdp_pc_ppp is income-like and right-skewed; log makes the coefficient
    # read as "per 1% change in GDP" and stops extreme values dominating.
    log_gdp = log(gdp_pc_ppp),

    # Standard IWISE bands (Young et al. 2019), useful if the linearity check
    # later suggests a threshold rather than a smooth gradient.
    iwise_cat = cut(iwisescore, breaks = c(-Inf, 2, 11, 23, Inf),
                    labels = c("No-to-marginal", "Low", "Moderate", "High"),
                    ordered_result = TRUE),

    # Numeric versions of the four other indices, for the Check 7 correlation
    # matrix only. They are not model predictors.
    across(all_of(v.indices), ~as.numeric(as.character(.x)),
           .names = "{.col}_num")
  )


# ---- SPLIT INTO TWO ANALYTIC SAMPLES ----------------------------------------
# Each outcome gets its own frame. This is the point: a respondent missing
# INCOME_5 should still count towards the FLI model, and vice versa.

# FLI sample: everyone with a computed index (3-item and 4-item responders).
d.fl <- d.iwise %>%
  drop_na(INDEX_FL_ord) %>%
  mutate(INDEX_FL_ord = droplevels(INDEX_FL_ord))
# droplevels() is a safety net only - no level should be empty here. It matters
# if you ever restrict the sample further, since polr() would otherwise try to
# estimate cut-points for categories containing nobody.

# Income sample.
d.inc <- d.iwise %>%
  drop_na(INCOME_5) %>%
  mutate(INCOME_5 = droplevels(INCOME_5))

cat("\n--- ANALYTIC SAMPLES ---\n")
cat("Full data:    ", nrow(d.iwise), "\n")
cat("FLI sample:   ", nrow(d.fl),
    sprintf("(%.1f%% of full)\n", 100 * nrow(d.fl) / nrow(d.iwise)))
cat("Income sample:", nrow(d.inc),
    sprintf("(%.1f%% of full)\n", 100 * nrow(d.inc) / nrow(d.iwise)))

cat("\nFLI outcome levels after droplevels():\n")
print(levels(d.fl$INDEX_FL_ord))
cat("Income outcome levels:\n")
print(levels(d.inc$INCOME_5))


# ---- VARIABLE LISTS PER MODEL -----------------------------------------------
# The FLI list includes the four items so partial responders are flagged as
# incomplete. (Redundant now that d.fl already filters on them, but it keeps
# the definition explicit and lets the diagnostics below report on them.)

# NOTE: the four FLI items are deliberately NOT in this list. Including them
# would drop 3-item responders, which is the choice we are not making.
v.model.fl  <- c("INDEX_FL_ord", v.main, "hhsize", "log_gdp", v.wgi,
                 v.categorical, v.cluster)
v.model.inc <- c("INCOME_5", v.main, "hhsize", "log_gdp", v.wgi,
                 v.categorical, v.cluster)


# =============================================================================
# CHECK 2: MISSING DATA
# -----------------------------------------------------------------------------
# WHY: the model deletes any row missing ANY variable in the formula. That is
# harmless only if the deleted rows are a random subset. If they differ
# systematically from the rows you keep, every coefficient is biased.
#
# NOTE: these figures now describe your ACTUAL analytic samples, so these are
# the numbers to report in the paper.
# =============================================================================

# 2a. Missingness per variable, within each analytic sample.
cat("\n=== FLI SAMPLE ===\n") #divide by total (79804)
d.fl %>%
  select(any_of(v.model.fl)) %>%
  naniar::miss_var_summary() %>%
  print(n = 30)

cat("\n=== INCOME SAMPLE ===\n")
d.inc %>%
  select(any_of(v.model.inc)) %>%
  naniar::miss_var_summary() %>%
  print(n = 30)

# 2b. Which variables go missing TOGETHER?
#     All bars are COUNTS: top bars = rows in that combination, left bars =
#     total missing per variable, dots = which variables are involved.
naniar::gg_miss_upset(d.fl  %>% select(any_of(v.model.fl)))
naniar::gg_miss_upset(d.inc %>% select(any_of(v.model.inc)))

# 2c. How many rows does each model lose to covariate missingness, on top of
#     the outcome filtering already applied?
d.fl  <- d.fl  %>% mutate(cc = complete.cases(select(., any_of(v.model.fl))))
d.inc <- d.inc %>% mutate(cc = complete.cases(select(., any_of(v.model.inc))))

cat(sprintf("\nFLI:    %d of %d retained (%.1f%% lost to covariates)\n",
            sum(d.fl$cc), nrow(d.fl), 100 * mean(!d.fl$cc)))
cat(sprintf("Income: %d of %d retained (%.1f%% lost to covariates)\n",
            sum(d.inc$cc), nrow(d.inc), 100 * mean(!d.inc$cc)))

# 2d. Do dropped rows DIFFER from kept rows?
#     Under MCAR the groups look alike. A large gap means listwise deletion is
#     not neutral, and you say so in the limitations.
bind_rows(
  d.fl %>% group_by(complete = cc) %>%
    summarise(model = "INDEX_FL_ord", n = n(),
              across(all_of(c("iwisescore", "hhsize", "gdp_pc_ppp")),
                     ~round(mean(.x, na.rm = TRUE), 2)), .groups = "drop"),
  d.inc %>% group_by(complete = cc) %>%
    summarise(model = "INCOME_5", n = n(),
              across(all_of(c("iwisescore", "hhsize", "gdp_pc_ppp")),
                     ~round(mean(.x, na.rm = TRUE), 2)), .groups = "drop")
) %>% relocate(model) %>% print()

# 2e. Which countries survive, and which are lost?
#     A country at 100% means the module was not fielded there - a COVERAGE
#     gap, not response bias. Imputation cannot help with those.
cat("\n--- Country coverage: full data vs analytic samples ---\n")
d.iwise %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(n_full = n(), .groups = "drop") %>%
  left_join(d.fl  %>% filter(cc) %>% count(across(all_of(v.cluster)),
                             name = "n_fl"),  by = v.cluster) %>%
  left_join(d.inc %>% filter(cc) %>% count(across(all_of(v.cluster)),
                             name = "n_inc"), by = v.cluster) %>%
  mutate(across(c(n_fl, n_inc), ~replace_na(.x, 0L)),
         pct_lost_fl  = round(100 * (1 - n_fl  / n_full), 1),
         pct_lost_inc = round(100 * (1 - n_inc / n_full), 1)) %>%
  arrange(desc(pct_lost_fl)) %>%
  print(n = Inf)

cat("\nCountries in FLI model:   ",
    d.fl %>% filter(cc) %>% distinct(across(all_of(v.cluster))) %>% nrow(), "\n")
cat("Countries in income model:",
    d.inc %>% filter(cc) %>% distinct(across(all_of(v.cluster))) %>% nrow(), "\n")

# 2f. For countries losing most of their cases, WHICH variable is responsible?
#     The answer determines what to do:
#       missing OUTCOME           -> accept, nothing to model
#       missing COUNTRY covariate -> fix the merge, or drop that covariate
#       missing INDIVIDUAL covar. -> consider multiple imputation
diagnose_drops <- function(data, cc_flag, v.model, threshold = 0.5) {
  bad <- data %>%
    group_by(across(all_of(v.cluster))) %>%
    summarise(pct = mean(!.data[[cc_flag]]), .groups = "drop") %>%
    filter(pct > threshold) %>%
    pull(!!sym(v.cluster))
  if (!length(bad)) { cat("No country above the threshold.\n"); return(invisible(NULL)) }
  data %>%
    filter(.data[[v.cluster]] %in% bad) %>%
    group_by(across(all_of(v.cluster))) %>%
    summarise(across(any_of(v.model), ~round(100 * mean(is.na(.x)), 1)),
              .groups = "drop") %>%
    print(width = Inf)
}

cat("\n--- FLI: countries losing >50% ---\n")
diagnose_drops(d.fl, "cc", v.model.fl)
cat("\n--- Income: countries losing >50% ---\n")
diagnose_drops(d.inc, "cc", v.model.inc)

# Countries entirely absent from a sample never appear above, because they
# have no rows left to summarise. Find them against the full data:
cat("\n--- Countries absent from the FLI sample entirely ---\n")
setdiff(unique(d.iwise[[v.cluster]]), unique(d.fl[[v.cluster]]))
cat("--- Countries absent from the income sample entirely ---\n")
setdiff(unique(d.iwise[[v.cluster]]), unique(d.inc[[v.cluster]]))

# For those, check what is missing in the FULL data:
d.iwise %>%
  filter(.data[[v.cluster]] %in%
           setdiff(unique(d.iwise[[v.cluster]]), unique(d.fl[[v.cluster]]))) %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(across(any_of(v.model.fl), ~round(100 * mean(is.na(.x)), 1)),
            .groups = "drop") %>%
  print(width = Inf)


# =============================================================================
# CHECK 3: OUTCOME DISTRIBUTION AND SPARSE CATEGORIES
# -----------------------------------------------------------------------------
# WHY: one cut-point is estimated per threshold (6 for a 7-level outcome, 4 for
# a 5-level one), each from respondents near that threshold. A level with very
# few cases gives an unstable cut-point.
#
# JUDGE ON COUNTS, NOT PERCENTAGES. With N in the tens of thousands, a level
# holding 3% of the sample still has well over a thousand respondents, which
# is ample. There is no canonical minimum percentage in the literature.
#
# AND NOTE: because the index is a deterministic function of four items, a
# level's frequency reflects how many response patterns map onto it - not
# sampling sparsity. A structurally rare score is not a data problem.
# =============================================================================

cat("\n--- INDEX_FL_ord (all respondents with a computed index) ---\n")
print(table(d.fl$INDEX_FL_ord, useNA = "ifany"))
print(round(100 * prop.table(table(d.fl$INDEX_FL_ord)), 1))

cat("\n--- INCOME_5 ---\n")
print(table(d.inc$INCOME_5, useNA = "ifany"))
print(round(100 * prop.table(table(d.inc$INCOME_5)), 1))

cat("\nSmallest FLI level:   ", min(table(d.fl$INDEX_FL_ord)), "cases\n")
cat("Smallest income level:", min(table(d.inc$INCOME_5)), "cases\n")
# In the low thousands = fine. In the low dozens = consider collapsing.

# Visual.
d.fl %>%
  count(INDEX_FL_ord) %>%
  ggplot(aes(INDEX_FL_ord, n)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.4, size = 3) +
  labs(title = "INDEX_FL_ord, all respondents with a computed index",
       subtitle = paste0("N = ", nrow(d.fl),
                         "; levels 25 and 75 arise only from 3-item responses"),
       x = "Financial Life Index", y = "n") +
  theme_minimal()

# 3b. SENSITIVITY: are 3-item and 4-item responders measured comparably?
# -----------------------------------------------------------------------------
# This is the check that justifies pooling them. Levels 25 and 75 exist ONLY in
# the 3-item group, so the two groups do not span the same set of possible
# scores. Three things to look at:

# (i) Which levels each group can actually reach.
cat("\n--- INDEX_FL_ord by number of items answered ---\n")
print(table(d.fl$INDEX_FL_ord, d.fl$n_fli_answered))

# (ii) Do the two groups differ on your key variables? If they look alike,
#      pooling is defensible and you say so in a footnote. If they differ,
#      report it as a limitation and run the 4-item-only model as a robustness
#      check (see 3c).
d.fl %>%
  group_by(fli_4item) %>%
  summarise(n = n(),
            mean_index = round(mean(as.numeric(as.character(INDEX_FL_ord))), 2),
            across(all_of(c("iwisescore", "hhsize", "gdp_pc_ppp")),
                   ~round(mean(.x, na.rm = TRUE), 2)),
            pct_female = round(100 * mean(female == "Female", na.rm = TRUE), 1),
            pct_urban  = round(100 * mean(urban_imp == "Peri-urban/urban",
                                          na.rm = TRUE), 1),
            .groups = "drop")

# (iii) Is the 3-item group concentrated in particular countries? If one item
#       was not fielded in certain surveys, this is a DESIGN feature, not
#       respondent refusal - which strengthens the case for pooling but also
#       means the 25/75 levels are country-specific artefacts.
d.fl %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(n = n(), pct_3item = round(100 * mean(!fli_4item), 1),
            .groups = "drop") %>%
  arrange(desc(pct_3item)) %>%
  print(n = Inf)

# (iv) Among 3-item responders, WHICH item is missing? If it is nearly always
#      the same one, that confirms a design/fielding explanation.
d.fl %>%
  filter(!fli_4item) %>%
  summarise(across(all_of(v.fli.items), ~sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "item", values_to = "n_missing") %>%
  arrange(desc(n_missing))

# 3c. ROBUSTNESS SAMPLE. Keep this defined so you can refit the final model on
#     4-item responders only and confirm your conclusions do not depend on the
#     pooling decision. Note it has 7 levels, not 9.
d.fl.4item <- d.fl %>%
  filter(fli_4item) %>%
  mutate(INDEX_FL_ord = droplevels(INDEX_FL_ord))

cat("\nRobustness sample (4-item only): N =", nrow(d.fl.4item),
    "| levels:", paste(levels(d.fl.4item$INDEX_FL_ord), collapse = ", "), "\n")


# ONLY IF a level turns out to be genuinely tiny. Collapse on substantive
# grounds and report that you did it.
# d.fl <- d.fl %>%
#   mutate(INDEX_FL_c = fct_collapse(INDEX_FL_ord,
#                                    "66.67+" = c("66.67", "83.33")))


# =============================================================================
# CHECK 4: PREDICTOR DISTRIBUTIONS AND RARE CATEGORIES
# -----------------------------------------------------------------------------
# WHY: a factor level with a handful of respondents gives a coefficient too
# imprecise to interpret. Implausible values are usually data errors. Heavy
# skew distorts estimates.
#
# Run on each analytic sample - the samples differ, so the distributions can too.
# =============================================================================

# 4a. Categorical predictors.
for (v in v.categorical) {
  cat("\n---", v, "(FLI sample) ---\n")
  print(table(d.fl[[v]], useNA = "ifany"))
}

# 4b. Continuous predictors.
d.fl %>%
  select(any_of(c("iwisescore", "hhsize", "gdp_pc_ppp", "log_gdp", v.wgi))) %>%
  summarise(across(everything(),
                   list(min  = ~min(.x, na.rm = TRUE),
                        med  = ~median(.x, na.rm = TRUE),
                        mean = ~mean(.x, na.rm = TRUE),
                        max  = ~max(.x, na.rm = TRUE),
                        sd   = ~sd(.x, na.rm = TRUE)))) %>%
  pivot_longer(everything(), names_to = c("var", "stat"),
               names_pattern = "^(.*)_(min|med|mean|max|sd)$") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  print(n = Inf)
# LOOKING FOR: impossible values; mean far from median (= skew).

# 4c. Histograms.
d.fl %>%
  select(any_of(c("iwisescore", "hhsize", "gdp_pc_ppp", "log_gdp"))) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 40) +
  facet_wrap(~name, scales = "free") +
  labs(title = "Continuous predictors, FLI sample") +
  theme_minimal()

# 4d. How is iwisescore distributed across the standard bands?
table(d.fl$iwise_cat, useNA = "ifany")


# =============================================================================
# CHECK 5: OUTCOME x CATEGORICAL PREDICTOR CROSS-TABS
# -----------------------------------------------------------------------------
# WHY: empty cells (zero respondents in a predictor-level x outcome-level
# combination) are the main cause of non-convergence and of separation.
# Re-run now, because the outcome has 7 levels instead of 9.
# =============================================================================

crosstab_check <- function(data, outcome, vars) {
  for (v in vars) {
    cat("\n===", outcome, "x", v, "===\n")
    tb <- table(data[[v]], data[[outcome]])
    print(tb)
    if (any(tb == 0))            cat(">>> WARNING: empty cell(s)\n")
    else if (any(tb < 5))        cat(">>> NOTE: cell(s) with fewer than 5 cases\n")
    else                         cat(">>> ok\n")
  }
}

crosstab_check(d.fl,  "INDEX_FL_ord", v.categorical)
crosstab_check(d.inc, "INCOME_5",     v.categorical)
# IF EMPTY CELLS: collapse the offending predictor category.


# =============================================================================
# CHECK 6: SAMPLE SIZE / EVENTS PER VARIABLE
# -----------------------------------------------------------------------------
# WHY: need enough cases in the SMALLEST outcome category to estimate all
# coefficients reliably. Rule of thumb: at least 10 per parameter.
# NOTE: a k-level factor costs k-1 parameters, not 1.
# Both N and the smallest category have changed, so this must be re-run.
# =============================================================================

rhs     <- paste("iwisescore + female + urban_imp + age_gp_profile +",
                 "maritalstatus + hhsize + employment + education")

# Contextual version. ONE governance term only - see Check 7.
# NB: this must match whichever WGI family v.wgi points at (line ~87).
rhs_ctx <- paste(rhs, "+ log_gdp +", v.wgi[6])

epv_check <- function(data, outcome, rhs) {
  d.cc <- data %>% drop_na(all_of(c(outcome, all.vars(as.formula(paste("~", rhs))))))
  n.p  <- ncol(model.matrix(as.formula(paste("~", rhs)), data = d.cc))
  tb   <- table(droplevels(d.cc[[outcome]]))
  cat(sprintf("\n%s\n  N = %d | parameters = %d | smallest category = %d | EPV = %.1f %s\n",
              outcome, nrow(d.cc), n.p, min(tb), min(tb) / n.p,
              ifelse(min(tb) / n.p < 10, "<- BELOW 10", "")))
}

epv_check(d.fl,  "INDEX_FL_ord", rhs_ctx)
epv_check(d.inc, "INCOME_5",     rhs_ctx)


# =============================================================================
# CHECK 7: MULTICOLLINEARITY
# -----------------------------------------------------------------------------
# WHY: if two predictors carry nearly the same information, the model cannot
# separate their effects - unstable coefficients, wide standard errors, signs
# that flip when you add or remove a variable.
#
# This is a property of the PREDICTORS ONLY. The outcome is irrelevant, which
# is why it can be checked before fitting. Run it on the FLI sample, since
# that is the smaller and more restricted of the two.
# =============================================================================

# 7a. Correlation matrix. Catches PAIRWISE overlap only.
d.fl %>%
  select(any_of(c("iwisescore", "hhsize", "log_gdp", v.wgi,
                  paste0(v.indices, "_num")))) %>%
  cor(use = "pairwise.complete.obs") %>%
  round(2) %>%
  print()
# LOOKING FOR: |r| above ~0.8. Expect it among the WGI indicators, and
# possibly between INDEX_YD and education (flagged in your codebook).

# 7b. VIF. Catches overlap involving THREE OR MORE variables, which pairwise
#     correlation misses. VIF = 4 means this standard error is twice as wide
#     as it would be with uncorrelated predictors.
#     A random outcome gives the identical answer - VIF ignores the outcome.
set.seed(2024)
m.vifprobe <- lm(as.formula(paste("rnorm(nrow(d.fl)) ~", rhs_ctx)), data = d.fl)
car::vif(m.vifprobe)
# For FACTORS this returns GVIF. Read GVIF^(1/(2*Df)) and SQUARE it before
# comparing to the usual thresholds of 5 or 10.

# 7c. All six WGI indicators together, to see how bad it is.
m.vifwgi <- lm(as.formula(paste("rnorm(nrow(d.fl)) ~",
                                paste(v.wgi, collapse = " + "))), data = d.fl)
car::vif(m.vifwgi)

# 7d. IF COLLINEAR: pick ONE on theoretical grounds, or reduce all six to a
#     single component. Run the PCA at COUNTRY level, not respondent level, or
#     the largest country dominates the components.
d.wgi.country <- d.iwise %>%
  select(all_of(c(v.cluster, v.wgi))) %>%
  distinct() %>%
  drop_na()

pca.wgi <- prcomp(d.wgi.country %>% select(all_of(v.wgi)), scale. = TRUE)
summary(pca.wgi)        # PC1 usually explains 80-90% of WGI variance
pca.wgi$rotation[, 1]   # loadings: all six should point the same direction

# If you use it, merge PC1 onto ALL THREE frames so they stay consistent.
d.pc1 <- d.wgi.country %>%
  mutate(wgi_pc1 = pca.wgi$x[, 1]) %>%
  select(all_of(v.cluster), wgi_pc1)

d.iwise <- left_join(d.iwise, d.pc1, by = v.cluster)
d.fl    <- left_join(d.fl,    d.pc1, by = v.cluster)
d.inc   <- left_join(d.inc,   d.pc1, by = v.cluster)


# =============================================================================
# CHECK 8: CLUSTERING / NON-INDEPENDENCE
# -----------------------------------------------------------------------------
# WHY: all regression assumes independent observations. Yours are nested in
# countries, and log_gdp / wgi_* take ONE value per country. Ignoring this
# makes standard errors far too small, so you find effects that are not there.
#
# The Moulton problem: tens of thousands of respondents, but only as many
# independent observations of the country-level variables as you have
# countries.
#
# The ICC uses an EMPTY model (no predictors), so it describes the data rather
# than your model - which is why it belongs here.
# =============================================================================

# 8a. How many clusters, and how big?
d.fl %>%
  count(across(all_of(v.cluster))) %>%
  summarise(n_countries = n(), min_n = min(n),
            median_n = median(n), max_n = max(n))
# Fewer than ~40 clusters means you need a small-sample correction (CR2) later.

# 8b. ICC: what share of outcome variance sits BETWEEN countries?
m.icc.fl <- lme4::lmer(
  as.numeric(as.character(INDEX_FL_ord)) ~ 1 + (1 | country), data = d.fl)
cat("\nICC, INDEX_FL_ord:\n"); print(performance::icc(m.icc.fl))

m.icc.inc <- lme4::lmer(
  as.numeric(as.character(INCOME_5)) ~ 1 + (1 | country), data = d.inc)
cat("\nICC, INCOME_5:\n"); print(performance::icc(m.icc.inc))

# LOOKING FOR: above ~0.05 means clustering is NOT ignorable.
# EXPECT INCOME_5 TO BE NEAR ZERO. Quintiles are constructed WITHIN country,
# so every country has ~20% in each - there is almost nothing between
# countries to explain, and country-level covariates cannot explain it.
# That is a substantive problem with the income model, not a bug.

# 8c. Visual: how much do country means differ?
d.fl %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(mean_fl = mean(as.numeric(as.character(INDEX_FL_ord)), na.rm = TRUE),
            .groups = "drop") %>%
  ggplot(aes(x = reorder(.data[[v.cluster]], mean_fl), y = mean_fl)) +
  geom_point() + coord_flip() +
  labs(title = "Country mean, INDEX_FL_ord", x = NULL, y = "Mean") +
  theme_minimal(base_size = 7)

# YOUR FIX (applied when you fit, not here):
#   Option A - cluster-robust SEs with a CR2 small-sample correction
#              (clubSandwich). Keeps coefficients, corrects inference.
#   Option B - multilevel ordinal model, ordinal::clmm() with (1 | country).
#              Preferred with country-level predictors.


# =============================================================================
# CHECK 9: SEPARATION
# -----------------------------------------------------------------------------
# WHY: separation means a predictor or combination perfectly predicts the
# outcome. Symptoms: coefficients above 15, standard errors in the thousands,
# convergence warnings. A DATA problem - more iterations never fix it.
#
# The ordinal model works by cutting the outcome at each threshold, so each
# cut can be tested in advance. 7 levels = 6 cuts (was 8).
# =============================================================================

check_separation <- function(data, outcome, rhs) {
  data <- data %>% drop_na(all_of(outcome))
  data[[outcome]] <- droplevels(data[[outcome]])
  lv <- levels(data[[outcome]])
  cat("\n--- Separation:", outcome, "(", length(lv) - 1, "cuts ) ---\n")
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

check_separation(d.fl,  "INDEX_FL_ord", rhs)
check_separation(d.inc, "INCOME_5",     rhs)

# IF DETECTED: collapse the offending category, drop the variable, or plan on
# Firth penalised likelihood (brglm2::brglmFit) at the modelling stage.


# =============================================================================
# WHAT COMES NEXT
# =============================================================================

 # 10. PROPORTIONAL ODDS - the assumption specific to this model. Tests whether
# each predictor's effect is the same across every cut of the outcome. It
# is frequently violated. Brant test, ordinal::nominal_test(), and a plot
# of coefficients across cuts. With your N the formal tests reject on
# trivial departures, so the plot is the most honest of the three.

 # 11. LINEARITY OF THE LOGIT - whether continuous predictors relate linearly
# to the log-odds. That scale only exists once a model defines it, and the
# relationship is PARTIAL, so it depends on the other covariates.

# 12. INFLUENTIAL OBSERVATIONS - Cook's distance, standardised residuals, and
# a sensitivity refit without the most influential 1%.

# 13. GOODNESS OF FIT - Lipsitz and Pulkstenis-Robinson tests.    