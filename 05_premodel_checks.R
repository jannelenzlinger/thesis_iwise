# =============================================================================
# 02 PRE-MODEL CHECKS  |  Ordinal logistic regression
# Outcomes: FLI_3item, INCOME_5 | Main predictor: iwisescore
#
# RULE: this script only READS data. It never modifies or saves d.fl / d.inc.
# If a check shows something needs changing, change 01_data_prep.R and re-run
# it before running this script again.
#
# STRUCTURE
# ---------
#   Section 0   Packages and data
#   Check 1     Variable types and level ordering (confirmation)
#   Check 2     Missing data (incl. model weights)
#   Check 3     Outcome distribution and sparse categories
#   Check 4     Predictor distributions and rare categories
#   Check 5     Outcome x categorical predictor cross-tabs
#   Check 6     Sample size / events per variable
#   Check 7     Multicollinearity
#   Check 8     Clustering / non-independence
#   Check 9     Separation
#   Next        Checks that require a fitted model
# =============================================================================


# ---- 0. PACKAGES AND DATA ---------------------------------------------------

pkgs <- c("tidyverse",        # data wrangling and plotting
          "naniar",           # missing data summaries and plots
          "car",              # vif(): multicollinearity
          "lme4",             # lmer(): empty model for the ICC
          "performance",      # icc()
          "detectseparation") # separation check

missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing_pkgs)) install.packages(missing_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

select <- dplyr::select   # MASS and others mask dplyr::select; be explicit

# Loads every data frame and variable list created in 01_data_prep.R
dat <- readRDS("data_prepared.rds")
list2env(dat, envir = .GlobalEnv)
rm(dat)


# =============================================================================
# CHECK 1: VARIABLE TYPES AND LEVEL ORDERING
# -----------------------------------------------------------------------------
# WHY: the model must know the outcome is ORDERED. Levels sorted alphabetically
# ("0","100","25") produce a model that runs fine and is completely wrong.
# The recode itself is in 01_data_prep.R (Section 4); this only confirms it.
# =============================================================================

class(d.fl$FLI_3item)
class(d.inc$INCOME_5)
levels(d.fl$FLI_3item)     # should be 0 < 25 < 50 < 75 < 100
levels(d.inc$INCOME_5)

str(d.iwise %>% select(any_of(c(v.outcomes, v.continuous, v.categorical))))


# ---- ANALYTIC SAMPLES -------------------------------------------------------

cat("\n--- ANALYTIC SAMPLES ---\n")
cat("Full data:    ", nrow(d.iwise), "\n")
cat("FLI sample:   ", nrow(d.fl),
    sprintf("(%.1f%% of full)\n", 100 * nrow(d.fl) / nrow(d.iwise)))
cat("Income sample:", nrow(d.inc),
    sprintf("(%.1f%% of full)\n", 100 * nrow(d.inc) / nrow(d.iwise)))

# FLI_3item should cover exactly the same respondents as Gallup's index.
cat("\nGallup index vs FLI_3item coverage:\n")
print(table(gallup = !is.na(d.iwise$INDEX_FL_ord),
            fli3   = !is.na(d.iwise$FLI_3item)))


# =============================================================================
# CHECK 2: MISSING DATA
# -----------------------------------------------------------------------------
# WHY: the model deletes any row missing ANY variable in the formula. Harmless
# only if the deleted rows are a random subset.
# These figures describe the ACTUAL analytic samples: report these numbers.
# bas_rate_pp is country-level, so where missing it is missing for a WHOLE
# country.
# =============================================================================

# 2a. Missingness per variable, within each analytic sample.
cat("\n=== FLI SAMPLE ===\n")
d.fl %>% select(any_of(v.model.fl)) %>%
  naniar::miss_var_summary() %>% print(n = 30)

cat("\n=== INCOME SAMPLE ===\n")
d.inc %>% select(any_of(v.model.inc)) %>%
  naniar::miss_var_summary() %>% print(n = 30)

# 2b. Which variables go missing TOGETHER?
naniar::gg_miss_upset(d.fl  %>% select(any_of(v.model.fl)))
naniar::gg_miss_upset(d.inc %>% select(any_of(v.model.inc)))

# 2c. Rows lost to covariate missingness (cc flag created in 01, Section 9).
cat(sprintf("\nFLI:    %d of %d retained (%.1f%% lost to covariates)\n",
            sum(d.fl$cc), nrow(d.fl), 100 * mean(!d.fl$cc)))
cat(sprintf("Income: %d of %d retained (%.1f%% lost to covariates)\n",
            sum(d.inc$cc), nrow(d.inc), 100 * mean(!d.inc$cc)))

# 2d. Do dropped rows DIFFER from kept rows?
bind_rows(
  d.fl %>% group_by(complete = cc) %>%
    summarise(model = "FLI_3item", n = n(),
              across(all_of(c("iwisescore", "hhsize", "gdp_pc_ppp", v.jmp)),
                     ~round(mean(.x, na.rm = TRUE), 2)), .groups = "drop"),
  d.inc %>% group_by(complete = cc) %>%
    summarise(model = "INCOME_5", n = n(),
              across(all_of(c("iwisescore", "hhsize", "gdp_pc_ppp", v.jmp)),
                     ~round(mean(.x, na.rm = TRUE), 2)), .groups = "drop")
) %>% relocate(model) %>% print()

# 2e. Which countries survive, and which are lost?
#     100% lost = module not fielded (coverage gap, not response bias).
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

# 2e-bis. Which countries have NO bas_rate_pp at all?
cat("\n--- Countries with no bas_rate_pp ---\n")
d.iwise %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(pct_missing_jmp = round(100 * mean(is.na(.data[[v.jmp]])), 1),
            .groups = "drop") %>%
  filter(pct_missing_jmp > 0) %>%
  print(n = Inf)

# 2f. For countries losing most of their cases, WHICH variable is responsible?
#       missing OUTCOME           -> accept, nothing to model
#       missing COUNTRY covariate -> fix the merge (in 01), or drop the covariate
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

# Countries entirely absent from a sample:
cat("\n--- Countries absent from the FLI sample entirely ---\n")
print(setdiff(unique(d.iwise[[v.cluster]]), unique(d.fl[[v.cluster]])))
cat("--- Countries absent from the income sample entirely ---\n")
print(setdiff(unique(d.iwise[[v.cluster]]), unique(d.inc[[v.cluster]])))

d.iwise %>%
  filter(.data[[v.cluster]] %in%
           setdiff(unique(d.iwise[[v.cluster]]), unique(d.fl[[v.cluster]]))) %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(across(any_of(c(v.model.fl, "WP2319", "WP30", "WP31")),
                   ~round(100 * mean(is.na(.x)), 1)),
            .groups = "drop") %>%
  print(width = Inf)

# 2g. MODEL WEIGHTS (w_fit, created in 01 Section 9)
#     w_fit should sum to the complete-case N. Extreme values mean a few
#     respondents carry a lot of influence: compare weighted and unweighted
#     models later if so.
check_weights <- function(data, label) {
  w <- data$w_fit[data$cc]
  cat(sprintf("\n%s: sum(w_fit) = %.0f | complete cases = %d | missing weight = %d\n",
              label, sum(w), sum(data$cc), sum(is.na(data[[v.weight]]))))
  print(round(quantile(w, c(0, .01, .5, .99, 1)), 2))
}
check_weights(d.fl,  "FLI")
check_weights(d.inc, "Income")


# =============================================================================
# CHECK 3: OUTCOME DISTRIBUTION AND SPARSE CATEGORIES
# -----------------------------------------------------------------------------
# WHY: one cut-point per threshold, each estimated from respondents near it.
# JUDGE ON COUNTS, NOT PERCENTAGES. A structurally rare score (few response
# patterns map onto it) is not a data problem.
# =============================================================================

cat("\n--- FLI_3item ---\n")
print(table(d.fl$FLI_3item, useNA = "ifany"))
print(round(100 * prop.table(table(d.fl$FLI_3item)), 1))

cat("\n--- INCOME_5 ---\n")
print(table(d.inc$INCOME_5, useNA = "ifany"))
print(round(100 * prop.table(table(d.inc$INCOME_5)), 1))

cat("\nSmallest FLI level:   ", min(table(d.fl$FLI_3item)), "cases\n")
cat("Smallest income level:", min(table(d.inc$INCOME_5)), "cases\n")
# Low thousands = fine. Low dozens = consider collapsing (in 01).

d.fl %>%
  count(FLI_3item) %>%
  ggplot(aes(FLI_3item, n)) +
  geom_col(fill = "#3D4222") +
  geom_text(aes(label = n), vjust = -0.4, size = 3) +
  labs(title = "Financial Life Index (3-item version)",
       subtitle = paste0("N = ", nrow(d.fl),
                         "; WP2319 + WP30 + WP31, uniform construction"),
       x = "Financial Life Index", y = "n") +
  theme_minimal()


# 3b. SENSITIVITY: does WP88 availability still matter?
# (i) Both groups should reach all five levels.
cat("\n--- FLI_3item by WP88 availability ---\n")
print(table(d.fl$FLI_3item, has_wp88 = d.fl$has_wp88))

# (ii) How the new index maps onto Gallup's.
cat("\n--- FLI_3item (rows) vs Gallup INDEX_FL_ord (cols) ---\n")
print(table(d.fl$FLI_3item, d.fl$INDEX_FL_ord))

cat("\nCorrelation between the two versions: ",
    round(cor(as.numeric(as.character(d.fl$FLI_3item)),
              as.numeric(as.character(d.fl$INDEX_FL_ord)),
              use = "pairwise.complete.obs"), 3), "\n")

# (iii) Do respondents with and without WP88 differ on key variables?
d.fl %>%
  group_by(has_wp88) %>%
  summarise(n = n(),
            mean_index = round(mean(as.numeric(as.character(FLI_3item))), 2),
            across(all_of(c("iwisescore", "hhsize", "gdp_pc_ppp", v.jmp)),
                   ~round(mean(.x, na.rm = TRUE), 2)),
            pct_female = round(100 * mean(female == "Female", na.rm = TRUE), 1),
            pct_urban  = round(100 * mean(urban_imp == "Peri-urban/urban",
                                          na.rm = TRUE), 1),
            .groups = "drop")

# (iv) Is WP88 non-administration a COUNTRY-level pattern?
d.fl %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(n = n(), pct_no_wp88 = round(100 * mean(!has_wp88), 1),
            .groups = "drop") %>%
  arrange(desc(pct_no_wp88)) %>%
  print(n = Inf)

# 3c. ROBUSTNESS SAMPLE (created in 01 Section 10).
cat("\nRobustness sample (Gallup 4-item, WP88 countries only): N =",
    nrow(d.fl.gallup), "| complete cases =", sum(d.fl.gallup$cc),
    "| levels:", paste(levels(d.fl.gallup$INDEX_FL_ord), collapse = ", "), "\n")


# =============================================================================
# CHECK 4: PREDICTOR DISTRIBUTIONS AND RARE CATEGORIES
# -----------------------------------------------------------------------------
# WHY: rare factor levels give imprecise coefficients; implausible values are
# usually data errors; heavy skew distorts estimates.
# =============================================================================

# 4a. Categorical predictors.
for (v in v.categorical) {
  cat("\n---", v, "(FLI sample) ---\n")
  print(table(d.fl[[v]], useNA = "ifany"))
}

# 4b. Continuous predictors.
d.fl %>%
  select(any_of(c("iwisescore", "hhsize", "gdp_pc_ppp", "log_gdp", v.wgi,
                  v.jmp))) %>%
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
  select(any_of(c("iwisescore", "hhsize", "gdp_pc_ppp", "log_gdp", v.jmp))) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 40, fill = "#3D4222") +
  facet_wrap(~name, scales = "free") +
  labs(title = "Continuous predictors, FLI sample") +
  theme_minimal()

# 4c-bis. How much does bas_rate_pp vary across COUNTRY-YEARS?
cy %>%
  summarise(n_country_years = sum(!is.na(bas_rate_pp)),
            median   = median(bas_rate_pp, na.rm = TRUE),
            iqr_lo   = quantile(bas_rate_pp, .25, na.rm = TRUE),
            iqr_hi   = quantile(bas_rate_pp, .75, na.rm = TRUE),
            pct_flat = round(100 * mean(abs(bas_rate_pp) < 0.05, na.rm = TRUE), 1))

# 4d. iwisescore across the standard bands.
table(d.fl$iwise_cat, useNA = "ifany")


# =============================================================================
# CHECK 5: OUTCOME x CATEGORICAL PREDICTOR CROSS-TABS
# -----------------------------------------------------------------------------
# WHY: empty cells cause non-convergence and separation.
# =============================================================================

crosstab_check <- function(data, outcome, vars) {
  for (v in vars) {
    cat("\n===", outcome, "x", v, "===\n")
    tb <- table(data[[v]], data[[outcome]])
    print(tb)
    if (any(tb == 0))     cat(">>> WARNING: empty cell(s)\n")
    else if (any(tb < 5)) cat(">>> NOTE: cell(s) with fewer than 5 cases\n")
    else                  cat(">>> ok\n")
  }
}

crosstab_check(d.fl,  "FLI_3item", v.categorical)
crosstab_check(d.inc, "INCOME_5",  v.categorical)
# IF EMPTY CELLS: collapse the offending category (in 01).


# =============================================================================
# CHECK 6: SAMPLE SIZE / EVENTS PER VARIABLE
# -----------------------------------------------------------------------------
# WHY: need enough cases in the SMALLEST outcome category. Rule of thumb: at
# least 10 per parameter. A k-level factor costs k-1 parameters.
# rhs and rhs_ctx are defined in 01 Section 8 (rhs_ctx uses GE only).
# =============================================================================

epv_check <- function(data, outcome, rhs) {
  d.cc <- data %>% drop_na(all_of(c(outcome, all.vars(as.formula(paste("~", rhs))))))
  n.p  <- ncol(model.matrix(as.formula(paste("~", rhs)), data = d.cc))
  tb   <- table(droplevels(d.cc[[outcome]]))
  cat(sprintf("\n%s\n  N = %d | parameters = %d | smallest category = %d | EPV = %.1f %s\n",
              outcome, nrow(d.cc), n.p, min(tb), min(tb) / n.p,
              ifelse(min(tb) / n.p < 10, "<- BELOW 10", "")))
}

epv_check(d.fl,  "FLI_3item", rhs_ctx)
epv_check(d.inc, "INCOME_5",  rhs_ctx)


# =============================================================================
# CHECK 7: MULTICOLLINEARITY
# -----------------------------------------------------------------------------
# WHY: overlapping predictors give unstable coefficients and wide SEs.
# A property of the PREDICTORS ONLY, so it can be checked before fitting.
# =============================================================================

# 7a. Correlation matrix (all six WGI, to document their overlap).
d.fl %>%
  select(any_of(c("iwisescore", "hhsize", "log_gdp", v.wgi, v.jmp,
                  paste0(v.indices, "_num")))) %>%
  cor(use = "pairwise.complete.obs") %>%
  round(2) %>%
  print()
# LOOKING FOR: |r| above ~0.8. Expect it among the WGI indicators, which is
# why the main model uses GE only.

# 7b. VIF on the main model (GE only). A random outcome gives the identical
#     answer - VIF ignores the outcome.
set.seed(2024)
m.vifprobe <- lm(as.formula(paste("rnorm(nrow(d.fl)) ~", rhs_ctx)), data = d.fl)
car::vif(m.vifprobe)
# For FACTORS read GVIF^(1/(2*Df)), SQUARED, against thresholds of 5 or 10.
# Earlier run with wgi_cc_sc: log_gdp 2.00, wgi_cc_sc 1.73, iwisescore 1.09.
# Re-check now that GE replaces CC.

# 7c. Country-year collinearity: IWISE vs JMP (cy created in 01 Section 6).
cy %>%
  select(-country_year) %>%
  cor(use = "pairwise.complete.obs", method = "spearman") %>%
  round(2)
# RESULT: bas_rate_pp x IWISE = 0.14. JMP levels -0.53 to -0.69 with IWISE and
# 0.79-0.96 with each other, so use at most ONE JMP level in any model.


# =============================================================================
# CHECK 8: CLUSTERING / NON-INDEPENDENCE
# -----------------------------------------------------------------------------
# WHY: respondents are nested in countries, and log_gdp / wgi / bas_rate_pp
# take ONE value per country-year. Ignoring this makes SEs far too small.
# The ICC uses an EMPTY model, so it describes the data, not your model.
# =============================================================================

# 8a. How many clusters, and how big?
d.fl %>%
  count(across(all_of(v.cluster))) %>%
  summarise(n_countries = n(), min_n = min(n),
            median_n = median(n), max_n = max(n))
# 74 countries, 500 to 3,503 respondents each.

d.fl %>% filter(cc) %>% count(across(all_of(v.cluster))) %>%
  summarise(n_countries_cc = n())

# 8b. ICC: share of outcome variance BETWEEN countries.
icc_formula <- function(outcome) {
  as.formula(paste0("as.numeric(as.character(", outcome, ")) ~ 1 + (1 | ",
                    v.cluster, ")"))
}

m.icc.fl <- lme4::lmer(icc_formula("FLI_3item"), data = d.fl)
cat("\nICC, FLI_3item:\n"); print(performance::icc(m.icc.fl))

m.icc.inc <- lme4::lmer(icc_formula("INCOME_5"), data = d.inc)
cat("\nICC, INCOME_5:\n"); print(performance::icc(m.icc.inc))
# Above ~0.05 = clustering NOT ignorable.
# EXPECT INCOME_5 NEAR ZERO: quintiles are constructed WITHIN country.

# 8c. Visual: country means.
d.fl %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(mean_fl = mean(as.numeric(as.character(FLI_3item)), na.rm = TRUE),
            .groups = "drop") %>%
  ggplot(aes(x = reorder(.data[[v.cluster]], mean_fl), y = mean_fl)) +
  geom_point() + coord_flip() +
  labs(title = "Country mean, FLI_3item", x = NULL, y = "Mean") +
  theme_minimal(base_size = 7)

# FIX (applied when fitting, not here): polr() with weights = w_fit and CR2
# cluster-robust SEs (clubSandwich).


# =============================================================================
# CHECK 9: SEPARATION
# -----------------------------------------------------------------------------
# WHY: a predictor that perfectly predicts the outcome gives huge coefficients
# and SEs. A DATA problem - more iterations never fix it.
# Each cut of the ordinal outcome is tested as a binary model.
# Uses the individual-level rhs; swap in rhs_ctx for the full model.
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

check_separation(d.fl,  "FLI_3item", rhs)
check_separation(d.inc, "INCOME_5",  rhs)
# IF DETECTED: collapse the category (in 01), drop the variable, or use Firth
# penalised likelihood (brglm2::brglmFit) at the modelling stage.


# =============================================================================
# WHAT COMES NEXT - checks that require a FITTED model (03_models.R)
# =============================================================================

# 10. PROPORTIONAL ODDS - Brant test, ordinal::nominal_test(), and a plot of
#     coefficients across cuts. With your N the formal tests reject on trivial
#     departures, so the plot is the most honest of the three.

# 11. LINEARITY OF THE LOGIT - iwisescore is zero-inflated (median 3, mean
#     6.6), so expect this to flag; iwise_cat is the alternative. Test
#     bas_rate_pp with splines (it can be negative).

# 12. INFLUENTIAL OBSERVATIONS - Cook's distance, standardised residuals, refit
#     without the most influential 1%. Also refit dropping one country at a
#     time and watch the country-level estimates.

# 13. GOODNESS OF FIT - Lipsitz and Pulkstenis-Robinson tests.