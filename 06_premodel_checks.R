# =============================================================================
# PRE-MODEL CHECKS  |  Ordinal logistic regression
# Outcomes: FLI_3item, INCOME_5 | Main predictor: iwisescore
#
# RULE: this script only READS data. It never modifies or saves d.fl / d.inc.
# If a check shows something needs changing, change 05_data_prep.R and re-run
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

missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))   # packages not yet installed
if (length(missing_pkgs)) install.packages(missing_pkgs)        # install only the missing ones
invisible(lapply(pkgs, library, character.only = TRUE))         # load all packages silently

select <- dplyr::select   # MASS and others mask dplyr::select; be explicit

# Loads every data frame and variable list created in 05_data_prep.R
dat <- readRDS("data/iwise_data_prepared.rds")                  # read saved list from 05
list2env(dat, envir = .GlobalEnv)                               # unpack list items into workspace
rm(dat)                                                         # remove the now-redundant list


# =============================================================================
# CHECK 1: VARIABLE TYPES AND LEVEL ORDERING
# -----------------------------------------------------------------------------
# WHY: the model must know the outcome is ORDERED. Levels sorted alphabetically
# ("0","100","25") produce a model that runs fine and is completely wrong.
# The recode itself is in 01_data_prep.R (Section 4); this only confirms it.
# =============================================================================

class(d.fl$FLI_3item)      # should be "ordered" "factor"
class(d.inc$INCOME_5)      # should be "ordered" "factor"
levels(d.fl$FLI_3item)     # should be 0 < 25 < 50 < 75 < 100
levels(d.inc$INCOME_5)     # should be quintiles in order, lowest first

str(d.iwise %>% select(any_of(c(v.outcomes, v.continuous, v.categorical))))  # type of every model variable


# ---- ANALYTIC SAMPLES -------------------------------------------------------

cat("\n--- ANALYTIC SAMPLES ---\n")                                # section header
cat("Full data:    ", nrow(d.iwise), "\n")                         # total respondents
cat("FLI sample:   ", nrow(d.fl),                                  # respondents with FLI_3item
    sprintf("(%.1f%% of full)\n", 100 * nrow(d.fl) / nrow(d.iwise)))   # as % of full data
cat("Income sample:", nrow(d.inc),                                 # respondents with INCOME_5
    sprintf("(%.1f%% of full)\n", 100 * nrow(d.inc) / nrow(d.iwise)))  # as % of full data

# FLI_3item should cover exactly the same respondents as Gallup's index.
cat("\nGallup index vs FLI_3item coverage:\n")                     # section header
print(table(gallup = !is.na(d.iwise$INDEX_FL_ord),                 # has Gallup index (TRUE/FALSE)
            fli3   = !is.na(d.iwise$FLI_3item)))                   # x has FLI_3item (TRUE/FALSE)


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
cat("\n=== FLI SAMPLE ===\n")                                      # section header
d.fl %>% select(any_of(v.model.fl)) %>%                            # FLI model variables only
  naniar::miss_var_summary() %>% print(n = 30)                     # n and % missing per variable

cat("\n=== INCOME SAMPLE ===\n")                                   # section header

d.inc %>% select(any_of(v.model.inc)) %>%                          # income model variables only
  naniar::miss_var_summary() %>% print(n = 30)                     # n and % missing per variable

# 2b. Which variables go missing TOGETHER?
naniar::gg_miss_upset(d.fl  %>% select(any_of(v.model.fl)))        # missingness patterns, FLI
naniar::gg_miss_upset(d.inc %>% select(any_of(v.model.inc)))       # missingness patterns, income

# 2c. Rows lost to covariate missingness (cc flag created in 05, Section 9).
cat(sprintf("\nFLI:    %d of %d retained (%.1f%% lost to covariates)\n",   # template for FLI line
            sum(d.fl$cc), nrow(d.fl), 100 * mean(!d.fl$cc)))       # kept n, total n, % lost
cat(sprintf("Income: %d of %d retained (%.1f%% lost to covariates)\n",     # template for income line
            sum(d.inc$cc), nrow(d.inc), 100 * mean(!d.inc$cc)))    # kept n, total n, % lost

# 2d. Do dropped rows DIFFER from kept rows?
bind_rows(                                                         # stack FLI and income results
  d.fl %>% group_by(complete = cc) %>%                             # FLI: split kept vs dropped
    summarise(model = "FLI_3item", n = n(),                        # label and group size
              across(all_of(c("iwisescore", "hhsize", "gdp_pc_ppp", v.jmp, v.gbd)),  # variables to compare
                     ~round(mean(.x, na.rm = TRUE), 2)), .groups = "drop"),          # mean per group
  d.inc %>% group_by(complete = cc) %>%                            # income: split kept vs dropped
    summarise(model = "INCOME_5", n = n(),                         # label and group size
              across(all_of(c("iwisescore", "hhsize", "gdp_pc_ppp", v.jmp, v.gbd)),  # variables to compare
                     ~round(mean(.x, na.rm = TRUE), 2)), .groups = "drop")           # mean per group
) %>% relocate(model) %>% print()                                  # model column first, then print

# 2e. Which countries survive, and which are lost?
#     100% lost = module not surveyed (coverage gap, not response bias).
cat("\n--- Country coverage: full data vs analytic samples ---\n")  # section header
d.iwise %>%                                                        # start from full data
  group_by(across(all_of(v.cluster))) %>%                          # one group per country-year
  summarise(n_full = n(), .groups = "drop") %>%                    # respondents in full data
  left_join(d.fl  %>% filter(cc) %>% count(across(all_of(v.cluster)),     # FLI complete cases per country-year
                                           name = "n_fl"),  by = v.cluster) %>%   # add as n_fl
  left_join(d.inc %>% filter(cc) %>% count(across(all_of(v.cluster)),     # income complete cases per country-year
                                           name = "n_inc"), by = v.cluster) %>%   # add as n_inc
  mutate(across(c(n_fl, n_inc), ~replace_na(.x, 0L)),              # absent country-year = 0 cases
         pct_lost_fl  = round(100 * (1 - n_fl  / n_full), 1),      # % lost from FLI model
         pct_lost_inc = round(100 * (1 - n_inc / n_full), 1)) %>%  # % lost from income model
  arrange(desc(pct_lost_fl)) %>%                                   # worst FLI losses first
  print(n = Inf)                                                   # print all rows

cat("\nCountries in FLI model:   ",                                # label
    d.fl %>% filter(cc) %>% distinct(across(all_of(v.cluster))) %>% nrow(), "\n")   # n country-years with complete cases
cat("Countries in income model:",                                  # label
    d.inc %>% filter(cc) %>% distinct(across(all_of(v.cluster))) %>% nrow(), "\n")  # n country-years with complete cases

# 2e-bis. Which country-years are missing any country-level covariate?
v.ctx <- c("log_gdp", v.wgi.main, v.jmp, v.gbd)                 # country-level covariates in rhs_ctx
d.iwise %>%
  group_by(across(all_of(v.cluster))) %>%                       # one row per country-year
  summarise(across(all_of(v.ctx), ~round(100 * mean(is.na(.x)), 1)),  # % missing per variable
            .groups = "drop") %>%
  filter(if_any(all_of(v.ctx), ~.x > 0)) %>%                    # keep country-years with any gap
  print(n = Inf)

# 2f. For countries losing most of their cases, WHICH variable is responsible?
#       missing OUTCOME
#       missing COUNTRY covariate 
#       missing INDIVIDUAL covar. 
diagnose_drops <- function(data, cc_flag, v.model, threshold = 0.5) {   # define helper function
  bad <- data %>%                                                  # find high-loss country-years
    group_by(across(all_of(v.cluster))) %>%                        # one group per country-year
    summarise(pct = mean(!.data[[cc_flag]]), .groups = "drop") %>%  # share of rows not complete
    filter(pct > threshold) %>%                                    # keep those above threshold (50%)
    pull(!!sym(v.cluster))                                         # return their IDs as a vector
  if (!length(bad)) { cat("No country above the threshold.\n"); return(invisible(NULL)) }  # stop if none
  data %>%                                                         # for the high-loss ones
    filter(.data[[v.cluster]] %in% bad) %>%                        # keep only those country-years
    group_by(across(all_of(v.cluster))) %>%                        # one group per country-year
    summarise(across(any_of(v.model), ~round(100 * mean(is.na(.x)), 1)),  # % missing per model variable
              .groups = "drop") %>%                                # remove grouping
    print(width = Inf)                                             # print all columns
}                                                                  # end of function

cat("\n--- FLI: countries losing >50% ---\n")                      # section header
diagnose_drops(d.fl, "cc", v.model.fl)                             # run for FLI sample
cat("\n--- Income: countries losing >50% ---\n")                   # section header
diagnose_drops(d.inc, "cc", v.model.inc)                           # run for income sample

# Countries entirely absent from a sample:
cat("\n--- Countries absent from the FLI sample entirely ---\n")   # section header
print(setdiff(unique(d.iwise[[v.cluster]]), unique(d.fl[[v.cluster]])))    # in full data but not FLI
cat("--- Countries absent from the income sample entirely ---\n")  # section header
print(setdiff(unique(d.iwise[[v.cluster]]), unique(d.inc[[v.cluster]])))   # in full data but not income

d.iwise %>%                                                        # start from full data
  filter(.data[[v.cluster]] %in%                                   # keep country-years that are...
           setdiff(unique(d.iwise[[v.cluster]]), unique(d.fl[[v.cluster]]))) %>%   # ...absent from FLI
  group_by(across(all_of(v.cluster))) %>%                          # one group per country-year
  summarise(across(any_of(c(v.model.fl, "WP2319", "WP30", "WP31")),         # model vars + FLI items
                   ~round(100 * mean(is.na(.x)), 1)),              # % missing each
            .groups = "drop") %>%                                  # remove grouping
  print(width = Inf)                                               # print all columns

# 2g. MODEL WEIGHTS (w_fit, created in 05 Section 9)
#     w_fit should sum to the complete-case N. Extreme values mean a few
#     respondents carry a lot of influence: compare weighted and unweighted
#     models later if so.
check_weights <- function(data, label) {                           # define helper function
  w <- data$w_fit[data$cc]                                         # weights of complete cases only
  cat(sprintf("\n%s: sum(w_fit) = %.0f | complete cases = %d | missing weight = %d\n",   # output template
              label, sum(w), sum(data$cc), sum(is.na(data[[v.weight]]))))  # sum of weights, n cc, n missing weight
  print(round(quantile(w, c(0, .01, .5, .99, 1)), 2))              # min, 1%, median, 99%, max
}                                                                  # end of function
check_weights(d.fl,  "FLI")                                        # run for FLI sample
check_weights(d.inc, "Income")                                     # run for income sample


# =============================================================================
# CHECK 3: OUTCOME DISTRIBUTION AND SPARSE CATEGORIES
# -----------------------------------------------------------------------------
# WHY: one cut-point per threshold, each estimated from respondents near it.
# JUDGE ON COUNTS, NOT PERCENTAGES. A structurally rare score (few response
# patterns map onto it) is not a data problem.
# =============================================================================

cat("\n--- FLI_3item ---\n")                                        # section header
print(table(d.fl$FLI_3item, useNA = "ifany"))                      # counts per FLI level
print(round(100 * prop.table(table(d.fl$FLI_3item)), 1))           # % per FLI level

cat("\n--- INCOME_5 ---\n")                                        # section header
print(table(d.inc$INCOME_5, useNA = "ifany"))                      # counts per quintile
print(round(100 * prop.table(table(d.inc$INCOME_5)), 1))           # % per quintile

cat("\nSmallest FLI level:   ", min(table(d.fl$FLI_3item)), "cases\n")   # smallest FLI category
cat("Smallest income level:", min(table(d.inc$INCOME_5)), "cases\n")     # smallest income category
# Low thousands = fine. Low dozens = consider collapsing (in 01).

d.fl %>%                                                           # FLI sample
  count(FLI_3item) %>%                                             # n per level
  ggplot(aes(FLI_3item, n)) +                                      # level on x, count on y
  geom_col(fill = "#3D4222") +                                     # bar chart
  geom_text(aes(label = n), vjust = -0.4, size = 3) +              # count above each bar
  labs(title = "Financial Life Index (3-item version)",            # plot title
       subtitle = paste0("N = ", nrow(d.fl),                       # subtitle with sample size
                         "; WP2319 + WP30 + WP31, uniform construction"),   # and item description
       x = "Financial Life Index", y = "n") +                      # axis labels
  theme_minimal()                                                  # clean theme


# 3b. SENSITIVITY: does WP88 availability still matter?
# (i) Both groups should reach all five levels.
cat("\n--- FLI_3item by WP88 availability ---\n")                  # section header
print(table(d.fl$FLI_3item, has_wp88 = d.fl$has_wp88))             # FLI level x WP88 fielded

# (ii) How the new index maps onto Gallup's.
cat("\n--- FLI_3item (rows) vs Gallup INDEX_FL_ord (cols) ---\n")  # section header
print(table(d.fl$FLI_3item, d.fl$INDEX_FL_ord))                    # cross-tab of the two indices

cat("\nCorrelation between the two versions: ",                    # label
    round(cor(as.numeric(as.character(d.fl$FLI_3item)),            # FLI_3item as numeric score
              as.numeric(as.character(d.fl$INDEX_FL_ord)),         # Gallup index as numeric score
              use = "pairwise.complete.obs"), 3), "\n")            # ignore missing pairs; 3 decimals

# (iii) Do respondents with and without WP88 differ on key variables?
d.fl %>%                                                           # FLI sample
  group_by(has_wp88) %>%                                           # split by WP88 fielded
  summarise(n = n(),                                               # group size
            mean_index = round(mean(as.numeric(as.character(FLI_3item))), 2),   # mean FLI score
            across(all_of(c("iwisescore", "hhsize", "gdp_pc_ppp", v.jmp)),     # key covariates
                   ~round(mean(.x, na.rm = TRUE), 2)),             # mean of each
            pct_female = round(100 * mean(female == "Female", na.rm = TRUE), 1),   # % female
            pct_urban  = round(100 * mean(urban_imp == "Peri-urban/urban",         # % urban
                                          na.rm = TRUE), 1),       # ignoring missing
            .groups = "drop")                                      # remove grouping

# (iv) Is WP88 non-administration a COUNTRY-level pattern?
d.fl %>%                                                           # FLI sample
  group_by(across(all_of(v.cluster))) %>%                          # one group per country-year
  summarise(n = n(), pct_no_wp88 = round(100 * mean(!has_wp88), 1),  # n and % without WP88
            .groups = "drop") %>%                                  # remove grouping
  arrange(desc(pct_no_wp88)) %>%                                   # highest % first
  print(n = Inf) #yes, if missing then 100% 

# 3c. ROBUSTNESS SAMPLE (created in 01 Section 10).
cat("\nRobustness sample (Gallup 4-item, WP88 countries only): N =",   # label
    nrow(d.fl.gallup), "| complete cases =", sum(d.fl.gallup$cc),      # total n and complete cases
    "| levels:", paste(levels(d.fl.gallup$INDEX_FL_ord), collapse = ", "), "\n")   # outcome levels


# =============================================================================
# CHECK 4: PREDICTOR DISTRIBUTIONS AND RARE CATEGORIES
# -----------------------------------------------------------------------------
# WHY: rare factor levels give imprecise coefficients; implausible values are
# usually data errors; heavy skew distorts estimates.
# =============================================================================

# 4a. Categorical predictors.
for (v in v.categorical) {                                         # loop over categorical variables
  cat("\n---", v, "(FLI sample) ---\n")                            # header with variable name
  print(table(d.fl[[v]], useNA = "ifany"))                         # counts per level, incl. NA
}                                                                  # end loop

# 4b. Continuous predictors.
d.fl %>%                                                           # FLI sample
  select(any_of(c("iwisescore", "iwise_mean", "hhsize", "gdp_pc_ppp",   # continuous variables...
                  "log_gdp", v.wgi, v.jmp, v.gbd))) %>%            # ...incl. WGI, JMP, DALY
  summarise(across(everything(),                                   # for every selected variable
                   list(min  = ~min(.x, na.rm = TRUE),             # minimum
                        med  = ~median(.x, na.rm = TRUE),          # median
                        mean = ~mean(.x, na.rm = TRUE),            # mean
                        max  = ~max(.x, na.rm = TRUE),             # maximum
                        sd   = ~sd(.x, na.rm = TRUE)))) %>%        # standard deviation
  pivot_longer(everything(), names_to = c("var", "stat"),          # reshape: one row per var-stat
               names_pattern = "^(.*)_(min|med|mean|max|sd)$") %>% # split name into var and stat
  pivot_wider(names_from = stat, values_from = value) %>%          # one row per var, stats as columns
  mutate(across(where(is.numeric), ~round(.x, 2))) %>%             # round to 2 decimals
  print(n = Inf)                                                   # print all rows
# LOOKING FOR: impossible values; mean far from median (= skew).

# 4c. Histograms.
d.fl %>%                                                           # FLI sample
  select(any_of(c("iwisescore", "hhsize", "gdp_pc_ppp", "log_gdp", v.jmp, v.gbd))) %>%   # variables to plot
  pivot_longer(everything()) %>%                                   # long format for faceting
  ggplot(aes(value)) +                                             # value on x
  geom_histogram(bins = 40, fill = "#3D4222") +                    # histogram, 40 bins
  facet_wrap(~name, scales = "free") +                             # one panel per variable
  labs(title = "Continuous predictors, FLI sample") +              # plot title
  theme_minimal()                                                  # clean theme

# 4c-bis. How much do the country-level covariates vary across COUNTRY-YEARS?

# bas_rate_pp: change in at-least-basic coverage (pp per year)
cy %>%                                                             # country-year data
  summarise(n_country_years = sum(!is.na(bas_rate_pp)),            # country-years with data
            median   = median(bas_rate_pp, na.rm = TRUE),          # median
            iqr_lo   = quantile(bas_rate_pp, .25, na.rm = TRUE),   # 25th percentile
            iqr_hi   = quantile(bas_rate_pp, .75, na.rm = TRUE),   # 75th percentile
            pct_flat = round(100 * mean(abs(bas_rate_pp) < 0.05, na.rm = TRUE), 1))   # % with ~no change

# gbd_cy: country-year DALY rate per 100,000 (checks only, not in models)
cy %>%                                                             # country-year data
  summarise(n_country_years = sum(!is.na(gbd_cy)),                 # country-years with data
            median = median(gbd_cy, na.rm = TRUE),                 # median
            iqr_lo = quantile(gbd_cy, .25, na.rm = TRUE),          # 25th percentile
            iqr_hi = quantile(gbd_cy, .75, na.rm = TRUE))          # 75th percentile

# 4d. iwisescore across the standard bands.
table(d.fl$iwise_cat, useNA = "ifany")                             # n per IWISE band, incl. NA


# =============================================================================
# CHECK 5: OUTCOME x CATEGORICAL PREDICTOR CROSS-TABS
# -----------------------------------------------------------------------------
# WHY: empty cells cause non-convergence and separation.
# =============================================================================

crosstab_check <- function(data, outcome, vars) {                  # define helper function
  for (v in vars) {                                                # loop over predictors
    cat("\n===", outcome, "x", v, "===\n")                         # header
    tb <- table(data[[v]], data[[outcome]])                        # predictor x outcome counts
    print(tb)                                                      # show table
    if (any(tb == 0))     cat(">>> WARNING: empty cell(s)\n")      # flag empty cells
    else if (any(tb < 5)) cat(">>> NOTE: cell(s) with fewer than 5 cases\n")   # flag sparse cells
    else                  cat(">>> ok\n")                          # otherwise fine
  }                                                                # end loop
}                                                                  # end of function

crosstab_check(d.fl,  "FLI_3item", v.categorical)                  # run for FLI
crosstab_check(d.inc, "INCOME_5",  v.categorical)                  # run for income
# IF EMPTY CELLS: collapse the offending category (in 01).


# =============================================================================
# CHECK 6: SAMPLE SIZE / EVENTS PER VARIABLE
# -----------------------------------------------------------------------------
# WHY: need enough cases in the SMALLEST outcome category. Rule of thumb: at
# least 10 per parameter. A k-level factor costs k-1 parameters.
# rhs and rhs_ctx are defined in 01 Section 8 (rhs_ctx uses GE only).
# =============================================================================

epv_check <- function(data, outcome, rhs) {                        # define helper function
  d.cc <- data %>% drop_na(all_of(c(outcome, all.vars(as.formula(paste("~", rhs))))))   # complete cases for this formula
  n.p  <- ncol(model.matrix(as.formula(paste("~", rhs)), data = d.cc))   # number of parameters (incl. dummies)
  tb   <- table(droplevels(d.cc[[outcome]]))                       # counts per outcome level
  cat(sprintf("\n%s\n  N = %d | parameters = %d | smallest category = %d | EPV = %.1f %s\n",   # output template
              outcome, nrow(d.cc), n.p, min(tb), min(tb) / n.p,    # N, parameters, smallest cat, EPV
              ifelse(min(tb) / n.p < 10, "<- BELOW 10", "")))      # flag if EPV < 10
}                                                                  # end of function

epv_check(d.fl,  "FLI_3item", rhs_ctx)                             # run for FLI, full model
epv_check(d.inc, "INCOME_5",  rhs_ctx)                             # run for income, full model


# =============================================================================
# CHECK 7: MULTICOLLINEARITY
# -----------------------------------------------------------------------------
# WHY: overlapping predictors give unstable coefficients and wide SEs.
# A property of the PREDICTORS ONLY, so it can be checked before fitting.
# Individual-level r for country-level variables is inflated by cluster size;
# report the country-year results (7c, 7e) for those.
# =============================================================================

d.cc.fl <- d.fl %>% filter(cc)                                 # FLI model sample

# 7a. Individual-level correlation matrices: Pearson AND Spearman
d.cor.fl <- d.cc.fl %>%                                        # FLI model sample
  mutate(age_gp_num = as.numeric(age_gp_profile)) %>%          # 1-4, so it can be correlated
  select(any_of(c("iwisescore", "age", "age_gp_num", "hhsize", "log_gdp",   # variables to correlate...
                  v.wgi, v.jmp, v.gbd, paste0(v.indices, "_num"))))         # ...incl. WGI, JMP, DALY, indices

for (m in c("pearson", "spearman")) {                          # run both methods
  cat("\n---", toupper(m), "---\n")                            # header with method name
  print(round(cor(d.cor.fl, use = "pairwise.complete.obs", method = m), 2))  # correlation matrix
}                                                              # end loop
# LOOKING FOR: |r| above ~0.8 (expected among the WGI).
# Large Pearson vs Spearman gaps = non-linearity or outliers.

# 7b. VIF on the main model (GE only). VIF ignores the outcome.
set.seed(2024)                                                 # reproducible random outcome
m.vifprobe <- lm(as.formula(paste("rnorm(nrow(d.fl)) ~", rhs_ctx)), data = d.fl)   # dummy model: random y, real predictors
car::vif(m.vifprobe)                                           # (G)VIF per predictor
# For FACTORS read GVIF^(1/(2*Df)), SQUARED, against thresholds of 5 or 10.
# Watch: iwise_mean, log_gdp, wgi_ge_sc, and the DALY rate vs age_gp_profile.

# 7c. Country-year correlation matrices (cy from 01 Section 6, incl. gbd_cy)
cy.num <- cy %>% select(-country_year)                         # numeric columns only
for (m in c("pearson", "spearman")) {                          # run both methods
  cat("\n--- COUNTRY-YEAR", toupper(m), "---\n")               # header with method name
  print(round(cor(cy.num, use = "pairwise.complete.obs", method = m), 2))    # correlation matrix
}                                                              # end loop
# Check iwise_mean against log_gdp, wgi_ge_sc, bas_rate_pp and gbd_cy.
# Earlier result: JMP levels 0.79-0.96 with each other -> max ONE per model.


# ---- 7d. DALY RATE vs AGE ---------------------------------------------------
# Within a country-year the DALY rate varies ONLY by age group, so it may
# overlap with age_gp_profile in the model.

# 7d-i. DALY rate by age group
d.cc.fl %>%                                                    # FLI model sample
  group_by(age_gp_profile) %>%                                 # one row per age group
  summarise(n      = n(),                                      # respondents per group
            median = median(.data[[v.gbd]]),                   # median DALY rate
            mean   = mean(.data[[v.gbd]]),                     # mean DALY rate
            min    = min(.data[[v.gbd]]),                      # lowest DALY rate
            max    = max(.data[[v.gbd]]),                      # highest DALY rate
            .groups = "drop") %>%                              # remove grouping
  mutate(across(where(is.numeric), ~round(.x, 1)))             # round to 1 decimal
# EXPECT: rates rise with age, steepest in 50+.

# 7d-ii. Correlation with age: overall and WITHIN country-year
d.cc.fl %>%                                                    # FLI model sample
  mutate(age_gp_num = as.numeric(age_gp_profile)) %>%          # 1-4
  group_by(across(all_of(v.cluster))) %>%                      # one group per country-year
  mutate(daly_within = .data[[v.gbd]] - mean(.data[[v.gbd]])) %>%  # remove country-year mean
  ungroup() %>%                                                # remove grouping
  summarise(                                                   # compute three correlations
    r_overall_agegp = cor(.data[[v.gbd]], age_gp_num, method = "spearman"),   # DALY vs age group, all data
    r_overall_age   = cor(.data[[v.gbd]], age, method = "spearman",           # DALY vs age in years
                          use = "complete.obs"),               # continuous age
    r_within_agegp  = cor(daly_within, age_gp_num, method = "spearman")       # DALY vs age group, within country-year
  ) %>%                                                        # end summarise
  mutate(across(everything(), ~round(.x, 2)))                  # round to 2 decimals
# r_within near 1 = within countries, the DALY rate is essentially an age proxy.

# 7d-iii. How much DALY variation is between country-years vs by age?
r2_cy  <- summary(lm(reformulate(paste0("factor(", v.cluster, ")"), v.gbd),   # DALY ~ country-year (as categories)
                     data = d.cc.fl))$r.squared                               # country-year only
r2_add <- summary(lm(reformulate(c(paste0("factor(", v.cluster, ")"),         # DALY ~ country-year (as categories)
                                   "age_gp_profile"), v.gbd),                 # + age group
                     data = d.cc.fl))$r.squared                               # + age group
# High r2_cy = DALY mainly varies BETWEEN countries.
# Large jump when adding age = age carries much of the variation.

cat(sprintf("\nDALY variance explained: country-year %.2f | + age group %.2f\n",
            r2_cy, r2_add))


# ---- 7e. GDP vs GOVERNMENT EFFECTIVENESS -----------------------------------
# Both are country-year level, so assess on cy.

cy %>%                                                         # country-year data
  summarise(pearson  = round(cor(log_gdp, .data[[v.wgi.main]],     # Pearson r
                                 use = "complete.obs"), 2),        # complete pairs only
            spearman = round(cor(log_gdp, .data[[v.wgi.main]],     # Spearman r
                                 method = "spearman", use = "complete.obs"), 2),   # complete pairs only
            n_cy     = sum(!is.na(log_gdp) & !is.na(.data[[v.wgi.main]])))         # country-years used
# |r| > 0.8 = keep one in the main model, the other as a sensitivity analysis.

cy %>%                                                         # country-year data
  ggplot(aes(log_gdp, .data[[v.wgi.main]])) +                  # GDP on x, GE on y
  geom_point() +                                               # one point per country-year
  geom_smooth(method = "lm", se = FALSE, colour = "#3D4222") + # linear trend line
  labs(title = "log GDP per capita vs Government Effectiveness",   # plot title
       x = "log GDP per capita (PPP)",                         # x-axis label
       y = "Government Effectiveness (0-100)") +               # y-axis label
  theme_minimal()                                              # clean theme

# =============================================================================
# CHECK 8: CLUSTERING / NON-INDEPENDENCE
# -----------------------------------------------------------------------------
# WHY: respondents are nested in countries, and log_gdp / wgi / bas_rate_pp
# take ONE value per country-year. Ignoring this makes SEs far too small.
# The ICC uses an EMPTY model, so it describes the data, not your model.
# DALY rate varies by age group within a country-year, unlike the other country-level variables.
# =============================================================================

# 8a. How many clusters, and how big?
d.fl %>%                                                           # FLI sample
  count(across(all_of(v.cluster))) %>%                             # respondents per country-year
  summarise(n_countries = n(), min_n = min(n),                     # n clusters, smallest cluster
            median_n = median(n), max_n = max(n))                  # median and largest cluster
# 74 countries, 500 to 3,503 respondents each.

d.fl %>% filter(cc) %>% count(across(all_of(v.cluster))) %>%       # complete cases per country-year
  summarise(n_countries_cc = n())                                  # n clusters in the model

# does Palestine 2025 drop out of the models?
# wgtnorm_using should be missing for 2025, so its rows should be cc = FALSE
# in both samples. If 2025 shows cc = TRUE, Palestine counts as TWO clusters.
for (d.name in c("d.fl", "d.inc")) {                               # loop over both samples
  cat("\n--- Palestine by country_year,", d.name, "---\n")         # header
  get(d.name) %>%                                                  # fetch the data frame by name
    filter(grepl("palestin", countrynew, ignore.case = TRUE)) %>%  # Palestine rows only
    count(country_year, wave, cc) %>%                              # n per country-year, wave, cc
    print()                                                        # show result
}                                                                  # end loop


# 8b. ICC: share of outcome variance BETWEEN countries.
icc_formula <- function(outcome) {                                 # define helper function
  as.formula(paste0("as.numeric(as.character(", outcome, ")) ~ 1 + (1 | ",   # outcome as number, intercept only
                    v.cluster, ")"))                               # random intercept per country-year
}                                                                  # end of function

m.icc.fl <- lme4::lmer(icc_formula("FLI_3item"), data = d.fl)      # empty model, FLI
cat("\nICC, FLI_3item:\n"); print(performance::icc(m.icc.fl))      # print ICC, FLI

m.icc.inc <- lme4::lmer(icc_formula("INCOME_5"), data = d.inc)     # empty model, income
cat("\nICC, INCOME_5:\n"); print(performance::icc(m.icc.inc))      # print ICC, income
# Above ~0.05 = clustering NOT ignorable.
# EXPECT INCOME_5 NEAR ZERO: quintiles are constructed WITHIN country.

# 8c. Visual: country means.
d.fl %>%                                                           # FLI sample
  group_by(across(all_of(v.cluster))) %>%                          # one group per country-year
  summarise(mean_fl = mean(as.numeric(as.character(FLI_3item)), na.rm = TRUE),   # mean FLI score
            .groups = "drop") %>%                                  # remove grouping
  ggplot(aes(x = reorder(.data[[v.cluster]], mean_fl), y = mean_fl)) +   # country-years sorted by mean
  geom_point() + coord_flip() +                                    # dot plot, horizontal
  labs(title = "Country mean, FLI_3item", x = NULL, y = "Mean") +  # labels
  theme_minimal(base_size = 7)                                     # small text for many countries

# >>> CHANGED: fix applied when fitting is now a random intercept per
# country-year, ordinal::clmm(... + (1 | country_year)), with iwisescore raw
# plus iwise_mean (was: polr() with CR2 cluster-robust SEs).
# <<< END CHANGED


# =============================================================================
# CHECK 9: SEPARATION
# -----------------------------------------------------------------------------
# WHY: a predictor that perfectly predicts the outcome gives huge coefficients
# and SEs. A DATA problem - more iterations never fix it.
# Each cut of the ordinal outcome is tested as a binary model.
# Uses the individual-level rhs; swap in rhs_ctx for the full model.
# =============================================================================

check_separation <- function(data, outcome, rhs) {                 # define helper function
  data <- data %>% drop_na(all_of(outcome))                        # drop rows missing the outcome
  data[[outcome]] <- droplevels(data[[outcome]])                   # remove empty outcome levels
  lv <- levels(data[[outcome]])                                    # outcome levels in order
  cat("\n--- Separation:", outcome, "(", length(lv) - 1, "cuts ) ---\n")   # header with n cuts
  for (k in seq_len(length(lv) - 1)) {                             # loop over each cut-point
    d.k <- data %>%                                                # build binary outcome
      mutate(y_bin = as.integer(as.numeric(.data[[outcome]]) > k)) # 1 = above cut k, 0 = at/below
    fit <- glm(as.formula(paste("y_bin ~", rhs)), data = d.k,      # logistic model for this cut
               family = binomial,                                  # binary outcome
               method = detectseparation::detect_separation)       # test for separation instead of fitting
    cat(sprintf("  cut above %-8s -> %s\n", lv[k],                 # print cut label
                ifelse(fit$outcome, "SEPARATION DETECTED", "ok"))) # and result
  }                                                                # end loop
}                                                                  # end of function

check_separation(d.fl,  "FLI_3item", rhs)                          # run for FLI
check_separation(d.inc, "INCOME_5",  rhs)                          # run for income
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