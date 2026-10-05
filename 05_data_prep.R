# =============================================================================
# 01 DATA PREPARATION  |  Ordinal logistic regression
# Dataset: iwise_all_contextual
# Outcomes: FLI_3item (Financial Life Index, 3-item), INCOME_5 (Income Quintiles)
# Main predictor: iwisescore (Water Insecurity Experiences, 0-36)
#
# RULE: this is the ONLY script that changes data. Every recode post-merge, derived
# variable, sample restriction and weight lives here. 06_premodel_checks.R and
# later scripts only READ data_prepared.rds.
#
# STRUCTURE
# ---------
#   0   Packages
#   1   Load data
#   2   Financial Life Index recalculation (recode, validate, build 3-item)
#   3   Variable groups
#   4   Variable types and level ordering
#   5   Derived variables
#   6   Country-year dataset (for Checks 4c-bis and 7c)
#   7   Split into analytic samples
#   8   Per-model variable lists and model formulas
#   9   Complete-case flags and model weights
#   10  Robustness sample (Gallup 4-item)
#   11  Save
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------

pkgs <- c("tidyverse")

missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing_pkgs)) install.packages(missing_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

select <- dplyr::select   # MASS and others mask dplyr::select; be explicit


# ---- 1. LOAD DATA -----------------------------------------------------------
# Load from FILE, not from an object already in memory, so the script runs
# standalone in a fresh session.

d.iwise <- readRDS("data/iwise_all_contextual.rds")


# ---- 2. FINANCIAL LIFE INDEX RECALCULATION ----------------------------------

# ---- 2.1 RECODE EACH ITEM TO 0/1 --------------------------------------------
# The positive answer is code 1 for ALL FOUR items. (Do not confuse these with
# the food/shelter items WP40/WP43, where the positive answer is "No" = 2.)

d.iwise <- d.iwise %>%
  mutate(
    # WP2319: 1 = "living comfortably on present income".
    # if_else() returns NA when the input is NA, so codes 2-6 become 0 while a
    # genuinely unanswered item stays NA.
    s_2319 = if_else(WP2319 == 1, 1, 0),

    # WP30: 1 = "Satisfied" with standard of living.
    s_30 = if_else(WP30 == 1, 1, 0),

    # WP31: 1 = standard of living "getting better". "The same" scores 0.
    s_31 = if_else(WP31 == 1, 1, 0),

    # WP88: 1 = local economic conditions "getting better".
    # Needed only for the validation in 2.2.
    s_88 = if_else(WP88 == 1, 1, 0)
  )

# Sanity check: each should be 0/1 with NAs only where the item was unanswered.
d.iwise %>%
  summarise(across(c(s_2319, s_30, s_31, s_88),
                   list(n_1  = ~sum(.x == 1, na.rm = TRUE),
                        n_0  = ~sum(.x == 0, na.rm = TRUE),
                        n_NA = ~sum(is.na(.x)))))

table(raw = d.iwise$WP2319, recoded = d.iwise$s_2319, useNA = "ifany")
table(raw = d.iwise$WP30,   recoded = d.iwise$s_30,   useNA = "ifany")
table(raw = d.iwise$WP31,   recoded = d.iwise$s_31,   useNA = "ifany")


# ---- 2.2 VALIDATE BY REPRODUCING GALLUP'S OWN 4-ITEM INDEX ------------------
# Do NOT skip this. If the recoding is wrong, the 3-item version inherits the
# error silently.

d.iwise <- d.iwise %>%
  mutate(
    # Number of the three non-WP2319 items answered.
    n_other = rowSums(!is.na(across(c(s_30, s_31, s_88)))),

    # Mean of whichever of the three are available. The varying denominator is
    # the source of Gallup's 25/75 levels.
    mean_other = rowMeans(across(c(s_30, s_31, s_88)), na.rm = TRUE),

    # Gallup's rule: WP2319 answered AND at least 2 of the other 3.
    FLI_check = if_else(
      !is.na(s_2319) & n_other >= 2,
      100 * (s_2319 + mean_other) / 2,
      NA_real_
    )
  )

# Agreement with Gallup's variable. Should be ~1.00.
d.iwise %>%
  filter(!is.na(INDEX_FL_ord)) %>%
  summarise(
    n = n(),
    agree = mean(abs(FLI_check -
                       as.numeric(as.character(INDEX_FL_ord))) < 0.01,
                 na.rm = TRUE)
  ) # agree = 1


# ---- 2.3 BUILD THE 3-ITEM INDEX ---------------------------------------------
# Same recipe, WP88 removed. Fixed denominator of 2 for everyone.

d.iwise <- d.iwise %>%
  mutate(
    # Eligibility: WP2319 plus BOTH remaining items.
    # Possible results: 0, 25, 50, 75, 100.
    FLI_3item_num = if_else(
      !is.na(s_2319) & !is.na(s_30) & !is.na(s_31),
      100 * (s_2319 + (s_30 + s_31) / 2) / 2,
      NA_real_
    ),

    # Ordered factor built from a NUMERIC vector, so levels sort numerically
    # (0 < 25 < 50 < 75 < 100) rather than alphabetically.
    FLI_3item = factor(FLI_3item_num, ordered = TRUE)
  )


# ---- 3. VARIABLE GROUPS -----------------------------------------------------

v.outcomes <- c("FLI_3item", "INCOME_5")
v.main     <- "iwisescore"      # main predictor, 0-36
v.ctrl     <- "iwise4"          # alternative 4-item version, 0-12

# The WGI come in two suffix families, one row per country-year:
#   _est = governance estimate, standard normal units, approx. -2.5 to 2.5
#   _sc  = the same estimate mapped onto a 0-100 scale
# NEVER use both for the same dimension in one model (perfect collinearity).
v.wgi.est  <- c("wgi_va_est", "wgi_pv_est", "wgi_ge_est",
                "wgi_rq_est", "wgi_rl_est", "wgi_cc_est")
v.wgi.sc   <- c("wgi_va_sc",  "wgi_pv_sc",  "wgi_ge_sc",
                "wgi_rq_sc",  "wgi_rl_sc",  "wgi_cc_sc")
v.wgi      <- v.wgi.sc     # all six, kept for Check 7 and sensitivity analyses
v.wgi.main <- "wgi_ge_sc"  # Government Effectiveness: the main model

# Pooled regression weight (READ_ME). Rescaled to w_fit in Section 9.
v.weight <- "wgtnorm_using"

# JMP covariate used in the models. Country-year level.
# Change of at-least-basic water service, percentage points per year.
v.jmp <- "bas_rate_pp"

# All JMP variables, for the country-year correlation matrix only (Check 7c).
v.jmp.all <- c("bas_rate_pp", "bas_level", "sm_level",
               "prem_level", "avail_level", "qual_level")


# GBD covariate. Country-year level.
# Enteric disease DALYs, matched, by age group(GBD 2023;
# 2024/2025 surveys matched to 2023, see gbd_carried_fwd).
v.gbd <- "daly_enteric_rate"   # RAW rate: descriptives and checks only

# >>> ADDED: model version of the DALY rate (natural log, created in Section 5).
# v.gbd stays the raw rate, because Section 4 converts v.continuous to numeric
# before Section 5 creates this variable.
v.gbd.model <- "daly_enteric_model"
# <<< END ADDED


v.continuous  <- c("iwisescore", "hhsize", "gdp_pc_ppp", v.wgi, v.jmp, v.gbd)
v.categorical <- c("female", "urban_imp", "age_gp_profile",
                   "maritalstatus", "employment", "education")

# Other Gallup indices. Missing for 30-73% of respondents (INDEX_YD 70-73%
# [excluded], INDEX_CA 40-47%, INDEX_PH 29-31%), so INDEX_PH and INDEX_CA are
# NOT predictors and are kept only for the correlation matrix in Check 7.
# INDEX_FS (Food & Shelter) is complete, so it IS a model predictor.
v.indices     <- c("INDEX_PH", "INDEX_FS", "INDEX_CA")
v.index.model <- "INDEX_FS"   # the one I want to use
# cluster = country_year (was iso3c), so the random intercept and
# the country mean (Section 5.1) are defined on the same unit. Palestine 2025
# drops out of the models via wgtnorm_using.
v.cluster <- "country_year"

v.fli.items <- c("WP2319", "WP30", "WP31", "WP88")


# ---- 4. VARIABLE TYPES AND LEVEL ORDERING -----------------------------------
# The FIRST level of each factor is the reference category. Idempotent: safe
# to re-run.

d.iwise <- d.iwise %>%
  mutate(
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


# ---- 5. DERIVED VARIABLES ---------------------------------------------------
# Create every derived variable on d.iwise BEFORE splitting. Anything mutated
# afterwards will not propagate into d.fl / d.inc.
#
# OUTCOME NOTE: FLI_3item (WP2319, WP30, WP31) is used from here on. WP88 was
# not administered to 38,704 respondents, which put Gallup's index on two
# different scales. FLI_3item has a fixed denominator and 5 levels, at no loss
# of sample. Gallup's INDEX_FL_ord is kept for the robustness sample (Section 10).

d.iwise <- d.iwise %>%
  mutate(
    # Flags WHO had WP88 fielded (sample-composition checks in Check 3b).
    n_fli_answered = rowSums(!is.na(across(all_of(v.fli.items)))),
    has_wp88       = !is.na(WP88),

    # >>> CHANGED (comment only): corrected interpretation
    # GDP is right-skewed; OR per 1-unit log = per 2.72-fold increase in GDP.
    # <<< END CHANGED
    log_gdp = log(gdp_pc_ppp),

    # >>> ADDED: DALY rate is right-skewed (mean 592, SD 699) and on a much
    # larger scale than the other predictors, which caused the clmm()
    # "very large eigenvalue" warning. Natural log, as for GDP.
    # OR per 1-unit log = per 2.72-fold increase in the DALY rate.
    daly_enteric_model = log(daly_enteric_rate),
    # <<< END ADDED

    # Standard IWISE bands (Young et al. 2019), in case linearity fails.
    iwise_cat = cut(iwisescore, breaks = c(-Inf, 2, 11, 23, Inf),
                    labels = c("No-to-marginal", "Low", "Moderate", "High"),
                    ordered_result = TRUE),

    # Numeric versions of the other indices, for Check 7 only.
    across(all_of(v.indices), ~as.numeric(as.character(.x)),
           .names = "{.col}_num")
  )

# >>> ADDED: the log is undefined for 0 or negative values, so check first.
stopifnot(min(d.iwise$daly_enteric_rate, na.rm = TRUE) > 0)
summary(d.iwise$daly_enteric_model)   # check: plausible range, no -Inf
# <<< END ADDED

# ---- 5.1 COUNTRY MEAN IWISE (within/between split, Mundlak) -----------------
# Model uses RAW iwisescore + country mean:
#   iwisescore coefficient = within-country effect
#   iwise_mean coefficient = contextual effect (between minus within)
#   between-country effect = sum of the two
# Computed here, BEFORE the split, so the mean uses every respondent with an
# IWISE score, not only each model's complete cases.
# Weight = wgt2: a within-country quantity (READ_ME: country-specific analysis,
# correct India weights). Same weight as the country-year data in Section 6.
 
d.iwise <- d.iwise %>%
  group_by(country_year) %>%
  mutate(
    .ok        = !is.na(iwisescore) & !is.na(wgt2),
    iwise_mean = if (any(.ok)) weighted.mean(iwisescore[.ok],
                                             wgt2[.ok])
                 else NA_real_,
    .ok        = NULL
  ) %>%
  ungroup()

 
# Ethiopia income classification
d.iwise$country_income_group[d.iwise$country_name == "Ethiopia"] <- "Low income"
 
 
# ---- 6. COUNTRY-YEAR DATASET ------------------------------------------------
# One row per country-year, for country-level checks (02: Checks 4c-bis, 7c).
# gbd_cy = average DALY rate across all adults in the country-year
# (wgt2-weighted mean of the age-group rates). Checks only, NOT a model variable.

cy <- d.iwise %>%
  group_by(country_year) %>%                                   # one group per country-year
  summarise(
    # Proportion (0-1) with moderate-to-high WI (observed iwisescore >= 12)
    iwise_mh = weighted.mean(iwisescore >= 12, wgt2, na.rm = TRUE),
    # Country-year DALY rate
    gbd_cy = {
      ok <- !is.na(.data[[v.gbd]]) & !is.na(wgt2)             # rows with rate AND weight
      if (any(ok)) weighted.mean(.data[[v.gbd]][ok], wgt2[ok]) else NA_real_
    },
    # Country-level variables: identical within a country-year, so take the first
    across(all_of(c("iwise_mean", v.jmp.all, "log_gdp", v.wgi)), first),
    .groups = "drop"                                           # remove grouping
  )

# Checks: one row per country-year, plausible DALY values
stopifnot(nrow(cy) == n_distinct(d.iwise$country_year))
summary(cy$gbd_cy)
 
 
# ---- 7. SPLIT INTO ANALYTIC SAMPLES -----------------------------------------
# Each outcome gets its own frame, so a respondent missing INCOME_5 still
# counts towards the FLI model, and vice versa.
# droplevels() stops polr() estimating cut-points for empty categories.
 
d.fl <- d.iwise %>%
  drop_na(FLI_3item) %>%
  mutate(FLI_3item = droplevels(FLI_3item))
 
d.inc <- d.iwise %>%
  drop_na(INCOME_5) %>%
  mutate(INCOME_5 = droplevels(INCOME_5))
 
 
# ---- 8. PER-MODEL VARIABLE LISTS AND FORMULAS -------------------------------
# These define the complete-case sample for each MAIN model.
# The FLI items are deliberately NOT included: FLI_3item already encodes its
# eligibility rule. Only the main WGI (GE) is included; sensitivity models
# with other WGIs define their own complete-case sample.
 
# >>> CHANGED: v.gbd -> v.gbd.model (log DALY rate) in both lists and rhs_ctx
v.model.fl  <- c("FLI_3item", v.main, "iwise_mean", "hhsize", v.index.model, "log_gdp",
                 v.wgi.main, v.jmp, v.gbd.model, v.categorical, v.cluster, v.weight)
v.model.inc <- c("INCOME_5",  v.main, "iwise_mean", "hhsize", v.index.model, "log_gdp",
                 v.wgi.main, v.jmp, v.gbd.model, v.categorical, v.cluster, v.weight)
# <<< END CHANGED
 
rhs     <- paste("iwisescore + female + urban_imp + age_gp_profile +",
                 "maritalstatus + hhsize + employment + education +",
                 v.index.model) #right hand side of the model
 
# Contextual version (with country level variables). ONE governance term only.
# >>> CHANGED: v.gbd -> v.gbd.model
rhs_ctx <- paste(rhs, "+ iwise_mean + log_gdp +", v.wgi.main, "+", v.jmp, "+", v.gbd.model)
# <<< END CHANGED
 
 
# ---- 9. COMPLETE-CASE FLAGS AND MODEL WEIGHTS -------------------------------
# cc = row has every variable the main model needs (including the weight).
#
# w_fit: polr() treats weights as frequency counts, so the weight is rescaled
# to sum to the COMPLETE-CASE N of each model. Rows are NOT dropped here, so
# Check 2 can still compare kept and dropped cases. w_fit is NA where cc is
# FALSE.
 
add_cc_weight <- function(data, v.model) {
  data %>%
    mutate(
      cc    = complete.cases(select(., any_of(v.model))),
      w_fit = if_else(cc, .data[[v.weight]] / mean(.data[[v.weight]][cc]),
                      NA_real_)
    )
}
 
d.fl  <- add_cc_weight(d.fl,  v.model.fl)
d.inc <- add_cc_weight(d.inc, v.model.inc)
 
 
# ---- 10. ROBUSTNESS SAMPLE --------------------------------------------------
# Gallup's own 4-item index, only where WP88 was fielded. 7 levels.
# Gets its own cc flag and weight, since its outcome differs.
 
d.fl.gallup <- d.fl %>%
  filter(has_wp88) %>%
  mutate(INDEX_FL_ord = droplevels(INDEX_FL_ord)) %>%
  add_cc_weight(replace(v.model.fl, v.model.fl == "FLI_3item", "INDEX_FL_ord"))
 
 
# ---- 11. SAVE ---------------------------------------------------------------
# Frames AND variable definitions saved together, so they stay in sync.
 
saveRDS(list(d.iwise       = d.iwise,
             d.fl          = d.fl,
             d.inc         = d.inc,
             d.fl.gallup   = d.fl.gallup,
             cy            = cy,
             rhs           = rhs,
             rhs_ctx       = rhs_ctx,
             v.outcomes    = v.outcomes,
             v.main        = v.main,
             v.ctrl        = v.ctrl,
             v.categorical = v.categorical,
             v.continuous  = v.continuous,
             v.wgi         = v.wgi,
             v.wgi.est     = v.wgi.est,
             v.wgi.sc      = v.wgi.sc,
             v.wgi.main    = v.wgi.main,
             v.weight      = v.weight,
             v.jmp         = v.jmp,
             v.jmp.all     = v.jmp.all,
             v.gbd         = v.gbd,
             v.gbd.model   = v.gbd.model,   # >>> ADDED
             v.indices     = v.indices,
             v.index.model = v.index.model,
             v.cluster     = v.cluster,
             v.fli.items   = v.fli.items,
             v.model.fl    = v.model.fl,
             v.model.inc   = v.model.inc),
        "data/iwise_data_prepared.rds")
 
file.exists("data/iwise_data_prepared.rds")
