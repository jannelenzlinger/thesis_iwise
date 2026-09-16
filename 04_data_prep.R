# =============================================================================
# 01 DATA PREPARATION  |  Ordinal logistic regression
# Dataset: iwise_analysis
# Outcomes: FLI_3item (Financial Life Index, 3-item), INCOME_5 (Income Quintiles)
# Main predictor: iwisescore (Water Insecurity Experiences, 0-36)
#
# RULE: this is the ONLY script that changes data. Every recode, derived
# variable, sample restriction and weight lives here. 02_premodel_checks.R and
# later scripts only READ data_prepared.rds.
# If a check makes you change something, change it HERE and re-run this script.
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

d.iwise <- readRDS("data/iwise_analysis.rds")


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

v.continuous  <- c("iwisescore", "hhsize", "gdp_pc_ppp", v.wgi, v.jmp)
v.categorical <- c("female", "urban_imp", "age_gp_profile",
                   "maritalstatus", "employment", "education")

# Other Gallup indices. Missing for 30-73% of respondents (INDEX_YD 70-73%
# [excluded], INDEX_CA 40-47%, INDEX_PH 29-31%), so INDEX_PH and INDEX_CA are
# NOT predictors and are kept only for the correlation matrix in Check 7.
# INDEX_FS (Food & Shelter) is complete, so it IS a model predictor.
v.indices     <- c("INDEX_PH", "INDEX_FS", "INDEX_CA")
v.index.model <- "INDEX_FS"   # the one I want to use
v.cluster <- "iso3c"

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

    # GDP is right-skewed; log gives "per 1% change" interpretation.
    log_gdp = log(gdp_pc_ppp),

    # Standard IWISE bands (Young et al. 2019), in case linearity fails.
    iwise_cat = cut(iwisescore, breaks = c(-Inf, 2, 11, 23, Inf),
                    labels = c("No-to-marginal", "Low", "Moderate", "High"),
                    ordered_result = TRUE),

    # Numeric versions of the other indices, for Check 7 only.
    across(all_of(v.indices), ~as.numeric(as.character(.x)),
           .names = "{.col}_num")
  )

# Ethiopia income classification
d.iwise$country_income_group[d.iwise$country_name == "Ethiopia"] <- "Low income"


# ---- 6. COUNTRY-YEAR DATASET ------------------------------------------------
# JMP and WGI vary only between country-years, so their correlation with IWISE
# is assessed at that level (Checks 4c-bis and 7c). wgt2 is used because these
# are within-country means (READ_ME: country-specific analysis).

cy <- d.iwise %>%
  group_by(country_year) %>%
  summarise(
    iwise_mean = weighted.mean(iwisescore,  wgt2, na.rm = TRUE),
    # Proportion (0-1) with moderate-to-high WI: observed iwisescore >= 12.
    # Built from iwisescore, NOT iwise12_imp, so it matches the model predictor.
    iwise_mh   = weighted.mean(iwisescore >= 12, wgt2, na.rm = TRUE),
    across(all_of(c(v.jmp.all, "log_gdp")), first),
    .groups = "drop"
  )


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

v.model.fl  <- c("FLI_3item", v.main, "hhsize", v.index.model, "log_gdp",
                 v.wgi.main, v.jmp, v.categorical, v.cluster, v.weight)
v.model.inc <- c("INCOME_5",  v.main, "hhsize", v.index.model, "log_gdp",
                 v.wgi.main, v.jmp, v.categorical, v.cluster, v.weight)

rhs     <- paste("iwisescore + female + urban_imp + age_gp_profile +",
                 "maritalstatus + hhsize + employment + education +",
                 v.index.model) #right hand side of the model

# Contextual version (with country level variables). ONE governance term only.
rhs_ctx <- paste(rhs, "+ log_gdp +", v.wgi.main, "+", v.jmp)


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
             v.indices     = v.indices,
             v.index.model = v.index.model,
             v.cluster     = v.cluster,
             v.fli.items   = v.fli.items,
             v.model.fl    = v.model.fl,
             v.model.inc   = v.model.inc),
        "data_prepared.rds")

cat("\nSaved data_prepared.rds\n")