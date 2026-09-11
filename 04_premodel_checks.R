# =============================================================================
# PRE-MODEL CHECKS  |  Ordinal logistic regression
# Dataset: iwise_analysis
# Outcomes: INDEX_FL_ord (Financial Life Index), INCOME_5 (Income Quintiles)
# Main predictor: iwisescore (Water Insecurity Experiences, 0-36)
# JMP covariate: bas_rate_pp (recent change in at-least-basic water, pp/year)
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Nine checks, in order, on the DATA. It STOPS BEFORE THE MODEL IS FITTED.
# The four checks that require a fitted model are listed at the end.
# See premodel_checks_explained.md for the reasoning behind each one.
#
# CHANGES FOR JMP: every line added or edited for bas_rate_pp is tagged # [JMP]
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


# Load from FILE, not from an object already in memory, so the script runs
# standalone in a fresh session tomorrow. Use whichever matches your file.
d.iwise <- readRDS("data/iwise_analysis.rds")

# --- SETUP: Financial Life Index Recalculation --------------------------------------
# ---- STEP 1: RECODE EACH ITEM TO 0/1 ----------------------------------------
# The positive answer is code 1 for ALL FOUR items. (Do not confuse these with
# the food/shelter items WP40/WP43, where the positive answer is "No" = 2.)

d.iwise <- d.iwise %>%
  mutate(

    # WP2319: 1 = "living comfortably on present income".
    # if_else() returns 1 when the condition is TRUE, 0 when FALSE, and - this
    # is the important part - NA when the input is NA. So codes 2-6 (getting
    # by, difficult, very difficult, DK, refused) all become 0, exactly as
    # Gallup specifies, while a genuinely unanswered item stays NA and is
    # excluded from the calculation.
    s_2319 = if_else(WP2319 == 1, 1, 0),

    # WP30: 1 = "Satisfied" with standard of living.
    # Codes 2 (dissatisfied), 3 (DK), 4 (refused) -> 0.
    s_30 = if_else(WP30 == 1, 1, 0),

    # WP31: 1 = standard of living "getting better".
    # Codes 2 (the same), 3 (getting worse), 4 (DK), 5 (refused) -> 0.
    # NOTE "the same" scores 0 - only improvement counts as positive.
    s_31 = if_else(WP31 == 1, 1, 0),

    # WP88: 1 = local economic conditions "getting better".
    # Same structure as WP31. Needed only for the validation in Step 2.
    s_88 = if_else(WP88 == 1, 1, 0)
  )

# Sanity check: each should be 0/1 with NAs only where the item was unanswered.
d.iwise %>%
  summarise(across(c(s_2319, s_30, s_31, s_88),
                   list(n_1 = ~sum(.x == 1, na.rm = TRUE),
                        n_0 = ~sum(.x == 0, na.rm = TRUE),
                        n_NA = ~sum(is.na(.x)))))


# ---- STEP 2: VALIDATE BY REPRODUCING GALLUP'S OWN 4-ITEM INDEX --------------
# Do NOT skip this. If your recoding is wrong, the 3-item version inherits the
# error silently. Reconstruct Gallup's index and check it matches theirs.

d.iwise <- d.iwise %>%
  mutate(

    # How many of the three non-WP2319 items does this respondent have?
    # across() selects the three columns; !is.na() gives TRUE where answered;
    # rowSums() counts the TRUEs per respondent (TRUE counts as 1).
    n_other = rowSums(!is.na(across(c(s_30, s_31, s_88)))),

    # The second component: the MEAN of whichever of the three are available.
    # na.rm = TRUE is what makes the denominator vary - with all three it
    # divides by 3, with two it divides by 2. That varying denominator is the
    # entire source of the 25/75 levels.
    mean_other = rowMeans(across(c(s_30, s_31, s_88)), na.rm = TRUE),

    # Gallup's eligibility rule: WP2319 answered AND at least 2 of the other 3.
    FLI_check = if_else(
      !is.na(s_2319) & n_other >= 2,
      100 * (s_2319 + mean_other) / 2,   # average the two components, x100
      NA_real_                            # otherwise no index
    )
  )

# Agreement with Gallup's variable. Should be ~1.00.
# The < 0.01 tolerance allows for floating-point rounding (33.33 vs 33.333...).
d.iwise %>%
  filter(!is.na(INDEX_FL_ord)) %>%
  summarise(
    n = n(),
    agree = mean(abs(FLI_check -
                     as.numeric(as.character(INDEX_FL_ord))) < 0.01,
                 na.rm = TRUE)
  ) #agree = 1


# ---- STEP 3: BUILD THE 3-ITEM INDEX -----------------------------------------
# Same recipe, WP88 removed. The second component is now the mean of WP30 and
# WP31 only - a FIXED denominator of 2 for everyone, which is the whole point.

d.iwise <- d.iwise %>%
  mutate(

    # Eligibility: WP2319 plus BOTH remaining items. Gallup's "at least 2 of 3"
    # becomes "2 of 2" here. Requiring both keeps the second component from
    # resting on a single question, and costs nothing - every respondent who
    # currently has an index already has WP30 and WP31.
    FLI_3item_num = if_else(
      !is.na(s_2319) & !is.na(s_30) & !is.na(s_31),

      # Component 1 = s_2319, which is 0 or 1.
      # Component 2 = (s_30 + s_31) / 2, which is 0, 0.5, or 1.
      # Average them and multiply by 100.
      # Possible results: 0, 25, 50, 75, 100 - five levels, and every
      # respondent's score is built the same way.
      100 * (s_2319 + (s_30 + s_31) / 2) / 2,

      NA_real_
    ),

    # Ordered factor, so polr()/clm() treat it as ordinal.
    # as.numeric() first so levels sort 0 < 25 < 50 < 75 < 100 rather than
    # alphabetically ("0","100","25","50","75") - the same trap as Check 1.
    FLI_3item = factor(FLI_3item_num, ordered = TRUE)
  )



# ---- 1. LOAD DATA AND DEFINE VARIABLE GROUPS --------------------------------
# Work on a copy so the raw object stays clean and you can always restart.

v.outcomes    <- c("FLI_3item", "INCOME_5")
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

v.wgi <- v.wgi.sc   # preference

# [JMP] JMP covariate used in the models. Country-year level, like log_gdp and
# the WGI. Change of at-least-basic water service, percentage points per year.
v.jmp <- "bas_rate_pp"

v.continuous  <- c("iwisescore", "hhsize", "gdp_pc_ppp", v.wgi, v.jmp)  # [JMP]
v.categorical <- c("female", "urban_imp", "age_gp_profile",
                   "maritalstatus", "employment", "education")
# The other Gallup indices are currently NOT used as predictors. They are missing
# for 30-73% of respondents (INDEX_YD 70-73% [excluded], INDEX_CA 40-47%, INDEX_PH 29-31%;
# INDEX_FS is complete), so including them would cut the analytic sample to
# roughly a quarter. Kept here only for the correlation matrix in Check 7.
v.indices     <- c("INDEX_PH", "INDEX_FS", "INDEX_CA")
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
class(d.iwise$FLI_3item)      # good
class(d.iwise$INCOME_5)
levels(d.iwise$FLI_3item)     # good
levels(d.iwise$INCOME_5)

# 1b. Recode. Safe to run even if already correct - it is idempotent.
#     as.character() -> the labels; as.numeric() -> real numbers;
#     factor() on a NUMERIC vector sorts levels numerically, which removes the
#     alphabetical hazard entirely.
d.iwise <- d.iwise %>%
  mutate(
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

    across(all_of(v.continuous), as.numeric)   # [JMP] now includes bas_rate_pp
  )

# 1c. Confirm the recode did what you expect.
str(d.iwise %>% select(any_of(c(v.outcomes, v.continuous, v.categorical))))


# ---- SETUP: DERIVED VARIABLES FIRST -----------------------------------------
# IMPORTANT: create every derived variable on d.iwise BEFORE splitting.
# Anything you mutate afterwards will not propagate into d.fl / d.inc.
#
# OUTCOME NOTE: the FLI used from here on is FLI_3item, your own index built
# from WP2319, WP30 and WP31 (recoding above, validated at agree = 1.00 against
# Gallup's version). WP88 is excluded because it was not administered to 38,704
# respondents, which forced Gallup's index onto two different scales - a
# three-item denominator for some respondents and a two-item denominator for
# others. That is what produced levels 25 and 75.
#
# FLI_3item has a FIXED denominator, so every respondent is measured the same
# way. It has 5 levels (0, 25, 50, 75, 100) instead of 9, at no loss of sample.
# Gallup's original INDEX_FL_ord is retained for the robustness check in 3c.

v.fli.items <- c("WP2319", "WP30", "WP31", "WP88")

d.iwise <- d.iwise %>%
  mutate(
    # How many FLI items did each respondent answer? With FLI_3item this no
    # longer affects the outcome scale - it only flags WHO had WP88 fielded,
    # which matters for the sample-composition checks in 3b.
    n_fli_answered = rowSums(!is.na(across(all_of(v.fli.items)))),
    has_wp88       = !is.na(WP88),

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

# not forget ethiopia income classification
d.iwise$country_income_group[d.iwise$country_name == "Ethiopia"] <- "Low income"



# ---- COUNTRY-YEAR COLLINEARITY: IWISE vs JMP -------------------------------
# JMP varies only between country-years, so the relevant correlation is at
# that level. At respondent level it would look artificially small.
# [JMP] Renamed from v.jmp to v.jmp.all, so it no longer overwrites the model
# covariate list defined in Section 1.

v.jmp.all <- c("bas_rate_pp", "bas_level", "sm_level",
               "prem_level", "avail_level", "qual_level")

cy <- d.iwise %>%
  group_by(country_year) %>%
  summarise(
    iwise_mean = weighted.mean(iwisescore,  wgt2, na.rm = TRUE),
    iwise_mh   = weighted.mean(iwise12_imp, wgt2, na.rm = TRUE),  # % mod-to-high WI
    across(all_of(c(v.jmp.all, "log_gdp")), first)
  )

cy %>%
  select(-country_year) %>%
  cor(use = "pairwise.complete.obs", method = "spearman") %>%
  round(2)
# RESULT: bas_rate_pp x IWISE = 0.14. JMP levels -0.53 to -0.69 with IWISE and
# 0.79-0.96 with each other, so use at most ONE JMP level in any model.





# ---- SPLIT INTO TWO ANALYTIC SAMPLES ----------------------------------------
# Each outcome gets its own frame. This is the point: a respondent missing
# INCOME_5 should still count towards the FLI model, and vice versa.

# FLI sample: everyone with a computed FLI_3item score.
d.fl <- d.iwise %>%
  drop_na(FLI_3item) %>%
  mutate(FLI_3item = droplevels(FLI_3item))
# droplevels() is a safety net - no level should be empty here. It matters if
# you ever restrict the sample further, since polr() would otherwise try to
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

cat("\nFLI_3item levels:\n"); print(levels(d.fl$FLI_3item))
cat("Income outcome levels:\n"); print(levels(d.inc$INCOME_5))

# Confirm the fix worked: FLI_3item should cover exactly the same respondents
# as Gallup's index, with no one gained or lost.
cat("\nGallup index vs FLI_3item coverage:\n")
print(table(gallup = !is.na(d.iwise$INDEX_FL_ord),
            fli3   = !is.na(d.iwise$FLI_3item))) #correct


# ---- VARIABLE LISTS PER MODEL -----------------------------------------------
# The four FLI items are deliberately NOT in these lists - FLI_3item already
# encodes its own eligibility rule, and adding the items would drop everyone
# without WP88 for no reason.
# [JMP] v.jmp added, so every missing-data check in Check 2 now counts it.

v.model.fl  <- c("FLI_3item", v.main, "hhsize", "log_gdp", v.wgi, v.jmp,
                 v.categorical, v.cluster)
v.model.inc <- c("INCOME_5", v.main, "hhsize", "log_gdp", v.wgi, v.jmp,
                 v.categorical, v.cluster)


# =============================================================================
# CHECK 2: MISSING DATA
# -----------------------------------------------------------------------------
# WHY: the model deletes any row missing ANY variable in the formula. That is
# harmless only if the deleted rows are a random subset. If they differ
# systematically from the rows you keep, every coefficient is biased.
#
# NOTE: these figures describe the ACTUAL analytic samples, so these are the
# numbers to report in the paper.
#
# [JMP] bas_rate_pp is country-level, so where it is missing it is missing for
# a WHOLE country. Expect N to fall versus the pre-JMP run; 2e/2f show where.
# =============================================================================

# 2a. Missingness per variable, within each analytic sample.
cat("\n=== FLI SAMPLE ===\n")
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
    summarise(model = "FLI_3item", n = n(),
              across(all_of(c("iwisescore", "hhsize", "gdp_pc_ppp", v.jmp)),  # [JMP]
                     ~round(mean(.x, na.rm = TRUE), 2)), .groups = "drop"),
  d.inc %>% group_by(complete = cc) %>%
    summarise(model = "INCOME_5", n = n(),
              across(all_of(c("iwisescore", "hhsize", "gdp_pc_ppp", v.jmp)),  # [JMP]
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

# [JMP] 2e-bis. Which countries have NO bas_rate_pp at all? These drop out of
#     both models entirely once the JMP covariate is included.
cat("\n--- Countries with no bas_rate_pp ---\n")
d.iwise %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(pct_missing_jmp = round(100 * mean(is.na(.data[[v.jmp]])), 1),
            .groups = "drop") %>%
  filter(pct_missing_jmp > 0) %>%
  print(n = Inf)

# 2f. For countries losing most of their cases, WHICH variable is responsible?
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
print(setdiff(unique(d.iwise[[v.cluster]]), unique(d.fl[[v.cluster]])))
cat("--- Countries absent from the income sample entirely ---\n")
print(setdiff(unique(d.iwise[[v.cluster]]), unique(d.inc[[v.cluster]])))

# For those, check what is missing in the FULL data:
d.iwise %>%
  filter(.data[[v.cluster]] %in%
           setdiff(unique(d.iwise[[v.cluster]]), unique(d.fl[[v.cluster]]))) %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(across(any_of(c(v.model.fl, "WP2319", "WP30", "WP31")),
                   ~round(100 * mean(is.na(.x)), 1)),
            .groups = "drop") %>%
  print(width = Inf)


# =============================================================================
# CHECK 3: OUTCOME DISTRIBUTION AND SPARSE CATEGORIES
# -----------------------------------------------------------------------------
# WHY: one cut-point is estimated per threshold - FOUR for the 5-level
# FLI_3item, four for INCOME_5 - each from respondents near that threshold. A
# level with very few cases gives an unstable cut-point.
#
# JUDGE ON COUNTS, NOT PERCENTAGES. With N in the tens of thousands, a level
# holding 3% of the sample still has well over a thousand respondents, which is
# ample. There is no canonical minimum percentage in the literature.
#
# AND NOTE: because the index is a deterministic function of three items, a
# level's frequency reflects how many response patterns map onto it - not
# sampling sparsity. A structurally rare score is not a data problem.
# =============================================================================

cat("\n--- FLI_3item ---\n")
print(table(d.fl$FLI_3item, useNA = "ifany"))
print(round(100 * prop.table(table(d.fl$FLI_3item)), 1))

cat("\n--- INCOME_5 ---\n")
print(table(d.inc$INCOME_5, useNA = "ifany"))
print(round(100 * prop.table(table(d.inc$INCOME_5)), 1))

cat("\nSmallest FLI level:   ", min(table(d.fl$FLI_3item)), "cases\n")
cat("Smallest income level:", min(table(d.inc$INCOME_5)), "cases\n")
# In the low thousands = fine. In the low dozens = consider collapsing.

# Visual.
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
# -----------------------------------------------------------------------------
# The MEASUREMENT problem is solved - every respondent now sits on the same
# 5-level scale regardless of whether WP88 was fielded. What remains is a
# SAMPLE COMPOSITION question: the countries where WP88 was not administered
# may differ from the rest, which matters for the 3c robustness check.

# (i) Both groups should now reach all five levels. This is the demonstration
#     that the recoding fixed the scale problem.
cat("\n--- FLI_3item by WP88 availability ---\n")
print(table(d.fl$FLI_3item, has_wp88 = d.fl$has_wp88)) #correct

# (ii) How the new index maps onto Gallup's. Shows exactly what dropping WP88
#      does to each respondent's score.
cat("\n--- FLI_3item (rows) vs Gallup INDEX_FL_ord (cols) ---\n")
print(table(d.fl$FLI_3item, d.fl$INDEX_FL_ord))

cat("\nCorrelation between the two versions: ",
    round(cor(as.numeric(as.character(d.fl$FLI_3item)),
              as.numeric(as.character(d.fl$INDEX_FL_ord)),
              use = "pairwise.complete.obs"), 3), "\n")

# (iii) Do respondents with and without WP88 differ on your key variables?
d.fl %>%
  group_by(has_wp88) %>%
  summarise(n = n(),
            mean_index = round(mean(as.numeric(as.character(FLI_3item))), 2),
            across(all_of(c("iwisescore", "hhsize", "gdp_pc_ppp", v.jmp)),  # [JMP]
                   ~round(mean(.x, na.rm = TRUE), 2)),
            pct_female = round(100 * mean(female == "Female", na.rm = TRUE), 1),
            pct_urban  = round(100 * mean(urban_imp == "Peri-urban/urban",
                                          na.rm = TRUE), 1),
            .groups = "drop")

# (iv) Is WP88 non-administration a COUNTRY-level pattern? Countries at 100%
#      confirm it was a fielding decision, not respondent refusal - which is
#      what justifies dropping the item rather than the respondents.
d.fl %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(n = n(), pct_no_wp88 = round(100 * mean(!has_wp88), 1),
            .groups = "drop") %>%
  arrange(desc(pct_no_wp88)) %>%
  print(n = Inf)


# 3c. ROBUSTNESS SAMPLE. Refit the final model using Gallup's own 4-item index
#     on the subsample where WP88 was fielded, to confirm your conclusions do
#     not depend on the decision to drop it. Note this has 7 levels, and covers
#     only the countries where WP88 was administered.
d.fl.gallup <- d.fl %>%
  filter(has_wp88) %>%
  mutate(INDEX_FL_ord = droplevels(INDEX_FL_ord))

cat("\nRobustness sample (Gallup 4-item, WP88 countries only): N =",
    nrow(d.fl.gallup),
    "| levels:", paste(levels(d.fl.gallup$INDEX_FL_ord), collapse = ", "), "\n")


# =============================================================================
# CHECK 4: PREDICTOR DISTRIBUTIONS AND RARE CATEGORIES
# -----------------------------------------------------------------------------
# WHY: a factor level with a handful of respondents gives a coefficient too
# imprecise to interpret. Implausible values are usually data errors. Heavy
# skew distorts estimates.
# =============================================================================

# 4a. Categorical predictors.
for (v in v.categorical) {
  cat("\n---", v, "(FLI sample) ---\n")
  print(table(d.fl[[v]], useNA = "ifany"))
}

# 4b. Continuous predictors.
d.fl %>%
  select(any_of(c("iwisescore", "hhsize", "gdp_pc_ppp", "log_gdp", v.wgi,
                  v.jmp))) %>%                                          # [JMP]
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
  select(any_of(c("iwisescore", "hhsize", "gdp_pc_ppp", "log_gdp",
                  v.jmp))) %>%                                          # [JMP]
  pivot_longer(everything()) %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 40, fill = "#3D4222") +
  facet_wrap(~name, scales = "free") +
  labs(title = "Continuous predictors, FLI sample") +
  theme_minimal()

# [JMP] 4c-bis. How much does bas_rate_pp actually vary? Country-level, so
#     judge it across COUNTRY-YEARS, not respondents. If most values sit near
#     zero, the covariate has little variation to explain anything with.
cy %>%
  summarise(n_country_years = sum(!is.na(bas_rate_pp)),
            median   = median(bas_rate_pp, na.rm = TRUE),
            iqr_lo   = quantile(bas_rate_pp, .25, na.rm = TRUE),
            iqr_hi   = quantile(bas_rate_pp, .75, na.rm = TRUE),
            pct_flat = round(100 * mean(abs(bas_rate_pp) < 0.05, na.rm = TRUE), 1))

# 4d. How is iwisescore distributed across the standard bands?
table(d.fl$iwise_cat, useNA = "ifany")


# =============================================================================
# CHECK 5: OUTCOME x CATEGORICAL PREDICTOR CROSS-TABS
# -----------------------------------------------------------------------------
# WHY: empty cells (zero respondents in a predictor-level x outcome-level
# combination) are the main cause of non-convergence and of separation.
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

crosstab_check(d.fl,  "FLI_3item", v.categorical)
crosstab_check(d.inc, "INCOME_5",  v.categorical)
# IF EMPTY CELLS: collapse the offending predictor category.


# =============================================================================
# CHECK 6: SAMPLE SIZE / EVENTS PER VARIABLE
# -----------------------------------------------------------------------------
# WHY: need enough cases in the SMALLEST outcome category to estimate all
# coefficients reliably. Rule of thumb: at least 10 per parameter.
# NOTE: a k-level factor costs k-1 parameters, not 1.
# =============================================================================

rhs     <- paste("iwisescore + female + urban_imp + age_gp_profile +",
                 "maritalstatus + hhsize + employment + education")

# Contextual version. ONE governance term only - see Check 7.
# NB: this must match whichever WGI family v.wgi points at.
# [JMP] bas_rate_pp added. Every check using rhs_ctx (6, 7b) now includes it.
rhs_ctx <- paste(rhs, "+ log_gdp +", v.wgi[6], "+", v.jmp)

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
# WHY: if two predictors carry nearly the same information, the model cannot
# separate their effects - unstable coefficients, wide standard errors, signs
# that flip when you add or remove a variable.
#
# This is a property of the PREDICTORS ONLY. The outcome is irrelevant, which
# is why it can be checked before fitting.
#
# [JMP] The country-year IWISE-vs-JMP matrix is in the setup section above.
# =============================================================================

# 7a. Correlation matrix. Catches PAIRWISE overlap only.
d.fl %>%
  select(any_of(c("iwisescore", "hhsize", "log_gdp", v.wgi, v.jmp,   # [JMP]
                  paste0(v.indices, "_num")))) %>%
  cor(use = "pairwise.complete.obs") %>%
  round(2) %>%
  print()
# LOOKING FOR: |r| above ~0.8. Expect it among the WGI indicators.

# 7b. VIF. Catches overlap involving THREE OR MORE variables, which pairwise
#     correlation misses. VIF = 4 means this standard error is twice as wide
#     as it would be with uncorrelated predictors.
#     A random outcome gives the identical answer - VIF ignores the outcome.
set.seed(2024)
m.vifprobe <- lm(as.formula(paste("rnorm(nrow(d.fl)) ~", rhs_ctx)), data = d.fl)
car::vif(m.vifprobe)
# For FACTORS this returns GVIF. Read GVIF^(1/(2*Df)) and SQUARE it before
# comparing to the usual thresholds of 5 or 10.
# [JMP] Before JMP: log_gdp 2.00, wgi_cc_sc 1.73, iwisescore 1.09. Compare the
# new values against these; expect log_gdp to rise a little (r = -0.53).

# =============================================================================
# CHECK 8: CLUSTERING / NON-INDEPENDENCE
# -----------------------------------------------------------------------------
# WHY: all regression assumes independent observations. Yours are nested in
# countries, and log_gdp / wgi_* / bas_rate_pp take ONE value per country-year.  # [JMP]
# Ignoring this makes standard errors far too small, so you find effects that
# are not there.
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
# 74 countries, 500 to 3,503 respondents each
# [JMP] Rerun on the complete-case sample to see how many survive with JMP:
d.fl %>% filter(cc) %>% count(across(all_of(v.cluster))) %>%
  summarise(n_countries_cc = n())

# 8b. ICC: what share of outcome variance sits BETWEEN countries?
m.icc.fl <- lme4::lmer(
  as.numeric(as.character(FLI_3item)) ~ 1 + (1 | country), data = d.fl)
cat("\nICC, FLI_3item:\n"); print(performance::icc(m.icc.fl))

m.icc.inc <- lme4::lmer(
  as.numeric(as.character(INCOME_5)) ~ 1 + (1 | country), data = d.inc)
cat("\nICC, INCOME_5:\n"); print(performance::icc(m.icc.inc))

# LOOKING FOR: above ~0.05 means clustering is NOT ignorable.
# EXPECT INCOME_5 TO BE NEAR ZERO. Quintiles are constructed WITHIN country,
# so every country has ~20% in each - there is almost nothing between countries
# to explain, and country-level covariates cannot explain it. That is a
# substantive problem with the income model, not a bug.

# 8c. Visual: how much do country means differ?
d.fl %>%
  group_by(across(all_of(v.cluster))) %>%
  summarise(mean_fl = mean(as.numeric(as.character(FLI_3item)), na.rm = TRUE),
            .groups = "drop") %>%
  ggplot(aes(x = reorder(.data[[v.cluster]], mean_fl), y = mean_fl)) +
  geom_point() + coord_flip() +
  labs(title = "Country mean, FLI_3item", x = NULL, y = "Mean") +
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
# The ordinal model works by cutting the outcome at each threshold, so each cut
# can be tested in advance. 5 levels = 4 cuts (was 8 with the 9-level version).
#
# [JMP] Unchanged: this check uses the individual-level rhs, so it does not
# include log_gdp, the WGI or bas_rate_pp. Swap rhs for rhs_ctx below if you
# want it on the full model.
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

# IF DETECTED: collapse the offending category, drop the variable, or plan on
# Firth penalised likelihood (brglm2::brglmFit) at the modelling stage.


# =============================================================================
# SAVE FOR THE MODELLING SCRIPT
# -----------------------------------------------------------------------------
# One saveRDS only. Saving the frames AND the variable definitions together
# keeps them in sync - you cannot accidentally load an old d.fl with a new rhs.
# Re-run this script whenever you change anything upstream, or the modelling
# script will silently analyse stale data.
# =============================================================================

saveRDS(list(d.iwise      = d.iwise,
             d.fl         = d.fl,
             d.inc        = d.inc,
             d.fl.gallup  = d.fl.gallup,   # robustness sample, Gallup 4-item
             rhs          = rhs,
             rhs_ctx      = rhs_ctx,       # [JMP] now includes bas_rate_pp
             v.main       = v.main,
             v.categorical = v.categorical,
             v.continuous = v.continuous,
             v.wgi        = v.wgi,
             v.jmp        = v.jmp,         # [JMP]
             v.indices    = v.indices,
             v.cluster    = v.cluster,
             v.model.fl   = v.model.fl,
             v.model.inc  = v.model.inc),
        "data_prepared.rds")

cat("\nSaved data_prepared.rds\n")
cat("Load in the modelling script with:\n")
cat('  dat <- readRDS("data_prepared.rds"); list2env(dat, envir = .GlobalEnv)\n')


# =============================================================================
# WHAT COMES NEXT - checks that require a FITTED model
# =============================================================================

# 10. PROPORTIONAL ODDS - the assumption specific to this model. Tests whether
#     each predictor's effect is the same across every cut of the outcome. It
#     is frequently violated. Brant test, ordinal::nominal_test(), and a plot
#     of coefficients across cuts. With your N the formal tests reject on
#     trivial departures, so the plot is the most honest of the three.

# 11. LINEARITY OF THE LOGIT - whether continuous predictors relate linearly to
#     the log-odds. That scale only exists once a model defines it, and the
#     relationship is PARTIAL, so it depends on the other covariates.
#     iwisescore is zero-inflated (median 3, mean 6.6), so expect this to flag;
#     iwise_cat is the ready-made alternative.
#     [JMP] Test bas_rate_pp here too. It can be negative, so use splines
#     rather than a log or Box-Tidwell term.

# 12. INFLUENTIAL OBSERVATIONS - Cook's distance, standardised residuals, and a
#     sensitivity refit without the most influential 1%.
#     [JMP] With a country-level covariate, also check influential COUNTRIES:
#     refit dropping one country at a time and watch the bas_rate_pp estimate.

# 13. GOODNESS OF FIT - Lipsitz and Pulkstenis-Robinson tests.