# =============================================================================
# ASSUMPTION CHECKS FOR LOGISTIC REGRESSION
# Outcomes: INDEX_FL_ord (Financial Life Index) and INCOME_5 (Income Quintiles)
# Main predictor: iwisescore (Water Insecurity Experiences, 0-36)
# =============================================================================
#
# READ THIS FIRST -------------------------------------------------------------
#
# Neither of your outcomes is binary:
#   INDEX_FL_ord has 9 ordered levels (0, 16.7, 25, 33.3, 50, 66.7, 75, 83.3, 100)
#   INCOME_5     has 5 ordered levels (quintiles)
#
# ORDINAL logistic regression (the proportional odds / cumulative logit model), 
# everything from binary logistic still
# applies, PLUS the proportional odds assumption (Section 8), which is the one
# that usually fails in practice.
#
# The script therefore does two things:
#   (a) runs every check on the ordinal models (the recommended analysis)
#   (b) gives you the binary versions of the same checks in case you decide to
#       dichotomise (e.g. bottom quintile vs rest), which is a defensible
#       fallback if proportional odds is badly violated.
#
# One more structural issue: your contextual covariates (gdp_pc_ppp, wgi_*) vary
# only at country level, and respondents are clustered in countries. That
# violates the independence-of-observations assumption. Section 7 quantifies it
# and shows the two standard fixes.
#
# =============================================================================

# ---- 0. PACKAGES ------------------------------------------------------------
# Install once, then comment out.
pkgs <- c(
  "tidyverse",        # data wrangling + ggplot
  "MASS",             # polr(): ordinal logistic regression
  "ordinal",          # clm() + nominal_test(): better proportional odds testing
  "brant",            # brant(): omnibus proportional odds test
  "car",              # vif(): multicollinearity
  "detectseparation", # complete/quasi-complete separation
  "performance",      # check_collinearity, check_outliers, model diagnostics
  "broom",            # tidy model output
  "generalhoslem",    # Lipsitz / Pulkstenis-Robinson GOF tests for ordinal models
  "ResourceSelection",# Hosmer-Lemeshow test (binary case)
  "lme4",             # glmer(): multilevel fallback for clustering
  "sandwich", "lmtest", "clubSandwich", # cluster-robust standard errors
  "naniar",           # missing data patterns
  "sjPlot"            # plotting model results
)
# install.packages(setdiff(pkgs, rownames(installed.packages())))
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(2024)


# ---- 1. LOAD AND PREPARE DATA ----------------------------------------------
# Replace with your own path. haven::read_dta() if it's a Stata .dta file.
# df <- haven::read_dta("your_data.dta")
# df <- readRDS("your_data.rds")

# Define your variable groups once, so the rest of the script is generic.
outcomes <- c("INDEX_FL_ord", "INCOME_5")

main_pred <- "iwisescore"          # 0-36 integer
ctrl_pred <- "iwise4"              # 0-12 integer, alternative specification

# Continuous / count covariates -> these are the ones that need a linearity check
continuous <- c("iwisescore", "hhsize", "gdp_pc_ppp",
                "wgi_va_xx", "wgi_pv_xx", "wgi_ge_xx",
                "wgi_rq_xx", "wgi_rl_xx", "wgi_cc_xx")

# Categorical covariates -> these need cell-count checks, not linearity checks
categorical <- c("female", "urban_imp", "age_gp_profile",
                 "maritalstatus", "employment", "education")

# Other indices you may want as covariates or outcomes
indices <- c("INDEX_PH", "INDEX_FS", "INDEX_CA", "INDEX_YD")

cluster_var <- "country"           # rename to match your data

df <- df %>%
  mutate(
    # Ordinal outcomes MUST be ordered factors for polr()/clm() to work.
    # as.numeric(as.character(.)) first, in case they arrive as labelled/haven
    # vectors, then order the levels numerically.
    INDEX_FL_ord = factor(as.numeric(as.character(INDEX_FL_ord)),
                          ordered = TRUE),
    INCOME_5     = factor(as.numeric(as.character(INCOME_5)),
                          ordered = TRUE),
    # Unordered factors for the categorical covariates
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
    across(all_of(continuous), as.numeric)
  )

# INCOME_5 is not available for Zimbabwe (per your codebook). Drop it from the
# income models only, otherwise you silently lose Zimbabwe from both models and
# the two samples stop being comparable.
df_income <- df %>% filter(!is.na(INCOME_5))

cat("\nN total:", nrow(df), " | N with INCOME_5:", nrow(df_income), "\n")


# ---- 2. MISSING DATA --------------------------------------------------------
# Not a model assumption as such, but complete-case analysis assumes data are
# Missing Completely At Random. If missingness is related to your predictors,
# listwise deletion biases the estimates. Check before you check anything else,
# because every diagnostic below is computed on the complete cases.

model_vars <- c(outcomes, main_pred, ctrl_pred, continuous, categorical,
                cluster_var)

df %>%
  select(any_of(model_vars)) %>%
  naniar::miss_var_summary() %>%
  print(n = 30)

naniar::gg_miss_upset(df %>% select(any_of(model_vars)))

# Compare complete cases vs dropped cases on the main predictor. A large,
# systematic difference is a warning sign.
df %>%
  mutate(complete = complete.cases(select(., any_of(model_vars)))) %>%
  group_by(complete) %>%
  summarise(n = n(), mean_iwise = mean(iwisescore, na.rm = TRUE))


# ---- 3. OUTCOME DISTRIBUTION AND SPARSE CELLS -------------------------------
# Ordinal models estimate one intercept (cut-point) per threshold. With 9 levels
# you need 8 cut-points; if some levels hold very few respondents those
# cut-points are unstable and the proportional odds test loses power.
# Rule of thumb: collapse any level with < ~2% of the sample.

table(df$INDEX_FL_ord, useNA = "ifany")
prop.table(table(df$INDEX_FL_ord)) %>% round(3)

table(df_income$INCOME_5, useNA = "ifany")

# Cross-tabulate the outcome against every categorical predictor. Zero cells are
# the main cause of separation (Section 6) and of non-convergence.
for (v in categorical) {
  cat("\n--- INDEX_FL_ord x", v, "---\n")
  print(table(df[[v]], df$INDEX_FL_ord))
}

# If INDEX_FL_ord is too sparse, a defensible collapse into 4-5 levels:
df <- df %>%
  mutate(INDEX_FL_c = cut(as.numeric(as.character(INDEX_FL_ord)),
                          breaks = c(-Inf, 25, 50, 75, Inf),
                          labels = c("Low", "Lower-mid", "Upper-mid", "High"),
                          ordered_result = TRUE))
table(df$INDEX_FL_c)


# ---- 4. SAMPLE SIZE / EVENTS PER VARIABLE -----------------------------------
# The classic rule is >= 10 events per estimated parameter (Peduzzi et al. 1996).
# For an ordinal model, apply it to the smallest outcome category, which is the
# binding constraint.

fml_full <- as.formula(
  paste("~", paste(c(main_pred, categorical, "hhsize", "gdp_pc_ppp",
                     "wgi_ge_xx", "wgi_cc_xx"), collapse = " + "))
)
n_params <- ncol(model.matrix(fml_full, data = df))
min_cell <- min(table(df$INDEX_FL_ord))

cat("\nParameters:", n_params,
    "| smallest outcome category:", min_cell,
    "| EPV:", round(min_cell / n_params, 1), "\n")
# EPV < 10 -> consider collapsing outcome levels, dropping covariates, or
# penalised estimation (logistf / brglm2).


# ---- 5. MULTICOLLINEARITY ---------------------------------------------------
# Logistic regression assumes predictors are not near-linear combinations of one
# another. Two specific risks in your setup:
#   (i)  the six WGI governance indicators are typically correlated at r > 0.8
#   (ii) your codebook flags INDEX_YD as possibly collinear with education
#
# VIF is computed on a linear model with the same right-hand side — the metric
# only concerns the predictor matrix, so the outcome's distribution is
# irrelevant. For factors, read the GVIF^(1/(2*df)) column and square it before
# comparing to the usual VIF thresholds (5 or 10).

# 5a. Correlation matrix of the continuous block
df %>%
  select(all_of(continuous), any_of(indices)) %>%
  cor(use = "pairwise.complete.obs") %>%
  round(2) %>%
  print()

# 5b. VIF
vif_probe <- lm(as.numeric(as.character(INDEX_FL_ord)) ~ .,
                data = df %>%
                  select(INDEX_FL_ord, all_of(main_pred),
                         all_of(categorical), hhsize, gdp_pc_ppp,
                         starts_with("wgi_")) %>%
                  drop_na())
car::vif(vif_probe)

# 5c. Same check directly on the fitted ordinal model
m_fl <- MASS::polr(
  INDEX_FL_ord ~ iwisescore + female + urban_imp + age_gp_profile +
    maritalstatus + hhsize + employment + education + gdp_pc_ppp,
  data = df, Hess = TRUE, method = "logistic"
)
performance::check_collinearity(m_fl)

# If the WGI indicators are collinear (they almost certainly are): pick ONE on
# theoretical grounds, or reduce all six to a single component:
wgi_pca <- prcomp(df %>% select(starts_with("wgi_")) %>% drop_na(),
                  scale. = TRUE)
summary(wgi_pca)   # PC1 usually explains 80-90% of WGI variance


# ---- 6. SEPARATION ----------------------------------------------------------
# Complete or quasi-complete separation means some predictor (or combination)
# perfectly predicts the outcome. Symptoms: enormous coefficients, standard
# errors in the thousands, convergence warnings. It is a data problem, not a
# model problem — more iterations will not fix it.
#
# detectseparation works on binary glm, so test each cumulative split of the
# ordinal outcome: this mirrors exactly what the cumulative logit model does.

check_separation <- function(data, outcome, rhs) {
  lv <- levels(data[[outcome]])
  for (k in seq_len(length(lv) - 1)) {
    d <- data %>%
      mutate(y_bin = as.integer(as.numeric(.data[[outcome]]) > k)) %>%
      drop_na(y_bin)
    fit <- glm(as.formula(paste("y_bin ~", rhs)), data = d,
               family = binomial,
               method = detectseparation::detect_separation)
    cat("Split above level", lv[k], "->",
        ifelse(fit$outcome, "SEPARATION DETECTED", "ok"), "\n")
  }
}

rhs <- "iwisescore + female + urban_imp + age_gp_profile + maritalstatus + hhsize + employment + education"
check_separation(df, "INDEX_FL_ord", rhs)
check_separation(df_income, "INCOME_5", rhs)

# Fix if detected: collapse the offending category, drop the variable, or use
# Firth penalised likelihood (brglm2::brglmFit).


# ---- 7. INDEPENDENCE OF OBSERVATIONS ----------------------------------------
# Your respondents are nested in countries, and gdp_pc_ppp / wgi_* are constant
# within country. Ordinary standard errors assume independent observations and
# will be far too small here — the classic Moulton problem: you have thousands
# of respondents but only as many independent observations of the contextual
# covariates as you have countries.

# 7a. How much outcome variance sits between countries? (ICC)
icc_probe <- lme4::lmer(
  as.numeric(as.character(INDEX_FL_ord)) ~ 1 + (1 | country), data = df
)
performance::icc(icc_probe)
# ICC > ~0.05 means clustering is not ignorable.

# 7b. Fix 1 — cluster-robust SEs. Keeps polr point estimates, corrects inference.
lmtest::coeftest(m_fl, vcov = sandwich::vcovCL(m_fl, cluster = df$country))
# With few clusters (< ~40 countries) use the small-sample correction instead:
clubSandwich::coef_test(m_fl, vcov = "CR2", cluster = df$country)

# 7c. Fix 2 — multilevel ordinal model with a country random intercept.
#     Preferred if you want to interpret between-country variation.
m_fl_mixed <- ordinal::clmm(
  INDEX_FL_ord ~ iwisescore + female + urban_imp + age_gp_profile +
    maritalstatus + hhsize + employment + education + gdp_pc_ppp +
    (1 | country),
  data = df, link = "logit"
)
summary(m_fl_mixed)

# 7d. If the Gallup World Poll design weights are in the file, use them.
#     Unweighted estimates are biased for population quantities.
# m_fl_w <- MASS::polr(..., weights = wgt, data = df)
# Better: survey::svyolr() with the full design, so SEs account for the design.


# ---- 8. PROPORTIONAL ODDS ASSUMPTION ----------------------------------------
# THE key assumption specific to ordinal logistic regression. It states that the
# effect of each predictor is the same across every cumulative split of the
# outcome — i.e. the odds ratio for iwisescore is identical whether you compare
# "0 vs everything above" or "83.3 vs 100". If it fails, a single coefficient
# per predictor is a misleading summary.

# 8a. Brant test (Brant 1990). Omnibus test plus one test per variable.
#     H0 = proportional odds holds. p < 0.05 = violation.
#     Caution: with n in the tens of thousands this test rejects on trivial
#     departures. Read it alongside 8b and 8c, don't treat p < 0.05 as fatal.
brant::brant(m_fl)

# 8b. Likelihood-ratio tests per variable, via ordinal::clm(). More informative
#     than Brant because it tells you exactly which predictors to relax.
m_fl_clm <- ordinal::clm(
  INDEX_FL_ord ~ iwisescore + female + urban_imp + age_gp_profile +
    maritalstatus + hhsize + employment + education + gdp_pc_ppp,
  data = df, link = "logit"
)
ordinal::nominal_test(m_fl_clm)   # significant row = that variable violates PO
ordinal::scale_test(m_fl_clm)     # significant = scale effects, consider clm scale=

# 8c. Graphical check — the most honest one. Fit a separate binary logit at each
#     cumulative split and plot the coefficients. Roughly flat lines with
#     overlapping CIs = proportional odds is reasonable; a clear trend or sign
#     flip = it is not.
po_plot_data <- function(data, outcome, rhs) {
  lv <- levels(data[[outcome]])
  map_dfr(seq_len(length(lv) - 1), function(k) {
    d <- data %>%
      mutate(y_bin = as.integer(as.numeric(.data[[outcome]]) > k))
    glm(as.formula(paste("y_bin ~", rhs)), data = d, family = binomial) %>%
      broom::tidy(conf.int = TRUE) %>%
      mutate(split = paste0("> ", lv[k]), split_n = k)
  })
}

po_dat <- po_plot_data(df, "INDEX_FL_ord", rhs)

po_dat %>%
  filter(term != "(Intercept)") %>%
  ggplot(aes(x = split_n, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  facet_wrap(~ term, scales = "free_y") +
  labs(x = "Cumulative split of the outcome", y = "Log-odds coefficient",
       title = "Proportional odds check: coefficients across cumulative splits",
       subtitle = "Flat = assumption plausible; trending = violated") +
  theme_minimal()

# 8d. If proportional odds fails, in increasing order of complexity:
#   - Partial proportional odds: relax it only for the offending variables
m_fl_ppo <- ordinal::clm(
  INDEX_FL_ord ~ iwisescore + female + urban_imp + hhsize + gdp_pc_ppp,
  nominal = ~ education,          # education gets its own coefficient per split
  data = df, link = "logit"
)
anova(m_fl_clm, m_fl_ppo)         # LRT: does relaxing it improve fit?
#   - Multinomial logistic regression (nnet::multinom), ignoring the ordering
#   - Continuation ratio or adjacent category models
#   - Dichotomise and run binary logistic (see Section 11)


# ---- 9. LINEARITY OF THE LOGIT ----------------------------------------------
# Logistic regression assumes each CONTINUOUS predictor is linearly related to
# the log-odds — not to the outcome itself. Only iwisescore, hhsize, gdp_pc_ppp
# and the wgi_* variables need this check; factors are exempt by construction.
# This is plausibly the assumption most at risk for you: water insecurity may
# have a threshold effect rather than a smooth gradient.

# 9a. Box-Tidwell test. Add x * log(x) to the model; a significant coefficient
#     on that term means the relationship is not linear on the logit scale.
#     x must be strictly positive, so shift variables containing zeros.
box_tidwell <- function(data, outcome, vars, rhs) {
  d <- data %>% drop_na(all_of(c(outcome, vars)))
  for (v in vars) {
    x <- d[[v]]
    if (min(x, na.rm = TRUE) <= 0) x <- x - min(x, na.rm = TRUE) + 0.5
    d$.bt <- x * log(x)
    f <- as.formula(paste(outcome, "~", rhs, "+ .bt"))
    fit <- MASS::polr(f, data = d, Hess = TRUE)
    ct <- coef(summary(fit))
    p  <- 2 * (1 - pnorm(abs(ct[".bt", "t value"])))
    cat(sprintf("%-14s Box-Tidwell p = %.4f %s\n", v, p,
                ifelse(p < 0.05, "<- NON-LINEAR", "")))
  }
}
box_tidwell(df, "INDEX_FL_ord", c("iwisescore", "hhsize", "gdp_pc_ppp"), rhs)

# 9b. Empirical logit plot — shows you the SHAPE of any non-linearity, which the
#     p-value does not. Bin the predictor, compute the observed log-odds of
#     being above the outcome median in each bin, and plot.
emp_logit_plot <- function(data, outcome, var, n_bins = 20) {
  med <- median(as.numeric(as.character(data[[outcome]])), na.rm = TRUE)
  data %>%
    mutate(y_bin = as.integer(as.numeric(as.character(.data[[outcome]])) > med),
           bin = ntile(.data[[var]], n_bins)) %>%
    drop_na(y_bin, bin) %>%
    group_by(bin) %>%
    summarise(x = mean(.data[[var]], na.rm = TRUE),
              n = n(),
              p = (sum(y_bin) + 0.5) / (n() + 1),   # Haldane correction
              .groups = "drop") %>%
    mutate(emp_logit = log(p / (1 - p))) %>%
    ggplot(aes(x = x, y = emp_logit)) +
    geom_point(aes(size = n)) +
    geom_smooth(method = "loess", se = TRUE, formula = y ~ x) +
    geom_smooth(method = "lm", se = FALSE, colour = "red",
                linetype = "dashed", formula = y ~ x) +
    labs(x = var, y = "Empirical logit",
         title = paste("Linearity of the logit:", var),
         subtitle = "Red dashed = linear fit. Systematic departure = violation.") +
    theme_minimal()
}
emp_logit_plot(df, "INDEX_FL_ord", "iwisescore")
emp_logit_plot(df, "INDEX_FL_ord", "hhsize")
emp_logit_plot(df, "INDEX_FL_ord", "gdp_pc_ppp")

# 9c. If non-linear, compare specifications formally rather than eyeballing:
m_lin  <- ordinal::clm(INDEX_FL_ord ~ iwisescore + female + urban_imp,
                       data = df)
m_quad <- ordinal::clm(INDEX_FL_ord ~ poly(iwisescore, 2) + female + urban_imp,
                       data = df)
m_spl  <- ordinal::clm(INDEX_FL_ord ~ splines::ns(iwisescore, 4) + female +
                         urban_imp, data = df)
AIC(m_lin, m_quad, m_spl)
anova(m_lin, m_quad, m_spl)

# A substantively attractive alternative for IWISE: use the standard cut-points
# (0-2 no-to-marginal, 3-11 low, 12-23 moderate, 24-36 high insecurity) and
# enter water insecurity as a factor. Then linearity stops being an assumption
# at all, and the coefficients are easier to interpret.
df <- df %>%
  mutate(iwise_cat = cut(iwisescore, breaks = c(-Inf, 2, 11, 23, Inf),
                         labels = c("No-to-marginal", "Low", "Moderate", "High"),
                         ordered_result = TRUE))


# ---- 10. INFLUENTIAL OBSERVATIONS AND OUTLIERS ------------------------------
# Not strictly an assumption, but a small number of extreme cases can drive the
# whole result. polr objects have no built-in influence measures, so run these
# on the equivalent binary splits, which is where influence actually shows up.

d_bin <- df %>%
  mutate(y_bin = as.integer(as.numeric(as.character(INDEX_FL_ord)) >
                              median(as.numeric(as.character(INDEX_FL_ord)),
                                     na.rm = TRUE))) %>%
  drop_na(y_bin, all_of(c(main_pred, categorical)), hhsize)

m_bin <- glm(as.formula(paste("y_bin ~", rhs)), data = d_bin, family = binomial)

# 10a. Cook's distance. Conventional cutoff 4/n, though with large n almost
#      nothing exceeds it — inspect the largest values relatively instead.
infl <- broom::augment(m_bin) %>% mutate(row = row_number())

infl %>%
  ggplot(aes(row, .cooksd)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 4 / nrow(infl), colour = "red", linetype = "dashed") +
  labs(title = "Cook's distance", y = "Cook's D", x = "Observation") +
  theme_minimal()

infl %>% slice_max(.cooksd, n = 20) %>% select(row, .cooksd, .std.resid)

# 10b. Standardised deviance residuals. |value| > 3 = poorly fitted case.
infl %>%
  filter(abs(.std.resid) > 3) %>%
  summarise(n_outliers = n(), pct = 100 * n() / nrow(infl))

infl %>%
  ggplot(aes(.fitted, .std.resid)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = c(-3, 3), colour = "red", linetype = "dashed") +
  labs(title = "Standardised residuals vs fitted",
       x = "Fitted probability", y = "Std. deviance residual") +
  theme_minimal()

# 10c. Leverage and a combined check
performance::check_outliers(m_bin, method = c("cook", "pareto"))

# 10d. Sensitivity: does dropping the top 1% of influence change conclusions?
drop_ids <- infl %>% slice_max(.cooksd, n = ceiling(0.01 * nrow(infl))) %>% pull(row)
m_bin_sens <- update(m_bin, data = d_bin[-drop_ids, ])
broom::tidy(m_bin) %>%
  select(term, est_full = estimate) %>%
  left_join(broom::tidy(m_bin_sens) %>% select(term, est_drop = estimate),
            by = "term") %>%
  mutate(pct_change = 100 * (est_drop - est_full) / abs(est_full))
# Changes above ~10-20% on your main predictor deserve a footnote in the paper.


# ---- 11. GOODNESS OF FIT ----------------------------------------------------
# Fit is not an assumption, but a badly fitting model makes assumption tests
# hard to interpret, so check it in the same pass.

# 11a. Ordinal models: Lipsitz and Pulkstenis-Robinson tests.
#      H0 = model fits. p < 0.05 = lack of fit.
generalhoslem::lipsitz.test(m_fl)
generalhoslem::pulkrob.chisq(m_fl, catvars = c("female", "urban_imp",
                                               "education"))

# 11b. Binary models: Hosmer-Lemeshow. Known to be sensitive to the number of
#      groups g; report a couple of values rather than cherry-picking one.
ResourceSelection::hoslem.test(m_bin$y, fitted(m_bin), g = 10)
ResourceSelection::hoslem.test(m_bin$y, fitted(m_bin), g = 20)

# 11c. Pseudo-R2 and classification (descriptive, not a test)
performance::r2_nagelkerke(m_fl)
performance::model_performance(m_bin)


# ---- 12. REPEAT EVERYTHING FOR INCOME_5 -------------------------------------
# Same checks, second outcome, Zimbabwe excluded.

m_inc <- MASS::polr(
  INCOME_5 ~ iwisescore + female + urban_imp + age_gp_profile +
    maritalstatus + hhsize + employment + education + gdp_pc_ppp,
  data = df_income, Hess = TRUE
)

m_inc_clm <- ordinal::clm(
  INCOME_5 ~ iwisescore + female + urban_imp + age_gp_profile +
    maritalstatus + hhsize + employment + education + gdp_pc_ppp,
  data = df_income, link = "logit"
)

performance::check_collinearity(m_inc)     # Section 5
brant::brant(m_inc)                        # Section 8
ordinal::nominal_test(m_inc_clm)           # Section 8
box_tidwell(df_income, "INCOME_5",
            c("iwisescore", "hhsize", "gdp_pc_ppp"), rhs)   # Section 9
emp_logit_plot(df_income, "INCOME_5", "iwisescore")
generalhoslem::lipsitz.test(m_inc)         # Section 11

# NOTE ON INCOME_5 SPECIFICALLY -----------------------------------------------
# Income quintiles are constructed WITHIN country. That means the quintile is a
# relative rank, not an absolute income level, so country-level covariates
# (gdp_pc_ppp, wgi_*) cannot explain variation in it by construction — every
# country has exactly 20% of respondents in each quintile. Including them in the
# INCOME_5 model is close to meaningless unless you have raw income instead.
# Also: hhsize is in your codebook for per-capita income, so check whether the
# quintiles are already per-capita adjusted before adding hhsize as a covariate,
# or you will be adjusting for it twice.


# ---- 13. SUMMARY TABLE ------------------------------------------------------
# Paste the results into a single table for your methods appendix.

assumption_summary <- tribble(
  ~Assumption,                    ~Test,                          ~Section,
  "Correct outcome scale",        "Level frequencies, sparse cells",       "3",
  "Adequate sample size",         "Events per variable >= 10",             "4",
  "No multicollinearity",         "VIF / GVIF^(1/2df) < 2.2",              "5",
  "No separation",                "detect_separation per split",           "6",
  "Independent observations",     "ICC, cluster-robust SE, clmm",          "7",
  "Proportional odds",            "Brant, nominal_test, coefficient plot", "8",
  "Linearity of the logit",       "Box-Tidwell, empirical logit plot",     "9",
  "No influential observations",  "Cook's D, std. residuals, sensitivity", "10",
  "Adequate fit",                 "Lipsitz, Pulkstenis-Robinson",          "11"
)
print(assumption_summary, n = Inf)

# =============================================================================
# END
# =============================================================================

