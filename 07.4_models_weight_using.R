# =============================================================================
# 08_models_wgt_using.R
# Main models refitted with wgt_using (= wgt2, Palestine 2025 missing).
#   Outcomes : INCOME_5 (d.inc), FLI_3item (d.fl)
#   IWISE    : continuous (iwisescore) and categorical (iwise_cat)
#   M0       : unadjusted          M1 : fully adjusted (rhs_ctx)
#   Random   : intercept + IWISE slope per country_year (as in 07.1 / 07.2)
#   Weights  : wgt_using, NOT rescaled (sums roughly to n).
# Saves: models_wgt_using_cont.rds, models_wgt_using_cat.rds
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
library(tidyverse)
library(ordinal)                                           # clmm()


# ---- 1. DATA ----------------------------------------------------------------
list2env(readRDS("data/iwise_data_prepared.rds"), envir = globalenv())  # d.fl, d.inc, rhs_ctx, v.*

prep <- function(data, y) {                                # analytic sample for one outcome
  data %>%
    filter(cc, !is.na(wgt_using)) %>%                      # complete cases; drops Palestine 2025
    mutate(country_year = factor(.data[[v.cluster]]),      # grouping factor
           iwise_cat    = factor(iwise_cat, ordered = FALSE),  # one dummy per category
           iwise_num    = as.integer(iwise_cat) - 1,       # 0,1,2,3: linear random slope only
           "{y}" := factor(.data[[y]], ordered = TRUE),    # clmm() needs an ordered outcome
           w_fit        = wgt_using)                       # not rescaled
}

d.inc.cc <- prep(d.inc, "INCOME_5")
d.fl.cc  <- prep(d.fl,  "FLI_3item")

levels(d.fl.cc$iwise_cat)                                  # first level = lowest WI (reference)
table(d.fl.cc$iwise_cat, d.fl.cc$iwise_num)                # diagonal, no sparse categories


# ---- 2. FORMULAS AND HELPERS ------------------------------------------------
rhs_ctx_cat <- str_replace(rhs_ctx, "\\biwisescore\\b", "iwise_cat")
rhs_ctx_cat                                                # check: iwise_cat in, iwisescore out

re_term      <- "(1 | country_year)"                       # random intercept only
re_slope     <- "(1 + iwisescore | country_year)"          # continuous slope
re_slope_cat <- "(1 + iwise_cat | country_year)"           # full: random slope per category
re_slope_num <- "(1 + iwise_num | country_year)"           # linear: one slope across categories

make_f <- function(y, rhs, re) as.formula(paste(y, "~", rhs, "+", re))

fit_clmm <- function(data, y, rhs, re) {                   # same settings as 07.1 / 07.2
  clmm(make_f(y, rhs, re), data = data, weights = w_fit,
       link = "logit", Hess = TRUE, nAGQ = 1)
}


# ---- 3. CONTINUOUS IWISE (SLOW) ---------------------------------------------
m.inc.0.wu <- fit_clmm(d.inc.cc, "INCOME_5",  "iwisescore", re_term)   # income, unadjusted (intercept only, 07.1)
m.inc.1.wu <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx,      re_slope)  # income, adjusted
m.fl.0.wu  <- fit_clmm(d.fl.cc,  "FLI_3item", "iwisescore", re_slope)  # FLI, unadjusted
m.fl.1.wu  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx,      re_slope)  # FLI, adjusted

saveRDS(mget(c("m.inc.0.wu", "m.inc.1.wu", "m.fl.0.wu", "m.fl.1.wu")),
        "models_wgt_using_cont.rds")


# ---- 4. CATEGORICAL IWISE (SLOW) --------------------------------------------
m.inc.0.cat.wu <- fit_clmm(d.inc.cc, "INCOME_5",  "iwise_cat", re_term)       # income, unadjusted (intercept only, 07.2)
m.inc.1.cat.wu <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx_cat, re_slope_cat)  # income, adjusted (full)
m.fl.0.cat.wu  <- fit_clmm(d.fl.cc,  "FLI_3item", "iwise_cat", re_slope_num)  # FLI, unadjusted (linear, 07.2)
m.fl.1.cat.wu  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx_cat, re_slope_cat)  # FLI, adjusted (full)

saveRDS(mget(c("m.inc.0.cat.wu", "m.inc.1.cat.wu", "m.fl.0.cat.wu", "m.fl.1.cat.wu")),
        "models_wgt_using_cat.rds")
# Reload: list2env(readRDS("models_wgt_using_cont.rds"), envir = globalenv())
#         list2env(readRDS("models_wgt_using_cat.rds"),  envir = globalenv())