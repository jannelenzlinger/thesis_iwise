# =============================================================================
# 07.6_models_weightless.R
# Main models refitted with wgt_using (= wgt2, Palestine 2025 missing).
#   Outcomes : INCOME_5 (d.inc), FLI_3item (d.fl)
#   IWISE    : continuous (iwisescore) and categorical (iwise_cat)
#   M0       : unadjusted          M1 : fully adjusted (rhs_ctx)
#   Random   : intercept + IWISE slope per country_year (as in 07.1 / 07.2)
#   Weights  : NO weight used
#   Saves: models_wgt_using_cont.rds, models_wgt_using_cat.rds
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
  clmm(make_f(y, rhs, re), data = data, link = "logit", Hess = TRUE, nAGQ = 1)
}


# ---- 3. CONTINUOUS IWISE (SLOW) ---------------------------------------------
m.inc.0.no.wgt <- fit_clmm(d.inc.cc, "INCOME_5",  "iwisescore", re_term)   # income, unadjusted (intercept only, 07.1)
m.inc.1.no.wgt <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx,      re_slope)  # income, adjusted
m.fl.0.wu.no.wgt  <- fit_clmm(d.fl.cc,  "FLI_3item", "iwisescore", re_slope)  # FLI, unadjusted
m.fl.1.no.wgt  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx,      re_slope)  # FLI, adjusted

saveRDS(mget(c("m.inc.0.wu", "m.inc.1.wu", "m.fl.0.wu", "m.fl.1.wu")),
        "models_wgt_using_cont.rds")


# ---- 4. CATEGORICAL IWISE (SLOW) --------------------------------------------
m.inc.0.cat.no.wgt <- fit_clmm(d.inc.cc, "INCOME_5",  "iwise_cat", re_term)       # income, unadjusted (intercept only, 07.2)
m.inc.1.cat.no.wgt <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_ctx_cat, re_slope_cat)  # income, adjusted (full)
m.fl.0.cat.no.wgt  <- fit_clmm(d.fl.cc,  "FLI_3item", "iwise_cat", re_slope_num)  # FLI, unadjusted (linear, 07.2)
m.fl.1.cat.no.wgt <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_ctx_cat, re_slope_cat)  # FLI, adjusted (full)

saveRDS(mget(c("m.inc.0.cat.no.wgt", "m.inc.1.cat.no.wgt", "m.fl.0.cat.no.wgt", "m.fl.1.cat.no.wgt")),
        "models_no_wgt_cat.rds")

list2env(readRDS("models_no_wgt_cont.rds"), envir = globalenv())
list2env(readRDS("models_no_wgt_cat.rds"),  envir = globalenv())









s <- summary(m.fl.1.wu)
s$condHess                                   # NaN or > 1e6 = problem
max(abs(m.fl.1.wu$gradient))                 # large = not converged
sum(is.na(s$coefficients[, 2]))              # must be 0
VarCorr(m.fl.1.wu)                           # SDs plausible? correlation not ±1?
summary(d.fl.cc$w_fit)                       # very large single weights?



# =============================================================================
# 07.6_models_weightless.R
# Main models refitted WITHOUT weights, on the same sample as the
# wgt_using models (Palestine 2025 dropped), so results are comparable.
#   Outcomes : INCOME_5 (d.inc), FLI_3item (d.fl)
#   IWISE    : continuous (iwisescore) and categorical (iwise_cat)
#   M0       : unadjusted          M1 : fully adjusted (rhs_ctx)
#   Random   : intercept + IWISE slope per country_year (as in 07.1 / 07.2)
#   Weights  : NONE
#   Fitting  : 8 models in parallel (future + furrr)
#   Saves: models_noweight_cont.rds, models_noweight_cat.rds
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
library(tidyverse)
library(ordinal)                                           # clmm()
install.packages("future")                                            # parallel backend
install.packages("furrr")                                             # future_map()
library(future)
library(furrr)

# ---- 1. DATA ----------------------------------------------------------------
list2env(readRDS("data/iwise_data_prepared.rds"), envir = globalenv())  # d.fl, d.inc, rhs_ctx, v.*

prep <- function(data, y) {                                # analytic sample for one outcome
  data %>%
    filter(cc, !is.na(wgt_using)) %>%                      # same sample as weighted models
    mutate(country_year = factor(.data[[v.cluster]]),      # grouping factor
           iwise_cat    = factor(iwise_cat, ordered = FALSE),  # one dummy per category
           iwise_num    = as.integer(iwise_cat) - 1,       # 0,1,2,3: linear random slope only
           "{y}" := factor(.data[[y]], ordered = TRUE))    # clmm() needs an ordered outcome
}

d.inc.cc <- prep(d.inc, "INCOME_5")
d.fl.cc  <- prep(d.fl,  "FLI_3item")

levels(d.fl.cc$iwise_cat)                                  # first level = lowest WI (reference)
table(d.fl.cc$iwise_cat, d.fl.cc$iwise_num)                # diagonal, no sparse categories


# ---- 2. FORMULAS AND HELPERS ------------------------------------------------
rhs_ctx_cat <- str_replace(rhs_ctx, "\\biwisescore\\b", "iwise_cat")
rhs_ctx_cat                                                # check: iwise_cat in, iwisescore out

re_slope     <- "(1 + iwisescore | country_year)"          # continuous slope
re_slope_cat <- "(1 + iwise_cat | country_year)"           # full: random slope per category
re_slope_num <- "(1 + iwise_num | country_year)"           # linear: one slope across categories


make_f <- function(y, rhs, re) as.formula(paste(y, "~", rhs, "+", re))

fit_clmm <- function(data, y, rhs, re) {                   # same settings, no weights
  clmm(make_f(y, rhs, re), data = data, link = "logit", Hess = TRUE, nAGQ = 1)
}


# ---- 3. MODEL SPECIFICATIONS ------------------------------------------------
# Each spec carries its own data, so a worker only receives the dataset it needs.
specs <- list(
  inc_cont_M0 = list(data = d.inc.cc, y = "INCOME_5",  rhs = "iwisescore", re = re_slope),
  inc_cont_M1 = list(data = d.inc.cc, y = "INCOME_5",  rhs = rhs_ctx,      re = re_slope),
  fl_cont_M0  = list(data = d.fl.cc,  y = "FLI_3item", rhs = "iwisescore", re = re_slope),
  fl_cont_M1  = list(data = d.fl.cc,  y = "FLI_3item", rhs = rhs_ctx,      re = re_slope),
  inc_cat_M0  = list(data = d.inc.cc, y = "INCOME_5",  rhs = "iwise_cat",  re = re_slope_cat),
  inc_cat_M1  = list(data = d.inc.cc, y = "INCOME_5",  rhs = rhs_ctx_cat,  re = re_slope_cat),
  fl_cat_M0   = list(data = d.fl.cc,  y = "FLI_3item", rhs = "iwise_cat",  re = re_slope_cat),
  fl_cat_M1   = list(data = d.fl.cc,  y = "FLI_3item", rhs = rhs_ctx_cat,  re = re_slope_cat)
)


# ---- 4. FIT IN PARALLEL -----------------------------------------------------
dir.create("models/noweight", recursive = TRUE, showWarnings = FALSE)

n_workers <- min(length(specs), future::availableCores() - 1)
plan(multisession, workers = n_workers)

t0 <- Sys.time()
fits <- future_map2(
  specs, names(specs),
  \(s, name) {
    fit <- tryCatch(fit_clmm(s$data, s$y, s$rhs, s$re),
                    error = \(e) e)                        # one failure doesn't stop the rest
    saveRDS(fit, file.path("models/noweight", paste0(name, ".rds")))  # backup per model
    fit
  },
  .options = furrr_options(packages = "ordinal",
                           seed = TRUE,
                           scheduling = Inf)               # one model per task: better load balance
)
Sys.time() - t0

plan(sequential)                                           # close workers


# ---- 5. CHECK ---------------------------------------------------------------
failed <- names(fits)[map_lgl(fits, \(m) inherits(m, "error"))]
failed                                                     # should be character(0)
map(fits[failed], conditionMessage)

map_dbl(fits[setdiff(names(fits), failed)], \(m) m$cond.H) # condition number of Hessian (< 1e6 ok)


# ---- 6. SAVE ----------------------------------------------------------------
models_cont <- fits[str_detect(names(fits), "_cont_")]
models_cat  <- fits[str_detect(names(fits), "_cat_")]

saveRDS(models_cont, "models_noweight_cont.rds")
saveRDS(models_cat,  "models_noweight_cat.rds")
