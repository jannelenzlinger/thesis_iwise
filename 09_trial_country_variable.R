# =============================================================================
# 07.7_models_slope_only.R
# Main models refitted with a random IWISE slope only (no random intercept).
#   Outcomes : INCOME_5 (d.inc), FLI_3item (d.fl)
#   IWISE    : continuous (iwisescore)
#   M0       : unadjusted          M1 : fully adjusted (rhs_ctx)
#   Fixed    : country_year added as fixed effect (toggle: fe_country)
#   Random   : IWISE slope per country_year, NO random intercept
#   Sample   : Palestine 2025 excluded explicitly
#   Weights  : fitted twice: wgt_using and unweighted
#   Fitting  : 8 models in parallel (future + furrr)
#   Saves: models_slope_only_wgt.rds, models_slope_only_noweight.rds
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
library(tidyverse)
library(ordinal)                                           # clmm()
library(future)                                            # parallel backend
library(furrr)                                             # future_map()


# ---- 1. DATA ----------------------------------------------------------------
list2env(readRDS("data/iwise_data_prepared.rds"), envir = globalenv())  # d.fl, d.inc, rhs_ctx, v.*

is_pal25 <- function(x) {                                  # flags Palestine 2025 country-waves
  x <- as.character(x)
  str_detect(x, regex("palestin", ignore_case = TRUE)) & str_detect(x, "2025")
}

prep <- function(data, y) {                                # analytic sample for one outcome
  data %>%
    mutate(country_year = factor(.data[[v.cluster]])) %>%  # grouping / fixed-effect factor
    filter(cc,
           !is_pal25(country_year),                        # explicit: drop Palestine 2025
           !is.na(wgt_using)) %>%                          # same sample for both weightings
    mutate(country_year = droplevels(country_year),        # remove empty levels
           "{y}"        := factor(.data[[y]], ordered = TRUE),  # clmm() needs an ordered outcome
           w_fit        = wgt_using)                       # not rescaled
}

d.inc.cc <- prep(d.inc, "INCOME_5")
d.fl.cc  <- prep(d.fl,  "FLI_3item")

# Check: Palestine 2022 still in, 2025 out
str_subset(levels(d.inc.cc$country_year), regex("palestin", ignore_case = TRUE))
str_subset(levels(d.fl.cc$country_year),  regex("palestin", ignore_case = TRUE))
stopifnot(!any(is_pal25(d.inc.cc$country_year)),
          !any(is_pal25(d.fl.cc$country_year)))


# ---- 2. FORMULAS AND HELPERS ------------------------------------------------
fe_country <- TRUE                                         # FALSE = country_year only as grouping

re_slope_only <- "(0 + iwisescore | country_year)"         # random slope, no random intercept

rhs_M0 <- "iwisescore"
rhs_M1 <- rhs_ctx                                          # contains iwisescore
if (fe_country) {
  rhs_M0 <- paste(rhs_M0, "+ country_year")
  rhs_M1 <- paste(rhs_M1, "+ country_year")
}
rhs_M0; rhs_M1                                             # check

make_f <- function(y, rhs, re) as.formula(paste(y, "~", rhs, "+", re))

fit_clmm <- function(data, y, rhs, re, weighted) {         # same settings as 07.1 / 07.2
  if (weighted) {
    clmm(make_f(y, rhs, re), data = data, weights = w_fit,
         link = "logit", Hess = TRUE, nAGQ = 1)
  } else {
    clmm(make_f(y, rhs, re), data = data,
         link = "logit", Hess = TRUE, nAGQ = 1)
  }
}


# ---- 3. MODEL SPECIFICATIONS ------------------------------------------------
# 2 outcomes x 2 adjustment levels x 2 weightings = 8 models
grid <- expand_grid(outcome  = c("inc", "fl"),
                    model    = c("M0", "M1"),
                    weighted = c(TRUE, FALSE))

specs <- pmap(grid, \(outcome, model, weighted) list(
  data     = if (outcome == "inc") d.inc.cc else d.fl.cc,
  y        = if (outcome == "inc") "INCOME_5" else "FLI_3item",
  rhs      = if (model == "M0") rhs_M0 else rhs_M1,
  re       = re_slope_only,
  weighted = weighted
))
names(specs) <- with(grid, paste(outcome, model, if_else(weighted, "wgt", "noweight"), sep = "_"))
names(specs)


# ---- 4. FIT IN PARALLEL -----------------------------------------------------
dir.create("models/slope_only", recursive = TRUE, showWarnings = FALSE)

n_workers <- min(length(specs), future::availableCores() - 1)  # lower if RAM runs out
plan(multisession, workers = n_workers)

t0 <- Sys.time()
fits <- future_map2(
  specs, names(specs),
  \(s, name) {
    fit <- tryCatch(fit_clmm(s$data, s$y, s$rhs, s$re, s$weighted),
                    error = \(e) e)                        # one failure doesn't stop the rest
    saveRDS(fit, file.path("models/slope_only", paste0(name, ".rds")))  # backup per model
    fit
  },
  .options = furrr_options(packages = "ordinal",
                           seed = TRUE,
                           scheduling = Inf)               # one model per task
)
Sys.time() - t0

plan(sequential)                                           # close workers


# ---- 5. CHECK ---------------------------------------------------------------
failed <- names(fits)[map_lgl(fits, \(m) inherits(m, "error"))]
failed                                                     # should be character(0)
map(fits[failed], conditionMessage)

ok <- setdiff(names(fits), failed)
map_dbl(fits[ok], \(m) m$cond.H)                           # condition number (< 1e6 ok)
map_dbl(fits[ok], \(m) coef(m)["iwisescore"])              # IWISE fixed effect, wgt vs noweight


vars <- all.vars(as.formula(paste("~", rhs_ctx)))

d.inc.cc %>%
  group_by(country_year) %>%
  summarise(across(all_of(vars), n_distinct)) %>%
  summarise(across(-country_year, max)) %>%
  pivot_longer(everything(), names_to = "var", values_to = "max_distinct") %>%
  filter(max_distinct == 1)          # these would be absorbed



# ---- 6. SAVE ----------------------------------------------------------------
saveRDS(fits[str_detect(names(fits), "_wgt$")],      "models_slope_only_wgt.rds")
saveRDS(fits[str_detect(names(fits), "_noweight$")], "models_slope_only_noweight.rds")