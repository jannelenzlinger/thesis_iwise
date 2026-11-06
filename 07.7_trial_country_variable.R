# =============================================================================
# 07.5_fl_model_options.R
# Sensitivity of the fully adjusted FLI model (M1) to the country-level part
#   Outcome:   FLI_3item        Predictor: iwisescore (continuous)
#   Weights:   wgt_using, rescaled to mean 1 (w_fit)
#
#   Models (all on the SAME complete-case sample, so logLik/AIC are comparable):
#     base    M1 as before:  all country covariates + (1 + iwisescore | country_year)
#     no_gdp  1. without log_gdp_c
#     no_wgi  2. without wgi_ge_sc_c
#     no_ri   3. no random intercept, GDP + WGI kept:     (0 + iwisescore | country_year)
#     fe_cy   4. no random intercept, no country covariates,
#                country_year as a fixed factor:          (0 + iwisescore | country_year)
#
#   The random IWISE slope is kept in all models (Heisig & Schaeffer 2019).
#   Set re = NULL in a spec to fit it without any random effect (clm()).
#
#   Runs the models in parallel (future + furrr), one model per CPU core.
#   Saves:  models_fl_options.rds   (named list of fits)
#           fl_options_compare.csv  (comparison table)
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
library(dplyr)
library(ordinal)
library(future)
library(furrr)


# ---- 1. DATA ----------------------------------------------------------------
# Load your data the same way as in 07.x. You need `d.fl.cc`: the FLI analysis
# sample, with country_year and wgt_using.
list2env(readRDS("data/iwise_data_prepared.rds"), envir = globalenv())


# ---- 2. VARIABLES -----------------------------------------------------------
# CHECK: these must match `rhs_ctx` in your M1 (names taken from your output).
y       <- "FLI_3item"
x_iwise <- "iwisescore"
x_indiv <- c("female", "urban_imp", "age_gp_profile", "maritalstatus",
             "hhsize_c", "employment", "education", "INDEX_FS")
x_ctry  <- c("log_gdp_c", "wgi_ge_sc_c", "bas_rate_pp", "daly_enteric_model_c")
g       <- "country_year"

re_full  <- "(1 + iwisescore | country_year)"   # random intercept + slope
re_slope <- "(0 + iwisescore | country_year)"   # random slope only


# ---- 3. ONE COMMON SAMPLE ---------------------------------------------------
# Complete cases over every variable used in any model, so all models are fitted
# on the same rows.
d.opt <- d.fl.cc %>%
  select(all_of(c(y, x_iwise, x_indiv, x_ctry, g, "wgt_using"))) %>%
  filter(if_all(everything(), ~ !is.na(.x))) %>%
  mutate(country_year = factor(country_year),
         w_fit = wgt_using / mean(wgt_using))     # mean-1 weights

cat("N =", nrow(d.opt), "| country-waves =", nlevels(d.opt$country_year), "\n")

# Variables that are constant within every country-wave. Model 4 has to drop
# these: the country_year dummies absorb them.
is_ctry_level <- function(v) {
  all(tapply(d.opt[[v]], d.opt[[g]], function(x) length(unique(x)) == 1))
}
ctry_level <- Filter(is_ctry_level, c(x_indiv, x_ctry))
cat("Constant within country_year (dropped in fe_cy):",
    paste(ctry_level, collapse = ", "), "\n")


# ---- 4. MODEL SPECIFICATIONS ------------------------------------------------
# ctry = country covariates in the model, re = random part (NULL = none),
# fe = country_year as a fixed factor.
specs <- list(
  base   = list(ctry = x_ctry,                          re = re_full,  fe = FALSE),
  no_gdp = list(ctry = setdiff(x_ctry, "log_gdp_c"),    re = re_full,  fe = FALSE),
  no_wgi = list(ctry = setdiff(x_ctry, "wgi_ge_sc_c"),  re = re_full,  fe = FALSE),
  no_ri  = list(ctry = x_ctry,                          re = re_slope, fe = FALSE),
  fe_cy  = list(ctry = character(0),                    re = re_slope, fe = TRUE)
)

make_formula <- function(s) {
  indiv <- if (s$fe) setdiff(x_indiv, ctry_level) else x_indiv
  rhs   <- c(x_iwise, indiv, s$ctry, if (s$fe) g, s$re)
  as.formula(paste(y, "~", paste(rhs, collapse = " + ")))
}

for (nm in names(specs)) cat(sprintf("%-7s %s\n", nm, deparse1(make_formula(specs[[nm]]))))


# ---- 5. FIT FUNCTION --------------------------------------------------------
ctrl <- clmm.control(method = "nlminb", maxIter = 200, maxLineIter = 200,
                     iter.max = 1000, eval.max = 2000)   # last two go to nlminb

fit_one <- function(s, data) {
  f  <- make_formula(s)
  t0 <- Sys.time()
  m  <- tryCatch(
    if (is.null(s$re)) {
      clm(f, data = data, weights = w_fit, link = "logit", Hess = TRUE)
    } else {
      clmm(f, data = data, weights = w_fit, link = "logit", Hess = TRUE,
           nAGQ = 1, control = ctrl)
    },
    error = function(e) e)                         # one failure doesn't stop the rest
  attr(m, "runtime_min") <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  m
}


# ---- 6. RUN IN PARALLEL -----------------------------------------------------
# Each worker gets its own copy of d.opt, so it needs its own share of RAM.
n_workers <- min(length(specs), parallelly::availableCores(omit = 1))  # leave 1 core free
cat("Fitting", length(specs), "models on", n_workers, "worker(s)\n")
plan(multisession, workers = n_workers)

fits <- future_map(specs, fit_one, data = d.opt,
                   .options = furrr_options(seed = TRUE, packages = "ordinal"))

plan(sequential)
saveRDS(fits, "models_fl_options.rds")
# fits <- readRDS("models_fl_options.rds")         # reload later without refitting
list2env(setNames(fits, paste0("m.fl.opt.", names(fits))), envir = globalenv())0
# gives m.fl.opt.base, m.fl.opt.no_gdp, m.fl.opt.no_wgi, m.fl.opt.no_ri, m.fl.opt.fe_cy

# ---- 7. CONVERGENCE CHECK ---------------------------------------------------
# rel = |Newton step| / SE per parameter. Below ~0.001 = converged.
conv_check <- function(m) {
  H <- m$Hessian; gr <- m$gradient
  rel <- tryCatch(abs(solve(H, gr)) / sqrt(diag(solve(H))),
                  error = function(e) rep(NA_real_, length(gr)))
  names(rel) <- names(gr)
  rel
}


# ---- 8. COMPARISON TABLE ----------------------------------------------------
summarise_fit <- function(m, name) {
  if (inherits(m, "error")) {
    return(tibble(model = name, error = conditionMessage(m)))
  }
  cf  <- summary(m)$coefficients
  b   <- cf[x_iwise, "Estimate"]; se <- cf[x_iwise, "Std. Error"]
  rel <- conv_check(m)
  sds <- if (inherits(m, "clmm")) attr(VarCorr(m)$country_year, "stddev") else c()
  ll  <- logLik(m)
  tibble(
    model        = name,
    n            = nobs(m),
    n_par        = attr(ll, "df"),
    logLik       = as.numeric(ll),
    AIC          = AIC(m),
    b_iwise      = b,
    se_iwise     = se,
    OR_iwise     = exp(b),
    OR_lo        = exp(b - 1.96 * se),
    OR_hi        = exp(b + 1.96 * se),
    sd_intercept = unname(sds["(Intercept)"]),     # NA = no random intercept
    sd_slope     = unname(sds[x_iwise]),
    condHess     = summary(m)$condHess,
    max_grad     = max(abs(m$gradient)),
    rel_iwise    = unname(rel[x_iwise]),
    rel_max      = max(rel, na.rm = TRUE),
    worst_par    = paste(names(sort(rel, decreasing = TRUE))[1:3], collapse = ", "),
    runtime_min  = attr(m, "runtime_min"),
    error        = NA_character_
  )
}

compare <- bind_rows(Map(summarise_fit, fits, names(fits))) %>%
  mutate(dAIC = AIC - min(AIC, na.rm = TRUE), .after = AIC)

print(as.data.frame(compare), digits = 4)
write.csv(compare, "fl_options_compare.csv", row.names = FALSE)

m.fl.1.wu$optRes$message          # old model
fits$base$optRes$message          # new model
m.fl.1.wu$control                 # limits the old model used
c(old = nobs(m.fl.1.wu), new = nobs(fits$base))   # same sample?

cov2cor(vcov(fits$base))[c("log_gdp_c", "wgi_ge_sc_c")]


# ---- 9. HOW TO READ THE TABLE -----------------------------------------------
# rel_max   < 0.001       converged. worst_par names what still moves.
# OR_iwise  similar in all models -> the IWISE result doesn't depend on how
#           country differences are handled. Report this as a robustness check.
# dAIC      0 = best fit; > 10 = clearly worse. Same rows in every model, so
#           comparable. fe_cy has many more parameters (n_par).
# sd_slope  should be similar across models; it is estimated in all of them.
# no_ri     without a random intercept, GDP/WGI/etc. must explain all country
#           differences in level. Expect a much worse logLik and too-small SEs.