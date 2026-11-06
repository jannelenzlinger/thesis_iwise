# =============================================================================
# 08.1_final_models.R
# Final cumulative link mixed models: IWISE (predictor) -> FLI and income quintile
#
#   m.fl.0   FLI_3item ~ iwisescore                                  (unadjusted)
#   m.fl.1   FLI_3item ~ iwisescore + individual covariates
#                        + GDP, JMP basic water, enteric DALYs        (fully adjusted, no WGI)
#   m.inc.0  INCOME_5  ~ iwisescore                                  (unadjusted)
#   m.inc.1  INCOME_5  ~ iwisescore + individual covariates + GDP, JMP, DALYs
#
#   Random part: random intercept + random IWISE slope per country_year,  # >>> CHANGED
#                except m.inc.0: random intercept only (in earlier runs  # >>> CHANGED
#                the random slope was not identified for this model).    # >>> CHANGED
#   Weights:     w_fit from 05 = wgtnorm_using rescaled to sum to the    # >>> CHANGED
#                complete-case N (each country-year weighted equally).   # >>> CHANGED
#   Link: logit. Laplace approximation (nAGQ = 1). Raised iteration
#   limits (ctrl, section 4).
#   Unadjusted and adjusted models of one outcome use the SAME rows.
#
#   Saves: models_final.rds        (named list of all fits)
#          fit_<name>.rds          (each model as soon as it finishes)
#          fit_log.txt             (progress log: open it while the script runs)
#          final_models_check.csv  (convergence + IWISE OR table)
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
library(dplyr)                                  # select(), filter(), mutate(), tibble(), bind_rows()
library(ordinal)                                # clmm(): cumulative link mixed models
library(future)                                 # parallel backend: plan(multisession)
library(furrr)                                  # future_imap(): purrr-style loop run in parallel


# ---- 1. DATA ----------------------------------------------------------------
# Loads every object saved by 05 into the workspace. Needed here: d.fl and
# d.inc, including the columns cc and w_fit (w_fit is NA where cc is FALSE).
list2env(readRDS("data/iwise_data_prepared.rds"), envir = globalenv())


# ---- 2. VARIABLES -----------------------------------------------------------
# CHECK: x_indiv must match the individual covariates of your fully adjusted model.
x_iwise <- "iwisescore"                         # predictor: continuous IWISE score (0-36)
x_indiv <- c("female", "urban_imp", "age_gp_profile", "maritalstatus",   # individual covariates
             "hhsize_c", "employment", "education", "INDEX_FS")          # (hhsize centred)
x_ctry  <- c("log_gdp_c", "bas_rate_pp", "daly_enteric_model_c")   # country-year covariates: log GDP
                                                                   # (centred), JMP basic water,
                                                                   # enteric DALYs (centred). WGI dropped
g       <- "country_year"                       # grouping variable for the random effects
w       <- "w_fit"                              # weight column, kept in the sample (section 3)

re_slope <- "(1 + iwisescore | country_year)"   # random intercept + random IWISE slope, correlated
re_int   <- "(1 | country_year)"                # random intercept only (used for m.inc.0)


# ---- 3. ANALYSIS SAMPLES ----------------------------------------------------
# The sample is defined in 05 (cc, w_fit). Here: keep it, fix the types.
prep <- function(data, y) {
  data %>%
    filter(cc) %>%                                          # 05 complete cases
    mutate(country_year = factor(country_year),             # grouping must be a factor
           across(all_of(y), ~ factor(.x, ordered = TRUE))) # ordered outcome; keeps level order
}
d.fl  <- prep(d.fl,  "FLI_3item")
d.inc <- prep(d.inc, "INCOME_5")

stopifnot(all.equal(sum(d.fl$w_fit),  nrow(d.fl)),          # 05 w_fit sums to N
          all.equal(sum(d.inc$w_fit), nrow(d.inc)))

# check that w_fit is the 05 version and sums to N.
# sum_w should equal n. If it doesn't, this sample differs from the 05
# complete cases: then uncomment the two rescaling lines below.
c(n = nrow(d.fl),  sum_w = sum(d.fl$w_fit))     # FLI
c(n = nrow(d.inc), sum_w = sum(d.inc$w_fit))    # income


# ---- 4. MODEL SPECIFICATIONS AND SETTINGS -----------------------------------
# One list entry per model: outcome (y), which dataset, fixed part (rhs), random part (re).
specs <- list(
  m.fl.0  = list(y = "FLI_3item", data = "d.fl",  rhs = x_iwise,                     re = re_slope),
  m.fl.1  = list(y = "FLI_3item", data = "d.fl",  rhs = c(x_iwise, x_indiv, x_ctry), re = re_slope),
  m.inc.0 = list(y = "INCOME_5",  data = "d.inc", rhs = x_iwise,                     re = re_int),   # random intercept only
  m.inc.1 = list(y = "INCOME_5",  data = "d.inc", rhs = c(x_iwise, x_indiv, x_ctry), re = re_slope)
)

# Builds "y ~ x1 + x2 + ... + (random part)" from a spec and turns the text into a formula.
make_formula <- function(s) as.formula(paste(s$y, "~", paste(c(s$rhs, s$re), collapse = " + ")))
# Prints every formula, so you can check them before the long fit.
for (nm in names(specs)) cat(sprintf("%-11s %s\n", nm, deparse1(make_formula(specs[[nm]]))))

# Optimiser settings. Outer loop (nlminb, all parameters): iter.max = max. steps,
# eval.max = max. likelihood evaluations (defaults 150 / 200 were too low).
# Inner loop (random effects per country_year): maxIter = max. Newton steps,
# maxLineIter = max. step halvings (defaults 50 / 50).
ctrl <- clmm.control(method = "nlminb", maxIter = 200, maxLineIter = 200,
                     iter.max = 1000, eval.max = 2000)


# ---- 5. FIT FUNCTION --------------------------------------------------------
# Fits one model, times it, saves it, writes a line to the log, returns the fit.
fit_one <- function(s, name, datasets) {
  t0 <- Sys.time()                                          # start time
  m  <- tryCatch(
    clmm(make_formula(s), data = datasets[[s$data]],        # formula + the right dataset
         weights = w_fit,                                   # sample size adjusted wgtnorm_using
                                                            # in the data after select() in section 3
         link = "logit",                                    # cumulative logit -> odds ratios
         Hess = TRUE,                                       # Hessian needed for SEs and conv_check()
         nAGQ = 1,                                          # Laplace approximation (needed for random slopes)
         control = ctrl),                                   # raised iteration limits from section 4
    error = function(e) e)                                  # one failure doesn't stop the rest
  attr(m, "runtime_min") <- as.numeric(difftime(Sys.time(), t0, units = "mins"))  # runtime in minutes
  saveRDS(m, paste0("fit_", name, ".rds"))                  # saved straight away (safe if R crashes)
  cat(format(Sys.time(), "%H:%M"), name, "done after", round(attr(m, "runtime_min")), "min",
      if (inherits(m, "error")) paste("- ERROR:", conditionMessage(m)) else "", "\n",
      file = "fit_log.txt", append = TRUE)                  # one log line per finished model
  m
}


# ---- 6. RUN IN PARALLEL -----------------------------------------------------
# One model per core. Each worker holds a copy of the data, so watch your RAM.
n_workers <- min(length(specs), parallelly::availableCores(omit = 1))   # 4 models, keep 1 core free
cat("Fitting", length(specs), "models on", n_workers, "worker(s). Progress: fit_log.txt\n")
plan(multisession, workers = n_workers)                     # start background R sessions

fits <- future_imap(specs, fit_one, datasets = list(d.fl = d.fl, d.inc = d.inc),  # fit_one(spec, name, datasets)
                    .options = furrr_options(seed = TRUE, packages = "ordinal"))  # load ordinal in each worker

plan(sequential)                                            # close the workers again
saveRDS(fits, "models_final.rds")                           # all four fits in one file
# fits <- readRDS("models_final.rds")                       # reload later without refitting
# list2env(fits, envir = globalenv())                       # -> m.fl.0, m.fl.1, m.inc.0, m.inc.1


# ---- 7. CONVERGENCE AND RESULTS ---------------------------------------------
# rel = |Newton step| / SE per parameter: how far the optimiser would still move
# each parameter, in units of its SE. Converged: rel_max < 0.001 and the
# optimiser message says "relative convergence" or "X-convergence".
conv_check <- function(m) {
  H <- m$Hessian; gr <- m$gradient                          # curvature and slope of the likelihood
  rel <- tryCatch(abs(solve(H, gr)) / sqrt(diag(solve(H))), # Newton step / SE
                  error = function(e) rep(NA_real_, length(gr)))  # NA if the Hessian can't be inverted
  setNames(rel, names(gr))                                  # name each value by its parameter
}

# One row of results and checks per model.
summarise_fit <- function(m, name) {
  if (inherits(m, "error")) return(tibble(model = name, error = conditionMessage(m)))  # failed fit
  cf  <- summary(m)$coefficients                            # estimates, SEs, z, p
  b   <- cf[x_iwise, "Estimate"]; se <- cf[x_iwise, "Std. Error"]   # IWISE log-odds and SE
  rel <- conv_check(m)
  sds <- attr(VarCorr(m)$country_year, "stddev")            # SDs of random intercept (and slope)
  tibble(
    model        = name,
    n_people     = nrow(m$model),                           # actual sample size
    sum_weights  = nobs(m),                                 # ordinal's nobs() = sum of weights; = n_people with w_fit
    b_iwise      = b,                                       # log-odds per IWISE point
    se_iwise     = se,
    OR_iwise     = exp(b),                                  # odds ratio per point
    OR_lo        = exp(b - 1.96 * se),                      # 95% CI (Wald)
    OR_hi        = exp(b + 1.96 * se),
    p_iwise      = cf[x_iwise, "Pr(>|z|)"],
    sd_intercept = unname(sds["(Intercept)"]),              # between-country-year SD of the intercept
    sd_slope     = unname(sds[x_iwise]),                    # SD of the IWISE slope; NA for m.inc.0  # >>> CHANGED
    logLik       = as.numeric(logLik(m)),                   # only comparable between models with the same rows + weights
    AIC          = AIC(m),
    optimiser    = m$optRes$message,                        # nlminb's stop message
    rel_iwise    = unname(rel[x_iwise]),                    # convergence of the IWISE estimate itself
    rel_max      = max(rel, na.rm = TRUE),                  # worst parameter
    worst_par    = paste(names(sort(rel, decreasing = TRUE))[1:3], collapse = ", "),  # 3 least converged
    condHess     = summary(m)$condHess,                     # Hessian condition number; < ~1e5 fine
    max_grad     = max(abs(m$gradient)),                    # largest gradient; close to 0 = converged
    converged    = rel_max < 0.001 & grepl("convergence", optimiser) & !grepl("false", optimiser),
                                                            # TRUE only if both criteria hold
    runtime_min  = attr(m, "runtime_min"),
    error        = NA_character_
  )
}

check <- bind_rows(Map(summarise_fit, fits, names(fits)))   # one row per model
print(as.data.frame(check), digits = 4)                     # print all columns
write.csv(check, "final_models_check.csv", row.names = FALSE)