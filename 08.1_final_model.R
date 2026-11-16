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
#   Random part: random intercept + random IWISE slope per country_year, 
#                except m.inc.0: random intercept only (in earlier runs 
#                the random slope was not identified for this model). 
#   Weights:     w_fit from 05 = wgtnorm_using rescaled to sum to the  
#                complete-case N (each country-year weighted equally). 
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

# ---- 8. PROPORTIONAL ODDS CHECK ---------------------------------------------
# clmm assumes ONE coefficient per predictor for all cuts of the outcome.
# Check: fit a separate binary model for each cut (Y > k) with the same fixed
# and random part, and compare its coefficients with the clmm estimate.
# Flat estimates around the red dashed line = PO holds; a trend or CIs that
# miss the line = violation. Signs match: positive = higher category in both.
#
# Starting here without refitting? Run sections 0-4, then:
# fits <- readRDS("models_final.rds")
 
library(lme4)                                   # glmer(): binary mixed logit per cut
library(purrr)                                  # map(), list_rbind()
library(ggplot2)                                # figure
 
po_models <- c("m.fl.1", "m.inc.1")             # fully adjusted models to check
cut_re    <- "same"                             # "same": glmer with the model's random part (comparable to clmm, slow)
                                                # "none": plain glm per cut (fast, rough first look)
ctrl_glmer <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
 
# One job per model x cut
jobs <- list_rbind(map(po_models, function(nm) {
  s <- specs[[nm]]
  K <- nlevels(get(s$data)[[s$y]])              # number of outcome categories
  tibble(model = nm, k = 1:(K - 1))
}))
jobs <- split(jobs, seq_len(nrow(jobs)))        # list of one-row tibbles for future_map()
 
# Fits the binary model for one cut, returns its coefficients
fit_cut <- function(job, specs, datasets, cut_re, ctrl_glmer) {
  s <- specs[[job$model]]
  d <- datasets[[s$data]]
  d$y_bin <- as.integer(as.integer(d[[s$y]]) > job$k)        # 1 = above category k
  fixed   <- paste("y_bin ~", paste(s$rhs, collapse = " + "))
  warn    <- character(0)
  m <- withCallingHandlers(
    tryCatch(
      if (cut_re == "same") {
        lme4::glmer(as.formula(paste(fixed, "+", s$re)), data = d, weights = w_fit,
                    family = binomial, nAGQ = 1, control = ctrl_glmer)
      } else {
        glm(as.formula(fixed), data = d, weights = w_fit, family = quasibinomial())
      },
      error = function(e) e),
    warning = function(w) {                     # non-integer weights warning is expected: drop it
      if (!grepl("non-integer", conditionMessage(w))) warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (inherits(m, "error"))
    return(tibble::tibble(model = job$model, cut = job$k, error = conditionMessage(m)))
  cf <- summary(m)$coefficients
  tibble::tibble(model = job$model, cut = job$k, term = rownames(cf),
                 est = cf[, 1], se = cf[, 2],
                 warning = if (length(warn)) paste(warn, collapse = "; ") else NA_character_,
                 error = NA_character_)
}
 
# Run in parallel (one cut per core; watch RAM as in section 6)
plan(multisession, workers = min(length(jobs), parallelly::availableCores(omit = 1)))
po_cuts <- future_map(jobs, fit_cut, specs = specs, datasets = list(d.fl = d.fl, d.inc = d.inc),
                      cut_re = cut_re, ctrl_glmer = ctrl_glmer,
                      .options = furrr_options(seed = TRUE, packages = c("lme4", "tibble")))
plan(sequential)
po_cuts <- list_rbind(po_cuts)
saveRDS(po_cuts, "po_cuts.rds")                 # reload later: po_cuts <- readRDS("po_cuts.rds")
 
# Failed fits or convergence warnings (check before reading the figure)
po_cuts %>% filter(!is.na(error) | !is.na(warning)) %>% distinct(model, cut, warning, error)
 
# clmm (PO) estimates: fixed effects only, thresholds dropped
po_ref <- list_rbind(map(po_models, function(nm) {
  b <- fits[[nm]]$beta
  tibble(model = nm, term = names(b), est_po = unname(b))
}))
 
po_cuts <- po_cuts %>%
  filter(is.na(error), term != "(Intercept)") %>%
  mutate(lo = est - 1.96 * se, hi = est + 1.96 * se) %>%
  left_join(po_ref, by = c("model", "term"))
 
# Figure per model
plot_po <- function(nm, title) {
  ggplot(filter(po_cuts, model == nm), aes(cut, est)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_hline(aes(yintercept = est_po), colour = "red", linetype = "dashed") +
    geom_pointrange(aes(ymin = lo, ymax = hi)) +
    facet_wrap(~ term, scales = "free_y") +
    scale_x_continuous(breaks = function(x) seq(1, floor(max(x)))) +
    labs(title = paste0(title, ": coefficients across cuts"),
         subtitle = "Points: binary model per cut (95% CI). Red dashed: clmm estimate.",
         x = "Cut (above category k)", y = "Log-odds (95% CI)") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(size = 8))
}
p_po_fl  <- plot_po("m.fl.1",  "FLI")
p_po_inc <- plot_po("m.inc.1", "Income")
p_po_fl
p_po_fl; p_po_inc
ggsave("po_check_fli.png",    p_po_fl,  width = 12, height = 8, dpi = 300)
ggsave("po_check_income.png", p_po_inc, width = 12, height = 8, dpi = 300)
 
# Summary per predictor: spread across cuts, cuts whose CI misses the clmm
# estimate, and sign changes. Sorted worst first.
po_summary <- po_cuts %>%
  group_by(model, term) %>%
  summarise(est_po      = first(est_po),
            min_cut     = min(est),
            max_cut     = max(est),
            spread      = max_cut - min_cut,
            n_cuts_off  = sum(est_po < lo | est_po > hi),
            sign_change = any(sign(est) != sign(first(est_po))),
            .groups = "drop") %>%
  arrange(model, desc(n_cuts_off), desc(spread))
print(as.data.frame(po_summary), digits = 3)
write.csv(po_summary, "po_check_summary.csv", row.names = FALSE)

