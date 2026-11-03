# =============================================================================
# 08_rf_comparison.R
# Robustness check: random forest (RF) vs multilevel ordinal logit (clmm)
#   Outcomes : FLI_3item, then INCOME_5 (as in 07)
#   Question : which model predicts the outcome better for people it has
#              NOT seen? Same predictors (rhs_ctx), same folds, same weights.
#   CV       : (A) random 5-fold          = new people, known country-years
#              (B) grouped 5-fold by       = new country-years
#                  country_year
#   Models   : null (category shares), clmm, RF, RF without IWISE
#   Metrics  : RPS (main), accuracy, mean absolute error (in categories)
# Needs iwise_data_prepared.rds from 05 (same as 07).
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
library(tidyverse)                   # data handling and ggplot2
library(ordinal)                     # clmm(): multilevel ordinal logit
library(ranger)                      # ranger(): fast random forest
library(furrr)

# ---- 1. SETTINGS (change here only) -----------------------------------------
v.wt    <- "wgt_using"           # weight for ALL models and metrics
k       <- 5                         # number of CV folds (3 = faster)
n.trees <- 500                       # trees per forest (500 is the usual default)
seed    <- 2026                      # one seed -> identical folds every run
re_cv   <- "(1 | country_year)"      # clmm random part: set to your final
                                     # structure, e.g. "(1 + iwisescore | country_year)"
n.workers <- max(1, parallel::detectCores() - 1)  # parallel R sessions for the CV (keep 1 core free)


# ---- 2. DATA ----------------------------------------------------------------
list2env(readRDS("iwise_data_prepared.rds"), envir = globalenv())  # restores d.fl, d.inc, rhs_ctx, v.* objects

x.vars      <- all.vars(as.formula(paste("~", rhs_ctx)))   # predictor names from rhs_ctx (same as clmm)
x.vars.noiw <- setdiff(x.vars, c("iwisescore", "iwise_mean"))  # same, without any IWISE information

prep_rf <- function(data, y) {                             # builds the analytic sample for one outcome
  data %>%
    filter(cc, !is.na(.data[[v.wt]])) %>%                  # complete cases + non-missing weight
    mutate(country_year = factor(.data[[v.cluster]]),      # grouping factor (clmm + grouped folds)
           "{y}" := factor(.data[[y]], ordered = TRUE),    # ordered outcome, lowest category first
           w = .data[[v.wt]]) %>%                          # weight under a short fixed name
    mutate(across(where(is.character), factor)) %>%        # ranger and clmm need factors, not text
    select(all_of(c(y, x.vars, "country_year", "w")))      # keep only what the models use
}

d.fl.rf  <- prep_rf(d.fl,  "FLI_3item")                    # FLI analytic sample
d.inc.rf <- prep_rf(d.inc, "INCOME_5")                     # income analytic sample

c(fl = nrow(d.fl.rf), inc = nrow(d.inc.rf))                # check: close to the N of the 07 models
summary(d.fl.rf$w)                                         # check: weights positive, no NA


# ---- 3. FOLDS ---------------------------------------------------------------
make_folds <- function(data, k, grouped) {                 # returns a fold number (1..k) per row
  set.seed(seed)                                           # same folds every run
  if (!grouped) {
    sample(rep(1:k, length.out = nrow(data)))              # (A) random: people spread over folds
  } else {
    g  <- levels(droplevels(data$country_year))            # all country-years
    fg <- setNames(sample(rep(1:k, length.out = length(g))), g)  # each country-year -> one fold
    unname(fg[as.character(data$country_year)])            # (B) a whole country-year is held out together
  }
}


# ---- 4. PREDICTION HELPERS --------------------------------------------------
# Every model returns an n x K matrix: one row per person, one column per
# outcome category, each row = predicted probabilities that sum to 1.

# 4a. clmm has no predict() method, so the probabilities are computed by hand:
#     P(Y <= j) = plogis(alpha_j - eta),   eta = X*beta + random effects
pred_clmm <- function(m, newdata, y, known_groups) {
  X <- model.matrix(as.formula(paste("~", rhs_ctx)),       # design matrix with the SAME coding
                    data = newdata,                        # as the fitted model (contrasts below)
                    contrasts.arg = m$contrasts)[, -1, drop = FALSE]  # drop intercept (thresholds replace it)
  b   <- m$beta                                            # fixed-effect coefficients
  eta <- drop(X[, names(b), drop = FALSE] %*% b)           # linear predictor, fixed part

  if (known_groups) {                                      # (A) country-year was in training data:
    re <- ranef(m)$country_year                            #     add its estimated random effects
    g  <- as.character(newdata$country_year)
    for (v in colnames(re)) {                              # intercept and (if used) IWISE slope
      xv  <- if (v == "(Intercept)") 1 else newdata[[v]]   # intercept multiplies 1, slope multiplies iwisescore
      u   <- re[g, v]                                      # deviation of each person's country-year
      u[is.na(u)] <- 0                                     # safety: unseen country-year -> 0
      eta <- eta + u * xv                                  # add the country-year deviation
    }
  }                                                        # (B) new country-year: random effect = 0 (its mean)

  cum <- plogis(outer(-eta, m$alpha, "+"))                 # cumulative probs P(Y <= j), j = 1..K-1
  p   <- cbind(cum, 1) - cbind(0, cum)                     # differences -> P(Y = j), j = 1..K
  colnames(p) <- levels(newdata[[y]])
  p
}

# 4b. RF predictions, columns put in category order
pred_rf <- function(rf, newdata, y) {
  p <- predict(rf, data = newdata)$predictions             # probability forest -> probability matrix
  p[, levels(newdata[[y]]), drop = FALSE]                  # same column order as the outcome levels
}

# 4c. Null model: everyone gets the (weighted) category shares of the training data
pred_null <- function(train, newdata, y) {
  sh <- tapply(train$w, train[[y]], sum)                   # weighted count per category
  sh[is.na(sh)] <- 0                                       # empty category -> 0
  matrix(sh / sum(sh), nrow(newdata), length(sh), byrow = TRUE,
         dimnames = list(NULL, levels(newdata[[y]])))
}


# ---- 5. METRICS -------------------------------------------------------------
# RPS = ranked probability score (0 = perfect, lower = better). It compares
#       the CUMULATIVE predicted and observed distributions, so a prediction
#       one category off is punished less than one four categories off.
# All metrics are weighted means over the test persons (weight w).
score <- function(p, obs, w) {
  K   <- ncol(p)
  o   <- as.integer(obs)                                   # observed category as 1..K
  O   <- outer(o, 1:K, "<=") * 1                           # observed cumulative (0/1)
  F   <- t(apply(p, 1, cumsum))                            # predicted cumulative
  rps <- rowSums((F[, -K, drop = FALSE] - O[, -K, drop = FALSE])^2) / (K - 1)  # RPS per person
  hat <- max.col(p, ties.method = "first")                 # most likely category
  tibble(rps = weighted.mean(rps, w),                      # main metric
         acc = weighted.mean(hat == o, w),                 # share predicted exactly right
         mae = weighted.mean(abs(hat - o), w))             # average distance in categories
}


# ---- 6. ONE FOLD ------------------------------------------------------------
run_fold <- function(data, y, test_id, grouped) {
  train <- data[!test_id, ]                                # model sees only these rows
  test  <- data[ test_id, ]                                # ...and is scored on these
  train$w_cv <- train$w / mean(train$w)                    # rescale to mean 1 (clmm reads weights as counts)

  # clmm: same model as 07, refitted on the training part only
  m <- clmm(as.formula(paste(y, "~", rhs_ctx, "+", re_cv)),
            data = train, weights = w_cv,                  # weight_using, rescaled
            link = "logit", Hess = FALSE, nAGQ = 1)        # no SEs needed -> Hess = FALSE (faster)

  # RF: probability forest, same predictors and weights
  fit_rf <- function(vars) {
    ranger(x = train[, vars], y = train[[y]],              # predictors and outcome
           probability  = TRUE,                            # predict probabilities, not just a class
           num.trees    = n.trees,
           case.weights = train$w_cv,                      # weighted sampling of persons
           seed         = seed,
           num.threads  = 1)                               # 1 core per forest: the folds already run in parallel
  }
  rf      <- fit_rf(x.vars)                                # all predictors
  rf.noiw <- fit_rf(x.vars.noiw)                           # without iwisescore and iwise_mean

  preds <- list(null       = pred_null(train, test, y),
                clmm       = pred_clmm(m, test, y, known_groups = !grouped),
                rf         = pred_rf(rf, test, y),
                rf_no_iwise = pred_rf(rf.noiw, test, y))

  imap_dfr(preds, ~ score(.x, test[[y]], test$w) %>% mutate(model = .y))  # one row per model
}


# ---- 7. RUN ALL FOLDS (IN PARALLEL) -----------------------------------------
# 2 outcomes x 2 schemes x k folds = 20 independent fits. Each one is a task;
# the tasks run at the same time in separate R sessions ("workers").
library(furrr)                                             # parallel version of purrr's map()

tasks <- expand_grid(y       = c("FLI_3item", "INCOME_5"), # one row per outcome x scheme x fold
                     grouped = c(FALSE, TRUE),
                     fold    = 1:k)

data_list  <- list(FLI_3item = d.fl.rf, INCOME_5 = d.inc.rf)  # data per outcome
folds_list <- list()                                       # folds made once here, so they are
for (y in names(data_list)) for (g in c(FALSE, TRUE))      # identical to the sequential version
  folds_list[[paste(y, g)]] <- make_folds(data_list[[y]], k, g)

plan(multisession, workers = n.workers)                    # start the workers (works on Windows, Mac, Linux)

cv_res <- future_pmap_dfr(tasks, function(y, grouped, fold) {
  test_id <- folds_list[[paste(y, grouped)]] == fold       # rows held out in this task
  run_fold(data_list[[y]], y, test_id, grouped) %>%
    mutate(fold    = fold,
           outcome = y,
           cv      = if (grouped) "B: new country-years" else "A: random")
}, .options = furrr_options(seed = TRUE,                   # reproducible random numbers in workers
                            packages = c("ordinal", "ranger", "dplyr", "purrr", "tibble")),
   .progress = TRUE)                                       # progress bar

plan(sequential)                                           # close the workers again

saveRDS(cv_res, "rf_cv_results.rds")                       # save so you don't rerun

cv_res %>%
  filter(model %in% c("clmm", "rf")) %>%
  select(outcome, cv, fold, model, rps) %>%
  pivot_wider(names_from = model, values_from = rps) %>%   # clmm and rf side by side per fold
  mutate(diff = clmm - rf) %>%                             # clmm minus rf, per fold
  arrange(outcome, cv, fold) %>%
  print(n = Inf)

# ---- 8. RESULTS TABLE -------------------------------------------------------
cv_tab <- cv_res %>%
  group_by(outcome, cv, model) %>%
  summarise(across(c(rps, acc, mae),                       # mean and SD over folds
                   list(m = mean, sd = sd)), .groups = "drop") %>%
  group_by(outcome, cv) %>%
  mutate(rpss = 1 - rps_m / rps_m[model == "null"]) %>%    # RPS skill score: share of null RPS removed
  ungroup() %>%
  mutate(model = factor(model, c("null", "clmm", "rf", "rf_no_iwise"))) %>%
  arrange(outcome, cv, model)

print(cv_tab, width = Inf)


# ---- 9. PLOT ----------------------------------------------------------------
ggplot(cv_res %>% mutate(model = factor(model, levels(cv_tab$model))),
       aes(model, rps)) +
  geom_point(alpha = .5) +                                 # one point per fold
  stat_summary(fun = mean, geom = "crossbar", width = .4) +  # bar = mean over folds
  facet_grid(outcome ~ cv, scales = "free_y") +
  labs(x = NULL, y = "Ranked probability score (lower = better)") +
  theme_bw()
ggsave("rf_cv_rps.png", width = 8, height = 5)


# ---- 10. IWISE IMPORTANCE (full-sample RF) ----------------------------------
# Permutation importance: shuffle one predictor, measure how much worse the
# forest predicts (out-of-bag). Larger = the RF relies on it more.
imp_rf <- function(data, y) {
  rf <- ranger(x = data[, x.vars], y = data[[y]],
               probability  = TRUE, num.trees = n.trees,
               case.weights = data$w / mean(data$w),
               importance   = "permutation",               # permutation, not impurity (less biased)
               seed = seed, num.threads = parallel::detectCores())
  tibble(outcome = y, term = names(rf$variable.importance),
         importance = rf$variable.importance) %>%
    arrange(desc(importance)) %>%
    mutate(rank = row_number())
}

imp <- bind_rows(imp_rf(d.fl.rf, "FLI_3item"),
                 imp_rf(d.inc.rf, "INCOME_5"))
saveRDS(imp, "rf_importance.rds")

imp %>% filter(term %in% c("iwisescore", "iwise_mean"))   # where does IWISE rank?
print(imp, n = Inf)
