# =============================================================================
# 07.4_models_wgt_using.R
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

list2env(readRDS("models_wgt_using_cont.rds"), envir = globalenv())
list2env(readRDS("models_wgt_using_cat.rds"),  envir = globalenv())


# ---- 5. LOAD MODELS (continuous IWISE only) ----------------------------------
# Run sections 0-2 first (data + helpers), skip 3-4 (slow refits).
list2env(readRDS("models_wgt_using_cont.rds"), envir = globalenv())

mods <- list(                                              # model + its analytic sample + outcome
  inc.0 = list(m = m.inc.0.wu, d = d.inc.cc, y = "INCOME_5"),
  inc.1 = list(m = m.inc.1.wu, d = d.inc.cc, y = "INCOME_5"),
  fl.0  = list(m = m.fl.0.wu,  d = d.fl.cc,  y = "FLI_3item"),
  fl.1  = list(m = m.fl.1.wu,  d = d.fl.cc,  y = "FLI_3item")
)


# ---- 6. CONVERGENCE / FIT CHECKS ---------------------------------------------
check_fit <- function(x) {
  s  <- summary(x$m)
  vc <- as.matrix(VarCorr(x$m)[[1]])
  cat("\n================", deparse(substitute(x)), "================\n")
  cat("condHess        :", s$condHess, "  (NaN or > 1e6 = problem)\n")
  cat("max |gradient|  :", max(abs(x$m$gradient)), "  (large = not converged)\n")
  cat("NA std. errors  :", sum(is.na(s$coefficients[, 2])), "  (must be 0)\n")
  cat("RE SDs          :", round(sqrt(diag(vc)), 4), "\n")
  if (ncol(vc) > 1)
    cat("RE correlation  :", round(cov2cor(vc)[1, 2], 3), "  (+-1 = problem)\n")
  cat("weights         : n =", nrow(x$d), " sum =", round(sum(x$d$w_fit)), "\n")
  print(summary(x$d$w_fit))                                # very large single weights?
}

for (nm in names(mods)) { cat("\n#####", nm); check_fit(mods[[nm]]) }


# ---- 7. IWISE EFFECT: beta, OR, OR per 12 points, ICC, country spread ---------
iwise_res <- function(m, k = 12) {
  ct <- summary(m)$coefficients
  b  <- ct["iwisescore", "Estimate"]
  se <- ct["iwisescore", "Std. Error"]
  vc <- as.matrix(VarCorr(m)[[1]])
  sd_sl <- if (ncol(vc) > 1) sqrt(vc[2, 2]) else NA       # no slope in income M0
  tibble(
    beta    = b, beta_lo = b - 1.96 * se, beta_hi = b + 1.96 * se,
    p       = ct["iwisescore", "Pr(>|z|)"],
    OR      = exp(b), OR_lo = exp(b - 1.96 * se), OR_hi = exp(b + 1.96 * se),
    OR12    = exp(k * b), OR12_lo = exp(k * (b - 1.96 * se)), OR12_hi = exp(k * (b + 1.96 * se)),
    ICC     = vc[1, 1] / (vc[1, 1] + pi^2 / 3),            # latent-scale ICC at IWISE = 0
    ORc_lo  = exp(b - 1.96 * sd_sl),                        # range of 95% of country ORs
    ORc_hi  = exp(b + 1.96 * sd_sl)
  )
}

res <- map_dfr(mods, ~ iwise_res(.x$m), .id = "model")
res %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% print(width = Inf)

# attenuation of beta1 after adjustment (%)
c(income = 100 * (1 - res$beta[res$model == "inc.1"] / res$beta[res$model == "inc.0"]),
  FLI    = 100 * (1 - res$beta[res$model == "fl.1"]  / res$beta[res$model == "fl.0"]))


# ---- 8. PREDICTED PROBABILITIES PER OUTCOME CATEGORY -------------------------
# Reference person: all covariates at reference level / centred mean, random effects = 0.
# P(Y <= j) = plogis(theta_j - beta * IWISE); category prob = difference of neighbours.
pred_probs <- function(x, iwise = c(0, 12, 24)) {
  th  <- x$m$alpha                                          # thresholds
  b   <- x$m$beta["iwisescore"]
  lev <- levels(x$d[[x$y]])
  map_dfr(iwise, function(v) {
    cum <- plogis(th - b * v)
    tibble(iwise = v, category = factor(lev, levels = lev), prob = diff(c(0, cum, 1)))
  })
}

probs <- map_dfr(mods, pred_probs, .id = "model")

probs %>%                                                   # table like slide 6
  mutate(prob = round(prob, 3)) %>%
  pivot_wider(names_from = category, values_from = prob) %>%
  print(n = Inf)

probs %>%                                                   # plot, adjusted models
  filter(model %in% c("inc.1", "fl.1")) %>%
  ggplot(aes(category, prob, fill = factor(iwise))) +
  geom_col(position = "dodge") +
  facet_wrap(~ model, scales = "free_x") +
  labs(x = NULL, y = "Predicted probability", fill = "IWISE score") +
  theme_minimal()

# ---- 9. PLOT: PREDICTED PROBABILITIES (100% stacked bars, as in pptx slide 6) ----
# ---- 9. PLOT: PREDICTED PROBABILITIES (100% stacked bars, as in pptx slide 6) ----
library(patchwork)                                         # combine the two panels

pal <- c("#E4EEDA", "#BCD6A4", "#86B36B", "#4F8A43", "#2A5530")   # light = lowest category

plot_probs <- function(df, mod, title, lab) {
  df %>%
    filter(model == mod) %>%
    mutate(iwise    = factor(paste("IWISE", iwise), levels = paste("IWISE", c(24, 12, 0))),
           category = as.character(category),               # map_dfr merged levels of both outcomes
           category = factor(category, levels = sort(unique(as.numeric(category))), labels = lab)) %>%
    ggplot(aes(x = prob, y = iwise, fill = category)) +
    geom_col(width = 0.75, colour = "white", linewidth = 0.3,
             position = position_stack(reverse = TRUE)) +          # lowest category on the left
    geom_text(aes(label = ifelse(prob >= 0.045, scales::percent(prob, accuracy = 1), ""),
                  colour = category),
              position = position_stack(vjust = 0.5, reverse = TRUE), size = 3.3) +
    scale_fill_manual(values = pal) +
    scale_colour_manual(values = c("#1F1F1F", "#1F1F1F", "#1F1F1F", "white", "white"), guide = "none") +
    scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(title = title, x = NULL, y = NULL, fill = NULL) +
    guides(fill = guide_legend(nrow = 1)) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(face = "bold", size = 11))
}

p.fl  <- plot_probs(probs, "fl.1",  "Financial Life Index category", c("0", "25", "50", "75", "100"))
p.inc <- plot_probs(probs, "inc.1", "Income quintile", c("Q1 (poorest)", "Q2", "Q3", "Q4", "Q5 (richest)"))

p.probs <- p.fl / p.inc +                               # stacked: one panel per row
  plot_annotation(caption = "Fully adjusted models, reference person (all covariates at reference level / centred mean, random effects = 0).")
p.probs

ggsave("fig_pred_probs_wu.png", p.probs, width = 9, height = 6, dpi = 300, bg = "white")



conv_check <- function(m) {
  H <- m$Hessian; g <- m$gradient
  step <- tryCatch(solve(H, g), error = function(e) NA)
  se   <- tryCatch(sqrt(diag(solve(H))), error = function(e) NA)
  c(max_grad = max(abs(g)),
    max_step_rel_se = max(abs(step) / se))   # < ~0.001: estimates have converged
}
conv_check(m.inc.1.wu)

H <- m.fl.1.wu$Hessian; g <- m.fl.1.wu$gradient
rel <- abs(solve(H, g)) / sqrt(diag(solve(H)))
sort(setNames(rel, names(g)), decreasing = TRUE)[1:5]



k <- mean(d.fl.cc$wgt_using)
k
rel / sqrt(k)          # step/SE on the mean-1 weight scale

rel[grep("iwise", names(rel))]

m.fl.1.wu$ST          # the raw L matrix
VarCorr(m.fl.1.wu)    # SDs and correlation, computed from L L'
