# =============================================================================
# 07.3_interactions_iwise_cont.R
# Effect modification of the CONTINUOUS IWISE association (adjusted M1)
#   Predictor: iwisescore (continuous)   Outcomes: FLI_3item, INCOME_5
#   Modifiers: gender (female), employment, urbanicity (urban_imp)
#
#   Step A  ONE AT A TIME (primary): M1 + iwisescore:X, one model per X.
#           Each model answers one pre-specified question: "does the IWISE
#           association differ by X?" (effect modification, Knol &
#           VanderWeele 2012). Holm correction over the 3 tests.
#   Step B  ALL TOGETHER (sensitivity): M1 + all three interactions.
#           First a global ("chunk") LR test of all interaction terms
#           (Harrell, RMS). Then each block adjusted for the others: an
#           included interaction can absorb an omitted one when modifiers
#           are correlated, e.g. gender and employment (Beiser-McGrath &
#           Beiser-McGrath 2020).
#   Report  Ratio of ORs (interaction terms) AND the IWISE OR within each
#           modifier group, with 95% CIs (Knol & VanderWeele 2012).
#           ORs are per `iwise_unit` points of the IWISE score.
#
#   Model:  clmm, logit link, weights = wgt_using,
#           random intercept + random IWISE slope per country(-wave):
#           (1 + iwisescore | country_year).
#           The base model M1 is fitted here with the same data, weights and
#           random part, so every interaction model is nested in it.
#           SLOW: 5 random-slope fits per outcome (10 in total).
#   Saves:  models_cont_interactions.rds (m.fl.int, m.inc.int: named lists)
#           interaction_results_cont.rds
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
library(tidyverse)
library(ordinal)


# ---- 1. DATA ----------------------------------------------------------------
list2env(readRDS("data/iwise_data_prepared.rds"), envir = globalenv())

prep <- function(data, y) {
  data %>%
    filter(cc) %>%
    mutate(.id          = row_number(),
           country_year = factor(.data[[v.cluster]]),
           "{y}" := factor(.data[[y]], ordered = TRUE))
}

d.fl.cc  <- prep(d.fl,  "FLI_3item")
d.inc.cc <- prep(d.inc, "INCOME_5")

# Random part: random intercept + random IWISE slope per country-wave
# (correlated). Modifiers are individual-level, so they get no random slope.
re <- "(1 + iwisescore | country_year)"

# Scale for reporting: OR per `iwise_unit` points of iwisescore (0-36).
# 1 = per point; e.g. set to sd(d.fl.cc$iwisescore) for per-SD ORs.
iwise_unit <- 1


# ---- 2. HELPERS -------------------------------------------------------------
make_f <- function(y, rhs, re) as.formula(paste(y, "~", rhs, "+", re))

fit_clmm <- function(data, y, rhs, re) {
  clmm(make_f(y, rhs, re), data = data, weights = wgt_using,
       link = "logit", Hess = TRUE, nAGQ = 1)
}

tidy_or <- function(m, label) {
  co <- summary(m)$coefficients[names(m$beta), , drop = FALSE]
  tibble(model = label, term = rownames(co),
         b = co[, 1], se = co[, 2], OR = exp(b),
         lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = co[, 4])
}

load_models <- function(file) {
  obj <- readRDS(file)
  list2env(obj, envir = globalenv())
  message("Loaded from ", file, ": ", paste(names(obj), collapse = ", "))
  invisible(names(obj))
}

fit_quality <- function(m) {
  c(condHess = suppressWarnings(summary(m)$condHess),
    max_grad = max(abs(m$gradient)),
    n_na_se  = sum(is.na(summary(m)$coefficients[names(m$beta), 2])))
}

# Interaction coefficients of iwisescore with one modifier term. R may name
# them "iwisescore:female" or "female:iwisescore", so both parts are
# identified by name, not by position.
int_terms <- function(m, modterm) {
  nm <- names(m$beta)[str_detect(names(m$beta), fixed(":"))]
  parts <- str_split(nm, fixed(":"))
  tibble(name = nm,
         iw   = map_lgl(parts, ~ any(.x == "iwisescore")),
         mod  = map_chr(parts, ~ .x[.x != "iwisescore"][1])) %>%
    filter(iw, !is.na(mod), startsWith(mod, modterm)) %>%
    select(-iw)
}

# IWISE OR (per iwise_unit points) WITHIN each group of the modifier:
#   reference group : exp(unit * b_iwise)
#   other group g   : exp(unit * (b_iwise + b_iwise:g)), SE from vcov.
# "reference" = the modifier's reference level (male, employed, rural).
# In the all-together model, the other modifiers are at their reference level.
iwise_by_group <- function(m, modterm, label, unit = iwise_unit) {
  b  <- m$beta
  V  <- vcov(m)[names(b), names(b)]
  it <- int_terms(m, modterm)
  map_dfr(c("reference", it$mod), function(g) {
    L <- setNames(numeric(length(b)), names(b))
    L["iwisescore"] <- 1
    if (g != "reference") L[it$name[it$mod == g]] <- 1
    est <- unit * sum(L * b)
    se  <- unit * sqrt(drop(t(L) %*% V %*% L))
    tibble(model = label, modifier = modterm, group = g,
           OR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
  })
}

lr_test <- function(m0, m1) {                              # nested LR test by hand
  lr <- 2 * as.numeric(logLik(m1) - logLik(m0))            # (anova() fails for fit_clmm models)
  df <- length(m1$beta) - length(m0$beta)                  # extra fixed effects only
  c(LR = lr, df = df, p = pchisq(lr, df, lower.tail = FALSE))
}

wald_block <- function(m, modterm) {                       # joint Wald test of one interaction
  bn <- int_terms(m, modterm)$name                         # block, inside a larger model
  b  <- m$beta[bn]
  V  <- vcov(m)[bn, bn, drop = FALSE]
  W  <- drop(t(b) %*% solve(V) %*% b)
  c(W = W, df = length(bn), p = pchisq(W, length(bn), lower.tail = FALSE))
}


# ---- 3. FORMULAS AND MODIFIERS ----------------------------------------------
# rhs_ctx = adjusted M1 right-hand side (from the prepared data), must contain
# the continuous iwisescore.
stopifnot(exists("rhs_ctx"), str_detect(rhs_ctx, "\\biwisescore\\b"))
rhs_ctx

# Each modifier is used exactly as it appears in rhs_ctx (e.g. "female" or
# "factor(employment)"), so main effect and interaction match.
rhs_terms <- str_trim(str_split(rhs_ctx, fixed("+"))[[1]])
term_of <- function(var) {
  hit <- rhs_terms[str_detect(rhs_terms, paste0("\\b", var, "\\b"))]
  if (length(hit) != 1) stop(var, " must appear in exactly one term of rhs_ctx (found ", length(hit), ")")
  hit
}
v.mod <- c(gender = "female", employment = "employment", urbanicity = "urban_imp")
t.mod <- map_chr(v.mod, term_of)
t.mod                                                      # check the terms

# Employment must be categorical, otherwise its interaction is linear (1 df)
sapply(v.mod, function(v) class(d.fl.cc[[v]])[1])
if (!is.factor(d.fl.cc$employment) && !str_detect(t.mod[["employment"]], fixed("factor(")))
  warning("employment is numeric in the data and in rhs_ctx: make it a factor first.")

# Formulas: base + one per modifier + all together
int_of  <- function(mods) paste0("iwisescore:", t.mod[mods], collapse = " + ")
rhs_int <- c(base = rhs_ctx,
             map_chr(set_names(names(t.mod)), ~ paste(rhs_ctx, "+", int_of(.x))),
             all  = paste(rhs_ctx, "+", int_of(names(t.mod))))
rhs_int                                                    # check: 5 formulas


# ---- 4. FIT OR LOAD ---------------------------------------------------------
# 5 models per outcome. Saved after every fit, so an interrupted run resumes.
# To force a refit: delete the .rds file and rm(m.fl.int, m.inc.int).
f.int <- "models_cont_interactions.rds"
if (file.exists(f.int)) load_models(f.int)
if (!exists("m.fl.int"))  m.fl.int  <- list()
if (!exists("m.inc.int")) m.inc.int <- list()

for (nm in names(rhs_int)) {
  if (is.null(m.fl.int[[nm]])) {
    message("FLI: ", nm);    m.fl.int[[nm]]  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_int[[nm]], re)
  }
  if (is.null(m.inc.int[[nm]])) {
    message("Income: ", nm); m.inc.int[[nm]] <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_int[[nm]], re)
  }
  saveRDS(list(m.fl.int = m.fl.int, m.inc.int = m.inc.int), f.int)
}

# Fit quality: all columns should look fine
round(sapply(m.fl.int,  fit_quality), 4)
round(sapply(m.inc.int, fit_quality), 4)

# Random intercept and slope variances (and their correlation) per outcome
map(list(FLI = m.fl.int$base, INC = m.inc.int$base), ~ VarCorr(.x))


# ---- 5. STEP A: ONE AT A TIME (PRIMARY) -------------------------------------
# LR test of each interaction block vs the base model; Holm over the 3 tests.
# With N ~ 80,000 small p-values are expected: judge the size of the ORs
# within groups (Section 7), not the p-value alone.
lr_tab <- function(ints, outcome) {
  imap_dfr(ints[names(ints) != "base"], function(m, nm) {
    r <- lr_test(ints$base, m)
    tibble(outcome = outcome, model = nm, LR = r[["LR"]], df = r[["df"]], p = r[["p"]])
  }) %>%
    mutate(p_holm = p.adjust(replace(p, model == "all", NA), "holm"))  # "all" = global test (Step B)
}
lr.int <- bind_rows(lr_tab(m.fl.int, "FLI"), lr_tab(m.inc.int, "INC"))
lr.int


# ---- 6. STEP B: ALL TOGETHER (SENSITIVITY) ----------------------------------
# 6a. Global chunk test = the "all" row of lr.int (all interaction terms at once).
lr.int %>% filter(model == "all")

# 6b. Each block within the all-together model (adjusted for the others)
round(sapply(t.mod, function(t) wald_block(m.fl.int$all,  t)), 4)   # FLI
round(sapply(t.mod, function(t) wald_block(m.inc.int$all, t)), 4)   # income

# 6c. Do the ratio of ORs change between one-at-a-time and all-together?
#     Large changes = the modifications overlap (e.g. gender vs employment).
int_or <- function(ints, outcome) {
  imap_dfr(ints[names(ints) != "base"], ~ tidy_or(.x, .y)) %>%
    filter(str_detect(term, fixed(":")), str_detect(term, "iwisescore")) %>%
    mutate(OR = exp(iwise_unit * b),                       # per iwise_unit points
           lo = exp(iwise_unit * (b - 1.96 * se)),
           hi = exp(iwise_unit * (b + 1.96 * se)),
           outcome = outcome,
           approach = if_else(model == "all", "all together", "one at a time"))
}
ror <- bind_rows(int_or(m.fl.int, "FLI"), int_or(m.inc.int, "INC"))
ror %>%                                                     # ratio of ORs, side by side
  select(outcome, term, approach, OR) %>%
  pivot_wider(names_from = approach, values_from = OR) %>%
  mutate(change_pct = 100 * (`all together` / `one at a time` - 1)) %>%
  print(n = Inf)


# ---- 7. IWISE ORs WITHIN MODIFIER GROUPS ------------------------------------
# Main output for the paper (Knol & VanderWeele 2012): the IWISE OR within
# each group, plus the ratio of ORs (ror above).
strat <- bind_rows(
  imap_dfr(t.mod, ~ iwise_by_group(m.fl.int[[.y]],  .x, "FLI, one at a time")),
  imap_dfr(t.mod, ~ iwise_by_group(m.inc.int[[.y]], .x, "INC, one at a time")),
  imap_dfr(t.mod, ~ iwise_by_group(m.fl.int$all,    .x, "FLI, all together")),
  imap_dfr(t.mod, ~ iwise_by_group(m.inc.int$all,   .x, "INC, all together")))
strat %>% print(n = Inf)

p.strat <- strat %>%
  filter(str_detect(model, "one at a time")) %>%
  ggplot(aes(x = OR, xmin = lo, xmax = hi, y = group)) +
  geom_vline(xintercept = 1, colour = "grey60", linetype = 2) +
  geom_pointrange(size = 0.3) +
  facet_grid(modifier ~ model, scales = "free_y", space = "free_y") +
  scale_x_log10() +
  labs(x = paste0("Adjusted OR per ", iwise_unit, "-point IWISE increase (95% CI, log scale)"),
       y = "Modifier group") +
  theme_minimal()
p.strat

dir.create("figures", showWarnings = FALSE)
ggsave("figures/interactions_iwise_cont.png", p.strat, width = 9, height = 7, dpi = 300, bg = "white")
saveRDS(list(lr.int = lr.int, ror = ror, strat = strat), "interaction_results_cont.rds")






