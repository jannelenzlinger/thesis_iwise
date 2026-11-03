# =============================================================================
# 07.3_interactions_iwise_cat.R
# Effect modification of the categorical IWISE association (adjusted M1)
#   Modifiers: gender (female), employment, urbanicity (urban_imp)
#
#   Step A  ONE AT A TIME (primary): M1 + iwise_cat:X, one model per X.
#           Each model answers one pre-specified question: "does the IWISE
#           association differ by X?" (effect modification, Knol &
#           VanderWeele 2012). Holm correction over the 3 tests.
#   Step B  ALL TOGETHER (sensitivity): M1 + all three interactions.
#           First a global ("chunk") LR test of all interaction terms
#           (Harrell, RMS). Then each block adjusted for the others: an
#           included interaction can absorb an omitted one when modifiers
#           are correlated, e.g. gender and employment (Beiser-McGrath &
#           Beiser-McGrath 2020).
#   Report  Ratio of ORs (interaction terms) AND the IWISE ORs within each
#           modifier group, with 95% CIs (Knol & VanderWeele 2012).
#
#   Random effects: random IWISE slope, per outcome exactly as chosen in 08
#           (re_cat["fl"], re_cat["inc"]), so each interaction model is
#           nested in its base model (m.fl.1.cat / m.inc.1.cat).
#           The modifiers are individual-level, so they get no random slope.
#           SLOW: 8 random-slope fits in total.
#   Needs:  models_cat_final.rds (08)
#   Saves:  models_cat_interactions.rds (m.fl.int, m.inc.int: named lists)
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
library(tidyverse)
library(ordinal)


# ---- 1. DATA ----------------------------------------------
list2env(readRDS("data/iwise_data_prepared.rds"), envir = globalenv())

prep <- function(data, y) {
  data %>%
    filter(cc) %>%
    mutate(.id          = row_number(),
           country_year = factor(.data[[v.cluster]]),
           iwise_cat    = factor(iwise_cat, ordered = FALSE),
           iwise_num    = as.integer(iwise_cat) - 1,       # needed if 08 chose the linear slope
           "{y}" := factor(.data[[y]], ordered = TRUE))
}

d.fl.cc  <- prep(d.fl,  "FLI_3item")
d.inc.cc <- prep(d.inc, "INCOME_5")
# If you releveled iwise_cat in 08, do the same here.


# ---- 2. HELPERS -------------------------------------------------------------
make_f <- function(y, rhs, re) as.formula(paste(y, "~", rhs, "+", re))

fit_clmm <- function(data, y, rhs, re) {
  clmm(make_f(y, rhs, re), data = data, weights = w_fit,
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

# Interaction coefficients that belong to one modifier term. R may name them
# "iwise_catX:femaleY" or "femaleY:iwise_catX" (order of first appearance in
# the formula), so both parts are identified by name, not by position.
int_terms <- function(m, modterm) {
  nm <- names(m$beta)[str_detect(names(m$beta), fixed(":"))]
  parts <- str_split(nm, fixed(":"))
  tibble(name = nm,
         cat  = map_chr(parts, ~ .x[str_starts(.x, "iwise_cat")][1]),
         mod  = map_chr(parts, ~ .x[!str_starts(.x, "iwise_cat")][1])) %>%
    filter(!is.na(cat), !is.na(mod), startsWith(mod, modterm))
}

# IWISE OR per category WITHIN each group of the modifier:
#   reference group : exp(b_cat)
#   other group g   : exp(b_cat + b_cat:g), SE from the covariance matrix.
# "reference" = the modifier's reference level (e.g. male, employed, rural).
# In the all-together model, the other modifiers are at their reference level.
iwise_by_group <- function(m, modterm, label) {
  b    <- m$beta
  V    <- vcov(m)[names(b), names(b)]
  main <- names(b)[str_starts(names(b), "iwise_cat") & !str_detect(names(b), fixed(":"))]
  it   <- int_terms(m, modterm)
  map_dfr(c("reference", unique(it$mod)), function(g) map_dfr(main, function(mt) {
    L <- setNames(numeric(length(b)), names(b))
    L[mt] <- 1
    if (g != "reference") L[it$name[it$cat == mt & it$mod == g]] <- 1
    est <- sum(L * b)
    se  <- sqrt(drop(t(L) %*% V %*% L))
    tibble(model = label, modifier = modterm, group = g,
           term = str_remove(mt, "^iwise_cat"),
           OR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
  }))
}

lr_test <- function(m0, m1) {                              # nested LR test by hand
  lr <- 2 * as.numeric(logLik(m1) - logLik(m0))            # (anova() fails for fit_clmm models)
  df <- length(m1$beta) - length(m0$beta)                  # extra fixed effects only
  c(LR = lr, df = df, p = pchisq(lr, df, lower.tail = FALSE))
}

wald_block <- function(m, modterm) {                       # joint Wald test of one interaction
  bn <- int_terms(m, modterm)$name                         # block, inside a larger model
  b  <- m$beta[bn]
  V  <- vcov(m)[bn, bn]
  W  <- drop(t(b) %*% solve(V) %*% b)
  c(W = W, df = length(bn), p = pchisq(W, length(bn), lower.tail = FALSE))
}


# ---- 3. BASE MODELS AND MODIFIERS -------------------------------------------
load_models("models_cat_final.rds")                        # m.*.cat, slope_type, re_cat from 08
slope_type                                                 # check: "full" or "linear" per model
re_cat                                                     # random part used below, per outcome

rhs_ctx_cat <- str_replace(rhs_ctx, "\\biwisescore\\b", "iwise_cat")

# Each modifier is used exactly as it appears in rhs_ctx (e.g. "female" or
# "factor(employment)"), so main effect and interaction match.
rhs_terms <- str_trim(str_split(rhs_ctx_cat, fixed("+"))[[1]])
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

# Formulas: one per modifier + all together
int_of  <- function(mods) paste0("iwise_cat:", t.mod[mods], collapse = " + ")
rhs_int <- c(map_chr(set_names(names(t.mod)), ~ paste(rhs_ctx_cat, "+", int_of(.x))),
             all = paste(rhs_ctx_cat, "+", int_of(names(t.mod))))
rhs_int                                                    # check: 4 formulas


# ---- 4. FIT OR LOAD ---------------------------------------------------------
# 4 models per outcome. Saved after every fit, so an interrupted run resumes.
# To force a refit: delete the .rds file and rm(m.fl.int, m.inc.int).
f.int <- "models_cat_interactions.rds"
if (file.exists(f.int)) load_models(f.int)
if (!exists("m.fl.int"))  m.fl.int  <- list()
if (!exists("m.inc.int")) m.inc.int <- list()

for (nm in names(rhs_int)) {
  if (is.null(m.fl.int[[nm]])) {
    message("FLI: ", nm);    m.fl.int[[nm]]  <- fit_clmm(d.fl.cc,  "FLI_3item", rhs_int[[nm]], re_cat[["fl"]])
  }
  if (is.null(m.inc.int[[nm]])) {
    message("Income: ", nm); m.inc.int[[nm]] <- fit_clmm(d.inc.cc, "INCOME_5",  rhs_int[[nm]], re_cat[["inc"]])
  }
  saveRDS(list(m.fl.int = m.fl.int, m.inc.int = m.inc.int), f.int)
}

# Fit quality: all columns should look fine
round(sapply(m.fl.int,  fit_quality), 4)
round(sapply(m.inc.int, fit_quality), 4)


# ---- 5. STEP A: ONE AT A TIME (PRIMARY) -------------------------------------
# LR test of each interaction block vs the base model; Holm over the 3 tests.
# With N ~ 80,000 small p-values are expected: judge the size of the ORs
# within groups (Section 7), not the p-value alone.
lr_tab <- function(base, ints, outcome) {
  imap_dfr(ints, function(m, nm) {
    r <- lr_test(base, m)
    tibble(outcome = outcome, model = nm, LR = r[["LR"]], df = r[["df"]], p = r[["p"]])
  }) %>%
    mutate(p_holm = p.adjust(replace(p, model == "all", NA), "holm"))  # "all" = global test (Step B)
}
lr.int <- bind_rows(lr_tab(m.fl.1.cat,  m.fl.int,  "FLI"),
                    lr_tab(m.inc.1.cat, m.inc.int, "INC"))
lr.int


# ---- 6. STEP B: ALL TOGETHER (SENSITIVITY) ----------------------------------
# 6a. Global chunk test = the "all" row of lr.int (all interaction terms at once).
lr.int %>% filter(model == "all")

# 6b. Each block within the all-together model (adjusted for the others)
round(sapply(t.mod, function(t) wald_block(m.fl.int$all,  t)), 4)   # FLI
round(sapply(t.mod, function(t) wald_block(m.inc.int$all, t)), 4)   # income

# 6c. Do the interaction ORs change between one-at-a-time and all-together?
#     Large changes = the modifications overlap (e.g. gender vs employment).
int_or <- function(ints, outcome) {
  imap_dfr(ints, ~ tidy_or(.x, .y)) %>%
    filter(str_detect(term, fixed(":")), str_detect(term, "iwise_cat")) %>%
    mutate(outcome = outcome, approach = if_else(model == "all", "all together", "one at a time"))
}
ror <- bind_rows(int_or(m.fl.int, "FLI"), int_or(m.inc.int, "INC"))
ror %>%                                                     # ratio of ORs, side by side
  select(outcome, term, approach, OR) %>%
  pivot_wider(names_from = approach, values_from = OR) %>%
  mutate(change_pct = 100 * (`all together` / `one at a time` - 1)) %>%
  print(n = Inf)


# ---- 7. IWISE ORs WITHIN MODIFIER GROUPS ------------------------------------
# Main output for the paper (Knol & VanderWeele 2012): the IWISE OR per
# category within each group, plus the ratio of ORs (ror above).
strat <- bind_rows(
  imap_dfr(t.mod, ~ iwise_by_group(m.fl.int[[.y]],  .x, "FLI, one at a time")),
  imap_dfr(t.mod, ~ iwise_by_group(m.inc.int[[.y]], .x, "INC, one at a time")),
  imap_dfr(t.mod, ~ iwise_by_group(m.fl.int$all,    .x, "FLI, all together")),
  imap_dfr(t.mod, ~ iwise_by_group(m.inc.int$all,   .x, "INC, all together")))
strat %>% print(n = Inf)

p.strat <- strat %>%
  filter(str_detect(model, "one at a time")) %>%
  ggplot(aes(term, OR, ymin = lo, ymax = hi, colour = group)) +
  geom_hline(yintercept = 1, colour = "grey60", linetype = 2) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.3) +
  facet_grid(modifier ~ model) +
  scale_y_log10() +
  labs(x = "WI category (vs reference)", y = "Adjusted OR (95% CI, log scale)",
       colour = "Modifier group") +
  theme_minimal()
p.strat

dir.create("figures", showWarnings = FALSE)
ggsave("figures/interactions_iwise_cat.png", p.strat, width = 9, height = 7, dpi = 300, bg = "white")
saveRDS(list(lr.int = lr.int, ror = ror, strat = strat), "interaction_results.rds")