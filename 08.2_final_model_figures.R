# =============================================================================
# 08.2_results_figures.R
# Results and figures from the final models of 08.1
#
#   1  Load data and the final models (models_final.rds)
#   2  Settings: colours, IWISE values, output folder
#   3  Figure 1: predicted probabilities at IWISE 0, 12, 24, 36 (reference person)
#   4  Figure 2: individual covariates, OR with 95% CI (adjusted models)
#   5  Figure 3: IWISE OR with 95% CI and the range across country-years
#   6  Figure 4: IWISE OR for each country-year (all labels: appendix)
#      6b Figure 4b: same, only the extremes labelled (slides)
#      6c Figure 4c: country-years by income group
#   7  Check: refit with IWISE in 12-point steps, compare with exp(12 * beta)
#
#   Models: m.fl.1 and m.inc.1 (fully adjusted, random intercept + IWISE slope)
#           for Figures 1-4. All four models for the check in Section 7.
#   Saves:  figures/fig1_pred_probs.png, fig2_covariates.png,
#           fig3_iwise_range.png, fig4_iwise_by_country.png,
#           fig4b_iwise_by_country_slim.png, fig4c_iwise_by_income_group.png
#           pred_probs.csv, iwise_by_country.csv, check_iwise12.csv
#           models_iwise12.rds (Section 7, slow: ~2 h per adjusted model)
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
library(dplyr)                                  # data handling
library(tidyr)                                  # pivot_wider()
library(ggplot2)                                # figures
library(patchwork)                              # combine figures (p1 / p2, p1 | p2)
library(ordinal)                                # clmm(), ranef(), condVar()
library(future)                                 # parallel refits (Section 7)
library(furrr)                                  # future_imap()


# ---- 1. DATA AND MODELS -----------------------------------------------------
list2env(readRDS("data/iwise_data_prepared.rds"), envir = globalenv())   # d.fl, d.inc from 05

# Same analysis samples as 08.1, Section 3 (needed for labels and the refit)
prep <- function(data, y) {
  data %>%
    filter(cc) %>%                                          # 05 complete cases
    mutate(country_year = factor(country_year),             # grouping must be a factor
           across(all_of(y), ~ factor(.x, ordered = TRUE))) # ordered outcome; keeps level order
}
d.fl  <- prep(d.fl,  "FLI_3item")
d.inc <- prep(d.inc, "INCOME_5")

fits <- readRDS("models_final.rds")                         # m.fl.0, m.fl.1, m.inc.0, m.inc.1

# Check: every model loaded, same rows as the data, optimiser message
stopifnot(all(sapply(fits, inherits, "clmm")))              # FALSE = a fit failed in 08.1
stopifnot(nrow(fits$m.fl.1$model)  == nrow(d.fl),
          nrow(fits$m.inc.1$model) == nrow(d.inc))
sapply(fits, function(m) m$optRes$message)                  # "relative convergence" = fine

mods <- list(FLI = fits$m.fl.1, Income = fits$m.inc.1)      # the models for Figures 1-4


# ---- 2. SETTINGS ------------------------------------------------------------
dir.create("figures", showWarnings = FALSE)

iwise_vals <- c(0, 12, 24, 36)                  # IWISE values for Figure 1 (36 = maximum)
or_unit    <- 1                                 # Figures 3 and 4: OR per 1 point (12 = per 12 points)
font       <- "Times New Roman"                                # "" = R default; e.g. "serif" or a font you have installed

# Reference person for Figure 1: every factor at its reference level and every
# centred covariate at 0 (= sample average). Non-zero values of uncentred
# covariates go here. bas_rate_pp = 0 means no change in SDG 6.1.
# INDEX_FS is NOT centred: 0 = the lowest Food & Shelter score. The slide used 0.
# For an average person use: c(INDEX_FS = mean(d.fl$INDEX_FS))
ref_values <- c(INDEX_FS = 0)


# MetBrewer "Archambault" (MetBrewer::MetPalettes$Archambault[[1]])
col_cat   <- c("#ab3329", "#ed968c", "#f9d14a", "#88a0dc", "#381a61")  # 5 outcome categories, lowest -> highest (brick red, salmon, yellow, light blue, dark purple)
col_title <- "#381a61"                          # dark purple: titles, points
col_head  <- "black"                          # plum: panel headings
col_band  <- "#e78429"                          # orange: range across country-years / average OR line
col_shade <- "#F2F2F2"                          # light grey: alternating row shading

theme_set(theme_minimal(base_size = 12, base_family = font) +
            theme(panel.grid = element_blank(),
                  plot.title = element_text(colour = col_head, face = "bold", size = 14),
                  axis.title = element_text(colour = "grey30", size = 10)))

pct <- function(x) paste0(round(100 * x), "%")  # 0.153 -> "15%"


# ---- 3. FIGURE 1: PREDICTED PROBABILITIES -----------------------------------
# clmm: logit P(Y <= j) = theta_j - eta, with eta = beta_IWISE * IWISE
# + beta * (reference values). Random effects = 0 (a typical country-year).
# P(Y = j) = P(Y <= j) - P(Y <= j-1).
# x and per: only for the refit in Section 7 (iwise12 = IWISE / 12).
pred_ref <- function(m, iwise = iwise_vals, ref = ref_values, x = "iwisescore", per = 1) {
  y   <- all.vars(formula(m))[1]                            # outcome name (first variable in the formula)
  lev <- levels(m$model[[y]])                               # outcome categories, lowest to highest
  ref <- ref[names(ref) %in% names(m$beta)]                 # keep only covariates the model has (unadjusted models: none)
  eta_ref <- sum(m$beta[names(ref)] * ref)                  # their contribution to eta: sum of beta * value
  bind_rows(lapply(iwise, function(v) {                     # repeat for every IWISE value in iwise_vals
    cum <- plogis(m$alpha - (m$beta[x] * v / per + eta_ref))   # P(Y <= j): inverse logit of (threshold - eta), one per threshold
    tibble(iwise = v, category = factor(lev, levels = lev), # one row per category
           prob = diff(c(0, cum, 1)))                       # P(Y = j) = difference of neighbouring cumulative probs; sums to 1
  }))
}

pp <- bind_rows(lapply(mods, pred_ref), .id = "outcome")   # predictions for both models, labelled FLI / Income

for (o in names(mods)) {                                    # one table per outcome
  pp %>% filter(outcome == o) %>%                           # rows of this outcome
    mutate(prob = round(100 * prob, 1)) %>%                 # proportion -> % with 1 decimal
    pivot_wider(names_from = category, values_from = prob) %>%   # one row per IWISE value, one column per category
    print()                                                 # show in the console
}
write.csv(pp, "pred_probs.csv", row.names = FALSE)          # save all predictions

# Text colour that is readable on each fill: white or black, whichever has the
# higher contrast with that fill (WCAG relative-luminance formula). Works for
# any palette, so it stays correct if col_cat changes.
text_on <- function(fill) {
  rgb <- grDevices::col2rgb(fill) / 255                     # hex -> red/green/blue on a 0-1 scale (one column per colour)
  lin <- ifelse(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)   # undo the sRGB gamma (linear light)
  L   <- colSums(lin * c(0.2126, 0.7152, 0.0722))           # relative luminance: weighted sum of R, G, B
  ifelse(1.05 / (L + 0.05) >= (L + 0.05) / 0.05, "white", "black")   # contrast with white vs with black: take the larger
}
text_on(col_cat)                                            # check: one text colour per category

plot_probs <- function(outcome, title, legend_labels, first_lab, last_lab) {
  d <- filter(pp, outcome == !!outcome) %>%                 # rows of this outcome (!! = the function argument, not the column)
    mutate(category = factor(category, levels = unique(category)))  # back to a factor in the right order (bind_rows() made it text)
  lo <- d$prob[as.integer(d$category) == 1]                         # lowest category, one value per IWISE value
  hi <- d$prob[as.integer(d$category) == nlevels(d$category)]       # highest category, one value per IWISE value
  n   <- length(iwise_vals)                                 # number of IWISE values (first = lowest, last = highest)
  cap <- paste0(first_lab, ": ", pct(lo[1]), " at IWISE ", iwise_vals[1],   # lowest category at the lowest IWISE ...
                " → ", pct(lo[n]), " at IWISE ", iwise_vals[n], "\n",       # ... and at the highest IWISE
                last_lab,  ": ", pct(hi[1]), " at IWISE ", iwise_vals[1],   # highest category at the lowest IWISE ...
                " → ", pct(hi[n]), " at IWISE ", iwise_vals[n])             # ... and at the highest IWISE
  ggplot(d, aes(prob, factor(iwise, levels = rev(iwise_vals)), fill = category)) +   # x = probability, y = IWISE value (first value on top)
    geom_col(width = 0.75, colour = "white", linewidth = 0.4,     # bars with thin white gaps between segments
             position = position_stack(reverse = TRUE)) +   # stack with the lowest category on the left
    geom_text(aes(label = ifelse(prob >= 0.045, pct(prob), "")),  # % label; none on segments under 4.5% (too thin)
              position = position_stack(vjust = 0.5, reverse = TRUE),   # centre of each segment, same order as the bars
              size = 3.4, family = font,                    # label size and font
              colour = text_on(col_cat)[as.integer(d$category)]) +   # readable text colour for each segment's fill
    scale_fill_manual(values = col_cat, labels = legend_labels) +    # fill colours and legend labels, by category order
    scale_x_continuous(expand = c(0, 0)) +                  # bars start and end exactly at the panel edges
    scale_y_discrete(labels = function(x) paste("IWISE", x)) +       # y labels: "IWISE 0", "IWISE 12", ...
    labs(title = title, x = NULL, y = NULL, fill = NULL, caption = cap) +   # title and caption; no axis titles
    theme(legend.position = "top", legend.justification = "left",   # legend above the bars, left-aligned
          axis.text.x = element_blank(), axis.text.y = element_text(size = 11),   # no x numbers (labels show the %)
          plot.caption = element_text(hjust = 0.5, face = "italic", size = 10)) +   # centred italic caption
    guides(fill = guide_legend(nrow = 1, keywidth = 0.6, keyheight = 0.6))  # legend in one row, small keys
}

p1 <- plot_probs("FLI",    "Financial Life Index category", c("0", "25", "50", "75", "100"),   # left panel
                 "Lowest category (0)", "Highest category (100)") |   # | = side by side (patchwork)
      plot_probs("Income", "Income quintile", paste0("Q", 1:5),                                 # right panel
                 "Poorest quintile", "Richest quintile")
p1                                                          # show the figure
ggsave("figures/fig1_pred_probs_colour.png", p1, width = 6.3, height = 3.5, dpi = 600, bg = "white",
       device = ragg::agg_png, scale = 1.6)          # scale > 1 = more room, smaller-looking text


# ---- 4. FIGURE 2: INDIVIDUAL COVARIATES -------------------------------------
# Names = coefficient names in the model, values = labels in the figure (top to bottom).
# INDEX_FS is left out, as on the slide: uncomment its line to show it.
cov_labels <- c(
  "femaleFemale"                   = "Female (vs male)",
  "urban_impPeri-urban/urban"      = "Peri-urban/urban (vs rural)",
  "age_gp_profile25-34"            = "Age 25-34 (vs 15-24)",
  "age_gp_profile35-49"            = "Age 35-49",
  "age_gp_profile50+"              = "Age 50+",
  "maritalstatusMarried/partnered" = "Married/partnered (vs single)",
  "maritalstatusSep/div/widowed"   = "Sep./div./widowed",
  "hhsize_c"                       = "Household size, +1 person",
  "employmentUnderemployed"        = "Underemployed (vs employed)",
  "employmentUnemployed"           = "Unemployed",
  "employmentOut of workforce"     = "Out of workforce",
  "educationSecondary"             = "Secondary (vs elementary)",
  "educationCollege"               = "College"
  # , "INDEX_FS"                   = "Food & Shelter index, +1 point"
)

cov_or <- bind_rows(lapply(mods, function(m) {
  cf <- coef(summary(m))                                    # estimates and SEs
  miss <- setdiff(names(cov_labels), rownames(cf))
  if (length(miss)) stop("Not in the model: ", paste(miss, collapse = ", "))   # check the names
  tibble(term = names(cov_labels), b = cf[term, "Estimate"], se = cf[term, "Std. Error"])
}), .id = "outcome") %>%
  mutate(OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),        # 95% Wald CI
         label = factor(cov_labels[term], levels = rev(cov_labels)),          # first label on top
         outcome = recode(outcome, FLI = "Financial Life Index", Income = "Income quintile"))

shade <- tibble(y = seq(length(cov_labels), 1, by = -2))    # every second row, from the top

p2 <- ggplot(cov_or, aes(OR, label)) +
  geom_blank() +                                            # sets up the axes before the shading
  annotate("rect", xmin = 1e-3, xmax = 1e3, ymin = shade$y - 0.5, ymax = shade$y + 0.5,
           fill = col_shade) +                              # alternating row shading
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(xmin = lo, xmax = hi), colour = col_title, linewidth = 0.6) +   # 95% CI
  geom_point(size = 2.4, colour = col_title) +
  geom_text(aes(x = Inf, label = sprintf("%.2f", OR)), hjust = 1.2, size = 3.6, family = font) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4), labels = c("0.5", "1", "2", "4")) +
  coord_cartesian(xlim = c(0.35, 10)) +                     # room on the right for the values
  facet_wrap(~ outcome) +
  labs(x = "Odds ratio (log scale), 95% CI", y = NULL) +
  theme(strip.text = element_text(colour = col_title, face = "bold", size = 13, hjust = 0),
        axis.text.y = element_text(size = 11, colour = "grey15"),
        axis.line.x = element_line(colour = "grey70"),
        panel.spacing = unit(1.5, "lines"))
p2
ggsave("figures/fig2_covariates.png", p2, width = 10, height = 6.5, dpi = 300, bg = "white")


# ---- 5. FIGURE 3: IWISE EFFECT ACROSS COUNTRY-YEARS -------------------------
# Dark point + bar: fixed IWISE OR with 95% CI (the average country-year).
# Light band: range that covers 95% of country-year ORs,
# exp(beta +/- 1.96 * SD of the random slope).
iwise_range <- bind_rows(lapply(mods, function(m) {
  cf <- coef(summary(m))["iwisescore", ]
  sd_slope <- attr(VarCorr(m)$country_year, "stddev")["iwisescore"]
  tibble(b = cf["Estimate"], se = cf["Std. Error"], sd_slope = unname(sd_slope))
}), .id = "outcome") %>%
  mutate(OR       = exp(or_unit * b),
         lo       = exp(or_unit * (b - 1.96 * se)),         # 95% CI of the average OR
         hi       = exp(or_unit * (b + 1.96 * se)),
         range_lo = exp(or_unit * (b - 1.96 * sd_slope)),   # 95% of country-year ORs
         range_hi = exp(or_unit * (b + 1.96 * sd_slope)),
         outcome  = factor(outcome, levels = rev(names(mods))))
print(iwise_range)

p3 <- ggplot(iwise_range, aes(y = outcome)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_segment(aes(x = range_lo, xend = range_hi, yend = outcome),
               linewidth = 9, colour = col_band, lineend = "round") +   # country-year range
  geom_segment(aes(x = lo, xend = hi, yend = outcome),
               linewidth = 1.6, colour = col_title) +                  # 95% CI
  geom_point(aes(x = OR), size = 3.5, shape = 21, fill = col_title, colour = "white", stroke = 1) +
  geom_text(aes(x = OR, label = sprintf("%.3f", OR)), nudge_y = 0.32,
            fontface = "bold", colour = col_title, size = 4, family = font) +
  geom_text(aes(x = range_lo, label = sprintf("%.3f", range_lo)), hjust = 1.35,
            size = 3.6, colour = "grey30", family = font) +
  geom_text(aes(x = range_hi, label = sprintf("%.3f", range_hi)), hjust = -0.35,
            size = 3.6, colour = "grey30", family = font) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  labs(x = if (or_unit == 1) "IWISE OR per point" else paste("IWISE OR per", or_unit, "points"),
       y = NULL) +
  theme(axis.text.y = element_text(size = 13, colour = "grey15"),
        axis.line.x = element_line(colour = "grey70"))
p3
ggsave("figures/fig3_iwise_range.png", p3, width = 6.5, height = 2.8, dpi = 300, bg = "white")


# ---- 6. FIGURE 4: IWISE EFFECT FOR EACH COUNTRY-YEAR ------------------------
# Country-year OR = exp(beta + u_j), u_j = predicted random slope (conditional
# mode). Approximate 95% interval: SE = sqrt(SE(beta)^2 + conditional variance
# of u_j). It ignores the covariance of beta and u_j, so read it as a guide.
# Labels: country name, plus the country_year code where a country has
# several waves. If you have a survey-year column, set v.year to its name.
v.country_name <- "country_name"
v.year         <- NULL                          # e.g. "YEAR_WAVE"

cy_labels <- function(d) {
  lab <- d %>% distinct(country_year, .keep_all = TRUE)
  name <- if (v.country_name %in% names(lab)) as.character(lab[[v.country_name]]) else as.character(lab$country_year)
  if (!is.null(v.year)) name <- paste(name, lab[[v.year]])
  dup <- name %in% name[duplicated(name)]                   # same name twice -> add the code
  name[dup] <- paste0(name[dup], " (", lab$country_year[dup], ")")
  tibble(country_year = as.character(lab$country_year), label = name)
}

by_cy <- bind_rows(lapply(names(mods), function(o) {
  m  <- mods[[o]]
  re <- ranef(m)$country_year                               # predicted random effects
  cv <- condVar(m)$country_year                             # their conditional variances
  cf <- coef(summary(m))["iwisescore", ]
  tibble(outcome      = o,
         country_year = rownames(re),
         u            = re$iwisescore,
         se           = sqrt(cf["Std. Error"]^2 + cv$iwisescore)) %>%
    mutate(OR = exp(or_unit * (cf["Estimate"] + u)),
           lo = exp(or_unit * (cf["Estimate"] + u - 1.96 * se)),
           hi = exp(or_unit * (cf["Estimate"] + u + 1.96 * se))) %>%
    left_join(cy_labels(if (o == "FLI") d.fl else d.inc), by = "country_year")
}))
write.csv(by_cy, "iwise_by_country.csv", row.names = FALSE)

by_cy %>% group_by(outcome) %>%                             # how many country-years below / above 1
  summarise(n = n(), OR_below_1 = sum(OR < 1), CI_below_1 = sum(hi < 1), CI_above_1 = sum(lo > 1))

plot_cy <- function(o, title) {
  d   <- filter(by_cy, outcome == o) %>% mutate(label = reorder(label, OR))   # sorted by OR
  avg <- filter(iwise_range, outcome == o)
  ggplot(d, aes(OR, label)) +
    geom_blank() +                                          # sets up the axes before the shading
    annotate("rect", xmin = avg$range_lo, xmax = avg$range_hi, ymin = -Inf, ymax = Inf,
             fill = col_shade) +                            # 95% range from Figure 3
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = avg$OR, colour = col_band, linewidth = 1) +       # average OR
    geom_linerange(aes(xmin = lo, xmax = hi), colour = "grey55") +            # 95% interval
    geom_point(size = 1.8, colour = col_title) +
    labs(title = title, x = if (or_unit == 1) "IWISE OR per point (95% interval)"
                            else paste("IWISE OR per", or_unit, "points (95% interval)"),
         y = NULL) +
    theme(axis.text.y = element_text(size = 7.5, colour = "grey15"),
          axis.line.x = element_line(colour = "grey70"))
}

p4 <- plot_cy("FLI", "Financial Life Index") | plot_cy("Income", "Income quintile")
p4
ggsave("figures/fig4_iwise_by_country.png", p4, width = 11, height = 13, dpi = 300, bg = "white")


# ---- 6b. FIGURE 4b: ALL COUNTRY-YEARS, ONLY THE EXTREMES LABELLED -----------
# For slides: same estimates as Figure 4, sorted from lowest (bottom) to highest (top) OR.
# Only the n_lab lowest and n_lab highest country-years are named.
n_lab <- 5                                                  # labels at each end

plot_cy_slim <- function(o, title) {
  d <- by_cy %>% filter(outcome == o) %>% arrange(OR) %>%
    mutate(rank = row_number(),                             # position after sorting
           lab  = ifelse(rank <= n_lab | rank > n() - n_lab, label, ""))   # only extremes named
  avg <- filter(iwise_range, outcome == o)
  ggplot(d, aes(OR, rank)) +
    annotate("rect", ymin = -Inf, ymax = Inf, xmin = avg$range_lo, xmax = avg$range_hi,
             fill = col_shade) +                            # 95% range from Figure 3
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +    # no association
    geom_vline(xintercept = avg$OR, colour = col_band, linewidth = 1) +     # average OR
    geom_linerange(aes(xmin = lo, xmax = hi), colour = "grey60") +          # 95% interval
    geom_point(size = 1.8, colour = col_title) +
    geom_text(aes(x = hi, label = lab), hjust = -0.1,                       # label right of the interval
              size = 3, colour = "grey20", family = font) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +           # room for the labels
    labs(title = title, y = "Country-years, sorted by OR",
         x = if (or_unit == 1) "IWISE OR per point" else paste("IWISE OR per", or_unit, "points")) +
    theme(axis.text.y = element_blank(),
          axis.line.x = element_line(colour = "grey70"))
}

p4b <- plot_cy_slim("FLI", "Financial Life Index") | plot_cy_slim("Income", "Income quintile")   # side by side
p4b
ggsave("figures/fig4b_iwise_by_country_slim.png", p4b, width = 11, height = 7, dpi = 300, bg = "white")


# ---- 6c. FIGURE 4c: COUNTRY-YEARS BY INCOME GROUP ---------------------------
# One dot per country-year, grouped by World Bank income group.
# CHECK the group labels first: table(d.fl$country_income_group).
# If they are written differently, adjust inc_lev (order: poorest first).
inc_lev <- c("Low income", "Lower middle income", "Upper middle income", "High income")
table(d.fl$country_income_group, useNA = "ifany")

grp <- bind_rows(FLI    = distinct(d.fl,  country_year, country_income_group),   # one row per country-year
                 Income = distinct(d.inc, country_year, country_income_group),
                 .id = "outcome") %>%
  mutate(country_year = as.character(country_year))         # same type as in by_cy

by_grp <- by_cy %>%
  left_join(grp, by = c("outcome", "country_year")) %>%
  mutate(country_income_group = factor(country_income_group,
                                       levels = rev(intersect(inc_lev, unique(country_income_group)))),
         outcome = recode(outcome, FLI = "Financial Life Index", Income = "Income quintile"))

by_grp %>% group_by(outcome, country_income_group) %>%     # numbers behind the figure
  summarise(n = n(), median_OR = median(OR), .groups = "drop") %>% print()

p4c <- ggplot(by_grp, aes(OR, country_income_group)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +      # no association
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = col_shade, colour = "grey60") +
  geom_jitter(height = 0.12, width = 0, size = 1.8, colour = col_title, alpha = 0.8) +  # one dot per country-year
  facet_wrap(~ outcome) +
  labs(x = if (or_unit == 1) "IWISE OR per point (one dot = one country-year)"
           else paste("IWISE OR per", or_unit, "points (one dot = one country-year)"),
       y = NULL) +
  theme(strip.text  = element_text(colour = col_title, face = "bold", size = 13, hjust = 0),
        axis.text.y = element_text(size = 11, colour = "grey15"),
        axis.line.x = element_line(colour = "grey70"),
        panel.spacing = unit(1.5, "lines"))
p4c
ggsave("figures/fig4c_iwise_by_income_group.png", p4c, width = 10, height = 4.5, dpi = 300, bg = "white")


# ---- 7. CHECK: REFIT WITH IWISE IN 12-POINT STEPS ---------------------------
# iwise12 = iwisescore / 12, so its OR is the OR per 12 points directly.
# The model is the same, only rescaled. Expected:
#   exp(beta_iwise12) = exp(12 * beta_iwisescore), same CI
#   SD of the random slope x 12, same logLik, same probabilities.
# Differences beyond rounding = an error in the calculation, or a model that
# hasn't converged. Same data, weights and settings as 08.1.
d.fl  <- mutate(d.fl,  iwise12 = iwisescore / 12)
d.inc <- mutate(d.inc, iwise12 = iwisescore / 12)

ctrl <- clmm.control(method = "nlminb", maxIter = 200, maxLineIter = 200,
                     iter.max = 1000, eval.max = 2000)       # same as 08.1

specs12 <- lapply(fits, function(m) {                       # same formula, iwisescore -> iwise12
  f <- gsub("\\biwisescore\\b", "iwise12", deparse1(formula(m)))
  list(formula = f, data = if (all.vars(formula(m))[1] == "FLI_3item") "d.fl" else "d.inc")
})
for (nm in names(specs12)) cat(sprintf("%-8s %s\n", nm, specs12[[nm]]$formula))   # check

fit12 <- function(s, name, datasets) {
  t0 <- Sys.time()
  m  <- tryCatch(clmm(as.formula(s$formula), data = datasets[[s$data]], weights = w_fit,
                      link = "logit", Hess = TRUE, nAGQ = 1, control = ctrl),
                 error = function(e) e)
  attr(m, "runtime_min") <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  saveRDS(m, paste0("fit12_", name, ".rds"))                # saved straight away
  m
}

if (file.exists("models_iwise12.rds")) {                    # already refitted: just load
  fits12 <- readRDS("models_iwise12.rds")
} else {                                                    # slow: one model per core
  plan(multisession, workers = min(length(specs12), parallelly::availableCores(omit = 1)))
  fits12 <- future_imap(specs12, fit12, datasets = list(d.fl = d.fl, d.inc = d.inc),
                        .options = furrr_options(seed = TRUE, packages = "ordinal"))
  plan(sequential)
  saveRDS(fits12, "models_iwise12.rds")
}

check12 <- bind_rows(lapply(names(fits), function(nm) {
  m1 <- fits[[nm]]; m12 <- fits12[[nm]]
  if (inherits(m12, "error")) return(tibble(model = nm, error = conditionMessage(m12)))
  c1  <- coef(summary(m1))["iwisescore", ]; c12 <- coef(summary(m12))["iwise12", ]
  sd1  <- attr(VarCorr(m1)$country_year,  "stddev")["iwisescore"]
  sd12 <- attr(VarCorr(m12)$country_year, "stddev")["iwise12"]
  p_diff <- max(abs(pred_ref(m1, iwise_vals)$prob -                         # probabilities at 0, 12, 24, 36
                    pred_ref(m12, iwise_vals, x = "iwise12", per = 12)$prob))
  tibble(model         = nm,
         OR12_calc     = exp(12 * c1["Estimate"]),                          # from the original model
         OR12_refit    = exp(c12["Estimate"]),                              # from the refit
         lo_calc       = exp(12 * (c1["Estimate"] - 1.96 * c1["Std. Error"])),
         lo_refit      = exp(c12["Estimate"] - 1.96 * c12["Std. Error"]),
         hi_calc       = exp(12 * (c1["Estimate"] + 1.96 * c1["Std. Error"])),
         hi_refit      = exp(c12["Estimate"] + 1.96 * c12["Std. Error"]),
         sd_slope_calc = 12 * unname(sd1),                                  # NA for m.inc.0 (no slope)
         sd_slope_refit = unname(sd12),
         logLik_orig   = as.numeric(logLik(m1)),
         logLik_refit  = as.numeric(logLik(m12)),
         max_prob_diff = p_diff,                                            # largest difference, 0-1 scale
         optimiser     = m12$optRes$message,
         error         = NA_character_)
}))
print(as.data.frame(check12), digits = 5)
write.csv(check12, "check_iwise12.csv", row.names = FALSE)
# Read: OR, CI and SD columns should agree to about 3 decimals, logLik to
# about 0.1, max_prob_diff < 0.001. Larger differences: check the optimiser
# message of both fits (Section 1 and here).