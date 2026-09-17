# =============================================================================
# Adding a 5-year SDG 6.1 progress rate to d.iwise_gov
# -----------------------------------------------------------------------------
# Input  : d.iwise_gdp_gov  -- individual-level, ~91,166 rows, 76 countries,
#                          keys COUNTRY_ISO3 (chr) and YEAR_CALENDAR (int)
# Output : d.iwise_jmp  -- same rows, plus JMP level + progress-rate columns
#
# Source : WHO/UNICEF Joint Monitoring Programme, Household World file
#          https://washdata.org/data/country/WLD/household/download
#          2025 update, series 2000-2024.
#
# Window : For collection year t, the rate is computed over [t-4, t] (5 annual
#          values). Where t > 2024 the window slides back to the last five
#          available years, i.e. 2020-2024, and jmp_lag records the gap.
#          Never extrapolated forward.
#
# NOTE   : JMP revises the ENTIRE back-series at each annual update. Record
#          your download date; a rerun next year will give different numbers.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. PACKAGES
# -----------------------------------------------------------------------------
# install.packages(c("readxl", "dplyr", "tidyr", "purrr", "stringr"))
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
d.iwise_gdp_gov <- readRDS("data/d_iwise_gdp_gov.rds")


 
# -----------------------------------------------------------------------------
# 2. THE FILE
#    Point this at the copy you already downloaded. No need to re-download.
#    If you ever do need to fetch it again, JMP is slow -- raise the timeout:
#      options(timeout = 600)
#      download.file("https://washdata.org/data/country/WLD/household/download",
#                    jmp_path, mode = "wb")
# -----------------------------------------------------------------------------
jmp_path <- "data/JMP_2025_WLD.xlsx"     # <-- EDIT to wherever your file is
stopifnot(file.exists(jmp_path))
# Record the vintage for your methods section.
jmp_download_date <- file.info(jmp_path)$mtime
jmp_download_date
 
 
# -----------------------------------------------------------------------------
# 3. READ THE RAW SHEET, NOT THE DISPLAY SHEET
#
#    The workbook has 9 sheets. "Water" is a FORMATTED DISPLAY sheet: every
#    cell is a formula that reads the raw sheet and converts it to a string --
#    "<1", ">99", "NA", "-". Parsing it means undoing that mangling.
#
#    "wat" is the underlying raw sheet: ONE header row, machine-readable
#    column names, plain numerics, NA for missing. Use it.
#
#    Naming convention in "wat":
#      _r = rural   _u = urban   _t = TOTAL (national)  <- you want _t
#      wat_basal_t      at least basic, national
#      wat_sm_t         safely managed, national
#      wat_imp_prem_t   improved + accessible on premises
#      wat_imp_av_t     improved + available when needed
#      wat_imp_qual_t   improved + free from contamination
#      wat_arc_sm_t     JMP's own annual rate of change (see section 3b)
# -----------------------------------------------------------------------------
excel_sheets(jmp_path)
 
jmp <- read_excel(jmp_path, sheet = "wat") |>
  select(iso3, country = name, year,
         bas   = wat_basal_t,
         sm    = wat_sm_t,
         prem  = wat_imp_prem_t,
         avail = wat_imp_av_t,
         qual  = wat_imp_qual_t,
         arc_sm  = wat_arc_sm_t,
         arc_bas = arc_wat_basal_t) |>
  mutate(iso3 = toupper(str_squish(iso3)),
         year = as.integer(year)) |>
  filter(!is.na(iso3), !is.na(year))
 
# CHECK A -- expected: 5,661 rows; 232 ISO3 codes; years 2000-2024.
dim(jmp)
n_distinct(jmp$iso3)
range(jmp$year)
summary(jmp[c("bas", "sm", "prem", "avail", "qual")]) #tuttobene
 
# CHECK A2 -- safely managed must not exceed any of its three components.
# NOTE: it is NOT always equal to the minimum. It is the joint proportion
# meeting all three criteria, so it is bounded above by the minimum and equals
# it in about 81% of rows, falling up to ~13pp below it elsewhere.
jmp |>
  filter(if_all(c(sm, prem, avail, qual), ~!is.na(.x))) |>
  mutate(excess = sm - pmin(prem, avail, qual)) |>
  summarise(n = n(), max_excess = max(excess),      # must be <= 0
            pct_equal = mean(abs(excess) < 0.01) * 100) #tuttobene
 
 
# -----------------------------------------------------------------------------
# 3b. WHY JMP'S OWN RATE COLUMN CANNOT BE USED
#     arc_sm and arc_bas are constant within country. Verified numerically
#     against this file: they equal (value_2024 - value_2000) / 24, i.e. a
#     full-period average over the whole 2000-2024 series.
#     That is not a recent rate and cannot be matched to a collection year.
#     Keep them only as a long-run comparator for your own 5-year figures.
# -----------------------------------------------------------------------------
jmp |> #check that proves arc_sm is a single constant per country rather than a year-by-year rate
  filter(!is.na(arc_sm)) |>
  group_by(iso3) |>
  summarise(n_distinct_arc = n_distinct(round(arc_sm, 6))) |>
  count(n_distinct_arc)                # expect all 1
 
 
# -----------------------------------------------------------------------------
# 4. CHECK YOUR OWN COUNTRIES ARE PRESENT
# -----------------------------------------------------------------------------
my_iso <- sort(unique(d.iwise_gdp_gov$iso3c))
length(my_iso)                                     # expect 76
 
# CHECK B -- countries of yours with no row in the JMP file at all.
setdiff(my_iso, unique(jmp$iso3)) #all there
 
# Kosovo and other non-UN entities are the usual absentees. If one of yours
# appears here under a code JMP writes differently, patch it rather than
# losing the country:
# jmp <- jmp |> mutate(iso3 = recode(iso3, "XKO" = "XKX"))
 
 
# -----------------------------------------------------------------------------
# 5. HOW MANY OF *YOUR* COUNTRIES HAVE A SAFELY-MANAGED SERIES?
#    Run this before writing any code that depends on it. This is the
#    decision point for the whole variable.
# -----------------------------------------------------------------------------
my_iso <- sort(unique(d.iwise_gdp_gov$iso3c))
length(my_iso)                                          # expect 76
 
coverage <- jmp |>
  filter(iso3 %in% my_iso, year >= 2020, year <= 2024) |>
  group_by(iso3) |>
  summarise(n_sm  = sum(!is.na(sm)),
            n_bas = sum(!is.na(bas)), .groups = "drop")
 
# CHECK D -- countries in your sample with NO safely-managed data at all.
# BENCHMARK from this file, across all 223 countries with a 2020-2024 window:
#   safely managed  159 complete, 61 with nothing
#   at least basic  216 complete,  2 with nothing
# Your 76 skew towards LMICs, so expect WORSE than 61/223 = 27% missing.
missing_sm <- coverage |> filter(n_sm == 0) |> pull(iso3)
length(missing_sm); missing_sm

# countries missing sm: arg, aus, bfa, bol, chn, cmr, com, gab, gin, jam, ken, 
# lbr, mli, mrt, mus, nam, ner, pan, slv, som, ury, ven, zmb

 
coverage |> summarise(full_sm = sum(n_sm == 5), none_sm = sum(n_sm == 0),
                      full_bas = sum(n_bas == 5), none_bas = sum(n_bas == 0))
 

coverage |> filter(n_sm > 0, n_sm < 5)      # all good (COG partial due to coding)
coverage |> filter(n_bas > 0, n_bas < 5)    # all good (COG partial due to coding)
coverage |> filter(n_bas == 0)              # No basic = ARG (Argentina)



# -----------------------------------------------------------------------------
# 6. THE RATE FUNCTION
#    Given a country's series and a target end-year, returns the window used
#    and three progress measures.
# -----------------------------------------------------------------------------
progress_rate <- function(df, target_year, value_col, window = 5) { #target year = GWP year
 
  s <- df |> #clean data, drop years where this particular column is missing, sort ascending
    filter(!is.na(.data[[value_col]])) |>
    arrange(year)
 
  if (nrow(s) == 0) return(tibble(win_start = NA_integer_, win_end = NA_integer_,
                                  n_yrs = 0L, level = NA_real_, rate_pp = NA_real_,
                                  rate_ols = NA_real_, gapclose = NA_real_)) #Returns the all-NA row rather than erroring (ARG)
 
  # End of window: the collection year, or the latest available year if the
  # series stops earlier. NEVER later than the collection year.
  end <- max(s$year[s$year <= target_year], na.rm = TRUE)
  if (is.infinite(end)) return(tibble(win_start = NA_integer_, win_end = NA_integer_,
                                      n_yrs = 0L, level = NA_real_, rate_pp = NA_real_,
                                      rate_ols = NA_real_, gapclose = NA_real_))
 
  
  
  start <- end - (window - 1)
  w <- s |> filter(year >= start, year <= end)
 
  if (nrow(w) < 2) return(tibble(win_start = min(w$year), win_end = end,
                                 n_yrs = nrow(w), level = w[[value_col]][nrow(w)],
                                 rate_pp = NA_real_, rate_ols = NA_real_,
                                 gapclose = NA_real_)) #bail out if data too short
 
  y0 <- w[[value_col]][1]; y1 <- w[[value_col]][nrow(w)]
  t0 <- w$year[1];         t1 <- w$year[nrow(w)]
  span <- t1 - t0
 
  # (a) Absolute change, percentage points per year.
  rate_pp <- (y1 - y0) / span
 
  # (b) OLS slope across all points in the window.
  rate_ols <- unname(coef(lm(w[[value_col]] ~ w$year))[2])
 
  # (c) Annual fractional closure of the gap to 100%. Ceiling-free.
  #     NA where the country was already at 100% at window start.
  gap0 <- 100 - y0
  gapclose <- if (gap0 <= 0) NA_real_ else 1 - ((100 - y1) / gap0)^(1 / span)
 
  tibble(win_start = t0, win_end = t1, n_yrs = nrow(w),
         level = y1, rate_pp = rate_pp, rate_ols = rate_ols, gapclose = gapclose)
}
 
 
# -----------------------------------------------------------------------------
# 7. BUILD THE COUNTRY x COLLECTION-YEAR LOOKUP
#    One row per (country, collection year) pair actually present in your data.
#    PSE contributes two rows, as before.
# -----------------------------------------------------------------------------
keys <- d.iwise_gdp_gov |>
  distinct(iso3c, year) |>
  arrange(iso3c, year)
 
nrow(keys)   # expect 77 if PSE is the only two-wave country
 
# Run the same rate function over all five series. sm and bas are the
# analysis variables; prem/avail/qual are diagnostics for interpreting sm.
series <- c("sm", "bas", "prem", "avail", "qual")
 
jmp_lookup <- keys |>
  mutate(res = map2(iso3c, year, function(iso, yr) {
    d <- jmp |> filter(iso3 == iso)
    map(series, ~progress_rate(d, yr, .x) |> rename_with(function(n) paste0(.x, "_", n))) |>
      bind_cols()
  })) |>
  unnest(res) |>
  mutate(
    jmp_lag = year - sm_win_end,        # 0 = contemporaneous
    # Lowest of the three components at the end of the window -- the binding
    # constraint. Safely managed is bounded above by this, and equals it about
    # 81% of the time. Descriptive only.
    sm_binding = case_when(
      is.na(prem_level) | is.na(avail_level) | is.na(qual_level) ~ NA_character_,
      pmin(prem_level, avail_level, qual_level) == prem_level  ~ "accessibility",
      pmin(prem_level, avail_level, qual_level) == avail_level ~ "availability",
      TRUE                                                      ~ "quality"
    )
  )
 
#CHECK
dim(jmp_lookup)          # 77 x 39
names(jmp_lookup)
jmp_lookup |> select(iso3c, year, sm_level, sm_rate_pp, bas_level, bas_rate_pp)

# CHECK E0 -- distribution of the binding constraint across your countries.
jmp_lookup |> count(sm_binding) #what is the weakest link (the lowest value)
 
# CHECK E -- every 2025 country should show lag 1 (series ends 2024).
jmp_lookup |> count(year, jmp_lag) #makessense
 
# CHECK F -- any window shorter than 5 years?
jmp_lookup |> filter(sm_n_yrs < 5) |> select(iso3c, year, sm_n_yrs) #24 years with 0
 
jmp_lookup |> filter(sm_n_yrs < 5) |> count(sm_n_yrs) #none with 3-4 years

jmp_lookup |> summarise(has_sm = sum(!is.na(sm_rate_pp)),
                        has_bas = sum(!is.na(bas_rate_pp)))

# CHECK G -- ceiling cases: countries at or near 100% throughout.
jmp_lookup |> filter(sm_level >= 99) |>
  select(iso3c, sm_level, sm_rate_pp, sm_gapclose) |> print(n = 40) #GBR and ISR have 99.(/99.5)
 
 
# -----------------------------------------------------------------------------
# 7b. BLANK RATES THAT ARE NOT MEASURED CHANGE
#     Two failure modes, both real in this sample:
#       ARG  no data after 2016 -- a 2021 survey gets a window ending 5 years
#            before fieldwork, and no safely managed series at all
#       COG  series ends 2021 and the final years are a straight-line
#            projection -- a 2-point rate of the model, not of the country
#     Rather than dropping these by name, apply a uniform rule. Each series
#     gets its OWN lag, because sm and bas can end in different years.
#
#     SET THESE DELIBERATELY -- they are analytic choices you must justify.
# -----------------------------------------------------------------------------
MAX_LAG  <- 2   # years between window end and fieldwork; 1 is normal for 2025
MIN_YRS  <- 3   # annual values needed before a slope means anything
 
for (p in series) {
  lag_col <- paste0(p, "_lag")
  jmp_lookup[[lag_col]] <- jmp_lookup$year - jmp_lookup[[paste0(p, "_win_end")]]
 
  bad <- (jmp_lookup[[lag_col]] > MAX_LAG) |
         (jmp_lookup[[paste0(p, "_n_yrs")]] < MIN_YRS)
  bad[is.na(bad)] <- FALSE          # already-NA rows need no blanking
 
  for (m in c("rate_pp", "rate_ols", "gapclose")) {
    jmp_lookup[[paste0(p, "_", m)]][bad] <- NA_real_
  }
}
 
# NOTE: _level, _win_start, _win_end and _n_yrs are deliberately NOT blanked.
# The level is still a valid (if stale) coverage figure, and the window
# columns are the audit trail showing why the rate was dropped.
 
# CHECK G2 -- which countries lost a rate, and why.
jmp_lookup |>
  filter(is.na(sm_rate_pp) | is.na(bas_rate_pp)) |>
  select(iso3c, year, sm_win_end, sm_n_yrs, sm_lag, sm_rate_pp,
         bas_win_end, bas_n_yrs, bas_lag, bas_rate_pp) |>
  arrange(desc(bas_lag)) |>
  print(n = 40) #24 are missing because no SM and ARG
 
# CHECK G3 -- final usable counts. This is the number your model runs on.
jmp_lookup |>
  summarise(n = n(),
            has_sm_rate  = sum(!is.na(sm_rate_pp)),
            has_bas_rate = sum(!is.na(bas_rate_pp)))
 
 
# -----------------------------------------------------------------------------
# 8. MERGE
# -----------------------------------------------------------------------------
d.iwise_gdp_gov_jmp <- d.iwise_gdp_gov |>
  left_join(jmp_lookup, by = c("iso3c", "year"),
            relationship = "many-to-one")

# CHECK H -- row count MUST be unchanged.
nrow(d.iwise_gdp_gov_jmp); nrow(d.iwise_gdp_gov_jmp) #good
stopifnot(nrow(d.iwise_gdp_gov) == nrow(d.iwise_gdp_gov_jmp)) #good
 
# CHECK I -- missingness at COUNTRY level, not respondent level.
d.iwise_gdp_gov_jmp |>
  distinct(iso3c, sm_level, sm_rate_pp, bas_level, bas_rate_pp) |>
  summarise(across(everything(), ~sum(is.na(.x)))) #correct 24 sm and one base rate (ARG)
 
# CHECK J -- do the two rate definitions agree on ranking?
with(distinct(d.iwise_gdp_gov_jmp, iso3c, sm_rate_pp, sm_gapclose),
     cor(sm_rate_pp, sm_gapclose, use = "complete.obs", method = "spearman")) #0.896

# High correlation means the ceiling isn't biting hard and either will do.
# Low correlation means the choice of definition changes your results --
# report both.
 
# CHECK K -- kinks. rate_pp uses only the endpoints; rate_ols uses all five.
# Measured on this file for the 2020-2024 window: the two agree almost
# perfectly (Spearman 0.998; only 1 of 159 countries differs by more than
# 0.15pp/yr). So kinks are NOT a practical concern here -- but run it, because
# the one country that does bend is worth knowing about.
jmp_lookup |>
  mutate(kink = abs(sm_rate_pp - sm_rate_ols)) |>
  filter(kink > 0.05) |>
  select(iso3c, year, sm_win_start, sm_win_end,
         sm_rate_pp, sm_rate_ols, kink, sm_binding) |>
  arrange(desc(kink)) |>
  print(n = 30) #gives 0, can be ignored
 
# CHECK L -- VARIANCE.
jmp_lookup |> #for sm
  summarise(n = sum(!is.na(sm_rate_pp)),
            median = median(sm_rate_pp, na.rm = TRUE),
            iqr_lo = quantile(sm_rate_pp, .25, na.rm = TRUE),
            iqr_hi = quantile(sm_rate_pp, .75, na.rm = TRUE),
            pct_flat = mean(abs(sm_rate_pp) < 0.05, na.rm = TRUE) * 100,
            pct_zero = mean(abs(sm_rate_pp) < 1e-9, na.rm = TRUE) * 100)
 
hist(jmp_lookup$sm_rate_pp, breaks = 30,
     main = "5-year safely managed rate, pp/year", xlab = "pp/year")
#Roughly unimodal around 0.2, right-skewed, with a handful of countries above 1.5 and a couple below
# Median: 0.230




jmp_lookup |> #for basic
  summarise(n = sum(!is.na(bas_rate_pp)),
            median = median(bas_rate_pp, na.rm = TRUE),
            iqr_lo = quantile(bas_rate_pp, .25, na.rm = TRUE),
            iqr_hi = quantile(bas_rate_pp, .75, na.rm = TRUE),
            pct_flat = mean(abs(bas_rate_pp) < 0.05, na.rm = TRUE) * 100,
            pct_zero = mean(abs(bas_rate_pp) < 1e-9, na.rm = TRUE) * 100)

hist(jmp_lookup$bas_rate_pp, breaks = 30,
     main = "5-year at least basic rate, pp/year", xlab = "pp/year")

# 
# Median: 0.338


# CHECK M -- your own 5-year rate against JMP's published 24-year average.
# They measure different periods, so they will not match; a country where they
# have OPPOSITE SIGNS has changed direction since 2000, which is exactly the
# kind of country your hypothesis is about.
#jmp_lookup |>
#  left_join(distinct(jmp, iso3, arc_sm), by = c("iso3c" = "iso3")) |>
#  filter(!is.na(sm_rate_pp), !is.na(arc_sm), sign(sm_rate_pp) != sign(arc_sm)) |>
#  select(iso3c, sm_rate_pp, arc_sm, sm_level) |>
#  arrange(sm_rate_pp) |>
#  print(n = 40)
 
# For any country worth a closer look, plot the raw series.
# plot_country <- function(iso, y0 = 2014) {
#   d <- jmp |> filter(iso3 == iso, year >= y0)
#   matplot(d$year, d[c("sm", "prem", "avail", "qual")], type = "l",
#           lty = 1, col = 1:4, xlab = "year", ylab = "%", main = iso)
#   legend("bottomright", c("safely managed", "on premises", "available",
#          "quality"), lty = 1, col = 1:4, bty = "n", cex = 0.8)
# }
# plot_country("XXX")
 

 
# =============================================================================
# 9. SAFELY MANAGED vs AT LEAST BASIC -- material for the supervisor meeting
#    Both variables already exist. This section does not create anything new;
#    it summarises them side by side so the choice can be argued rather than
#    assumed. Everything works on the 76 country-year rows, NOT on the
#    91,166 respondent rows -- the covariate varies only between countries.
# =============================================================================
 
# --- 9a. Headline comparison ------------------------------------------------
compare <- jmp_lookup |>
  summarise(
    across(c(sm_rate_pp, bas_rate_pp),
           list(n        = ~sum(!is.na(.x)),
                missing  = ~sum(is.na(.x)),
                median   = ~median(.x, na.rm = TRUE),
                iqr_lo   = ~quantile(.x, .25, na.rm = TRUE),
                iqr_hi   = ~quantile(.x, .75, na.rm = TRUE),
                min      = ~min(.x, na.rm = TRUE),
                max      = ~max(.x, na.rm = TRUE),
                sd       = ~sd(.x, na.rm = TRUE),
                pct_flat = ~mean(abs(.x) < 0.05, na.rm = TRUE) * 100,
                pct_zero = ~mean(abs(.x) < 1e-9, na.rm = TRUE) * 100,
                pct_neg  = ~mean(.x < 0, na.rm = TRUE) * 100))) |>
  pivot_longer(everything(),
               names_to = c("variable", "stat"),
               names_pattern = "^(sm|bas)_rate_pp_(.*)$") |>
  pivot_wider(names_from = variable, values_from = value) |>
  mutate(across(c(sm, bas), ~round(.x, 3)))
 
print(compare, n = 20)
 
# --- 9b. Ceiling exposure ---------------------------------------------------
# How many countries sit so high that the rate cannot move much.
jmp_lookup |>
  summarise(sm_ge99  = sum(sm_level  >= 99, na.rm = TRUE),
            sm_ge95  = sum(sm_level  >= 95, na.rm = TRUE),
            bas_ge99 = sum(bas_level >= 99, na.rm = TRUE),
            bas_ge95 = sum(bas_level >= 95, na.rm = TRUE))
 
# --- 9c. Do the two rank countries the same way? ----------------------------
# High correlation: the choice is presentational. Low: it is substantive.
with(jmp_lookup, cor(sm_rate_pp, bas_rate_pp,
                     use = "complete.obs", method = "spearman"))
 
with(jmp_lookup, cor(sm_rate_pp, bas_rate_pp,
                     use = "complete.obs", method = "pearson"))
 
plot(jmp_lookup$bas_rate_pp, jmp_lookup$sm_rate_pp,
     xlab = "at least basic, pp/year", ylab = "safely managed, pp/year",
     main = "Do the two rungs agree?"); abline(0, 1, lty = 2)
 
# --- 9d. Who is missing from which -----------------------------------------
jmp_lookup |>
  mutate(status = case_when(
    !is.na(sm_rate_pp) &  !is.na(bas_rate_pp) ~ "both",
     is.na(sm_rate_pp) &  !is.na(bas_rate_pp) ~ "basic only",
    !is.na(sm_rate_pp) &   is.na(bas_rate_pp) ~ "safely managed only",
    TRUE                                      ~ "neither")) |>
  count(status)
 
# --- 9e. Country table to take to the meeting -------------------------------
regions <- read_excel(jmp_path, sheet = "wat") |>
  distinct(iso3, country = name, region_sdg, region_income)


overview_table_sdg <- jmp_lookup |>
  left_join(regions, by = c("iso3c" = "iso3")) |>
  select(iso3c, country, region_sdg, region_income, year,
         sm_level, sm_rate_pp, sm_win_start, sm_win_end, sm_n_yrs, sm_lag,
         bas_level, bas_rate_pp, bas_win_start, bas_win_end, bas_n_yrs, bas_lag,
         sm_binding) |>
  arrange(is.na(sm_rate_pp), sm_rate_pp)
 
print(overview_table_sdg, n = 80)
library(readr)
write_csv(overview_table_sdg, "data/jmp_rate_comparison.csv")
 
# NOTE: region_sdg and region_income need to be carried through section 3.
# If the join above fails, add them to the select() there:
#   select(iso3, country = name, year, ..., region_sdg, region_income)






#delete the ones I don't need after discussion with Karin

library(dplyr)
d.iwise_gdp_gov_jmp <- d.iwise_gdp_gov_jmp |>
  dplyr::select(-dplyr::any_of(c("sm_rate_pp", "sm_rate_ols", "sm_gapclose",
                                 "sm_win_start", "sm_win_end", "sm_n_yrs", "sm_lag",
                                 "bas_rate_ols", "bas_gapclose")),
                -dplyr::matches("^(prem|avail|qual)_(rate_pp|rate_ols|gapclose|win_start|win_end|n_yrs|lag)$"))
 
#double-check
names(d.iwise_gdp_gov_jmp) |> grep(pattern = "^(sm|bas|prem|avail|qual|jmp)_", value = TRUE)
nrow(d.iwise_gdp_gov_jmp)   # must be 91,166
d.iwise_gdp_gov_jmp |>
  dplyr::distinct(iso3c, bas_level, bas_rate_pp, sm_level) |>
  dplyr::summarise(dplyr::across(everything(), ~sum(is.na(.x)))) #values unchanged
summary(d.iwise_gdp_gov_jmp$bas_rate_pp) #sanity check


saveRDS(d.iwise_gdp_gov_jmp, "data/iwise_gdp_gov_jmp.rds")
