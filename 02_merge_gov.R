# =============================================================================
# Merging the six Worldwide Governance Indicators into d.iwise_gdp
# -----------------------------------------------------------------------------
# Input  : d.iwise_gdp  -- individual-level, ~91,166 rows, 76 countries,
#                          keys `iso3c` (chr) and `year` (int), 2021-2025
# Output : d.iwise_gdp_gov  -- same rows, plus 12 WGI columns and match diagnostics
#
# Vintage: WGI 2025 revision (API last updated 2026-03-18), covering 1996-2024.
#          The 2025 revision retired the old VA.EST / GE.EST / ... codes and
#          replaced them with GOV_WGI_ prefixed IDs. It also recalculated the
#          full history back to 1996, so these values are NOT comparable to
#          WGI editions published before December 2025.
#
# Scales : Two are published side by side, both from the revised series.
#            _est  governance estimate, approx. -2.5 to +2.5
#            _sc   governance score, absolute 0-100, anchored to fixed
#                  benchmark countries
#
# NOTE   : WGI ends in 2024. Every 2025 collection year carries forward from
#          2024. This is expected, not an error.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. PACKAGES
#    The WDI package cannot be used here: it queries source 2 only, and the
#    WGI live in source 3. We call the API directly.
# -----------------------------------------------------------------------------
library(jsonlite)
library(purrr)
library(dplyr)
library(tidyr)

packageVersion("dplyr")   # need >= 1.1.0 for relationship = "many-to-many", 1.2.1

d.iwise_gdp <- readRDS("data/d_iwise_gdp.rds")

# -----------------------------------------------------------------------------
# 2. PULL THE TWELVE SERIES
#    Six dimensions x two scales. Also available per dimension, if you want
#    them later: .SE (standard error), .SC_LB / .SC_UB (90% CI bounds on the
#    score), .SR (number of sources).
# -----------------------------------------------------------------------------
fetch_wgi <- function(code, varname) {
  url <- paste0("https://api.worldbank.org/v2/country/all/indicator/", code,
                "?source=3&format=json&date=2018:2025&per_page=20000")
  res <- fromJSON(url, flatten = TRUE)

  if (is.null(res[[2]]) || length(res[[2]]) == 0) {
    stop("No data returned for ", code, ". Test the URL in a browser.")
  } #same request has to be made twelve times with only the indicator code changing.

  res[[2]] |>
    transmute(iso3c    = countryiso3code,
              wgi_year = as.integer(date),
              !!varname := as.numeric(value)) |>
    filter(iso3c != "")        # blank iso3 marks aggregates on this endpoint
}

codes <- c(
  # standard statistical unit, approx -2.5 to +2.5
  wgi_va_est = "GOV_WGI_VA.EST", wgi_pv_est = "GOV_WGI_PV.EST",
  wgi_ge_est = "GOV_WGI_GE.EST", wgi_rq_est = "GOV_WGI_RQ.EST",
  wgi_rl_est = "GOV_WGI_RL.EST", wgi_cc_est = "GOV_WGI_CC.EST",
  # absolute 0-100 score
  wgi_va_sc  = "GOV_WGI_VA.SC",  wgi_pv_sc  = "GOV_WGI_PV.SC",
  wgi_ge_sc  = "GOV_WGI_GE.SC",  wgi_rq_sc  = "GOV_WGI_RQ.SC",
  wgi_rl_sc  = "GOV_WGI_RL.SC",  wgi_cc_sc  = "GOV_WGI_CC.SC"
)

wgi <- imap(codes, \(code, varname) fetch_wgi(code, varname)) |>
  reduce(full_join, by = c("iso3c", "wgi_year")) |>
  mutate(iso3c = trimws(toupper(iso3c)))


# -----------------------------------------------------------------------------
# 3. CHECK A -- DID ALL TWELVE ARRIVE, AND ON THE RIGHT SCALES?
# -----------------------------------------------------------------------------

wgi |> #overview table of number of observations, min/max and mean
  select(starts_with("wgi_")) |>
  pivot_longer(everything(), names_to = "variable", values_to = "v") |>
  group_by(variable) |>
  summarise(n_obs = sum(!is.na(v)),
            min   = min(v, na.rm = TRUE),
            max   = max(v, na.rm = TRUE),
            mean  = round(mean(v, na.rm = TRUE), 2)) |>
  print(n = 12)

dim(wgi) # ~210 countries x 7 years

# Latest year with data. Expect 2024. The API omits missing years entirely
# rather than returning NA rows, so 2025 will simply be absent.
wgi |> filter(!is.na(wgi_ge_est)) |> summarise(max_year = max(wgi_year))


# -----------------------------------------------------------------------------
# 4. CHECK B -- COVERAGE OF YOUR 76 COUNTRIES
#    WGI is perception-based and does not depend on national accounts, so
#    Venezuela, which had no GDP, should have governance scores.
# -----------------------------------------------------------------------------
my_iso <- d.iwise_gdp |> distinct(iso3c) |> pull(iso3c)
length(my_iso)                                   # expect 76

setdiff(my_iso, unique(wgi$iso3c))               # expect character(0)
#correct


wgi |>
  filter(iso3c %in% c("VEN", "PSE", "LBN"), wgi_year >= 2021) |>
  select(iso3c, wgi_year, wgi_ge_est, wgi_ge_sc) |>
  arrange(iso3c, wgi_year) #ensuring no issues with the countries from gdp problems


# -----------------------------------------------------------------------------
# 5. MATCH ON COUNTRY + COLLECTION YEAR, WITH CARRY-FORWARD
#    For each country-year in your data, take the most recent WGI year at or
#    before the collection year. Never carry backward.
# -----------------------------------------------------------------------------

# 5a. The distinct country-years you need. Small table.
keys <- d.iwise_gdp |> distinct(iso3c, year)
nrow(keys) # 77 (76 + PSE wave 2)

# 5b. Best available WGI year per key. All twelve series are published as a
#     block, so presence of wgi_ge_est stands in for presence of the rest --
#     Check C below verifies that assumption held.
wgi_matched <- keys |> #ach of your 77 country-years matches every WGI year for that country — about 7 rows each
  left_join(wgi, by = "iso3c", relationship = "many-to-many") |> #tells dplyr "yes, I meant to do this"
  filter(wgi_year <= year, !is.na(wgi_ge_est)) |> #then filter for the one we need
  group_by(iso3c, year) |>
  slice_max(wgi_year, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    wgi_lag          = year - wgi_year,
    wgi_match_status = if_else(wgi_lag == 0, "exact year match", "carried forward")
  )

nrow(wgi_matched) == nrow(keys)   # TRUE means nothing was dropped

wgi_lookup <- keys |>
  left_join(wgi_matched, by = c("iso3c", "year")) |>
  mutate(wgi_match_status = replace_na(wgi_match_status, "no data published"))


# -----------------------------------------------------------------------------
# 6. CHECK C -- READ THIS BEFORE MERGING
#    Every 2025 collection year should show lag 1, everything else lag 0.
#    A lag of 2+ means that country is missing from recent rounds.
# -----------------------------------------------------------------------------
wgi_lookup |> count(wgi_match_status, wgi_lag) #1 for 39, 0 for 38

wgi_lookup |>
  filter(wgi_lag > 0 | is.na(wgi_lag)) |> #all the ones from 2025 have carried forwars
  arrange(desc(wgi_lag), iso3c)

# Confirms the block assumption in 5b: any row with an estimate should also
# have all eleven other values. Expect 0.
wgi_lookup |>
  filter(!is.na(wgi_ge_est)) |>
  filter(if_any(starts_with("wgi_"), is.na)) |>
  nrow() #0


# -----------------------------------------------------------------------------
# 7. THE MERGE
# -----------------------------------------------------------------------------
d.iwise_gdp_gov <- d.iwise_gdp |>
  left_join(wgi_lookup, by = c("iso3c", "year"))


# -----------------------------------------------------------------------------
# 8. CHECK D -- POST-MERGE VALIDATION
# -----------------------------------------------------------------------------

# D1. Row count unchanged. Must be TRUE. If FALSE, the join duplicated rows.
nrow(d.iwise_gdp_gov) == nrow(d.iwise_gdp)

# D2. What got added. Expect 12 WGI columns plus wgi_year, wgi_lag,
#     wgi_match_status.
setdiff(names(d.iwise_gdp_gov), names(d.iwise_gdp)) #true #correct

# D3. Missingness at country level, not respondent level.
d.iwise_gdp_gov |>
  distinct(iso3c, year, across(starts_with("wgi_"))) |>
  summarise(across(starts_with("wgi_"), \(x) sum(is.na(x)))) #zero

# D4. Which countries are missing, if any.
d.iwise_gdp_gov |> filter(is.na(wgi_ge_est)) |> distinct(iso3c, year) #zero

# D5. Eyeballing the values. Sort by a dimension whose ordering you already know
#     and confirm it looks plausible.
d.iwise_gdp_gov |>
  distinct(iso3c, year, wgi_year, wgi_ge_est, wgi_ge_sc, wgi_rl_est, wgi_cc_est) |>
  arrange(wgi_ge_est)


# -----------------------------------------------------------------------------
# 9. CHECK E -- HOW THE TWO SCALES RELATE
#  Look at Pearson's correlation. The 0-100 score is anchored to fixed benchmark countries, 
#  so within the revised series it should be close to a monotone transform of the estimate.
#  Expect r near 1. If so, treat them as one variable in two units and pick
#  ONE for any given model. Using both is collinearity by construction.
# -----------------------------------------------------------------------------
d.iwise_gdp_gov |>
  distinct(iso3c, year, wgi_ge_est, wgi_ge_sc) |>
  summarise(r = cor(wgi_ge_est, wgi_ge_sc, use = "complete.obs")) #r=1!

# Visual: If the relationship is linear, the two are interchangeable and
# the choice is purely presentational. If it bends, the score is a nonlinear
# rescaling.
plot(d.iwise_gdp_gov$wgi_ge_est, d.iwise_gdp_gov$wgi_ge_sc,
     xlab = "estimate", ylab = "0-100 score", pch = 16, cex = 0.4)
#perfectly linear

# -----------------------------------------------------------------------------
# 10. CHECK F -- COLLINEARITY ACROSS THE SIX DIMENSIONS
#     Source variables map to up to two dimensions, so the six correlate very
#     strongly by construction. Pairwise r above 0.9 is normal. Do not put all
#     six in one regression.
# -----------------------------------------------------------------------------
d.iwise_gdp_gov |>
  distinct(iso3c, year, across(ends_with("_est"))) |>
  select(ends_with("_est")) |>
  cor(use = "pairwise.complete.obs") |>
  round(2) #often very high, lowest is between rq & pv with 0.65

# -----------------------------------------------------------------------------
# 11. SAVE
# -----------------------------------------------------------------------------
saveRDS(d.iwise_gdp_gov, "data/d_iwise_gov.rds")
