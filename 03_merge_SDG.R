# =============================================================================
# Adding a 5-year SDG 6.1 progress rate to d.iwise_gov
# -----------------------------------------------------------------------------
# Input  : d.iwise_gov  -- individual-level, ~91,166 rows, 76 countries,
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
# 2. GET THE FILE
#    Try the programmatic download first. If it fails (the site occasionally
#    blocks non-browser requests), open the URL in a browser, save the .xlsx
#    into data/, and skip to section 3.
# -----------------------------------------------------------------------------
dir.create("data", showWarnings = FALSE)
jmp_path <- "data/JMP_2025_WLD_household.xlsx"

if (!file.exists(jmp_path)) {
  download.file(
    url      = "https://washdata.org/data/country/WLD/household/download",
    destfile = jmp_path,
    mode     = "wb"      # ESSENTIAL on Windows -- without it the xlsx corrupts
  )
}

# Record the vintage for your methods section.
jmp_download_date <- file.info(jmp_path)$mtime
jmp_download_date


# -----------------------------------------------------------------------------
# 3. INSPECT BEFORE PARSING
#    The water sheet has a MULTI-ROW header and three repeated column blocks
#    (national, then rural, then urban). Look at it before trusting any code.
# -----------------------------------------------------------------------------
excel_sheets(jmp_path)

# Print the first 6 rows raw, no header interpretation.
peek <- read_excel(jmp_path, sheet = "Water", col_names = FALSE, n_max = 6)
print(as.data.frame(peek)[, 1:14])

# You are looking for the row that contains "ISO3". Everything above it is
# title/blocking; that row (or the one below) holds the real column names.


# -----------------------------------------------------------------------------
# 4. READ THE WATER SHEET
#    Auto-detects the header row rather than hard-coding a skip, because the
#    JMP shifts the layout between updates.
# -----------------------------------------------------------------------------
raw_all <- read_excel(jmp_path, sheet = "Water", col_names = FALSE,
                      .name_repair = "minimal")

hdr_row <- which(apply(raw_all, 1, function(r)
  any(str_detect(toupper(as.character(r)), "^ISO3$"), na.rm = TRUE)))[1]

stopifnot(!is.na(hdr_row))
hdr_row   # sanity-check this against the peek above

jmp_raw <- read_excel(jmp_path, sheet = "Water",
                      skip = hdr_row - 1, .name_repair = "unique")

names(jmp_raw)[1:20]   # CHECK A -- read these before continuing


# -----------------------------------------------------------------------------
# 5. LOCATE THE COLUMNS
#    Column names repeat across the national / rural / urban blocks. The
#    NATIONAL block comes first, so the FIRST match is the one you want.
#    Verify with Check B before relying on it.
# -----------------------------------------------------------------------------
nm <- names(jmp_raw)

col_iso  <- which(str_detect(toupper(nm), "^ISO3"))[1]
col_year <- which(str_detect(toupper(nm), "^YEAR"))[1]
col_sm   <- which(str_detect(tolower(nm), "safely managed"))[1]
col_bas  <- which(str_detect(tolower(nm), "at least basic"))[1]

# The three components that safely managed is built from. JMP regresses each
# separately and takes the MINIMUM of the three, so these identify which
# constraint is binding in each country and whether that constraint changed.
col_prem  <- which(str_detect(tolower(nm), "premises"))[1]        # accessibility
col_avail <- which(str_detect(tolower(nm), "available"))[1]       # availability
col_qual  <- which(str_detect(tolower(nm), "contamination"))[1]   # quality

# CHECK B -- confirm you grabbed the NATIONAL block, not rural/urban, and that
# all seven indices are non-NA. If any component is NA the column name has
# changed in this vintage; look at names(jmp_raw) and set it by hand.
tibble(role = c("iso3", "year", "safely_managed", "at_least_basic",
                "on_premises", "available", "free_from_contamination"),
       index = c(col_iso, col_year, col_sm, col_bas, col_prem, col_avail, col_qual),
       name  = nm[c(col_iso, col_year, col_sm, col_bas, col_prem, col_avail, col_qual)])

jmp <- jmp_raw |>
  select(iso3  = all_of(col_iso),
         year  = all_of(col_year),
         sm    = all_of(col_sm),
         bas   = all_of(col_bas),
         prem  = all_of(col_prem),
         avail = all_of(col_avail),
         qual  = all_of(col_qual)) |>
  mutate(
    iso3 = as.character(iso3),
    year = suppressWarnings(as.integer(year)),
    # JMP writes "-", "<1", ">99" and blanks. Strip the inequality signs, then
    # coerce. ">99" becomes 99 and "<1" becomes 1: defensible, but note it.
    across(c(sm, bas, prem, avail, qual),
           ~suppressWarnings(as.numeric(str_remove_all(as.character(.x), "[<>]"))))
  ) |>
  filter(!is.na(iso3), !is.na(year))

# CHECK C -- range and coverage of the series.
range(jmp$year, na.rm = TRUE)              # expect 2000-2024
n_distinct(jmp$iso3)                       # expect ~234
summary(jmp[c("sm", "bas", "prem", "avail", "qual")])   # all should sit in [0, 100]

# CHECK C2 -- sm should equal the minimum of the three components wherever all
# three exist. If it doesn't, you have grabbed the wrong columns.
jmp |>
  filter(!is.na(sm), !is.na(prem), !is.na(avail), !is.na(qual)) |>
  mutate(diff = sm - pmin(prem, avail, qual)) |>
  summarise(n = n(), max_abs_diff = max(abs(diff)))     # expect diff ~ 0


# -----------------------------------------------------------------------------
# 6. HOW MANY OF *YOUR* COUNTRIES HAVE A SAFELY-MANAGED SERIES?
#    Run this before writing any code that depends on it. This is the
#    decision point for the whole variable.
# -----------------------------------------------------------------------------
my_iso <- sort(unique(d.iwise_gov$COUNTRY_ISO3))
length(my_iso)                                          # expect 76

coverage <- jmp |>
  filter(iso3 %in% my_iso, year >= 2016) |>
  group_by(iso3) |>
  summarise(n_sm  = sum(!is.na(sm)),
            n_bas = sum(!is.na(bas)), .groups = "drop")

# CHECK D -- countries in your sample with NO safely-managed data at all.
missing_sm <- coverage |> filter(n_sm == 0) |> pull(iso3)
length(missing_sm); missing_sm

# And countries absent from the JMP file entirely.
setdiff(my_iso, unique(jmp$iso3))


# -----------------------------------------------------------------------------
# 7. THE RATE FUNCTION
#    Given a country's series and a target end-year, returns the window used
#    and three progress measures.
# -----------------------------------------------------------------------------
progress_rate <- function(df, target_year, value_col, window = 5) {

  s <- df |>
    filter(!is.na(.data[[value_col]])) |>
    arrange(year)

  if (nrow(s) == 0) return(tibble(win_start = NA_integer_, win_end = NA_integer_,
                                  n_yrs = 0L, level = NA_real_, rate_pp = NA_real_,
                                  rate_ols = NA_real_, gapclose = NA_real_))

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
                                 gapclose = NA_real_))

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
# 8. BUILD THE COUNTRY x COLLECTION-YEAR LOOKUP
#    One row per (country, collection year) pair actually present in your data.
#    PSE contributes two rows, as before.
# -----------------------------------------------------------------------------
keys <- d.iwise_gov |>
  distinct(COUNTRY_ISO3, YEAR_CALENDAR) |>
  arrange(COUNTRY_ISO3, YEAR_CALENDAR)

nrow(keys)   # expect 77 if PSE is the only two-wave country

# Run the same rate function over all five series. sm and bas are the
# analysis variables; prem/avail/qual are diagnostics for interpreting sm.
series <- c("sm", "bas", "prem", "avail", "qual")

jmp_lookup <- keys |>
  mutate(res = map2(COUNTRY_ISO3, YEAR_CALENDAR, function(iso, yr) {
    d <- jmp |> filter(iso3 == iso)
    map(series, ~progress_rate(d, yr, .x) |> rename_with(function(n) paste0(.x, "_", n))) |>
      bind_cols()
  })) |>
  unnest(res) |>
  mutate(
    jmp_lag = YEAR_CALENDAR - sm_win_end,        # 0 = contemporaneous
    # Which of the three components is holding safely managed down at the end
    # of the window. NA where any component is missing.
    sm_binding = case_when(
      is.na(prem_level) | is.na(avail_level) | is.na(qual_level) ~ NA_character_,
      pmin(prem_level, avail_level, qual_level) == prem_level  ~ "accessibility",
      pmin(prem_level, avail_level, qual_level) == avail_level ~ "availability",
      TRUE                                                      ~ "quality"
    )
  )

# CHECK E0 -- distribution of the binding constraint across your countries.
jmp_lookup |> count(sm_binding)

# CHECK E -- every 2025 country should show lag 1 (series ends 2024).
jmp_lookup |> count(YEAR_CALENDAR, jmp_lag)

# CHECK F -- any window shorter than 5 years?
jmp_lookup |> filter(sm_n_yrs < 5) |> select(COUNTRY_ISO3, YEAR_CALENDAR, sm_n_yrs)

# CHECK G -- ceiling cases: countries at or near 100% throughout.
jmp_lookup |> filter(sm_level >= 99) |>
  select(COUNTRY_ISO3, sm_level, sm_rate_pp, sm_gapclose) |> print(n = 40)


# -----------------------------------------------------------------------------
# 9. MERGE
# -----------------------------------------------------------------------------
d.iwise_jmp <- d.iwise_gov |>
  left_join(jmp_lookup, by = c("COUNTRY_ISO3", "YEAR_CALENDAR"),
            relationship = "many-to-one")

# CHECK H -- row count MUST be unchanged.
nrow(d.iwise_gov); nrow(d.iwise_jmp)
stopifnot(nrow(d.iwise_gov) == nrow(d.iwise_jmp))

# CHECK I -- missingness at COUNTRY level, not respondent level.
d.iwise_jmp |>
  distinct(COUNTRY_ISO3, sm_level, sm_rate_pp, bas_level, bas_rate_pp) |>
  summarise(across(everything(), ~sum(is.na(.x))))

# CHECK J -- do the two rate definitions agree on ranking?
with(distinct(d.iwise_jmp, COUNTRY_ISO3, sm_rate_pp, sm_gapclose),
     cor(sm_rate_pp, sm_gapclose, use = "complete.obs", method = "spearman"))
# High correlation means the ceiling isn't biting hard and either will do.
# Low correlation means the choice of definition changes your results --
# report both.

# CHECK K -- kinks. rate_pp uses only the endpoints; rate_ols uses all five
# points. For a straight series they are identical. A gap means the series
# bends inside the window, which for safely managed usually means the binding
# component switched. Those countries need a look before you trust the rate.
jmp_lookup |>
  mutate(kink = abs(sm_rate_pp - sm_rate_ols)) |>
  filter(kink > 0.15) |>
  select(COUNTRY_ISO3, YEAR_CALENDAR, sm_win_start, sm_win_end,
         sm_rate_pp, sm_rate_ols, kink, sm_binding) |>
  arrange(desc(kink)) |>
  print(n = 30)

# For any country flagged above, plot the raw series and look at it.
# plot_country <- function(iso, y0 = 2014) {
#   d <- jmp |> filter(iso3 == iso, year >= y0)
#   matplot(d$year, d[c("sm", "prem", "avail", "qual")], type = "l",
#           lty = 1, col = 1:4, xlab = "year", ylab = "%", main = iso)
#   legend("bottomright", c("safely managed", "on premises", "available",
#          "quality"), lty = 1, col = 1:4, bty = "n", cex = 0.8)
# }
# plot_country("XXX")

saveRDS(d.iwise_jmp, "data/d_iwise_jmp.rds")


# =============================================================================
# METHODS NOTE (adapt and paste)
# -----------------------------------------------------------------------------
# Progress towards SDG 6.1 was measured using national estimates of the
# population using safely managed drinking water services (indicator 6.1.1),
# taken from the WHO/UNICEF Joint Monitoring Programme household World file
# (2025 update, series 2000-2024; downloaded [DATE] from washdata.org).
#
# For each country, a progress rate was computed over the five years ending in
# the year of data collection. Countries surveyed in 2025 (n = X) were assigned
# the 2020-2024 window, since the series ends in 2024. West Bank and Gaza
# contributed two waves, each with its own window.
#
# The primary measure is the annual change in coverage in percentage points,
# (y_end - y_start) / span, which approximates the share of the population
# whose service level changed per year and so corresponds most closely to the
# hypothesised mechanism. The OLS slope across all five values and the annual
# fractional closure of the gap to universal coverage are reported as
# sensitivity analyses. Safely managed estimates were unavailable for n = X
# countries, for which the at-least-basic series was used [OR: which were
# excluded]; the two are not interchangeable and were not pooled.
#
# LIMITATIONS
# 1. Ecological. The covariate is a national aggregate standing in for
#    individual experience; it indicates whether the country changed, not
#    whether the respondent's household did.
# 2. Modelled, not observed. JMP estimates are a least-squares line fitted to
#    available survey points, so the rate recovers the JMP model's own trend
#    rather than five independent observations. No standard error is attached.
# 3. Partly extrapolated. JMP projects the fitted line up to two years beyond
#    the last survey (further under some conditions), and the most recent
#    survey is typically two to six years old. For some countries the rate
#    therefore cannot reflect conditions around fieldwork.
# 4. Non-linear by construction. Safely managed is the minimum of three
#    separately fitted component series, so the trend bends where the binding
#    component switches. n = X countries showed such a kink within their
#    window.
# 5. Ceiling. Countries near universal coverage cannot record a large
#    percentage-point gain; a low rate there means no headroom, not stagnation.
# 6. Vintage-dependent. JMP revises the full back-series at each annual
#    update; figures correspond to the 2025 vintage, downloaded [DATE].
# 7. Direction unspecified. Rapid improvement could either lower reported
#    insecurity or raise expectations faster than conditions improve; the
#    covariate is entered without a directional prediction.
# =============================================================================

