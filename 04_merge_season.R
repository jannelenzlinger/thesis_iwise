# =============================================================================
# 05_add_season.R
#
# Adds season-at-interview and recent-rainfall variables from CHIRPS to
#   d.iwise_gdp_gov_jmp   ->   d.iwise_analysis_complete
#
# Produces, per country x calendar month:
#   - climatological season position (continuous + categorical)
#   - a homogeneity diagnostic saying whether a country-level season label
#     is defensible at all
# and per country x year x month of interview:
#   - actual rainfall, percent of normal, standardised anomaly
#
# CHIRPS v2.0, monthly, 0.05 deg, 50S-50N.
# Countries outside 50S-50N get NA and are handled separately at the end.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. CONFIGURATION  -- edit this block, nothing else should need changing
# -----------------------------------------------------------------------------

# Where to keep the downloaded rasters (a few GB; put this OUTSIDE the git repo
# or add it to .gitignore)
CHIRPS_DIR   <- "data/chirps_monthly"
POP_DIR      <- "data/population"
OUT_DIR      <- "data/derived"

# Years used to define the "normal" seasonal cycle
CLIM_YEARS   <- 1991:2020

# Three different windows, doing three different jobs. With a 12-month IWISE
# recall period these must be kept separate:
#
#   SEASON_WINDOW  climatological season position at interview. Where in the
#                  normal annual cycle the interview fell. Fixed within a
#                  country, varies only by calendar month.
#   RECENT_WINDOW  actual conditions shortly before the interview. This is the
#                  hypothesised source of recall bias.
#   RECALL_WINDOW  actual conditions over the period respondents were asked
#                  about. This is the true exposure and must be controlled for,
#                  otherwise "recent conditions" partly proxies for it.
SEASON_WINDOW <- 3
RECENT_WINDOW <- 3
RECALL_WINDOW <- 12

# Names of the identifying columns in d.iwise_gdp_gov_jmp.
# Adjust to whatever they are actually called in your merged frame.
VAR_ISO3     <- "iso3"        # ISO3 country code
VAR_YEAR     <- "year"        # year of interview
VAR_MONTH    <- "month"       # month of interview, 1-12, as integer

# Minimum share of a country's population that must share the national
# seasonal cycle for the country-level label to be considered usable.
CONCORD_MIN  <- 0.80

# Countries drier than this (mm/yr, population-weighted) get flagged: a
# "wet season" label is not meaningful when there is essentially no rain.
ARID_CUTOFF  <- 250

dir.create(CHIRPS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(POP_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR,    recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 1. PACKAGES
# -----------------------------------------------------------------------------
# install.packages(c("terra","sf","exactextractr","geodata","rnaturalearth",
#                    "rnaturalearthdata","countrycode","tidyverse"))

library(tidyverse)
library(terra)
library(sf)
library(exactextractr)
library(countrycode)


# -----------------------------------------------------------------------------
# 2. WHICH MONTHS DO WE NEED?
# -----------------------------------------------------------------------------
# The climatology years, plus every year present in the IWISE data (and the
# year before, so that a January interview has a full 3-month lookback).

stopifnot(exists("d.iwise_gdp_gov_jmp"))

survey_years <- d.iwise_gdp_gov_jmp[[VAR_YEAR]] |> unique() |> sort()
survey_years <- survey_years[!is.na(survey_years)]
need_years   <- sort(unique(c(CLIM_YEARS, survey_years, survey_years - 1)))

message("CHIRPS years needed: ", min(need_years), "-", max(need_years),
        " (", length(need_years) * 12, " monthly rasters)")


# -----------------------------------------------------------------------------
# 3. DOWNLOAD CHIRPS MONTHLY GEOTIFFS
# -----------------------------------------------------------------------------
# ~2 MB per file gzipped, so roughly 1 GB for 1991-2023. Runs once; the loop
# skips anything already on disk. CHC prefer FTP over HTTPS for scripted
# downloads.
#
# If you would rather not host the rasters yourself, the rOpenSci `chirps`
# package (chirps::get_chirps) pulls values for points or polygons from an API
# instead. It is much slower for a global job like this one, but fine if you
# only end up keeping 20-30 countries.
#
# Note: CHIRPS v3.0 now exists at
#   ftp://ftp.chc.ucsb.edu/pub/org/chc/products/CHIRPS/v3.0/monthly/global/tifs/
# v2.0 is used here because it is what the published literature uses.

chirps_url <- function(y, m) {
  sprintf(paste0("ftp://ftp.chc.ucsb.edu/pub/org/chc/products/CHIRPS-2.0/",
                 "global_monthly/tifs/chirps-v2.0.%d.%02d.tif.gz"), y, m)
}

for (y in need_years) {
  for (m in 1:12) {
    dest <- file.path(CHIRPS_DIR, sprintf("chirps-v2.0.%d.%02d.tif.gz", y, m))
    if (file.exists(dest) && file.size(dest) > 1e5) next
    try(download.file(chirps_url(y, m), dest, mode = "wb", quiet = TRUE),
        silent = TRUE)
  }
  message("downloaded ", y)
}

# terra reads gzipped tifs directly through GDAL's virtual filesystem,
# so there is no need to decompress.
chirps_path <- function(y, m) {
  file.path("/vsigzip", normalizePath(CHIRPS_DIR, mustWork = TRUE),
            sprintf("chirps-v2.0.%d.%02d.tif.gz", y, m))
}


# -----------------------------------------------------------------------------
# 4. COUNTRY POLYGONS AND POPULATION WEIGHTS
# -----------------------------------------------------------------------------

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
  transmute(
    iso3 = countrycode(admin, "country.name", "iso3c", warn = FALSE),
    name = admin
  ) |>
  filter(!is.na(iso3)) |>
  st_make_valid()

# keep only countries that appear in the IWISE data
iwise_iso3 <- unique(d.iwise_gdp_gov_jmp[[VAR_ISO3]])
world      <- world |> filter(iso3 %in% iwise_iso3)

missing_geom <- setdiff(iwise_iso3, world$iso3)
if (length(missing_geom))
  warning("No polygon matched for: ", paste(missing_geom, collapse = ", "))

# Population raster (GPW, ~4.6 km), resampled onto the CHIRPS grid so that
# exact_extract can use it as a weight layer.
template <- rast(chirps_path(CLIM_YEARS[1], 1))
NAflag(template) <- -9999

pop_file <- file.path(POP_DIR, "pop_on_chirps_grid.tif")
if (!file.exists(pop_file)) {
  pop_raw <- geodata::population(year = 2020, res = 2.5, path = POP_DIR)
  pop     <- resample(pop_raw, template, method = "bilinear")
  pop[is.na(pop)] <- 0
  writeRaster(pop, pop_file, overwrite = TRUE)
}
pop <- rast(pop_file)


# -----------------------------------------------------------------------------
# 5. MONTHLY CLIMATOLOGY, POPULATION-WEIGHTED, PER COUNTRY
# -----------------------------------------------------------------------------
# 12 layers: mean rainfall in each calendar month over CLIM_YEARS.

clim_file <- file.path(OUT_DIR, "chirps_climatology_1991_2020.tif")

if (!file.exists(clim_file)) {
  clim_layers <- vector("list", 12)
  for (m in 1:12) {
    stk <- rast(vapply(CLIM_YEARS, chirps_path, character(1), m = m))
    NAflag(stk) <- -9999
    clim_layers[[m]] <- app(stk, mean, na.rm = TRUE)
    message("climatology month ", m)
  }
  clim <- rast(clim_layers)
  names(clim) <- month.abb
  writeRaster(clim, clim_file, overwrite = TRUE)
}
clim <- rast(clim_file)
names(clim) <- month.abb

# Population-weighted country means, one row per country x month
clim_country <- exact_extract(clim, world, "weighted_mean", weights = pop,
                              progress = TRUE) |>
  as_tibble() |>
  setNames(month.abb) |>
  mutate(iso3 = world$iso3) |>
  pivot_longer(all_of(month.abb), names_to = "month_abb",
               values_to = "rain_clim_mm") |>
  mutate(month = match(month_abb, month.abb)) |>
  select(iso3, month, rain_clim_mm) |>
  arrange(iso3, month)


# -----------------------------------------------------------------------------
# 6. HOMOGENEITY DIAGNOSTIC
#    Is a single country-level season label defensible for this country?
# -----------------------------------------------------------------------------
# For every CHIRPS pixel we find the peak of the WINDOW-month running total
# (circular across the year). We then ask what share of the country's
# population lives in a pixel whose peak is within +/- 1 month of the national
# peak. High share -> one national season label describes nearly everyone.
# Low share -> the country contains genuinely different seasonal calendars and
# a country-level label is misleading, not just noisy.

circ_diff <- function(a, b) pmin((a - b) %% 12, (b - a) %% 12)

roll_circ <- function(v, k) {                       # circular running sum
  idx <- outer(seq_len(12), seq_len(k) - 1, function(i, j) ((i - j - 1) %% 12) + 1)
  rowSums(matrix(v[idx], nrow = 12))
}

homog <- map_dfr(seq_len(nrow(world)), function(i) {

  poly <- world[i, ]
  cl   <- crop(clim, poly) |> mask(vect(poly))
  pw   <- crop(pop,  poly) |> mask(vect(poly))

  vals <- as.data.frame(cl, na.rm = FALSE) |> as.matrix()
  w    <- as.data.frame(pw, na.rm = FALSE)[, 1]

  ok <- stats::complete.cases(vals) & !is.na(w) & w > 0
  if (sum(ok) < 5) {
    return(tibble(iso3 = poly$iso3, n_px = sum(ok), concordance = NA_real_,
                  agree_wetdry = NA_real_, n_wet_runs = NA_integer_))
  }
  vals <- vals[ok, , drop = FALSE]
  w    <- w[ok]

  # pixel-level peak of the running SEASON_WINDOW-month total
  roll_px   <- t(apply(vals, 1, roll_circ, k = SEASON_WINDOW))
  peak_px   <- max.col(roll_px, ties.method = "first")

  # national peak, population-weighted
  nat_clim  <- colSums(vals * w) / sum(w)
  peak_nat  <- which.max(roll_circ(nat_clim, SEASON_WINDOW))

  concordance <- sum(w[circ_diff(peak_px, peak_nat) <= 1]) / sum(w)

  # agreement on which months are "wetter than average" (a stricter test)
  wet_px  <- vals > (rowSums(vals) / 12)
  share_m <- colSums(wet_px * w) / sum(w)
  agree   <- mean(pmax(share_m, 1 - share_m))

  # number of separate wet spells nationally -> unimodal vs bimodal
  wet_nat <- nat_clim > mean(nat_clim)
  runs    <- sum(wet_nat & !wet_nat[c(12, 1:11)])

  tibble(iso3 = poly$iso3, n_px = sum(ok), concordance = concordance,
         agree_wetdry = agree, n_wet_runs = as.integer(runs))
}, .progress = TRUE)


# -----------------------------------------------------------------------------
# 7. SEASON CLASSIFICATION AT COUNTRY LEVEL
# -----------------------------------------------------------------------------

season_country <- clim_country |>
  group_by(iso3) |>
  arrange(month, .by_group = TRUE) |>
  mutate(
    rain_annual_mm   = sum(rain_clim_mm),

    # running SEASON_WINDOW-month climatological total ending at this month
    rain_clim_win_mm = roll_circ(rain_clim_mm, SEASON_WINDOW),

    # PRIMARY CONTINUOUS VARIABLE ------------------------------------------
    # share of the country's annual rainfall falling in the window ending
    # at this month. Scale-free, comparable across climates, no threshold.
    season_share     = rain_clim_win_mm / rain_annual_mm,

    # rank version, 1 = driest window of the year, 12 = wettest
    season_rank      = rank(rain_clim_win_mm, ties.method = "first"),

    # CATEGORICAL COMPANION ------------------------------------------------
    # terciles of the 12 windows: 4 months each. Use for tables and plots.
    season_cat = cut(season_rank, breaks = c(0, 4, 8, 12),
                     labels = c("dry", "intermediate", "wet")),

    # ROBUSTNESS DEFINITION ------------------------------------------------
    # classic "rainy season" = the fewest consecutive-ranked months that
    # together carry 70% of annual rainfall
    wet_70 = {
      o   <- order(rain_clim_mm, decreasing = TRUE)
      cum <- cumsum(rain_clim_mm[o]) / sum(rain_clim_mm)
      keep <- o[seq_len(which(cum >= 0.70)[1])]
      month %in% month[keep]
    },

    # simple above-average-month binary
    wet_binary = rain_clim_mm > (rain_annual_mm / 12)
  ) |>
  ungroup() |>
  left_join(homog, by = "iso3") |>
  mutate(
    arid_flag     = rain_annual_mm < ARID_CUTOFF,
    bimodal_flag  = n_wet_runs >= 2,
    # the inclusion rule: use this to subset, or as an interaction term
    season_usable = !is.na(concordance) & concordance >= CONCORD_MIN & !arid_flag
  )

write_csv(season_country, file.path(OUT_DIR, "season_country_month.csv"))

# Look at this before deciding anything about country inclusion.
season_country |>
  distinct(iso3, rain_annual_mm, concordance, agree_wetdry, n_wet_runs,
           arid_flag, bimodal_flag, season_usable) |>
  arrange(concordance) |>
  print(n = 100)


# -----------------------------------------------------------------------------
# 8. ACTUAL CONDITIONS AT INTERVIEW (not climatology)
# -----------------------------------------------------------------------------
# What the climatology cannot tell you: whether that particular window was
# unusually wet or dry. This is the variable that speaks most directly to
# "conditions right now" and it aggregates over space far better than a
# categorical label, because anomalies are already locally normalised.

# Country x year x month actual rainfall, population-weighted
actual_file <- file.path(OUT_DIR, "chirps_country_monthly.csv")

if (!file.exists(actual_file)) {
  actual <- map_dfr(need_years, function(y) {
    stk <- rast(vapply(1:12, function(m) chirps_path(y, m), character(1)))
    NAflag(stk) <- -9999
    ex <- exact_extract(stk, world, "weighted_mean", weights = pop,
                        progress = FALSE)
    tibble(iso3 = rep(world$iso3, each = 12),
           year = y,
           month = rep(1:12, times = nrow(world)),
           rain_mm = as.numeric(t(as.matrix(ex))))
  }, .progress = TRUE)
  write_csv(actual, actual_file)
}
actual <- read_csv(actual_file, show_col_types = FALSE)

# Two accumulations ending at each month, each standardised against the same
# calendar window across CLIM_YEARS:
#
#   recent  (RECENT_WINDOW months)  -> the hypothesised bias driver
#   recall  (RECALL_WINDOW months)  -> the actual reference period respondents
#                                      were asked about; the true exposure
#
# Note these are normalised anomalies (z-scores), not gamma-fitted SPI. For a
# proper SPI use SPEI::spi() on the country monthly series; the difference
# matters mainly in very dry, zero-inflated months.

accumulate <- function(df, k, suffix) {
  out <- df |>
    arrange(iso3, year, month) |>
    group_by(iso3) |>
    mutate(win_mm = slider::slide_dbl(rain_mm, sum, .before = k - 1,
                                      .complete = TRUE)) |>
    ungroup()

  ref <- out |>
    filter(year %in% CLIM_YEARS) |>
    group_by(iso3, month) |>
    summarise(win_mean = mean(win_mm, na.rm = TRUE),
              win_sd   = sd(win_mm,   na.rm = TRUE), .groups = "drop")

  out |>
    left_join(ref, by = c("iso3", "month")) |>
    mutate(z   = (win_mm - win_mean) / win_sd,
           pct = 100 * win_mm / win_mean) |>
    select(iso3, year, month, win_mm, z, pct) |>
    rename_with(~ paste0(c("rain_mm", "rain_z", "rain_pct_normal"), "_", suffix),
                c("win_mm", "z", "pct"))
}

anom <- accumulate(actual, RECENT_WINDOW, "recent") |>
  left_join(accumulate(actual, RECALL_WINDOW, "recall"),
            by = c("iso3", "year", "month"))


# -----------------------------------------------------------------------------
# 9. MERGE
# -----------------------------------------------------------------------------

season_join <- season_country |>
  select(iso3, month,
         rain_clim_mm, rain_clim_win_mm, rain_annual_mm,
         season_share, season_rank, season_cat, wet_70, wet_binary,
         concordance, agree_wetdry, n_wet_runs,
         arid_flag, bimodal_flag, season_usable)

d.iwise_analysis_complete <- d.iwise_gdp_gov_jmp |>
  mutate(
    .iso3  = .data[[VAR_ISO3]],
    .year  = as.integer(.data[[VAR_YEAR]]),
    .month = as.integer(.data[[VAR_MONTH]])
  ) |>
  left_join(season_join, by = c(".iso3" = "iso3", ".month" = "month")) |>
  left_join(anom,        by = c(".iso3" = "iso3", ".year" = "year",
                                ".month" = "month")) |>
  select(-.iso3, -.year, -.month)


# -----------------------------------------------------------------------------
# 10. CHECKS -- run these, do not skip
# -----------------------------------------------------------------------------

# (a) Does interview month actually vary within country-year? If most
#     country-years have n_months == 1, season is collinear with country-year
#     and cannot be separated from a country-year fixed effect.
d.iwise_analysis_complete |>
  count(.data[[VAR_ISO3]], .data[[VAR_YEAR]], .data[[VAR_MONTH]]) |>
  group_by(.data[[VAR_ISO3]], .data[[VAR_YEAR]]) |>
  summarise(n_months = n(), .groups = "drop") |>
  count(n_months)

# (b) How much of the sample survives the homogeneity rule?
d.iwise_analysis_complete |>
  count(season_usable) |>
  mutate(pct = round(100 * n / sum(n), 1))

# (c) Who failed to merge, and why
d.iwise_analysis_complete |>
  filter(is.na(season_share)) |>
  count(.data[[VAR_ISO3]], .data[[VAR_MONTH]]) |>
  print(n = 50)
# Expect: countries wholly north of 50N (no CHIRPS coverage), countries with
# no polygon match, and rows with a missing interview month.

# (d) Sanity check the season variable against something you know.
#     Bangladesh should peak Jun-Aug, Ghana should show two peaks, Chile
#     should peak Jun-Aug (southern hemisphere winter rainfall).
season_country |>
  filter(iso3 %in% c("BGD", "GHA", "CHL", "KEN")) |>
  select(iso3, month, rain_clim_mm, season_cat) |>
  pivot_wider(names_from = month, values_from = c(rain_clim_mm, season_cat))

saveRDS(d.iwise_analysis_complete,
        file.path(OUT_DIR, "d.iwise_analysis_complete.rds"))


# -----------------------------------------------------------------------------
# WHAT YOU NOW HAVE
# -----------------------------------------------------------------------------
# season_share            continuous, 0-1. Share of annual rainfall in the
#                         SEASON_WINDOW ending at the interview month.
#                         Primary season variable.
# season_rank             1-12 ordinal version.
# season_cat              dry / intermediate / wet. For descriptives.
# wet_70, wet_binary      alternative binary definitions, for robustness.
# rain_z_recent           standardised anomaly over RECENT_WINDOW months
#                         before the interview. The hypothesised bias driver.
# rain_z_recall           standardised anomaly over the full RECALL_WINDOW
#                         (12 months) the respondent was asked about. This is
#                         the true exposure -- a control, not a predictor of
#                         interest.
# rain_pct_normal_recent / _recall   same quantities in interpretable units.
# concordance             0-1. Share of population whose local seasonal cycle
#                         matches the national one. Your inclusion criterion.
# bimodal_flag            TRUE where the country has two rainy seasons.
# arid_flag               TRUE where annual rainfall is too low for "season"
#                         to mean much.
# season_usable           concordance >= CONCORD_MIN and not arid.
#
# THE IDENTIFYING LOGIC, given a 12-month recall period
# -----------------------------------------------------
# Every respondent reports on one full annual cycle, so the true amount of
# seasonal variation captured by the outcome is roughly constant regardless of
# when in the year the interview happened. Season at interview should
# therefore have no effect on the true quantity. If it predicts IWISE anyway,
# that is a measurement effect rather than a real one.
#
# The one thing that breaks this: a January and a July interview cover
# different (overlapping) 12-month periods, so genuine differences in realised
# conditions could masquerade as recall bias. rain_z_recall is what closes
# that gap.
#
# Suggested specification:
#   IWISE ~ season_share + rain_z_recent + rain_z_recall + <covariates>
#           + country-year fixed effects
# The coefficient of interest is on season_share (or rain_z_recent) after
# rain_z_recall is held constant. Restrict the main analysis to
# season_usable == TRUE, use the full sample as a sensitivity check, and
# report concordance in the methods -- it is the honest answer to "you only
# had country-level location".