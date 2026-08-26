# 01 Packages and loading data ----

install.packages(c("WDI", "countrycode", "dplyr", "readr"))

library(WDI)          # World Bank API client
library(countrycode)  # ISO3 <-> country name utilities
library(dplyr)
library(readr)

source("00_load.R")

#Create a copy of the dataset
mydata <- d.iwise |>
  rename(iso3c = COUNTRY_ISO3,
         year  = YEAR_CALENDAR)


# 02  - Validate the merge keys ---- 
#    Only iso3c and year are used to join. The country name is never used for
#    matching -- but it is useful for spotting rows where the name and the code
#    disagree, which means one of them is wrong.
#
#    NOTE: this is individual-level data (~91k rows, 76 countries), so every
#    diagnostic below ends in distinct() or count() to collapse to country
#    level. If any output says "omitted N rows", it printed at row level and
#    you are reading a truncated slice -- add distinct() before trusting it.
# -----------------------------------------------------------------------------
 
# Strip stray whitespace and force uppercase -- both silently break joins
mydata <- mydata |>
  mutate(iso3c = toupper(trimws(iso3c)))
 
## --- CHECK A: is every ISO3 code well-formed and recognised? -----------------

# countrycode() returns NA for codes it cannot resolve. Anything listed here is
# malformed (wrong length, typo) or is a valid World Bank code that ISO does
# not recognise -- XKX (Kosovo), CHI (Channel Islands). The latter are fine and
# will still merge; do not "fix" them.
mydata |>
  mutate(check = countrycode(iso3c, "iso3c", "country.name",
                             custom_match = c("XKX" = "Kosovo"),
                             warn = FALSE)) |>
  filter(is.na(check) | nchar(iso3c) != 3) |>
  distinct(iso3c, countrynew)

## --- CHECK B: do your codes and your names agree? ----------------------------

# Derives the name from the code and compares it to your own column. Names
# differ in spelling for legitimate reasons ("Czechia" vs "Czech Republic"),
# so read the output rather than treating every row as an error -- you are
# looking for rows where the two are genuinely DIFFERENT COUNTRIES, which
# means a mis-keyed code.
mydata |>
  mutate(name_from_code = countrycode(iso3c, "iso3c", "country.name",
                                      custom_match = c("XKX" = "Kosovo"),
                                      warn = FALSE)) |>
  filter(is.na(name_from_code) |
           tolower(countrynew) != tolower(name_from_code)) |>
  distinct(iso3c, countrynew, name_from_code) #all good
 

## --- CHECK C: watch for these specific codes ---------------------------------

# COD = Congo, Dem. Rep. (Kinshasa) | COG = Congo, Rep. (Brazzaville)
#   -- if you surveyed both, BOTH must appear below, on separate lines
# ROU = Romania (not ROM) | TLS = Timor-Leste  -- older datasets use dead codes
# PSE = West Bank and Gaza | TWN = Taiwan, NOT in the WDI at any code
mydata |>
  filter(iso3c %in% c("COD","COG","ROM","ROU","TLS","TMP","TWN","XKX","PSE")) |>
  distinct(iso3c, countrynew, year) #all fine
 
## --- CHECK D: does every row have exactly one year in range? -----------------
range(mydata$year, na.rm = TRUE)   # 2020 - 2025
sum(is.na(mydata$year))            # 0
 
### Drop rows with no country code ------------------------------------------
# Trailing blank rows are common from CSV/Excel imports. Check the count that
# disappears -- it should match what Check A flagged, and nothing more.
n_before <- nrow(mydata)
mydata   <- mydata |> filter(!is.na(iso3c), nchar(iso3c) == 3)
n_before - nrow(mydata)            # rows removed
 
## --- CHECK E: exactly ONE collection year per country? -----------------------
# Data is individual-level: many respondents per country, but by design
# each country was surveyed in a single year. This should return ZERO rows.
# Anything listed here has respondents split across two or more years, which
# breaks the one-year-per-country assumption.
mydata |>
  group_by(iso3c) |>
  summarise(n_years = n_distinct(year), years = paste(sort(unique(year)),
                                                      collapse = ", ")) |>
  filter(n_years > 1) #only PSE with 2022 and 2025
 
# --- CHECK F: expected number of countries, and respondents per country? -----
n_distinct(mydata$iso3c)           # 76
 


# --- 4. DOWNLOAD WORLD BANK DATA ------ 
#    start = 2019 gives slack for the carry-forward rule in section 6.
# -----------------------------------------------------------------------------
gdp_raw <- WDI(
  country   = "all",
  indicator = c(gdp_pc_ppp = "NY.GDP.PCAP.PP.KD"),
  start     = 2019,
  end       = 2025,
  extra     = TRUE          # adds iso3c, region, income group
)
 
# Save the raw pull IMMEDIATELY. The World Bank revises data between editions,
# so this file plus its date is what makes your results reproducible.
saveRDS(gdp_raw, paste0("data/wdi_raw_", Sys.Date(), ".rds")) #downloaded on 26.08.2026
 

# --- 5. CLEAN ----
#    The default download contains ~50 non-countries (World, Euro area,
#    "Low income", regions). Real countries have a non-Aggregates region.
# -----------------------------------------------------------------------------
gdp <- gdp_raw |>
  filter(region != "Aggregates") |>
  select(iso3c, year, gdp_pc_ppp) |>
  mutate(year = as.integer(year)) |>
  filter(!is.na(gdp_pc_ppp))

## CHECK G: coverage by year. The PPP series lags, so recent years are ------
# thinner. 
gdp |> count(year, name = "n_countries") |> arrange(year) #all good

# --- 6. REVIEW ----
# countries flagged by Check E as having more than one wave --------
# Two waves are mechanically fine -- each is matched to its own year's GDP.
# But inspect the GDP series for these countries before trusting the result,
# because a country surveyed in two distant years may have had a large real
# shock in between (war, currency collapse, recovery from one).
# Edit the vector to match whatever Check E returned.
multi_wave <- c("PSE")
 
gdp |>
  filter(iso3c %in% multi_wave) |>
  arrange(iso3c, year)
 
gdp_raw |>
  filter(iso3c %in% multi_wave) |>
  select(iso3c, year, gdp_pc_ppp) |>
  arrange(iso3c, year)
 
# use the corresponding waves


# --- 7. MERGE ----
#    Rule: match on collection year. If that year is unpublished for a country,
#    fall back to its most recent EARLIER year, and flag it.
#    The rule is fixed in advance and applied uniformly --  not decided
#    after seeing which countries drop out.
 
# 7a. Exact match on country + year
merged <- left_join(mydata, gdp, by = c("iso3c", "year"))
 
# 7b. Build the fallback for rows that found nothing
fallback <- merged |>
  filter(is.na(gdp_pc_ppp)) |>
  distinct(iso3c, year_collect = year) |>
  left_join(gdp |> rename(gdp_year = year, gdp_cf = gdp_pc_ppp),
            by = "iso3c", relationship = "many-to-many") |>
  filter(gdp_year <= year_collect) |>          # only ever look backwards
  group_by(iso3c, year_collect) |>
  slice_max(gdp_year, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(iso3c, year = year_collect, gdp_cf, gdp_year)
 
# 7c. Apply it, keeping a record of what was substituted
merged <- merged |>
  left_join(fallback, by = c("iso3c", "year")) |>
  mutate(
    gdp_carried_fwd = is.na(gdp_pc_ppp) & !is.na(gdp_cf),
    gdp_year_used   = if_else(is.na(gdp_pc_ppp), gdp_year, year),
    gdp_pc_ppp      = coalesce(gdp_pc_ppp, gdp_cf)
  ) |>
  select(-gdp_cf, -gdp_year)
 
 
# 7d. Document HOW each row got its value, so an NA is self-explaining later.
#     This is a base dataset -- no rows are dropped. Downstream models will
#     drop the NAs listwise; anything not using GDP keeps the full sample.
merged <- merged |>
  mutate(gdp_match_status = case_when(
    is.na(gdp_pc_ppp)          ~ "no data published",
    gdp_carried_fwd            ~ "carried forward",
    TRUE                       ~ "exact year match"
  ))
 
# Expected: 75 countries "exact year match" (PSE twice, once per wave),
# Lebanon "carried forward" from 2024, Venezuela "no data published".
merged |> distinct(iso3c, countrynew, year, gdp_match_status) |>
  count(gdp_match_status) #correct
 
merged |> distinct(iso3c, countrynew, year, gdp_year_used, gdp_match_status) |>
  filter(gdp_match_status != "exact year match") #correct


# --- 8. DIAGNOSTICS ----

## --- CHECK H: row count unchanged. Must be TRUE. -----------------------------
nrow(merged) == nrow(mydata)
 
##--- CHECK I: what is still missing, and why? --------------------------------
# Read the output three ways:
#   - a country you expected to be there  -> wrong/missing ISO3 code
#   - Taiwan or North Korea               -> not in the WDI at all, ever
#   - Venezuela, Syria, Eritrea, Somalia  -> genuine gaps in the source data
merged |>
  filter(is.na(gdp_pc_ppp)) |>
  count(countrynew, iso3c, year, name = "n_respondents_lost") #correct, only venezuela
 
## --- CHECK J: how much was carried forward, and by how many years? -----------
merged |>
  filter(gdp_carried_fwd) |>
  mutate(lag_years = year - gdp_year_used) |>
  distinct(countrynew, year, gdp_year_used, lag_years) |>
  arrange(desc(lag_years))
 
## --- CHECK K: eyeball the merge at COUNTRY level, not row level --------------
# each country shows the year you expect from your fieldwork records.
merged |>
  distinct(countrynew, iso3c, year, gdp_year_used, gdp_pc_ppp) |>
  arrange(year, countrynew) |>
  print(n = Inf) #correct, only Lebanon
 
# Each country should appear ONCE -- except PSE, which has two waves (2022 and
# 2025) and so legitimately appears twice with different GDP values. Any OTHER
# country appearing twice means the join attached two values within one
# country-year, which is an error.
merged |> distinct(countrynew, iso3c, year, gdp_pc_ppp) |>
  count(iso3c) |> filter(n > 1) #correct
 
## --- CHECK L: sanity-check the values against known figures ------------------
# e.g. a high-income country should land around 50,000-70,000; a low-income
# one around 1,000-3,000. Wildly different means the wrong series or a units
# problem.
summary(merged$gdp_pc_ppp)


# --- RENAME and save -----
d.iwise_gdp <- merged
rm(merged)
saveRDS(d.iwise_gdp, "data/d_iwise_gdp.rds")
