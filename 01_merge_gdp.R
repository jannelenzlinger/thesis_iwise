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
 
# --- CHECK A: is every ISO3 code well-formed and recognised? -----------------
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
 
# --- CHECK B: do your codes and your names agree? ----------------------------
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
  distinct(iso3c, countrynew, name_from_code)
 
# --- CHECK C: watch for these specific codes ---------------------------------
# COD = Congo, Dem. Rep. (Kinshasa) | COG = Congo, Rep. (Brazzaville)
#   -- if you surveyed both, BOTH must appear below, on separate lines
# ROU = Romania (not ROM) | TLS = Timor-Leste  -- older datasets use dead codes
# PSE = West Bank and Gaza | TWN = Taiwan, NOT in the WDI at any code
mydata |>
  filter(iso3c %in% c("COD","COG","ROM","ROU","TLS","TMP","TWN","XKX","PSE")) |>
  distinct(iso3c, countrynew, year)
 
# --- CHECK D: does every row have exactly one year in range? -----------------
range(mydata$year, na.rm = TRUE)   # expect 2021 - 2025
sum(is.na(mydata$year))            # expect 0
 
# --- Drop rows with no country code ------------------------------------------
# Trailing blank rows are common from CSV/Excel imports. Check the count that
# disappears -- it should match what Check A flagged, and nothing more.
n_before <- nrow(mydata)
mydata   <- mydata |> filter(!is.na(iso3c), nchar(iso3c) == 3)
n_before - nrow(mydata)            # rows removed
 
# --- CHECK E: exactly ONE collection year per country? -----------------------
# Your data is individual-level: many respondents per country, but by design
# each country was surveyed in a single year. This should return ZERO rows.
# Anything listed here has respondents split across two or more years, which
# breaks the one-year-per-country assumption.
mydata |>
  group_by(iso3c) |>
  summarise(n_years = n_distinct(year), years = paste(sort(unique(year)),
                                                      collapse = ", ")) |>
  filter(n_years > 1)
 
# --- CHECK F: expected number of countries, and respondents per country? -----
n_distinct(mydata$iso3c)           # expect 76
 
# Eyeball the sample sizes -- a country with 3 respondents is worth knowing
# about before it silently drives a result.
mydata |> count(iso3c, countrynew, year, sort = TRUE) |> print(n = Inf)