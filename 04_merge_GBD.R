library(dplyr)        # data manipulation
library(tidyr)        # pivot_wider
library(readr)        # read_csv
library(countrycode)  # country names -> ISO3

# 1. Load ---------------------------------------------------------------
d.iwise_gdp_gov_jmp <- readRDS("data/iwise_gdp_gov_jmp.rds")   # survey data
gbd_raw <- read_csv("Q:/Abteilungsprojekte/Sandec/7_PHIC/02-Projects/P07-GLO-IWISE/P07-02-Input-data/Data_GBD/GBD_enteric_data_29092026.csv")                                  # GBD: 17 age bands, Number + Rate

# Check A: exact labels before filtering
gbd_raw |> distinct(population_group_name, measure_name, sex_name,
                    cause_name, metric_name)                    # confirm filter strings
gbd_raw |> distinct(age_name) |> print(n = Inf)                 # expect the 17 bands (15-19 ... 95+)
range(gbd_raw$year)                                             # expect max 2023

# 2. Age lookup ---------------------------------------------------------
age_map <- tibble(
  age_name = c("15-19 years", "20-24 years",                    # -> 0: 15-24
               "25-29 years", "30-34 years",                    # -> 1: 25-34
               "35-39 years", "40-44 years", "45-49 years",     # -> 2: 35-49
               "50-54 years", "55-59 years", "60-64 years",     # -> 3: 50+
               "65-69 years", "70-74 years", "75-79 years",
               "80-84 years", "85-89 years", "90-94 years",
               "95+ years"),
  age_gp_profile = c(0, 0, 1, 1, 2, 2, 2, rep(3, 10))          # group code per band
)
bands_expected <- c(2, 2, 3, 10)                                # bands per group (0-3)

# 3. Filter -------------------------------------------------------------
gbd_long <- gbd_raw |>
  filter(cause_name   == "Enteric infections",                  # enteric infections
         measure_name == "DALYs (Disability-Adjusted Life Years)",  # DALYs only
         sex_name     == "Both",                                # both sexes
         metric_name  %in% c("Number", "Rate")) |>              # counts and rates per 100k
  inner_join(age_map, by = "age_name")                          # keep 17 bands, add group code

nrow(gbd_long) > 0                                              # FALSE = a label didn't match

# Check B: no duplicate rows (duplicates would break pivot_wider)
gbd_long |> count(location_name, year, age_name, metric_name) |>
  filter(n > 1)                                                 # expect 0 rows

# 4. Reshape and back-calculate population ------------------------------
gbd_bands <- gbd_long |>
  pivot_wider(id_cols = c(location_name, year, age_name, age_gp_profile),  # one row per band (unique combination of these becomes a row)
              names_from = metric_name, values_from = val) |>   # Number and Rate as columns
  mutate(pop = if_else(Rate > 0, Number / Rate * 1e5, NA_real_))  # band population; NA if rate = 0 (Number / rate * 100'000)

gbd_bands |> filter(is.na(Number) | is.na(Rate) | is.na(pop))   # expect 0 rows

# 5. Country codes ------------------------------------------------------
gbd_bands <- gbd_bands |>
  mutate(iso3c = countrycode(location_name, "country.name", "iso3c"))  # match to iso3c

# Check C: unmatched GBD locations (should be regions/global only)
gbd_bands |> filter(is.na(iso3c)) |> distinct(location_name)
gbd_bands |> filter(!is.na(iso3c)) |>
  summarise(same = n_distinct(location_name) == n_distinct(iso3c))    # expect TRUE

# 6. Aggregate to age groups ---------------------------------------
gbd <- gbd_bands |>
  filter(!is.na(iso3c)) |>                                      # drop potential regions
  group_by(iso3c, gbd_year = year, age_gp_profile) |>           # one row per country-year-group
  summarise(daly_diarr_rate = sum(Number) / sum(pop) * 1e5,     # crude rate per 100k
            n_bands = n(),                                      # bands used (for checking)
            .groups = "drop")

# Check D: each group built from the right number of bands
gbd |> filter(n_bands != bands_expected[age_gp_profile + 1])    # expect 0 rows
gbd <- gbd |> select(-n_bands)                                  # drop helper column

# 7. Pre-merge key checks -----------------------------------------------
my_iso <- d.iwise_gdp_gov_jmp |> distinct(iso3c) |> pull()      # your countries
length(my_iso)                                                  # expect 76

setdiff(my_iso, gbd$iso3c)                                      # Check E: expect none missing

# Check F: all 4 age groups present for every country-year you need
needed <- d.iwise_gdp_gov_jmp |>
  distinct(iso3c, year) |>                                        # survey country-years
  mutate(gbd_year = pmin(year, max(gbd$gbd_year))) |>             # 2024/25 -> 2023
  distinct(iso3c, gbd_year) |>                                    # GBD years needed
  cross_join(tibble(age_gp_profile = 0:3))                        # x 4 age groups

needed |>
  anti_join(gbd, by = c("iso3c", "gbd_year", "age_gp_profile"))   # expect 0 rows

# 8. Match year and merge -----------------------------------------------
max_gbd <- max(gbd$gbd_year)                                    # latest GBD year (2023)

d.iwise_all_contextual <- d.iwise_gdp_gov_jmp |>
  mutate(gbd_year_used   = pmin(year, max_gbd),                 # use 2023 for 2024/25 waves, pmin takes smallest of the two values
         gbd_carried_fwd = year > max_gbd,                      # flag carried-forward rows
         age_key = as.numeric(age_gp_profile)) |>               # strip Stata labels for join
  left_join(gbd, by = c("iso3c",
                        "gbd_year_used" = "gbd_year",
                        "age_key" = "age_gp_profile"),          # match country, year, age group
            relationship = "many-to-one")                       # many people -> one rate

# Check G: row count unchanged
stopifnot(nrow(d.iwise_all_contextual) == nrow(d.iwise_gdp_gov_jmp))

# Check H: missing rates should only be people with missing age
d.iwise_all_contextual |>
  filter(is.na(daly_diarr_rate)) |>
  count(iso3c, year, age_missing = is.na(age_key))              # expect age_missing = TRUE only

# Check I: plausibility by age group (rates usually highest in 50+)
gbd |> filter(iso3c %in% my_iso) |>
  group_by(age_gp_profile) |>
  summarise(min = min(daly_diarr_rate),
            median = median(daly_diarr_rate),
            max = max(daly_diarr_rate))

# Check J: carry-forward and spot-check
d.iwise_all_contextual |> distinct(year, gbd_year_used, gbd_carried_fwd)
d.iwise_all_contextual |>
  filter(iso3c %in% c("PSE", "VEN", "LBN")) |>
  distinct(iso3c, year, age_key, gbd_year_used, daly_diarr_rate) |>
  arrange(iso3c, year, age_key)

# 9. Save ---------------------------------------------------------------
saveRDS(d.iwise_all_contextual, "data/iwise_all_contextual.rds")
file.exists("data/iwise_all_contextual.rds")
