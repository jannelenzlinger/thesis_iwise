library(dplyr)
library(readr)
library(countrycode)

# 1. Load ---------------------------------------------------------------
d.iwise_gdp_gov_jmp <- readRDS("data/iwise_gdp_gov_jmp.rds")
gbd_raw <- read_csv("Q:/Abteilungsprojekte/Sandec/7_PHIC/02-Projects/P07-GLO-IWISE/P07-02-Input-data/GDB_diarrhea_data.csv") 

# Check A: labels in the raw file (confirm exact strings before filtering)
gbd_raw |> distinct(population_group_name, measure_name, sex_name,
                    age_name, cause_name, metric_name)
range(gbd_raw$year)                                     # expect max 2023

# 2. Filter -------------------------------------------------------------
gbd <- gbd_raw |>
  filter(cause_name   == "Diarrheal diseases",
         measure_name == "DALYs (Disability-Adjusted Life Years)",
         sex_name     == "Both",
         age_name     == "Age-standardized",
         metric_name  == "Rate")

nrow(gbd) > 0                                           # FALSE = a label didn't match

# 3. Country codes ------------------------------------------------------
gbd <- gbd |>
  mutate(iso3c = countrycode(location_name, "country.name", "iso3c"))

# Check B: GBD locations that failed to get a code (should be regions/global only)
gbd |> filter(is.na(iso3c)) |> distinct(location_name)
n_distinct(gbd$location_name) == n_distinct(gbd$iso3c)   # expect TRUE
gbd <- gbd |>
  filter(!is.na(iso3c)) |>
  select(iso3c, gbd_year = year,
         daly_diarr_rate  = val,
         daly_diarr_lower = lower,
         daly_diarr_upper = upper)

# Check C: one row per country-year (duplicates = subnational or population_group clash)
gbd |> count(iso3c, gbd_year) |> filter(n > 1)

# 4. Pre-merge key check ------------------------------------------------
my_iso <- d.iwise_gdp_gov_jmp |> distinct(iso3c) |> pull()
length(my_iso)                                          # expect 76

# Check D: your countries missing from GBD
setdiff(my_iso, gbd$iso3c) #zero

# Check E: every country has 2023 (needed for 2024/2025 waves)
setdiff(my_iso, gbd |> filter(gbd_year == 2023) |> pull(iso3c))
gbd |> filter(gbd_year == 2023, is.na(daly_diarr_rate))   # expect 0 rows

# 5. Match year and merge -----------------------------------------------
max_gbd <- max(gbd$gbd_year) #finds the latest year in the GBD data (2023) and stores i

d.iwise_analysis <- d.iwise_gdp_gov_jmp |>
  mutate(gbd_year_used    = pmin(year, max_gbd),        
         gbd_carried_fwd  = year > max_gbd) |>
  left_join(gbd, by = c("iso3c", "gbd_year_used" = "gbd_year"),
            relationship = "many-to-one")

# Check F: row count unchanged
stopifnot(nrow(d.iwise_analysis) == nrow(d.iwise_gdp_gov_jmp))

# Check G: missing values by country-year
d.iwise_analysis |>
  filter(is.na(daly_diarr_rate)) |>
  distinct(iso3c, year)                               

# Check H: plausibility (range; val within uncertainty interval)
summary(gbd$daly_diarr_rate)
gbd |> filter(daly_diarr_lower > daly_diarr_rate |
              daly_diarr_rate  > daly_diarr_upper)      # expect 0 rows

# Check I: carry-forward and spot-check
d.iwise_analysis |> distinct(year, gbd_year_used, gbd_carried_fwd)   
d.iwise_analysis |>
  filter(iso3c %in% c("PSE", "VEN", "LBN")) |>
  distinct(iso3c, year, gbd_year_used, daly_diarr_rate) |>       
  arrange(iso3c, year)

# 6. Save ---------------------------------------------------------------
saveRDS(d.iwise_analysis, "data/iwise_analysis.rds")
file.exists("data/iwise_analysis.rds")


summary(gbd$daly_diarr_rate)
