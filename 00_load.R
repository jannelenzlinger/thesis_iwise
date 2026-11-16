library(tidyverse)

# Load dataset ----
setwd("Q:/Abteilungsprojekte/Sandec/7_PHIC/02-Projects/P07-GLO-IWISE/P07-03-Analysis")
d.iwise <- read.csv("IWISE_all.csv", header = TRUE, sep = ",")
d.iwise <- d.iwise |>
  filter(!is.na(COUNTRY_ISO3), nchar(COUNTRY_ISO3) == 3) |>
  mutate(COUNTRY_ISO3 = toupper(trimws(COUNTRY_ISO3)))

summary(d.iwise)
str(d.iwise)

# Recode Variables ----

#IWISE
str(d.iwise$iwisescore)
summary(d.iwise$iwisescore)
#Income quintiles
d.iwise$INCOME_5 <- factor(d.iwise$INCOME_5,
                                        levels = c("1", "2", "3", "4", "5"),
                                        ordered = TRUE) #income quintiles as ordered categories
#FLI
d.iwise$INDEX_FL_ord <- factor(round(d.iwise$INDEX_FL, 4),
                               levels  = sort(unique(round(d.iwise$INDEX_FL, 4))),
                               labels  = c("0", "16.67", "25", "33.33", "50",
                                           "66.67", "75", "83.33", "100"),
                               ordered = TRUE) #Financial life index as ordered categories
#Food and Shelter Index
d.iwise$INDEX_FS <- factor(d.iwise$INDEX_FS,
                                levels = c("0", "50", "100"),
                                ordered = TRUE)
#Community Attachment Index
d.iwise$INDEX_CA <- factor(round(d.iwise$INDEX_CA, 4),
                            levels  = sort(unique(round(d.iwise$INDEX_CA, 4))),
                            labels  = c("0", "33.33", "50.00", "66.67", "100"),
                            ordered = TRUE) #Community Attachment Index as ordered categories
#Personal Health Index
d.iwise$INDEX_PH <- factor(round(d.iwise$INDEX_PH, 4),
                            levels  = sort(unique(round(d.iwise$INDEX_PH, 4))),
                            labels  = c("0", "20.00", "25.00", "33.33", "40.00", "50.00", "60.00", "66.67", "75", "80", "100"),
                            ordered = TRUE) #Personal Health Index as ordered categories

#Youth Development Index
d.iwise$INDEX_YD <- factor(round(d.iwise$INDEX_YD, 4),
                            levels  = sort(unique(round(d.iwise$INDEX_YD, 4))),
                            labels  = c("0", "33.33", "66.67", "100"),
                            ordered = TRUE) #Personal Health Index as ordered categories


#check for iwise items
library(dplyr)

items <- c("iwise_worry", "iwise_interrupt", "iwise_clothes", "iwise_plans", "iwise_food", "iwise_hands",
           "iwise_body", "iwise_drink", "iwise_angry", "iwise_sleep", "iwise_none", "iwise_shame")

setdiff(items, names(d.iwise))   # should return character(0)

check <- d.iwise %>%
  mutate(
    across(all_of(items), ~ ifelse(.x %in% 0:3, .x, NA)),  # keep only valid 0-3
    n_answered = rowSums(!is.na(across(all_of(items)))),
    has_score  = !is.na(iwisescore)                        # replace with your score name
  )

table(check$n_answered, check$has_score, useNA = "ifany")

check %>%
  filter(n_answered == 12, has_score) %>%
  summarise(mismatch = sum(rowSums(across(all_of(items))) != iwisescore))


library(dplyr)

# 1. Overall: response codes incl. NA
table(d.iwise$WP88, useNA = "ifany")

# 2. Total missing (NA = not asked)
sum(is.na(d.iwise$WP88))
mean(is.na(d.iwise$WP88)) * 100          # percent

# 3. Missingness per country
wp88_country <- d.iwise %>%
  group_by(country_year) %>%
  summarise(n = n(),
            n_missing = sum(is.na(WP88)),
            pct_missing = mean(is.na(WP88)) * 100) %>%
  arrange(desc(pct_missing))

wp88_country

# 4. How many countries have WP88 missing completely vs. partly
sum(wp88_country$pct_missing == 100)                                    # completely
sum(wp88_country$pct_missing > 0 & wp88_country$pct_missing < 100)      # partly

# 5. Which countries are completely missing
wp88_country %>% filter(pct_missing == 100) %>% pull(country_year)
wp88_country %>%
  filter(pct_missing > 0 & pct_missing < 100)
d.iwise %>%
  filter(country_year == 11) %>%
  distinct(countrynew, country_year)


# 6. Same check by country-wave
d.iwise %>%
  group_by(country_year) %>%
  summarise(pct_missing = mean(is.na(WP88)) * 100) %>%
  summarise(waves_complete_missing = sum(pct_missing == 100),
            waves_partly_missing   = sum(pct_missing > 0 & pct_missing < 100))