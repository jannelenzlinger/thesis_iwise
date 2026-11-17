# Descriptive analysis

# income quintiles: income_5
# %% Install packages
install.packages("ggplot2")
library(ggplot2)
install.packages(c("jsonlite", "curl", "sysfonts", "showtextdb", "showtext"))
install.packages("showtext")
library(showtext)

font_add_google("Source Serif 4", "source_serif")
showtext_auto()


# %% Set working directory and load data
setwd("Q:/Abteilungsprojekte/Sandec/7_PHIC/02-Projects/P07-GLO-IWISE/P07-03-Analysis")
d.iwise <- read.csv("IWISE_all.csv", header = TRUE, sep = ",")
summary(d.iwise)


#%% IWISE score

str(d.iwise$iwisescore_imp) #iwise
summary(d.iwise$iwisescore)
sum(d.iwise$iwisescore >= 12, na.rm = TRUE)
sum(d.iwise$iwisescore == 0, na.rm=TRUE)
sum(is.na(d.iwise$iwisescore))

ggplot(d.iwise, aes(x = d.iwise$iwisescore)) +
  geom_histogram(binwidth = 1, fill = "#3F4624", color = "white") +
  labs(title = "Distribution of IWISE Score", x = "IWISE Score", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5, size = 28),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 28))

ggsave("iwise_histogram.png", width = 8, height = 5, dpi = 300, scale = 0.7)


#%%
ggplot(d.iwise[d.iwise$iwisescore > 0, ], aes(x = iwisescore)) +
  geom_histogram(binwidth = 1, fill = "#3F4624", color = "white") +
  labs(title = "Distribution of IWISE Score, Excluding Zeros", x = NULL, y = NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 28),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 28))

ggsave("iwise_histogram_no_zeros.png", width = 8, height = 5, dpi = 300, scale = 0.7)


#%% Income Quintiles
str(d.iwise$INCOME_5) #income quintiles
summary(d.iwise$INCOME_5)
#redefine quintiles as ordered factors
d.iwise$INCOME_5 <- factor(d.iwise$INCOME_5, 
                                        levels = c("1", "2", "3", "4", "5"), 
                                        ordered = TRUE) 
#check distribution within country
table(d.iwise$INCOME_5, d.iwise$countrynew)
prop.table(table(d.iwise$INCOME_5, d.iwise$countrynew), margin = 2) * 100

#%%
quintile_pct <- as.data.frame(prop.table(table(d.iwise$countrynew, d.iwise$INCOME_5), margin = 1) * 100)
names(quintile_pct) <- c("country", "quintile", "percent")
p <- ggplot(quintile_pct, aes(x = quintile, y = country, fill = percent)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#3F4624") +
  labs(title = "Income Quintile Distribution by Country", x = "Quintile", y = NULL, fill = "%") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 24))
p

ggsave("Q:/Abteilungsprojekte/Sandec/7_PHIC/02-Projects/P07-GLO-IWISE/P07-03-Analysis/iwise_quintile_heatmap.png", plot = p, width = 10, height = 12, dpi = 300)
message("Save attempted")
list.files(getwd(), pattern = "quintile_heatmap")




#%% Financial Life Index
str(d.iwise$INDEX_FL) #FLI
summary(d.iwise$INDEX_FL)
sort(unique(d.iwise$INDEX_FL))
sum(d.iwise$INDEX_FL == 0, na.rm=TRUE)
sum(d.iwise$INDEX_FL == 100, na.rm=TRUE)
sum(is.na(d.iwise$INDEX_FL))

ggplot(d.iwise[!d.iwise$INDEX_FL %in% c(25, 75) & !is.na(d.iwise$INDEX_FL), ], 
       aes(y = factor(round(INDEX_FL, 1)))) +
  geom_bar(fill = "#3F4624") +
  labs(title = "Distribution of Financial Life Index", 
       x = "Count", y = "Financial Life Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 28),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 28))

ggsave("distribution_FLI.png", width = 8, height = 5, dpi = 300, scale = 0.7)

#%% Country overview

#country income group
str(d.iwise$country_income_group)
summary(d.iwise$country_income_group)
sort(unique(d.iwise$country_income_group))
sum(is.na(d.iwise$country_income_level))


lvls <- c("", "Low income", "Lower middle income",
          "Upper middle income", "High income")

d.iwise$country_income_group <- factor(d.iwise$country_income_group,
                                       levels = lvls,
                                       labels = c("Not classified", "Low income",
                                                  "Lower middle income",
                                                  "Upper middle income", "High income"),
                                       ordered = TRUE)

table(d.iwise$country_income_group, useNA = "ifany")

str(d.iwise$country_income_group)
sum(is.na(d.iwise$country_income_level))

country_lvl <- unique(d.iwise[, c("countrynew", "country_income_group")])
country_lvl <- country_lvl[country_lvl$country_income_group != "Not classified" |
                           is.na(country_lvl$country_income_group), ]
country_lvl <- droplevels(country_lvl)

levels(country_lvl$country_income_group)

ggplot(country_lvl, aes(x = country_income_group)) +
  geom_bar(fill = "#3F4624", color = "white") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 3.5) +
  labs(title = "Number of Countries by Income Group",
       x = "Income Group", y = "Number of Countries") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 28))

ggsave("distribution_countries_income.png", width = 8, height = 5, dpi = 300, scale = 0.7)



class(country_lvl$country_income_group)
levels(country_lvl$country_income_group)





iwise_items <- d.iwise[, paste0("WP", ...)]   # your 12 item columns
n_miss <- rowSums(is.na(iwise_items))
table(n_miss, is.na(d.iwise$iwise_score))




# ============================================================================
# Exclusion flow for the sample description (Results 3.1)
# Counts each exclusion step in order: a respondent is counted only at the
# FIRST step that removes them, so the steps add up to the final N.
# Needs: d.iwise (all respondents), v.model.fl, v.model.inc (from 05 script)
# ============================================================================

library(dplyr)

# Variable names in the model lists must exist in d.iwise; if this stops,
# swap in the uncentred names (e.g. hhsize instead of hhsize_c)
stopifnot(all(v.model.fl  %in% names(d.iwise)),
          all(v.model.inc %in% names(d.iwise)))

# ---- 1. Full dataset --------------------------------------------------------
d.iwise %>% summarise(n = n(),
                      countries = n_distinct(country_name),
                      country_years = n_distinct(country_year))

# Check that country_year 19 is Palestine 2025 (from earlier output)
d.iwise %>% filter(country_year == 19) %>% count(country_name)

# ---- 2. Flow function -------------------------------------------------------
flow <- function(data, steps) {
  keep <- rep(TRUE, nrow(data))
  out  <- tibble(step = "Full dataset", n_excluded = NA_integer_,
                 n_remaining = nrow(data),
                 countries   = n_distinct(data$country_name))
  for (s in names(steps)) {
    drop <- keep & steps[[s]]                 # newly excluded at this step
    keep <- keep & !steps[[s]]
    out  <- add_row(out, step = s, n_excluded = sum(drop),
                    n_remaining = sum(keep),
                    countries   = n_distinct(data$country_name[keep]))
  }
  list(table = out, keep = keep)
}

# ---- 3. Shared exclusions ---------------------------------------------------
shared <- with(d.iwise, list(
  "Palestine 2025 (no weights)" = country_year %in% 19,
  "Venezuela (no GDP)"          = country_name %in% "Venezuela, RB",
  "Argentina (outdated JMP)"    = country_name %in% "Argentina",
  "Missing IWISE item(s)"       = is.na(iwisescore)
))

# ---- 4. FLI sample ----------------------------------------------------------
fl <- flow(d.iwise, c(shared, with(d.iwise, list(
  "Egypt (FLI not fielded)"          = country_name %in% "Egypt, Arab Rep.",
  "Morocco (FLI not fielded)"        = country_name %in% "Morocco",
  "Missing FLI"                      = is.na(FLI_3item)
)), list(
  "Missing covariates" = !complete.cases(select(d.iwise, all_of(v.model.fl)))
)))
fl$table

# ---- 5. Income sample -------------------------------------------------------
inc <- flow(d.iwise, c(shared, with(d.iwise, list(
  "Zimbabwe (no income quintiles)" = country_name %in% "Zimbabwe",
  "Missing income quintile"        = is.na(INCOME_5)
)), list(
  "Missing covariates" = !complete.cases(select(d.iwise, all_of(v.model.inc)))
)))
inc$table

# ---- 6. Raw country sizes (for the text) ------------------------------------
d.iwise %>%
  filter(country_name %in% c("Argentina", "Venezuela, RB", "Zimbabwe",
                             "Egypt, Arab Rep.", "Morocco", "India") |
         country_year %in% 19) %>%
  count(country_name, country_year)

# ---- 7. India: FLI loss -----------------------------------------------------
d.iwise %>%
  filter(country_name == "India", !is.na(iwisescore)) %>%
  summarise(n = n(),
            missing_fli = sum(is.na(FLI_3item)),
            pct_missing = round(100 * mean(is.na(FLI_3item)), 1))

# ---- 8. Check against the model samples -------------------------------------
# Should match 73,109 / 72 and 83,075 / 73 (or your current numbers)
c(fli = sum(fl$keep), income = sum(inc$keep))
c(fli = sum(d.fl$cc), income = sum(d.inc$cc))

# ---- 9. IWISE: kept vs dropped (respondents with an IWISE score) -----------
iwise_gap <- function(keep) {
  d.iwise %>%
    mutate(kept = keep) %>%
    filter(!is.na(iwisescore)) %>%
    group_by(kept) %>%
    summarise(n = n(), mean_iwise = round(mean(iwisescore), 2), .groups = "drop")
}
iwise_gap(fl$keep)
iwise_gap(inc$keep)
