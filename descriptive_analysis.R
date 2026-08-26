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







### Testing the assumptions for the logistic regression

independence of observations

no multicollinearity

linearity of predictors with the log-odds


adequate events per variable.

