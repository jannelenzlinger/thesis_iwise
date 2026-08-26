# Loading the data

# Data prep Gallup ---- 
## Set working directory and load data
# %% Set working directory and load data
install.packages("tidyverse")
library(tidyverse)
setwd("Q:/Abteilungsprojekte/Sandec/7_PHIC/02-Projects/P07-GLO-IWISE/P07-03-Analysis")
d.iwise <- read.csv("IWISE_all.csv", header = TRUE, sep = ",")
summary(d.iwise)
str(d.iwise)

## Recode Variables ----

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


# Data prep external covariates ----
sum(d.iwise$year_calendar != d.iwise$year_wave, na.rm = TRUE) #checking that calendar/wave are the same (=0) --> use calendar for matching

## Progress towards SDG 6.1 ----


## CHIRPS seasonality ----


## GDP per Capita (WB) ----


## Governance & Corruption (WB) ----