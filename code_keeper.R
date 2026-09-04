# 7c. All six WGI indicators together, to show why only one can be used.
m.vifwgi <- lm(as.formula(paste("rnorm(nrow(d.fl)) ~",
                                paste(v.wgi, collapse = " + "))), data = d.fl)
car::vif(m.vifwgi)

# 7d. IF COLLINEAR: pick ONE on theoretical grounds, or reduce all six to a
#     single component. Run the PCA at COUNTRY level, not respondent level, or
#     the largest country dominates the components.
d.wgi.country <- d.iwise %>%
  select(all_of(c(v.cluster, v.wgi))) %>%
  distinct() %>%
  drop_na()

pca.wgi <- prcomp(d.wgi.country %>% select(all_of(v.wgi)), scale. = TRUE)
summary(pca.wgi)        # PC1 usually explains 80-90% of WGI variance
pca.wgi$rotation[, 1]   # loadings: all six should point the same direction -> TRUE

# If you use it, merge PC1 onto ALL THREE frames so they stay consistent.
d.pc1 <- d.wgi.country %>%
  mutate(wgi_pc1 = pca.wgi$x[, 1]) %>%
  select(all_of(v.cluster), wgi_pc1)

d.iwise <- left_join(d.iwise, d.pc1, by = v.cluster)
d.fl    <- left_join(d.fl,    d.pc1, by = v.cluster)
d.inc   <- left_join(d.inc,   d.pc1, by = v.cluster)

# --------------------------------------------------------
# --------------------------------------------------------
d.iwise %>%
  summarise(
    n_total   = n(),
    n_na      = sum(is.na(INCOME_4)),
    pct_na    = round(100 * mean(is.na(INCOME_4)), 2),
    n_zero    = sum(INCOME_4 == 0, na.rm = TRUE),
    pct_zero  = round(100 * mean(INCOME_4 == 0, na.rm = TRUE), 2),
    n_valid   = sum(!is.na(INCOME_4) & INCOME_4 > 0),
    n_neg     = sum(INCOME_4 < 0, na.rm = TRUE)
  ) %>%
  glimpse()

# Distribution of non-missing values, in daily international dollars
summary(d.iwise$INCOME_4 / 365)
quantile(d.iwise$INCOME_4 / 365, seq(0, 1, 0.1), na.rm = TRUE)

# Zeros and NAs by country — coverage gap or item non-response?
d.iwise %>%
  group_by(iso3c) %>%
  summarise(n = n(),
            pct_na   = round(100 * mean(is.na(INCOME_4)), 1),
            pct_zero = round(100 * mean(INCOME_4 == 0, na.rm = TRUE), 1),
            .groups = "drop") %>%
  filter(pct_na > 0 | pct_zero > 0) %>%
  arrange(desc(pct_na)) %>%
  print(n = Inf)

# How much of INCOME_4 is imputed rather than reported?
table(d.iwise$INCOME_7, useNA = "ifany")