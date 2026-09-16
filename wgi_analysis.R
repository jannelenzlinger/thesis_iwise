library(dplyr)
library(car)     # vif()
install.packages("psych")   # corr.test(), KMO()
install.packages("domir")# dominance analysis
library(psych)

gov <- c("wgi_va_est", "wgi_pv_est", "wgi_ge_est",
         "wgi_rq_est", "wgi_rl_est", "wgi_cc_est")

#==============================================================
# Country-level data (one row per country_year)
#==============================================================
cl <- d.iwise %>%
  filter(!is.na(iwise12_imp), !is.na(projectionwgt_svy_yr)) %>%
  group_by(country_year) %>%
  summarise(
    across(all_of(gov), ~ first(.x)),   # constant within country_year
    wi_prev = 100 * weighted.mean(iwise12_imp, projectionwgt_svy_yr),
    pop     = sum(projectionwgt_svy_yr),
    .groups = "drop"
  ) %>%
  filter(if_all(all_of(gov), ~ !is.na(.x)))

#==============================================================
# STEP 1: Collinearity and PCA
#==============================================================
corr.test(cl[gov])
vif(lm(reformulate(gov, "wi_prev"), data = cl))

pca <- prcomp(cl[gov], scale. = TRUE)
summary(pca)                                        # variance explained
pca$sdev^2                                          # eigenvalues (Kaiser > 1)
screeplot(pca, type = "lines"); abline(h = 1, lty = 2)
pca$rotation                                        # loadings
KMO(cor(cl[gov]))

cl$gov_pc1 <- pca$x[, 1]
if (mean(pca$rotation[, 1]) < 0) cl$gov_pc1 <- -cl$gov_pc1   # higher = better governance

#==============================================================
# STEP 2: Country-level preview
#==============================================================
vars <- c("wi_prev", gov, "gov_pc1")

corr.test(cl[vars])                                      # unweighted
round(cov.wt(cl[vars], wt = cl$pop, cor = TRUE)$cor, 3)  # population-weighted

# Bivariate OLS: R² and slope per indicator
sapply(c(gov, "gov_pc1"), function(g) {
  m <- lm(reformulate(g, "wi_prev"), data = cl)
  c(R2 = summary(m)$r.squared, b = unname(coef(m)[2]))
}) |> t() |> round(3)

# Dominance analysis (share of R² per indicator)
domir(reformulate(gov, "wi_prev"),
      function(fml) summary(lm(fml, data = cl))$r.squared)