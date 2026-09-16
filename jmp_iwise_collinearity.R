dat <- readRDS("data_prepared.rds"); list2env(dat, envir = .GlobalEnv)

# ---- COLLINEARITY: IWISE vs at-least-basic LEVEL ----------------------------
# bas_level is the coverage series bas_rate_pp is derived from, so check both
# the IWISE overlap and the level-vs-rate overlap.

# (a) Country-year correlations, with confidence intervals.
#     Spearman for the headline (no linearity assumption), Pearson for the CI.
pairs.bas <- list(c("iwise_mean", "bas_level"),
                  c("iwise_mean", "bas_rate_pp"),
                  c("bas_level",  "bas_rate_pp"),
                  c("bas_level",  "log_gdp"))

map_dfr(pairs.bas, function(p) {
  d <- cy %>% select(all_of(p)) %>% drop_na()
  ct <- cor.test(d[[1]], d[[2]])                       # Pearson, for the CI
  tibble(x = p[1], y = p[2], n = nrow(d),
         spearman = round(cor(d[[1]], d[[2]], method = "spearman"), 2),
         pearson  = round(ct$estimate, 2),
         ci_lo    = round(ct$conf.int[1], 2),
         ci_hi    = round(ct$conf.int[2], 2),
         r2_pct   = round(100 * ct$estimate^2))         # % shared variance
})

# (b) Is the relationship linear? A curve would mean Pearson understates it.
#     Expect a floor effect: above ~95% coverage, IWISE still varies a lot.
ggplot(cy, aes(bas_level, iwise_mean)) +
  geom_point(alpha = .7) +
  geom_smooth(method = "loess", se = FALSE, colour = "#3D4222") +
  labs(title = "IWISE vs at-least-basic water coverage, country-years",
       x = "At least basic (%)", y = "Mean IWISE score") +
  theme_minimal()

# (c) Individual level - what the model actually sees. Much lower than (a),
#     because most IWISE variation is WITHIN countries.
d.fl %>%
  select(iwisescore, bas_level, bas_rate_pp, log_gdp) %>%
  cor(use = "pairwise.complete.obs") %>%
  round(2)

# (d) VIF with bas_level added, in three versions. If bas_level and log_gdp
#     cannot coexist (they correlated 0.84 at country-year level), the
#     comparison shows which one to keep.
rhs.lvl     <- paste(rhs_ctx, "+ bas_level")              # everything
rhs.lvl.nog <- paste(rhs, "+", v.wgi[6], "+", v.jmp, "+ bas_level")  # no GDP

set.seed(2024)
for (f in list(base = rhs_ctx, plus_level = rhs.lvl, no_gdp = rhs.lvl.nog)) {
  m <- lm(as.formula(paste("rnorm(nrow(d.fl)) ~", f)), data = d.fl)
  cat("\n--- ", names(which(sapply(list(rhs_ctx, rhs.lvl, rhs.lvl.nog),
                                   identical, f))), "---\n")
  print(round(car::vif(m)[, 1], 2))   # column 1 = GVIF
}
# READ: for factors, square GVIF^(1/(2*Df)) before comparing to 5 or 10.