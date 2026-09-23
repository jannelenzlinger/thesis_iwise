# First modelling try

# 10. PROPORTIONAL ODDS - Brant test, ordinal::nominal_test(), and a plot of
#     coefficients across cuts. With your N the formal tests reject on trivial
#     departures, so the plot is the most honest of the three.

# 11. LINEARITY OF THE LOGIT - iwisescore is zero-inflated (median 3, mean
#     6.6), so expect this to flag; iwise_cat is the alternative. Test
#     bas_rate_pp with splines (it can be negative).

# 12. INFLUENTIAL OBSERVATIONS - Cook's distance, standardised residuals, refit
#     without the most influential 1%. Also refit dropping one country at a
#     time and watch the country-level estimates.

# 13. GOODNESS OF FIT - Lipsitz and Pulkstenis-Robinson tests.