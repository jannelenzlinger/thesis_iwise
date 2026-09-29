d.iwise <- readRDS("data/iwise_analysis.rds")
library(tidyverse)


v.fli.items <- c("WP2319", "WP30", "WP31", "WP88")

# Count answered items per respondent
d.iwise <- d.iwise %>%
  mutate(n_fli_answered = rowSums(!is.na(across(all_of(v.fli.items)))))

# How many answered how many items, and does an index exist for them?
d.iwise %>%
  count(n_fli_answered, index_present = !is.na(INDEX_FL_ord)) %>%
  mutate(pct_of_total = round(100 * n / sum(n), 1))

# Restricted to respondents who have a computed index
d.iwise %>%
  filter(!is.na(INDEX_FL_ord)) %>%
  count(n_fli_answered) %>%
  mutate(pct_of_analytic = round(100 * n / sum(n), 1))

# Which item is missing among 3-item responders?
d.iwise %>%
  filter(n_fli_answered == 3) %>%
  summarise(across(all_of(v.fli.items), ~sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "item", values_to = "n_missing") %>%
  arrange(desc(n_missing))

# Which index levels does each group reach?
table(d.iwise$INDEX_FL_ord, d.iwise$n_fli_answered)


table(d.iwise$WP88, useNA = "ifany")
table(d.iwise$WP2319, useNA = "ifany")


d.iwise %>%
  group_by(iso3c) %>%
  summarise(n = n(), pct_wp88_na = round(100 * mean(is.na(WP88)), 1)) %>%
  arrange(desc(pct_wp88_na)) %>% print(n = Inf)

# =============================================================================
# RECALCULATING THE FINANCIAL LIFE INDEX WITHOUT WP88
# =============================================================================
#
# WHY
# ---
# WP88 ("are economic conditions in your city getting better or worse?") was
# not administered to 38,704 of your 91,166 respondents. Gallup's rule is that
# an index is computed if WP2319 is answered plus at least TWO of the other
# three items - so respondents without WP88 still get an index, but it is built
# from a two-item second component instead of a three-item one.
#
# That is what creates levels 25 and 75: they are only reachable when the
# second component averages to exactly 0.5, which requires an EVEN number of
# items. So your 9-level outcome is really two different instruments pooled
# into one variable.
#
# Dropping WP88 for EVERYONE puts all respondents on one scale (5 levels:
# 0, 25, 50, 75, 100) at no loss of sample - anyone who has WP30 and WP31 gets
# a score, which is all 79,804 who currently have an index.
#
# TRADE-OFF: WP88 is the only community-level item. Without it the index
# measures personal financial circumstances only, which is a narrower construct
# than Gallup's. Name it differently in your write-up and state the deviation.
#
# GALLUP'S CODING RULES (World Poll Methodology, Feb 2026, "Financial Life
# Index", Index Construction):
#   - WP2319: "living comfortably on present income" = 1, everything else = 0
#   - The other three: positive answer = 1, everything else = 0
#   - "Everything else" EXPLICITLY INCLUDES don't-know and refused - they are
#     zeros, not missing
#   - Only a record with NO ANSWER AT ALL is ineligible for that item
#   - Score = 100 * (mean of WP2319 + mean of the other items) / 2
# =============================================================================


# ---- STEP 0: VERIFY THE RESPONSE CODES BEFORE RECODING ----------------------
# Never recode without looking at the actual values first. The codebook says
# what the values SHOULD be; the data says what they ARE.

table(d.iwise$WP2319, useNA = "ifany")
# EXPECTED (codebook): 1 Living comfortably | 2 Getting by | 3 Difficult
#                      4 Very difficult | 5 (DK) | 6 (Refused)
# Confirmed in your data: no NA at all, all 91,166 answered.

table(d.iwise$WP30, useNA = "ifany")
# EXPECTED: 1 Satisfied | 2 Dissatisfied | 3 (DK) | 4 (Refused)

table(d.iwise$WP31, useNA = "ifany")
# EXPECTED: 1 Getting better | 2 (The same) | 3 Getting worse | 4 (DK) | 5 (Refused)

table(d.iwise$WP88, useNA = "ifany")
# EXPECTED: 1 Getting better | 2 (The same) | 3 Getting worse | 4 (DK) | 5 (Refused)
# Confirmed in your data: codes 4 and 5 are PRESENT (1,382 and 73), so DK and
# Refused were NOT stripped to NA. The 38,704 NAs are genuine non-administration
# - a fielding decision, not respondent behaviour. That is what makes dropping
# WP88 defensible rather than arbitrary.


# ---- STEP 1: RECODE EACH ITEM TO 0/1 ----------------------------------------
# The positive answer is code 1 for ALL FOUR items. (Do not confuse these with
# the food/shelter items WP40/WP43, where the positive answer is "No" = 2.)

d.iwise <- d.iwise %>%
  mutate(

    # WP2319: 1 = "living comfortably on present income".
    # if_else() returns 1 when the condition is TRUE, 0 when FALSE, and - this
    # is the important part - NA when the input is NA. So codes 2-6 (getting
    # by, difficult, very difficult, DK, refused) all become 0, exactly as
    # Gallup specifies, while a genuinely unanswered item stays NA and is
    # excluded from the calculation.
    s_2319 = if_else(WP2319 == 1, 1, 0),

    # WP30: 1 = "Satisfied" with standard of living.
    # Codes 2 (dissatisfied), 3 (DK), 4 (refused) -> 0.
    s_30 = if_else(WP30 == 1, 1, 0),

    # WP31: 1 = standard of living "getting better".
    # Codes 2 (the same), 3 (getting worse), 4 (DK), 5 (refused) -> 0.
    # NOTE "the same" scores 0 - only improvement counts as positive.
    s_31 = if_else(WP31 == 1, 1, 0),

    # WP88: 1 = local economic conditions "getting better".
    # Same structure as WP31. Needed only for the validation in Step 2.
    s_88 = if_else(WP88 == 1, 1, 0)
  )

# Sanity check: each should be 0/1 with NAs only where the item was unanswered.
d.iwise %>%
  summarise(across(c(s_2319, s_30, s_31, s_88),
                   list(n_1 = ~sum(.x == 1, na.rm = TRUE),
                        n_0 = ~sum(.x == 0, na.rm = TRUE),
                        n_NA = ~sum(is.na(.x)))))


# ---- STEP 2: VALIDATE BY REPRODUCING GALLUP'S OWN 4-ITEM INDEX --------------
# Do NOT skip this. If your recoding is wrong, the 3-item version inherits the
# error silently. Reconstruct Gallup's index and check it matches theirs.

d.iwise <- d.iwise %>%
  mutate(

    # How many of the three non-WP2319 items does this respondent have?
    # across() selects the three columns; !is.na() gives TRUE where answered;
    # rowSums() counts the TRUEs per respondent (TRUE counts as 1).
    n_other = rowSums(!is.na(across(c(s_30, s_31, s_88)))),

    # The second component: the MEAN of whichever of the three are available.
    # na.rm = TRUE is what makes the denominator vary - with all three it
    # divides by 3, with two it divides by 2. That varying denominator is the
    # entire source of the 25/75 levels.
    mean_other = rowMeans(across(c(s_30, s_31, s_88)), na.rm = TRUE),

    # Gallup's eligibility rule: WP2319 answered AND at least 2 of the other 3.
    FLI_check = if_else(
      !is.na(s_2319) & n_other >= 2,
      100 * (s_2319 + mean_other) / 2,   # average the two components, x100
      NA_real_                            # otherwise no index
    )
  )

# Agreement with Gallup's variable. Should be ~1.00.
# The < 0.01 tolerance allows for floating-point rounding (33.33 vs 33.333...).
d.iwise %>%
  filter(!is.na(INDEX_FL_ord)) %>%
  summarise(
    n = n(),
    agree = mean(abs(FLI_check -
                     as.numeric(as.character(INDEX_FL_ord))) < 0.01,
                 na.rm = TRUE)
  ) #agree = 1


# ---- STEP 3: BUILD THE 3-ITEM INDEX -----------------------------------------
# Same recipe, WP88 removed. The second component is now the mean of WP30 and
# WP31 only - a FIXED denominator of 2 for everyone, which is the whole point.

d.iwise <- d.iwise %>%
  mutate(

    # Eligibility: WP2319 plus BOTH remaining items. Gallup's "at least 2 of 3"
    # becomes "2 of 2" here. Requiring both keeps the second component from
    # resting on a single question, and costs nothing - every respondent who
    # currently has an index already has WP30 and WP31.
    FLI_3item_num = if_else(
      !is.na(s_2319) & !is.na(s_30) & !is.na(s_31),

      # Component 1 = s_2319, which is 0 or 1.
      # Component 2 = (s_30 + s_31) / 2, which is 0, 0.5, or 1.
      # Average them and multiply by 100.
      # Possible results: 0, 25, 50, 75, 100 - five levels, and every
      # respondent's score is built the same way.
      100 * (s_2319 + (s_30 + s_31) / 2) / 2,

      NA_real_
    ),

    # Ordered factor, so polr()/clm() treat it as ordinal.
    # as.numeric() first so levels sort 0 < 25 < 50 < 75 < 100 rather than
    # alphabetically ("0","100","25","50","75") - the same trap as Check 1.
    FLI_3item = factor(FLI_3item_num, ordered = TRUE)
  )


# ---- STEP 4: CHECK THE RESULT -----------------------------------------------

# Distribution. Expect exactly 5 levels and N = 79,804. --> TRUE
table(d.iwise$FLI_3item, useNA = "ifany")
round(100 * prop.table(table(d.iwise$FLI_3item)), 1)

# Did anyone gain or lose an index relative to Gallup's version?
# Rows in "TRUE / FALSE" gained one; "FALSE / TRUE" lost one (should be zero) --> TRUE
table(has_gallup = !is.na(d.iwise$INDEX_FL_ord),
      has_3item  = !is.na(d.iwise$FLI_3item))

# How do the two versions map onto each other? Shows what dropping WP88 does
# to each respondent's score.
table(new = d.iwise$FLI_3item, old = d.iwise$INDEX_FL_ord)

# Are the two versions still measuring much the same thing? A high correlation
# supports treating the 3-item version as a valid substitute.
cor(as.numeric(as.character(d.iwise$FLI_3item)),
    as.numeric(as.character(d.iwise$INDEX_FL_ord)),
    use = "pairwise.complete.obs") #cor = 0.975

# Does the new index still behave sensibly against your main predictor?
# Expect mean iwisescore to fall as financial life improves. 
d.iwise %>%
  group_by(FLI_3item) %>%
  summarise(n = n(),
            mean_iwise = round(mean(iwisescore, na.rm = TRUE), 2),
            .groups = "drop") #somewhat true

# Confirm WP88 non-administration is a COUNTRY-level pattern, which is the
# evidence that this is a fielding decision rather than respondent refusal.
d.iwise %>%
  group_by(iso3c) %>%
  summarise(n = n(),
            pct_wp88_missing = round(100 * mean(is.na(WP88)), 1),
            .groups = "drop") %>%
  arrange(desc(pct_wp88_missing)) %>%
  print(n = Inf)


# ---- STEP 5: SWAP IT INTO THE PRE-MODEL SCRIPT ------------------------------
# In 01_premodel_checks.R, replace INDEX_FL_ord with FLI_3item in:
#   - v.outcomes
#   - the d.fl definition
#   - v.model.fl
#   - Checks 3, 5, 6, 8, 9
# Then re-run. Sample size is unchanged (79,804), so the missingness findings
# carry over, but the outcome now has 5 levels instead of 9 - which changes the
# cut-point count, the EPV denominator, and the number of separation cuts.

# =============================================================================