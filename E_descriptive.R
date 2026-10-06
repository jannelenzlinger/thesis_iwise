n_distinct(d.iwise$iso3c)                                    # total countries: expect 76
n_distinct(d.iwise$country_year)                             # total country-years: expect 77
d.iwise %>% distinct(country_year, country_name) %>%         # country-years missing from the FLI model
  filter(!country_year %in% d.fl$country_year[d.fl$cc])
d.fl %>%
  filter(country_name == "Argentina") %>%                    # Argentina only
  summarise(n = n(),                                         # respondents
            across(all_of(v.model.fl), ~ sum(is.na(.x)))) %>%  # NAs per model variable
  select(n, where(~ any(.x > 0)))                            # show only variables with NAs

d.iwise %>%
  filter(country_name == "Argentina") %>%                    # Argentina only
  distinct(across(all_of(v.jmp.all)))                        # all JMP variables, one row

cy %>%
  filter(bas_level >= 98) %>%                                # near-universal basic water
  select(country_year, bas_level, bas_rate_pp) %>%           # level and rate
  arrange(desc(bas_level))