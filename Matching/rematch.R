# Rematching script — TIP graduates vs. non-graduates
# Original match did not include enough pre-TIP earnings/employment, causing
# pre-trend violations in the event study. This rematch includes them.

library(tidyverse)
library(MatchIt)

# Load and prepare
df <- read_csv("/Users/camerondrayton/Downloads/TIP/pre_match_all_tip_data.csv",
               show_col_types = FALSE)

# Keep only graduates and non-graduates
df <- df %>%
  filter(Status %in% c("Graduated", "Non Graduated")) %>%
  mutate(
    graduated       = as.integer(Status == "Graduated"),
    log_earn_1yr    = log(earnings_1yr_preTIP + 1),
    log_earn_2yr    = log(earnings_2yr_preTIP + 1),
    age_num         = as.numeric(whole_Age_at_start),
    cohort_2022     = as.integer(cohort_2022),
    job_pre         = as.integer(job_right_before_tip),
    any_arrest      = as.integer(any_arrest_before_start),
    any_sentence    = as.integer(any_sentence_before_start),
    max_ogs         = as.numeric(max_ogs_before_start),
    num_arrests     = as.numeric(num_arrests_before_start),
    gender          = factor(replace_na(gender, "Unknown")),
    race_raw        = as.character(race),
    race            = factor(ifelse(is.na(race_raw) | race_raw %in% c("0", "", "NA"), "Unknown", race_raw)),
    age_band        = factor(case_when(
      age_num < 25             ~ "18-24",
      age_num >= 25 & age_num < 35 ~ "25-34",
      age_num >= 35            ~ "35+",
      TRUE                     ~ NA_character_
    ))
  ) %>%
  filter(!is.na(age_num), !is.na(max_ogs)) %>%
  select(-race_raw)

cat("Before matching — Graduated:", sum(df$graduated),
    "| Non-Graduated:", sum(1 - df$graduated), "\n")
cat("Race distribution (NAs =", sum(is.na(df$race)), "):\n")
print(table(df$race, useNA = "always"))

# Check pre-match balance
cat("\n=== Pre-match means by group ===\n")
df %>%
  group_by(graduated) %>%
  summarise(
    earn_1yr  = mean(earnings_1yr_preTIP),
    earn_2yr  = mean(earnings_2yr_preTIP),
    job_pre   = mean(job_pre),
    any_arrest = mean(any_arrest),
    any_sentence = mean(any_sentence),
    max_ogs   = mean(max_ogs, na.rm = TRUE),
    age       = mean(age_num, na.rm = TRUE),
    cohort_22 = mean(cohort_2022)
  ) %>%
  print()

# Propensity score matching
m_out <- matchit(
  graduated ~ log_earn_1yr + log_earn_2yr + job_pre +
    max_ogs + num_arrests,
  data    = df,
  method  = "nearest",
  ratio   = 1,
  replace = FALSE,
  exact   = ~age_band + cohort_2022 + gender + race,   # force exact match within age group, cohort, gender, and race
  caliper = 0.1,                        # reject matches with propensity score distance > 0.1 SD
  std.caliper = TRUE
)

print(summary(m_out, un = FALSE))

# Post-match balance check
# Standardized mean differences < 0.1 indicate good balance
print(summary(m_out)$sum.matched)

# Extract matched sample
matched_df <- match.data(m_out)

cat("\nMatched sample — Graduated:", sum(matched_df$graduated),
    "| Non-Graduated:", sum(1 - matched_df$graduated), "\n")

# Save the matched tip_ids so the panel data can be filtered to this new sample
matched_ids <- matched_df %>% select(tip_id, graduated, subclass, weights)
write_csv(matched_ids,
          "/Users/camerondrayton/Downloads/TIP/matched_ids_rematched.csv")

# Rebuild panel with rematched IDs
panel_raw <- read_csv("/Users/camerondrayton/Downloads/TIP/matched_panel_data_grad_non_grad (1).csv",
                      show_col_types = FALSE)

panel_rematched <- panel_raw %>%
  inner_join(matched_ids %>% select(tip_id), by = "tip_id")

cat("Panel rows (rematched):", nrow(panel_rematched),
    "| Individuals:", n_distinct(panel_rematched$tip_id), "\n")

write_csv(panel_rematched,
          "/Users/camerondrayton/Downloads/TIP/panel_rematched.csv")
