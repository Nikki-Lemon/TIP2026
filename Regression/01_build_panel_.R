# Build from Source Files

library(tidyverse)
library(lubridate)

setwd("~/Downloads/TIP")
RAW <- "~/Downloads"

# Reference CPI: Dec 2024 = 315.605. All earnings adjusted to Dec 2024 dollars.
REF_CPI <- 315.605

# Study window: TIP start dates must be on or after MIN_DATE; exit dates on or before MAX_DATE.
MIN_DATE <- as.Date("2018-01-01")
MAX_DATE <- as.Date("2024-09-30")

# Minimum days enrolled to count as a valid non-graduate comparison.
MIN_DAYS <- 30


# Status field values:
#   Graduated:     "Graduated-Employed", "Graduated-Unemployed", "Graduated-Unknown"
#   Non Graduated: "DNF: Asked to Leave", "DNF: Dropped Out"
#   All others:    dropped (NA)
students <- read_csv(file.path(RAW, "TIP_Clean.csv"), show_col_types = FALSE) %>%
  rename(
    tip_id = participant_id,
    DOB    = DoB
  ) %>%
  mutate(
    StartDate       = as.Date(StartDate),
    EndDate         = as.Date(EndDate),
    InterviewedDate = as.Date(InterviewedDate),
    DOB             = as.Date(DOB),

    status = case_when(
      Status %in% c("Graduated-Employed", "Graduated-Unemployed", "Graduated-Unknown") ~ "Graduated",
      Status %in% c("DNF: Asked to Leave", "DNF: Dropped Out")                         ~ "Non Graduated",
      TRUE ~ NA_character_
    ),

    tip_start     = coalesce(StartDate, InterviewedDate),
    tip_end       = EndDate,
    days_enrolled = as.numeric(difftime(tip_end, tip_start, units = "days")),

    treatment   = as.integer(status == "Graduated"),
    cohort_2022 = as.integer(!is.na(tip_start) & tip_start >= as.Date("2022-01-01")),

    age_at_start = coalesce(
      as.numeric(Age_at_start),
      as.numeric(difftime(tip_start, DOB, units = "days")) / 365.25
    ),
    age_band = case_when(
      age_at_start < 25                        ~ "18-24",
      age_at_start >= 25 & age_at_start < 35  ~ "25-34",
      age_at_start >= 35                       ~ "35+",
      TRUE                                     ~ NA_character_
    ),

    end_quarter = if_else(
      !is.na(tip_end),
      as.integer(paste0(year(tip_end), ceiling(month(tip_end) / 3))),
      NA_integer_
    ),

    # How much post-exit follow-up is available before the data cutoff?
    valid_years_post      = as.numeric(difftime(MAX_DATE, tip_end, units = "days")) / 365.25,
    n_quarters_post_avail = pmin(6L, pmax(0L, as.integer(floor(valid_years_post * 4))))
  ) %>%
  filter(
    !is.na(status),
    !is.na(tip_start),
    !is.na(tip_end),
    tip_start >= MIN_DATE,   
    tip_end   <= MAX_DATE    
  )

print(count(students, status, cohort_2022))


# Load demographics and merge
demos <- read_csv(file.path(RAW, "TIP_demographics_rg.csv"), show_col_types = FALSE) %>%
  select(tip_id, gender, race) %>%
  mutate(
    gender = replace_na(as.character(gender), "Unknown"),
    race   = replace_na(as.character(race),   "Unknown")
  )

students <- students %>%
  left_join(demos, by = "tip_id") %>%
  filter(!is.na(age_at_start), !is.na(race))


# Build CPI lookup: map each year-quarter to its inflation index
# Use mid-quarter month: Q1=Feb, Q2=May, Q3=Aug, Q4=Nov
cpi_raw <- read_csv(file.path(RAW, "inflation_data.csv"), show_col_types = FALSE)

cpi_lookup <- cpi_raw %>%
  select(Year, Feb, May, Aug, Nov) %>%
  rename(Q1 = Feb, Q2 = May, Q3 = Aug, Q4 = Nov) %>%
  pivot_longer(cols = Q1:Q4, names_to = "qtr_label", values_to = "CPI") %>%
  mutate(
    qtr     = as.integer(str_extract(qtr_label, "\\d")),
    quarter = as.integer(paste0(Year, qtr))
  ) %>%
  select(quarter, CPI)


# Load and inflation-adjust earnings
# tip_cohort_uiearnings.csv is employer-quarter level — sum to person-quarter first
earnings_raw <- read_csv(file.path(RAW, "tip_cohort_uiearnings.csv"), show_col_types = FALSE) %>%
  rename_with(tolower) %>%
  filter(!is.na(tip_id), !is.na(quarter), !is.na(earnings)) %>%
  mutate(
    tip_id  = as.integer(tip_id),
    quarter = as.integer(quarter),
    earnings = as.numeric(earnings)
  ) %>%
  filter(tip_id %in% students$tip_id) %>%
  left_join(cpi_lookup, by = "quarter") %>%
  mutate(
    # Inflation-adjust: convert to Dec 2024 dollars
    earnings_adj = if_else(!is.na(CPI), earnings * (REF_CPI / CPI), earnings)
  ) %>%
  group_by(tip_id, quarter) %>%
  summarise(
    earnings_q     = sum(earnings,     na.rm = TRUE),
    earnings_q_adj = sum(earnings_adj, na.rm = TRUE),
    n_employers_q  = n_distinct(employer_legal_name, na.rm = TRUE),
    .groups = "drop"
  )


# Compute pre-TIP earnings features for matching
# Earnings in the 1 and 2 years before tip_start
earn_features <- students %>%
  select(tip_id, tip_start) %>%
  left_join(earnings_raw, by = "tip_id") %>%
  mutate(
    qtr_year = as.integer(substr(as.character(quarter), 1, 4)),
    qtr_num  = as.integer(substr(as.character(quarter), 5, 5)),
    qtr_date = as.Date(paste0(qtr_year, "-", (qtr_num - 1) * 3 + 1, "-01"))
  ) %>%
  filter(!is.na(tip_start)) %>%
  group_by(tip_id) %>%
  summarise(
    earnings_1yr_preTIP = sum(earnings_q_adj[qtr_date >= tip_start - days(365) &
                                               qtr_date < tip_start], na.rm = TRUE),
    earnings_2yr_preTIP = sum(earnings_q_adj[qtr_date >= tip_start - days(730) &
                                               qtr_date < tip_start], na.rm = TRUE),
    job_pre_1yr         = as.integer(earnings_1yr_preTIP > 0),
    .groups = "drop"
  )

students <- students %>% left_join(earn_features, by = "tip_id") %>%
  mutate(
    earnings_1yr_preTIP = replace_na(earnings_1yr_preTIP, 0),
    earnings_2yr_preTIP = replace_na(earnings_2yr_preTIP, 0),
    job_pre_1yr         = replace_na(job_pre_1yr, 0L),
    log_earn_1yr        = log(earnings_1yr_preTIP + 1),
    log_earn_2yr        = log(earnings_2yr_preTIP + 1)
  )


# Load convictions and compute pre-TIP criminal history
convictions <- read_csv(file.path(RAW, "TIP_MergedSentencingData.csv"), show_col_types = FALSE) %>%
  mutate(
    tip_id  = as.integer(tip_id),
    dof     = as.Date(DOF),
    ogs     = suppressWarnings(as.numeric(OGS))
  ) %>%
  filter(tip_id %in% students$tip_id, !is.na(dof))

# Pre-TIP offense history: any conviction before tip_start, max OGS, count
offense_features <- students %>%
  select(tip_id, tip_start) %>%
  left_join(convictions, by = "tip_id") %>%
  filter(!is.na(tip_start), dof < tip_start) %>%
  group_by(tip_id) %>%
  summarise(
    any_sentence_before  = 1L,
    num_sentences_before = n(),
    max_ogs_before       = max(ogs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(max_ogs_before = if_else(is.infinite(max_ogs_before), NA_real_, max_ogs_before))

students <- students %>%
  left_join(offense_features, by = "tip_id") %>%
  mutate(
    any_sentence_before  = replace_na(any_sentence_before, 0L),
    num_sentences_before = replace_na(num_sentences_before, 0L),
    max_ogs_before       = replace_na(max_ogs_before, 0)
  )


# Build quarterly panel spine: all quarters from MIN_DATE through MAX_DATE
all_quarters <- tibble(
  date = seq.Date(
    from = floor_date(MIN_DATE, unit = "quarter"),
    to   = floor_date(MAX_DATE, unit = "quarter"),
    by   = "quarter"
  )
) %>%
  mutate(quarter = as.integer(paste0(year(date), quarter(date)))) %>%
  pull(quarter)

valid_ids <- students %>% filter(!is.na(end_quarter)) %>% pull(tip_id)

spine <- expand_grid(tip_id = valid_ids, quarter = all_quarters)

cat("\nSpine:", nrow(spine), "person-quarters |", n_distinct(spine$tip_id), "individuals\n")

# Offset from year(MIN_DATE) so qtr_seq = 1 for the first quarter in the spine
qtr_seq <- tibble(quarter = all_quarters) %>%
  mutate(
    qtr_year = as.integer(substr(as.character(quarter), 1, 4)),
    qtr_num  = as.integer(substr(as.character(quarter), 5, 5)),
    qtr_seq  = (qtr_year - year(MIN_DATE)) * 4 + qtr_num
  )

# Load arrest data
arrests_raw <- read_csv(file.path(RAW, "arrest_data.csv"), show_col_types = FALSE) %>%
  rename_with(tolower) %>%
  mutate(tip_id = as.integer(tip_id))

# Find whichever column holds the arrest date and coerce to Date
arrest_date_col <- intersect(names(arrests_raw), c("arrestdate", "arrest_date", "offensedate", "offense_date"))
if (length(arrest_date_col) == 0) stop("No arrest date column found in arrest_data.csv")

arrests <- arrests_raw %>%
  mutate(
    arrest_date = as.Date(.data[[arrest_date_col[1]]]),
    arr_year    = year(arrest_date),
    arr_qtr     = ceiling(month(arrest_date) / 3),
    quarter     = as.integer(paste0(arr_year, arr_qtr))
  ) %>%
  filter(!is.na(quarter), tip_id %in% students$tip_id) %>%
  group_by(tip_id, quarter) %>%
  summarise(arrested_q = 1L, .groups = "drop")

# Quarterly offense indicator from convictions (date of offense)
convictions_q <- convictions %>%
  mutate(
    off_year = year(dof),
    off_qtr  = ceiling(month(dof) / 3),
    quarter  = as.integer(paste0(off_year, off_qtr))
  ) %>%
  filter(!is.na(quarter)) %>%
  group_by(tip_id, quarter) %>%
  summarise(
    offended_q      = 1L,
    offense_count_q = n(),
    max_ogs_q       = suppressWarnings(max(ogs, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(max_ogs_q = if_else(is.infinite(max_ogs_q), NA_real_, max_ogs_q))


# Assemble quarterly panel
panel_q <- spine %>%
  left_join(students %>% select(tip_id, status, treatment, cohort_2022, end_quarter,
                                tip_start, tip_end, age_band, race, gender,
                                any_sentence_before), by = "tip_id") %>%
  left_join(earnings_raw, by = c("tip_id", "quarter")) %>%
  left_join(convictions_q, by = c("tip_id", "quarter")) %>%
  left_join(arrests, by = c("tip_id", "quarter")) %>%
  left_join(qtr_seq, by = "quarter") %>%
  mutate(
    # Fill NAs (0 = not employed / no offense that quarter)
    earnings_q     = replace_na(earnings_q,     0),
    earnings_q_adj = replace_na(earnings_q_adj, 0),
    n_employers_q  = replace_na(n_employers_q,  0L),
    offended_q     = replace_na(offended_q,     0L),
    offense_count_q = replace_na(offense_count_q, 0L),
    arrested_q     = replace_na(arrested_q,     0L),

    # Employment thresholds (nominal earnings)
    employed_q         = as.integer(earnings_q > 0),
    employed_full_q    = as.integer(earnings_q >= 2250),   # ~full-time at min wage
    employed_partial_q = as.integer(earnings_q > 0 & earnings_q < 2250),

    # graduated_tip: 1 from exit quarter onward for graduates only
    graduated_tip = as.integer(treatment == 1 & quarter >= end_quarter),

    # Relative quarter to exit (0 = exit quarter, negative = before)
    end_qtr_seq   = (as.integer(substr(as.character(end_quarter), 1, 4)) - year(MIN_DATE)) * 4 +
                     as.integer(substr(as.character(end_quarter), 5, 5)),
    rel_quarter   = qtr_seq - end_qtr_seq,

    year_quarter  = as.character(quarter),
    year          = as.integer(substr(as.character(quarter), 1, 4))
  ) %>%
  select(-qtr_year, -qtr_num, -qtr_seq, -end_qtr_seq)


# Aggregate to annual panel
# employed_full_year:    had earnings in evert quarter of that calendar year
# employed_partial_year: had earnings in only some but not all quarters of that year
panel_a <- panel_q %>%
  group_by(tip_id, year, treatment, cohort_2022, status,
           end_quarter, age_band, race, gender, any_sentence_before) %>%
  summarise(
    employed_year         = as.integer(sum(employed_q) > 0),
    employed_full_year    = as.integer(sum(employed_q) == 4L),            # nonzero earnings in all 4 quarters
    employed_partial_year = as.integer(sum(employed_q) > 0 & sum(employed_q) < 4L),  # nonzero in 1–3 quarters
    earnings_year         = sum(earnings_q,     na.rm = TRUE),
    earnings_year_adj     = sum(earnings_q_adj, na.rm = TRUE),
    offended_year         = as.integer(sum(offended_q)  > 0),
    arrested_year         = as.integer(sum(arrested_q)  > 0),
    .groups = "drop"
  ) %>%
  mutate(
    end_year      = as.integer(substr(as.character(end_quarter), 1, 4)),
    graduated_tip = as.integer(treatment == 1 & year >= end_year),
    rel_year      = as.integer(year - end_year),
    rel_year_bin  = as.integer(case_when(
      rel_year <= -4 ~ -4L,
      rel_year >=  4 ~  4L,
      TRUE           ~ as.integer(rel_year)
    ))
  )

cat("\nAnnual panel:", nrow(panel_a), "rows |",
    n_distinct(panel_a$tip_id), "individuals\n")


# Build person-wide dataset (one row per person, all features)
outcomes_18m <- panel_q %>%
  filter(rel_quarter >= 0, rel_quarter <= 5) %>%   # quarters 0–5 = 18 months
  group_by(tip_id) %>%
  summarise(
    employed_post_18m      = as.integer(sum(employed_q) > 0),
    employed_full_post_18m = as.integer(sum(employed_q) == n()),          # all available quarters
    employed_part_post_18m = as.integer(sum(employed_q) > 0 & sum(employed_q) < n()),
    any_offense_18m        = as.integer(sum(offended_q) > 0),
    any_arrest_18m         = as.integer(sum(arrested_q) > 0),
    .groups = "drop"
  )

person_wide <- students %>%
  select(tip_id, status, treatment, cohort_2022, age_at_start, age_band,
         gender, race, tip_start, tip_end, end_quarter,
         earnings_1yr_preTIP, earnings_2yr_preTIP, log_earn_1yr, log_earn_2yr,
         job_pre_1yr, any_sentence_before, num_sentences_before, max_ogs_before,
         valid_years_post, n_quarters_post_avail) %>%
  left_join(outcomes_18m, by = "tip_id") %>%
  mutate(
    across(starts_with("employed_post") | starts_with("any_"), ~ replace_na(.x, 0L)),
    # Flag: 18m outcomes are only reliable for people with a full 6-quarter window.
    # People who exited TIP less than 1.5 years before MAX_DATE have truncated outcomes.
    full_18m_window = as.integer(n_quarters_post_avail >= 6)
  )


# Load social services (annual, person x year)
#services <- read_csv(file.path(RAW, "human_services.csv"), show_col_types = FALSE) %>%
  #rename_with(tolower) %>%
  #mutate(tip_id = as.integer(tip_id)) %>%
  #filter(tip_id %in% students$tip_id)


# Save outputs
write_csv(panel_q,    "panel_quarterly.csv")
write_csv(panel_a,    "panel_annual.csv")
write_csv(person_wide, "person_wide.csv")
#write_csv(services,    "services_annual.csv")
