# =========================================================
# Panel regressions in R
# - OLS for earnings_year
# - clogit (survival package) for employed_year,
#   has_offense_year, recidivism_year
# =========================================================

# -----------------------------
# 0. Packages
# -----------------------------
packages <- c(
  "readr",
  "dplyr",
  "fixest",
  "survival",
  "broom"
)

to_install <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(to_install) > 0) install.packages(to_install)

library(readr)
library(dplyr)
library(fixest)
library(survival)
library(broom)

# -----------------------------
# 1. Read data
# -----------------------------
setwd("D:/Python Directory/Capstone/R")
panel_df <- read_csv("panel_df.csv", show_col_types = FALSE)
recid_panel_df <- read_csv("recid_panel_df.csv", show_col_types = FALSE)

# -----------------------------
# 2. Basic cleaning
# -----------------------------
panel_df <- panel_df %>%
  mutate(
    tip_id = as.numeric(tip_id),
    year = as.integer(year),
    treatment = as.numeric(treatment),
    post = as.numeric(post),
    cohort_2022 = as.numeric(cohort_2022),
    earnings_year = as.numeric(earnings_year),
    employed_year = as.numeric(employed_year),
    has_offense_year = as.numeric(has_offense_year),
    has_offense_by_year = as.numeric(has_offense_by_year),
    has_prior_offense_by_year = as.numeric(has_prior_offense_by_year),
    treatment_post = treatment * post,
    treatment_post_cohort = treatment * post * cohort_2022,
    log_earnings_year = log(earnings_year + 1)
  )

recid_panel_df <- recid_panel_df %>%
  mutate(
    tip_id = as.numeric(tip_id),
    year = as.integer(year),
    treatment = as.numeric(treatment),
    post = as.numeric(post),
    cohort_2022 = as.numeric(cohort_2022),
    recidivism_year = as.numeric(recidivism_year),
    has_prior_offense_by_year = as.numeric(has_prior_offense_by_year),
    treatment_post = treatment * post,
    treatment_post_cohort = treatment * post * cohort_2022
  )

cat("panel_df rows:", nrow(panel_df), "\n")
cat("recid_panel_df rows:", nrow(recid_panel_df), "\n")

# -----------------------------
# 3. OLS for earnings_year
# year FE, clustered by tip_id
# -----------------------------
model_earn_simple <- feols(
  earnings_year ~ treatment_post, 
  data = panel_df,
  cluster = ~tip_id
)

model_earn_ols <- feols(
  earnings_year ~ treatment_post + treatment_post_cohort + has_offense_by_year | factor(year),
  data = panel_df,
  cluster = ~tip_id
)

cat("\n==============================\n")
cat("OLS: earnings_year\n")
cat("==============================\n")
print(summary(model_earn_simple))
print(summary(model_earn_ols))

# optional log earnings version
model_log_earn_ols <- feols(
  log_earnings_year ~ treatment_post + treatment_post_cohort + has_offense_by_year | factor(year),
  data = panel_df,
  cluster = ~tip_id
)

cat("\n==============================\n")
cat("OLS: log_earnings_year\n")
cat("==============================\n")
print(summary(model_log_earn_ols))

# -----------------------------
# 4. Helper: keep only strata with outcome variation
# -----------------------------
keep_varying_strata <- function(df, outcome, strata_var = "tip_id") {
  df %>%
    group_by(.data[[strata_var]]) %>%
    mutate(.n_unique_y = n_distinct(.data[[outcome]], na.rm = TRUE)) %>%
    ungroup() %>%
    filter(.n_unique_y > 1) %>%
    select(-.n_unique_y)
}

# -----------------------------
# 5. clogit for employed_year
# -----------------------------
emp_df <- panel_df %>%
  select(
    tip_id, year, employed_year,
    treatment_post, treatment_post_cohort, has_offense_by_year
  ) 

cat("\nEmployment clogit rows used:", nrow(emp_df), "\n")
cat("Employment clogit unique tip_id used:", n_distinct(emp_df$tip_id), "\n")

model_emp_clogit <- clogit(
  employed_year ~ treatment_post + treatment_post_cohort + has_offense_by_year + factor(year) + strata(tip_id),
  data = emp_df,
  method = "efron"
)

cat("\n==============================\n")
cat("clogit: employed_year\n")
cat("==============================\n")
print(summary(model_emp_clogit))

# -----------------------------
# 6. clogit for has_offense_year
# -----------------------------
off_df <- panel_df %>%
  select(
    tip_id, year, has_offense_year,
    treatment_post, treatment_post_cohort, has_prior_offense_by_year
  )

cat("\nOffense clogit rows used:", nrow(off_df), "\n")
cat("Offense clogit unique tip_id used:", n_distinct(off_df$tip_id), "\n")

model_off_clogit <- clogit(
  has_offense_year ~ treatment_post + treatment_post_cohort + has_prior_offense_by_year + factor(year) + strata(tip_id),
  data = off_df,
  method = "efron"
)

cat("\n==============================\n")
cat("clogit: has_offense_year\n")
cat("==============================\n")
print(summary(model_off_clogit))

# -----------------------------
# 7. clogit for recidivism_year
# use recid_panel_df
# -----------------------------
recid_df <- recid_panel_df %>%
  select(
    tip_id, year, recidivism_year,
    treatment_post, treatment_post_cohort, has_prior_offense_by_year
  )

cat("\nRecidivism clogit rows used:", nrow(recid_df), "\n")
cat("Recidivism clogit unique tip_id used:", n_distinct(recid_df$tip_id), "\n")

model_recid_clogit <- clogit(
  recidivism_year ~ treatment_post + treatment_post_cohort + has_prior_offense_by_year + factor(year) + strata(tip_id),
  data = recid_df,
  method = "efron"
)

cat("\n==============================\n")
cat("clogit: recidivism_year\n")
cat("==============================\n")
print(summary(model_recid_clogit))

# -----------------------------
# 8. Odds ratios for clogit models
# -----------------------------
cat("\n==============================\n")
cat("Odds Ratios: employed_year clogit\n")
cat("==============================\n")
print(exp(coef(model_emp_clogit)))

cat("\n==============================\n")
cat("Odds Ratios: has_offense_year clogit\n")
cat("==============================\n")
print(exp(coef(model_off_clogit)))

cat("\n==============================\n")
cat("Odds Ratios: recidivism_year clogit\n")
cat("==============================\n")
print(exp(coef(model_recid_clogit)))

# -----------------------------
# 9. Tidy result tables
# -----------------------------
tidy_earn <- broom::tidy(model_earn_ols, conf.int = TRUE) %>%
  mutate(model = "earnings_year_ols")

tidy_log_earn <- broom::tidy(model_log_earn_ols, conf.int = TRUE) %>%
  mutate(model = "log_earnings_year_ols")

tidy_emp <- broom::tidy(model_emp_clogit, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(model = "employed_year_clogit")

tidy_off <- broom::tidy(model_off_clogit, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(model = "has_offense_year_clogit")

tidy_recid <- broom::tidy(model_recid_clogit, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(model = "recidivism_year_clogit")

all_results <- bind_rows(
  tidy_earn,
  tidy_log_earn,
  tidy_emp,
  tidy_off,
  tidy_recid
)

cat("\n==============================\n")
cat("Combined results table\n")
cat("==============================\n")
print(all_results)

# optional export
write_csv(all_results, "panel_regression_results_clogit.csv")