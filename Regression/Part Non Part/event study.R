# =========================================================
# Pre-trend / event-study analysis in R
# using panel_df.csv
# =========================================================

# -----------------------------
# 0. Packages
# -----------------------------
packages <- c(
  "readr",
  "dplyr",
  "fixest",
  "ggplot2",
  "purrr"
)

to_install <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(to_install) > 0) install.packages(to_install)

library(readr)
library(dplyr)
library(fixest)
library(ggplot2)
library(purrr)

# -----------------------------
# 1. Read data
# -----------------------------
setwd("D:/Python Directory/Capstone/R")
panel_df <- read_csv("panel_df_5.csv", show_col_types = FALSE)

# -----------------------------
# 2. Basic cleaning
# -----------------------------
panel_df <- panel_df %>%
  mutate(
    unit_id = as.character(unit_id),
    tip_id = as.numeric(tip_id),
    year = as.integer(year),
    year_offset = as.integer(year_offset),
    treatment = as.numeric(treatment),
    cohort_2022 = as.numeric(cohort_2022),
    earnings_year = as.numeric(earnings_year),
    employed_year = as.numeric(employed_year),
    has_offense_year = as.numeric(has_offense_year),
    recidivism_year = as.numeric(recidivism_year),
    log_earnings_year = log(earnings_year + 1)
  )

# quick check
cat("Rows:", nrow(panel_df), "\n")
cat("Unique unit_id:", dplyr::n_distinct(panel_df$unit_id), "\n")
cat("Event times:", sort(unique(panel_df$year_offset)), "\n")

# -----------------------------
# 3. Set event-study parameters
# -----------------------------
OMIT_YR <- -1

event_times <- sort(unique(panel_df$year_offset))
event_times <- event_times[!is.na(event_times)]

# keep a mapping for labels
event_name <- function(k) {
  if (k < 0) {
    paste0("event_m", abs(k))
  } else {
    paste0("event_p", k)
  }
}

event_map <- setNames(event_times, sapply(event_times, event_name))

# create event dummies
for (k in event_times) {
  nm <- event_name(k)
  panel_df[[nm]] <- as.numeric(panel_df$year_offset == k)
}

# omit -1
event_dummies <- sapply(event_times[event_times != OMIT_YR], event_name)

# treatment interactions
for (nm in event_dummies) {
  panel_df[[paste0(nm, "_tr")]] <- panel_df[[nm]] * panel_df$treatment
}

event_interactions <- paste0(event_dummies, "_tr")

cat("Omitted event year:", OMIT_YR, "\n")
cat("Event interactions used:\n")
print(event_interactions)

# -----------------------------
# 4. Outcomes
# -----------------------------
outcomes <- c(
  "log_earnings_year",
  "employed_year",
  "has_offense_year",
  "recidivism_year"
)

# -----------------------------
# 5. Helper to run event-study FE regression
#    unit FE + year FE, clustered by tip_id
# -----------------------------
run_event_study <- function(df, outcome, event_interactions) {
  
  needed_cols <- c("unit_id", "tip_id", "year", outcome, event_interactions)
  
  temp <- df %>%
    select(all_of(needed_cols)) %>%
    mutate(across(all_of(c(outcome, event_interactions)), as.numeric)) %>%
    mutate(across(all_of(c(outcome, event_interactions)), ~ ifelse(is.infinite(.x), NA, .x))) %>%
    filter(!is.na(.data[[outcome]]), !is.na(unit_id), !is.na(tip_id), !is.na(year))
  
  # keep only interactions with variation
  usable_event_interactions <- event_interactions[
    sapply(event_interactions, function(v) dplyr::n_distinct(temp[[v]], na.rm = TRUE) > 1)
  ]
  
  cat("\n==============================\n")
  cat("Outcome:", outcome, "\n")
  cat("Rows after cleaning:", nrow(temp), "\n")
  cat("Usable event interactions:\n")
  print(usable_event_interactions)
  
  if (length(usable_event_interactions) == 0) {
    stop(paste("No usable event interactions remain for", outcome))
  }
  
  if (dplyr::n_distinct(temp[[outcome]], na.rm = TRUE) <= 1) {
    stop(paste("Outcome has no variation for", outcome))
  }
  
  rhs <- paste(usable_event_interactions, collapse = " + ")
  
  # fixest formula with FE
  fml <- as.formula(
    paste0(outcome, " ~ ", rhs, " | unit_id + year")
  )
  
  model <- feols(
    fml,
    data = temp,
    cluster = ~tip_id
  )
  
  # extract only event-interaction terms
  coef_vec <- coef(model)
  se_vec <- se(model)
  p_vec <- pvalue(model)
  
  keep_terms <- names(coef_vec)[grepl("_tr$", names(coef_vec))]
  
  res <- data.frame(
    term = keep_terms,
    coef = unname(coef_vec[keep_terms]),
    se = unname(se_vec[keep_terms]),
    p_value = unname(p_vec[keep_terms]),
    stringsAsFactors = FALSE
  )
  
  # map term back to event time
  res$event_time <- sapply(
    gsub("_tr$", "", res$term),
    function(x) event_map[[x]]
  )
  
  res <- res %>%
    arrange(event_time) %>%
    mutate(
      ci_lower = coef - 1.96 * se,
      ci_upper = coef + 1.96 * se,
      outcome = outcome
    )
  
  list(
    model = model,
    results = res
  )
}

# -----------------------------
# 6. Run event-study for all outcomes
# -----------------------------
event_models <- list()
event_results <- list()

for (outcome in outcomes) {
  out <- run_event_study(panel_df, outcome, event_interactions)
  event_models[[outcome]] <- out$model
  event_results[[outcome]] <- out$results
}

# -----------------------------
# 7. Combine tidy results
# -----------------------------
event_results_table <- bind_rows(event_results) %>%
  arrange(outcome, event_time)

print(event_results_table)

# optional export
write_csv(event_results_table, "event_study_results.csv")

# -----------------------------
# 8. Add omitted period (-1) as zero for plotting
# -----------------------------
event_plot_table <- map_dfr(names(event_results), function(outcome) {
  temp <- event_results[[outcome]]
  
  omitted_row <- data.frame(
    term = "omitted_ref",
    coef = 0,
    se = 0,
    p_value = NA_real_,
    event_time = OMIT_YR,
    ci_lower = 0,
    ci_upper = 0,
    outcome = outcome,
    stringsAsFactors = FALSE
  )
  
  bind_rows(temp, omitted_row) %>%
    arrange(event_time)
})

print(event_plot_table)

# optional export
write_csv(event_plot_table, "event_study_plot_table.csv")

# -----------------------------
# 9. Pre-trend table only
# -----------------------------
pretrend_table <- event_plot_table %>%
  filter(event_time < 0) %>%
  arrange(outcome, event_time)

cat("\n==============================\n")
cat("Pre-trend coefficients\n")
cat("==============================\n")
print(pretrend_table)

write_csv(pretrend_table, "pretrend_table.csv")

# -----------------------------
# 10. Plot
# -----------------------------
title_map <- c(
  "log_earnings_year" = "Log Earnings",
  "employed_year" = "Employment",
  "has_offense_year" = "Any Offense",
  "recidivism_year" = "Recidivism"
)

ggplot(event_plot_table, aes(x = event_time, y = coef)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = OMIT_YR, linetype = "dashed") +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) +
  facet_wrap(~ outcome, scales = "free_y", labeller = as_labeller(title_map)) +
  labs(
    x = "Event Time (year_offset)",
    y = "Participant vs Control Effect",
    title = "Event-Study / Pre-Trend Analysis"
  ) +
  theme_minimal()

# optional save
ggsave("event_study_plot.png", width = 12, height = 6, dpi = 300)