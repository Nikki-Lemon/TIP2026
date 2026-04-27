# Matching, EDA, DiD, Event Studies, Conditional Logit

library(tidyverse)
library(MatchIt)
library(fixest)
library(survival)
library(broom)
library(car)

setwd("~/Downloads/TIP")


# Load data
person_wide <- read_csv("person_wide.csv", show_col_types = FALSE)
panel_q     <- read_csv("panel_quarterly.csv", show_col_types = FALSE) %>%
  mutate(year_quarter = factor(year_quarter))
panel_a     <- read_csv("panel_annual.csv", show_col_types = FALSE)

cat("Loaded:", nrow(person_wide), "individuals |",
    nrow(panel_q), "person-quarters |",
    nrow(panel_a), "person-years\n\n")


# Match graduates to non-graduates

cat("Before match — Graduates:", sum(person_wide$treatment == 1),
    "| Non-graduates:", sum(person_wide$treatment == 0), "\n")

match_data <- person_wide %>%
  filter(!is.na(age_at_start), !is.na(max_ogs_before)) %>%
  mutate(
    race   = factor(replace_na(as.character(race),   "Unknown")),
    gender = factor(replace_na(as.character(gender), "Unknown")),
    age_band    = factor(age_band),
    cohort_2022 = as.integer(cohort_2022)
  )

m_out <- matchit(
  treatment ~ log_earn_1yr + log_earn_2yr + job_pre_1yr +
    max_ogs_before + num_sentences_before,
  data        = match_data,
  method      = "nearest",
  ratio       = 1,
  replace     = FALSE,
  exact       = ~age_band + cohort_2022 + gender + race,
  caliper     = 0.1,
  std.caliper = TRUE
)

cat("\nPost-match balance:\n")
match_summary <- summary(m_out, un = FALSE)
bal_table <- match_summary$sum.matched
print(bal_table[, intersect(c("Means Treated", "Means Control", "Std. Mean Diff.",
                               "Means.Treated", "Means.Control"), colnames(bal_table))])

matched_df  <- match.data(m_out)
matched_ids <- matched_df %>% select(tip_id, treatment, subclass, weights)

cat("\nMatched sample — Graduates:", sum(matched_df$treatment == 1),
    "| Non-graduates:", sum(matched_df$treatment == 0), "\n")
cat("Pairs:", sum(matched_df$treatment == 1), "\n")
cat("Unmatched graduates (dropped):",
    sum(person_wide$treatment == 1) - sum(matched_df$treatment == 1), "\n\n")

write_csv(matched_ids, "matched_ids.csv")


# Filter panels to matched sample
panel_q_m <- panel_q %>%
  inner_join(matched_ids %>% select(tip_id, subclass), by = "tip_id") %>%
  mutate(
    treatment   = as.integer(treatment),
    cohort_2022 = as.integer(cohort_2022),
    graduated_tip = as.integer(graduated_tip)
  )

panel_a_m <- panel_a %>%
  inner_join(matched_ids %>% select(tip_id, subclass), by = "tip_id") %>%
  mutate(
    treatment   = as.integer(treatment),
    cohort_2022 = as.integer(cohort_2022),
    graduated_tip = as.integer(graduated_tip)
  )

person_m <- matched_df

cat("\nQuarterly panel (all periods):\n")
panel_q_m %>%
  group_by(treatment, cohort_2022) %>%
  summarise(
    employed_full    = round(mean(employed_full_q),    3),
    employed_partial = round(mean(employed_partial_q), 3),
    offended         = round(mean(offended_q),         3),
    n_persons        = n_distinct(tip_id),
    .groups = "drop"
  ) %>%
  print()

cat("\nPost-exit 18m outcomes (person-level):\n")
person_m %>%
  group_by(treatment, cohort_2022) %>%
  summarise(
    employed_full_18m = round(mean(employed_full_post_18m, na.rm = TRUE), 3),
    any_offense_18m   = round(mean(any_offense_18m, na.rm = TRUE), 3),
    n = n(),
    .groups = "drop"
  ) %>%
  print()


# Linear DiD — individual + year-quarter fixed effects, windowed <= 8 quarters post-exit

dw <- panel_q_m %>% filter(rel_quarter <= 6)
cat("Windowed N:", nrow(dw), "| Individuals:", n_distinct(dw$tip_id), "\n\n")

lpm_full <- feols(employed_full_q ~ graduated_tip + graduated_tip:cohort_2022 |
                    tip_id + year_quarter, cluster = ~tip_id, data = dw)

lpm_part <- feols(employed_partial_q ~ graduated_tip + graduated_tip:cohort_2022 |
                    tip_id + year_quarter, cluster = ~tip_id, data = dw)

lpm_off  <- feols(offended_q ~ graduated_tip + graduated_tip:cohort_2022 |
                    tip_id + year_quarter, cluster = ~tip_id, data = dw)

ols_earn <- feols(log(earnings_q_adj + 1) ~ graduated_tip + graduated_tip:cohort_2022 |
                    tip_id + year_quarter, cluster = ~tip_id, data = dw)

cat("LPM — Employed (full-time):\n");  print(summary(lpm_full))
cat("\nLPM — Employed (partial):\n");  print(summary(lpm_part))
cat("\nLPM — Offense:\n");             print(summary(lpm_off))
cat("\nLog-OLS — Earnings (adj):\n");  print(summary(ols_earn))


# Annual windowed LPM / log-OLS reference (same FE structure, aggregated to year)
# 1.5 years falls within annual year 1; rel_year <= 1 covers exit year + first full year
daw <- panel_a_m %>% filter(rel_year <= 1)
lpm_full_a <- feols(employed_full_year ~ graduated_tip + graduated_tip:cohort_2022 |
                      tip_id + year, cluster = ~tip_id, data = daw)
lpm_part_a <- feols(employed_partial_year ~ graduated_tip + graduated_tip:cohort_2022 |
                      tip_id + year, cluster = ~tip_id, data = daw)
lpm_off_a  <- feols(offended_year ~ graduated_tip + graduated_tip:cohort_2022 |
                      tip_id + year, cluster = ~tip_id, data = daw)
ols_earn_a <- feols(log(earnings_year_adj + 1) ~ graduated_tip + graduated_tip:cohort_2022 |
                      tip_id + year, cluster = ~tip_id, data = daw)

cat("\n=== Annual LPM / log-OLS DiD (windowed <= 1yr, aligns with 1.5-quarter window) ===\n")
cat("LPM — Employed (full):\n");     print(summary(lpm_full_a))
cat("\nLPM — Employed (partial):\n"); print(summary(lpm_part_a))
cat("\nLPM — Offense:\n");            print(summary(lpm_off_a))
cat("\nLog-OLS — Earnings:\n");       print(summary(ols_earn_a))


# Event studies — quarterly

panel_q_m <- panel_q_m %>%
  mutate(rel_quarter_bin = as.integer(case_when(
    rel_quarter <= -8 ~ -8L,
    rel_quarter >=  8 ~  8L,
    TRUE              ~ as.integer(rel_quarter)
  )))

dq_pre  <- panel_q_m %>% filter(cohort_2022 == 0)
dq_post <- panel_q_m %>% filter(cohort_2022 == 1)

esq_full  <- feols(employed_full_q    ~ i(rel_quarter_bin, treatment, ref = -1) | tip_id + year_quarter,
                   cluster = ~tip_id, data = panel_q_m)
esq_part  <- feols(employed_partial_q ~ i(rel_quarter_bin, treatment, ref = -1) | tip_id + year_quarter,
                   cluster = ~tip_id, data = panel_q_m)
esq_off   <- feols(offended_q         ~ i(rel_quarter_bin, treatment, ref = -1) | tip_id + year_quarter,
                   cluster = ~tip_id, data = panel_q_m)
esq_earn  <- feols(log(earnings_q_adj + 1) ~ i(rel_quarter_bin, treatment, ref = -1) | tip_id + year_quarter,
                   cluster = ~tip_id, data = panel_q_m)   # log-OLS, not LPM

esq_full_pre  <- feols(employed_full_q    ~ i(rel_quarter_bin, treatment, ref = -1) | tip_id + year_quarter,
                       cluster = ~tip_id, data = dq_pre)
esq_part_pre  <- feols(employed_partial_q ~ i(rel_quarter_bin, treatment, ref = -1) | tip_id + year_quarter,
                       cluster = ~tip_id, data = dq_pre)
esq_off_pre   <- feols(offended_q         ~ i(rel_quarter_bin, treatment, ref = -1) | tip_id + year_quarter,
                       cluster = ~tip_id, data = dq_pre)

esq_full_post <- feols(employed_full_q    ~ i(rel_quarter_bin, treatment, ref = -1) | tip_id + year_quarter,
                       cluster = ~tip_id, data = dq_post)
esq_part_post <- feols(employed_partial_q ~ i(rel_quarter_bin, treatment, ref = -1) | tip_id + year_quarter,
                       cluster = ~tip_id, data = dq_post)
esq_off_post  <- feols(offended_q         ~ i(rel_quarter_bin, treatment, ref = -1) | tip_id + year_quarter,
                       cluster = ~tip_id, data = dq_post)

png("event_study_quarterly_combined.png", width = 2400, height = 600, res = 150)
par(mfrow = c(1, 4))
iplot(esq_full,  main = "Employment (full-time)", xlab = "Quarters relative to exit")
iplot(esq_part,  main = "Employment (partial)",   xlab = "Quarters relative to exit")
iplot(esq_earn,  main = "Log earnings (adj)",     xlab = "Quarters relative to exit")
iplot(esq_off,   main = "Offense",                xlab = "Quarters relative to exit")
par(mfrow = c(1, 1))
dev.off()

png("event_study_quarterly_cohort_split.png", width = 3200, height = 1200, res = 150)
par(mfrow = c(2, 3))
iplot(esq_full_pre,  main = "Employ full — Pre-2022",    xlab = "Quarters relative to exit")
iplot(esq_part_pre,  main = "Employ partial — Pre-2022", xlab = "Quarters relative to exit")
iplot(esq_off_pre,   main = "Offense — Pre-2022",        xlab = "Quarters relative to exit")
iplot(esq_full_post, main = "Employ full — 2022+",       xlab = "Quarters relative to exit")
iplot(esq_part_post, main = "Employ partial — 2022+",    xlab = "Quarters relative to exit")
iplot(esq_off_post,  main = "Offense — 2022+",           xlab = "Quarters relative to exit")
par(mfrow = c(1, 1))
dev.off()


# Event studies — annual
da_pre  <- panel_a_m %>% filter(cohort_2022 == 0)
da_post <- panel_a_m %>% filter(cohort_2022 == 1)

esa_full  <- feols(employed_full_year    ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                   cluster = ~tip_id, data = panel_a_m)
esa_part  <- feols(employed_partial_year ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                   cluster = ~tip_id, data = panel_a_m)
esa_off   <- feols(offended_year         ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                   cluster = ~tip_id, data = panel_a_m)
esa_earn  <- feols(log(earnings_year_adj + 1) ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                   cluster = ~tip_id, data = panel_a_m)   # log-OLS, not LPM

esa_full_pre  <- feols(employed_full_year    ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year, cluster = ~tip_id, data = da_pre)
esa_part_pre  <- feols(employed_partial_year ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year, cluster = ~tip_id, data = da_pre)
esa_off_pre   <- feols(offended_year         ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year, cluster = ~tip_id, data = da_pre)
esa_full_post <- feols(employed_full_year    ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year, cluster = ~tip_id, data = da_post)
esa_part_post <- feols(employed_partial_year ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year, cluster = ~tip_id, data = da_post)
esa_off_post  <- feols(offended_year         ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year, cluster = ~tip_id, data = da_post)

png("event_study_annual_combined.png", width = 2400, height = 600, res = 150)
par(mfrow = c(1, 4))
iplot(esa_full,  main = "Employment (full-time)", xlab = "Years relative to exit")
iplot(esa_part,  main = "Employment (partial)",   xlab = "Years relative to exit")
iplot(esa_earn,  main = "Log earnings (adj)",     xlab = "Years relative to exit")
iplot(esa_off,   main = "Offense",                xlab = "Years relative to exit")
par(mfrow = c(1, 1))
dev.off()

png("event_study_annual_cohort_split.png", width = 3200, height = 1200, res = 150)
par(mfrow = c(2, 3))
iplot(esa_full_pre,  main = "Employ full — Pre-2022",    xlab = "Years relative to exit")
iplot(esa_part_pre,  main = "Employ partial — Pre-2022", xlab = "Years relative to exit")
iplot(esa_off_pre,   main = "Offense — Pre-2022",        xlab = "Years relative to exit")
iplot(esa_full_post, main = "Employ full — 2022+",       xlab = "Years relative to exit")
iplot(esa_part_post, main = "Employ partial — 2022+",    xlab = "Years relative to exit")
iplot(esa_off_post,  main = "Offense — 2022+",           xlab = "Years relative to exit")
par(mfrow = c(1, 1))
dev.off()


# Pre-trend tests-- are all pre-period event study coefficients are jointly zero.


pretrend_test <- function(model, time_var, label) {
  cn  <- names(coef(model))
  pre <- grep(paste0(time_var, "::-[2-9]|", time_var, "::-[1][0-9]"), cn, value = TRUE)
  if (length(pre) == 0) { cat(label, ": no pre-period terms\n"); return(invisible(NULL)) }
  test   <- linearHypothesis(model, paste(pre, "= 0"))
  pval   <- test$`Pr(>Chisq)`[2]
  result <- if (pval > 0.05) "PASS" else if (pval > 0.01) "MARGINAL" else "FAIL"
  cat(sprintf("%-55s  df=%d  p=%.3f  %s\n", label, length(pre), pval, result))
  invisible(test)
}

pretrend_test(esq_full,      "rel_quarter_bin", "Employ full (combined)")
pretrend_test(esq_part,      "rel_quarter_bin", "Employ partial (combined)")
pretrend_test(esq_earn,      "rel_quarter_bin", "Log earnings (combined)")
pretrend_test(esq_off,       "rel_quarter_bin", "Offense (combined)")
pretrend_test(esq_full_pre,  "rel_quarter_bin", "Employ full (pre-2022)")
pretrend_test(esq_part_pre,  "rel_quarter_bin", "Employ partial (pre-2022)")
pretrend_test(esq_off_pre,   "rel_quarter_bin", "Offense (pre-2022)")
pretrend_test(esq_full_post, "rel_quarter_bin", "Employ full (2022+)")
pretrend_test(esq_part_post, "rel_quarter_bin", "Employ partial (2022+)")
pretrend_test(esq_off_post,  "rel_quarter_bin", "Offense (2022+)")

pretrend_test(esa_full,      "rel_year_bin", "Employ full (combined)")
pretrend_test(esa_part,      "rel_year_bin", "Employ partial (combined)")
pretrend_test(esa_earn,      "rel_year_bin", "Log earnings (combined)")
pretrend_test(esa_off,       "rel_year_bin", "Offense (combined)")
pretrend_test(esa_full_pre,  "rel_year_bin", "Employ full (pre-2022)")
pretrend_test(esa_part_pre,  "rel_year_bin", "Employ partial (pre-2022)")
pretrend_test(esa_off_pre,   "rel_year_bin", "Offense (pre-2022)")
pretrend_test(esa_full_post, "rel_year_bin", "Employ full (2022+)")
pretrend_test(esa_part_post, "rel_year_bin", "Employ partial (2022+)")
pretrend_test(esa_off_post,  "rel_year_bin", "Offense (2022+)")


# Conditional logit

tidy_clogit <- function(model, label) {
  tidy(model, conf.int = TRUE) %>%
    filter(grepl("graduated_tip|treatment", term)) %>%
    select(term, estimate, conf.low, conf.high, p.value) %>%
    mutate(label = label, OR = round(exp(estimate), 3),
           OR_lo = round(exp(conf.low), 3), OR_hi = round(exp(conf.high), 3),
           .before = term)
}

cl_full <- clogit(
  employed_full_q ~ graduated_tip + graduated_tip:cohort_2022 +
    factor(year_quarter) + strata(tip_id),
  data = dw, method = "efron"
)

cl_part <- clogit(
  employed_partial_q ~ graduated_tip + graduated_tip:cohort_2022 +
    factor(year_quarter) + strata(tip_id),
  data = dw, method = "efron"
)

cl_off  <- clogit(
  offended_q ~ graduated_tip + graduated_tip:cohort_2022 +
    factor(year_quarter) + strata(tip_id),
  data = dw, method = "efron"
)

cat("\nResults (OR scale):\n")
print(bind_rows(
  tidy_clogit(cl_full, "Employed (full-time)"),
  tidy_clogit(cl_part, "Employed (partial)"),
  tidy_clogit(cl_off,  "Offense")
), n = Inf)

for (nm in c("cl_full", "cl_part", "cl_off")) {
  cat(nm, ":\n")
  print(linearHypothesis(get(nm), "graduated_tip + graduated_tip:cohort_2022 = 0"))
}

person_m2 <- person_m %>%
  filter(full_18m_window == 1) %>%
  mutate(
    cohort_2022         = as.integer(cohort_2022),
    treatment           = as.integer(treatment),
    any_sentence_before = as.integer(any_sentence_before)
  )


mp_full <- clogit(employed_full_post_18m ~ treatment + treatment:cohort_2022 +
                    any_sentence_before + strata(subclass),
                  data = person_m2, method = "efron")
mp_part <- clogit(employed_part_post_18m ~ treatment + treatment:cohort_2022 +
                    any_sentence_before + strata(subclass),
                  data = person_m2, method = "efron")
mp_off  <- clogit(any_offense_18m ~ treatment + treatment:cohort_2022 +
                    any_sentence_before + strata(subclass),
                  data = person_m2, method = "efron")

cat("\nResults (OR scale):\n")
print(bind_rows(
  tidy_clogit(mp_full, "Employed (full, 18m)"),
  tidy_clogit(mp_part, "Employed (partial, 18m)"),
  tidy_clogit(mp_off,  "Any offense (18m)")
), n = Inf)

for (nm in c("mp_full", "mp_part", "mp_off")) {
  cat(nm, ":\n")
  tryCatch(
    print(linearHypothesis(get(nm), "treatment + treatment:cohort_2022 = 0")),
    error = function(e) cat(
      "  Skipped — aliased coefficient (interaction not estimable):",
      conditionMessage(e), "\n"
    )
  )
}

# Sensitivity test — relax exact matching to cohort_2022 only
# Main spec: exact on age_band + cohort_2022 + gender + race, demographics out of PS.
# Sensitivity: exact on cohort_2022 only; age_band + gender + race enter the PS instead.

m_sens <- matchit(
  treatment ~ log_earn_1yr + log_earn_2yr + job_pre_1yr +
    max_ogs_before + num_sentences_before +
    age_band + gender + race,          # demographics now in PS, not exact cells
  data        = match_data,
  method      = "nearest",
  ratio       = 1,
  replace     = FALSE,
  exact       = ~cohort_2022,          # keep only temporal separation
  caliper     = 0.1,
  std.caliper = TRUE
)

sens_summary <- summary(m_sens, un = FALSE)
sens_bal <- sens_summary$sum.matched
print(sens_bal[, intersect(c("Means Treated", "Means Control", "Std. Mean Diff.",
                              "Means.Treated", "Means.Control"), colnames(sens_bal))])

matched_sens   <- match.data(m_sens)
matched_ids_s  <- matched_sens %>% select(tip_id, treatment, subclass, weights)

cat("\nSensitivity pairs:", sum(matched_sens$treatment == 1),
    "| Main spec pairs:", sum(matched_df$treatment == 1),
    "| Additional pairs:", sum(matched_sens$treatment == 1) - sum(matched_df$treatment == 1), "\n\n")

# Filter panels to sensitivity matched sample
panel_q_s <- panel_q %>%
  inner_join(matched_ids_s %>% select(tip_id, subclass), by = "tip_id") %>%
  mutate(
    treatment     = as.integer(treatment),
    cohort_2022   = as.integer(cohort_2022),
    graduated_tip = as.integer(graduated_tip),
    rel_quarter_bin = as.integer(case_when(
      rel_quarter <= -8 ~ -8L,
      rel_quarter >=  8 ~  8L,
      TRUE              ~ as.integer(rel_quarter)
    ))
  )

dw_s <- panel_q_s %>% filter(rel_quarter <= 8)

# Same LPM DiD on sensitivity sample
slpm_full <- feols(employed_full_q ~ graduated_tip + graduated_tip:cohort_2022 |
                     tip_id + year_quarter, cluster = ~tip_id, data = dw_s)
slpm_part <- feols(employed_partial_q ~ graduated_tip + graduated_tip:cohort_2022 |
                     tip_id + year_quarter, cluster = ~tip_id, data = dw_s)
slpm_off  <- feols(offended_q ~ graduated_tip + graduated_tip:cohort_2022 |
                     tip_id + year_quarter, cluster = ~tip_id, data = dw_s)
sols_earn <- feols(log(earnings_q_adj + 1) ~ graduated_tip + graduated_tip:cohort_2022 |
                     tip_id + year_quarter, cluster = ~tip_id, data = dw_s)

# Spec 2 matched-pair clogit on sensitivity sample
person_s <- matched_sens %>%
  filter(full_18m_window == 1) %>%
  mutate(
    cohort_2022         = as.integer(cohort_2022),
    treatment           = as.integer(treatment),
    any_sentence_before = as.integer(any_sentence_before)
  )

smp_full <- clogit(employed_full_post_18m ~ treatment + treatment:cohort_2022 +
                     any_sentence_before + strata(subclass),
                   data = person_s, method = "efron")
smp_part <- clogit(employed_part_post_18m ~ treatment + treatment:cohort_2022 +
                     any_sentence_before + strata(subclass),
                   data = person_s, method = "efron")
smp_off  <- clogit(any_offense_18m ~ treatment + treatment:cohort_2022 +
                     any_sentence_before + strata(subclass),
                   data = person_s, method = "efron")

# main spec vs. sensitivity

get_fe <- function(model, term) {
  # feols objects use $coeftable; lm objects use coef(summary())
  ct <- if (inherits(model, "fixest")) summary(model)$coeftable
        else coef(summary(model))
  if (is.null(ct) || !term %in% rownames(ct)) return(c(NA_real_, NA_real_))
  c(ct[term, "Estimate"], ct[term, "Pr(>|t|)"])
}

for (row in list(
  list("Employ full — graduated_tip",    lpm_full,  slpm_full),
  list("Employ full — ×2022+",           lpm_full,  slpm_full),
  list("Employ partial — graduated_tip", lpm_part,  slpm_part),
  list("Employ partial — ×2022+",        lpm_part,  slpm_part),
  list("Offense — graduated_tip",        lpm_off,   slpm_off),
  list("Log earn — graduated_tip",       ols_earn,  sols_earn)
)) {
  term <- if (grepl("×2022", row[[1]])) "graduated_tip:cohort_2022" else "graduated_tip"
  main <- get_fe(row[[2]], term)
  sens <- get_fe(row[[3]], term)
  cat(sprintf("%-30s  %9.3f %6.3f    %9.3f %6.3f\n",
              row[[1]], main[1], main[2], sens[1], sens[2]))
}

get_cl_or <- function(model, term) {
  sm <- summary(model)$coefficients
  if (!term %in% rownames(sm)) return(c(NA_real_, NA_real_))
  c(exp(sm[term, "coef"]), sm[term, "Pr(>|z|)"])
}

for (row in list(
  list("Employ full (18m) — treatment",    mp_full,  smp_full),
  list("Employ partial (18m) — treatment", mp_part,  smp_part),
  list("Offense (18m) — treatment",        mp_off,   smp_off)
)) {
  main <- get_cl_or(row[[2]], "treatment")
  sens <- get_cl_or(row[[3]], "treatment")
  cat(sprintf("%-30s  %9.3f %6.3f    %9.3f %6.3f\n",
              row[[1]],
              main[1], ifelse(is.na(main[2]), NA_real_, main[2]),
              sens[1], ifelse(is.na(sens[2]), NA_real_, sens[2])))
}

# Conditional earnings — intensive margin

# (a) Conditional on any employment
dw_emp  <- dw %>% filter(employed_q == 1)

ols_earn_cond <- feols(
  log(earnings_q_adj) ~ graduated_tip + graduated_tip:cohort_2022 |
    tip_id + year_quarter,
  cluster = ~tip_id, data = dw_emp
)
cat("Log earnings | employed_q == 1:\n")
print(summary(ols_earn_cond))

# (b) Conditional on full-time employment
dw_full <- dw %>% filter(employed_full_q == 1)
cat("\n(b) Full-time employment: N quarters =", nrow(dw_full),
    "| N individuals =", n_distinct(dw_full$tip_id), "\n")

ols_earn_full_cond <- feols(
  log(earnings_q_adj) ~ graduated_tip + graduated_tip:cohort_2022 |
    tip_id + year_quarter,
  cluster = ~tip_id, data = dw_full
)
print(summary(ols_earn_full_cond))

# (c) Descriptive: mean quarterly earnings by group × pre/post, conditional on employment

dw %>%
  filter(employed_q == 1) %>%
  mutate(period = if_else(rel_quarter >= 0, "post-exit", "pre-exit")) %>%
  group_by(treatment, cohort_2022, period) %>%
  summarise(
    mean_earn   = round(mean(earnings_q_adj), 0),
    median_earn = round(median(earnings_q_adj), 0),
    p25         = round(quantile(earnings_q_adj, 0.25), 0),
    p75         = round(quantile(earnings_q_adj, 0.75), 0),
    n_qtrs      = n(),
    n_persons   = n_distinct(tip_id),
    .groups = "drop"
  ) %>%
  arrange(cohort_2022, treatment, period) %>%
  print(n = Inf)

# (d) Comparison: unconditional vs. conditional earnings effect

for (row in list(
  list("Unconditional log-OLS (incl. zeros)",    ols_earn,           "graduated_tip"),
  list("Conditional: employed_q == 1",           ols_earn_cond,      "graduated_tip"),
  list("Conditional: employed_full_q == 1",      ols_earn_full_cond, "graduated_tip"),
  list("Uncond ×2022+ interaction",              ols_earn,           "graduated_tip:cohort_2022"),
  list("Cond (any empl) ×2022+ interaction",     ols_earn_cond,      "graduated_tip:cohort_2022"),
  list("Cond (full-time) ×2022+ interaction",    ols_earn_full_cond, "graduated_tip:cohort_2022")
)) {
  r <- get_fe(row[[2]], row[[3]])
  cat(sprintf("%-40s  %8.3f %6.3f\n", row[[1]], r[1], r[2]))
}
