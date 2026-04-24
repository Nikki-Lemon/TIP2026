# Conditional Logit - Grad vs. Non-Grad, Annual Panel
# Annual panel built from the quarterly panel (860 matched individuals)
# Outcomes: employed_full_year, employed_partial_year, offended_year (binary); earnings_year (OLS)

library(tidyverse)
library(fixest)
library(survival)
library(broom)
library(car)

setwd("/Users/camerondrayton/Downloads/TIP")

# Build annual panel from quarterly data
dq <- read_csv("panel_rematched_quarterly.csv", show_col_types = FALSE) %>%
  mutate(
    treatment       = as.integer(treatment),
    cohort_2022     = as.integer(cohort_2022),
    graduated_tip   = as.integer(graduated_tip),
    employed_q      = as.integer(employed_q),
    employed_full_q = as.integer(employed_full_q),
    offended_q      = as.integer(offended_q),
    year            = as.integer(substr(as.character(quarter), 1, 4)),
    end_year        = as.integer(substr(as.character(end_quarter), 1, 4))
  )

matched_ids <- read_csv("matched_ids_rematched.csv", show_col_types = FALSE)

# employed_partial_year = worked at least one quarter but never hit full-time threshold
da <- dq %>%
  group_by(tip_id, year, treatment, cohort_2022, end_year) %>%
  summarise(
    employed_year         = as.integer(sum(employed_q)      > 0),
    employed_full_year    = as.integer(sum(employed_full_q) > 0),
    employed_partial_year = as.integer(sum(employed_q) > 0 & sum(employed_full_q) == 0),
    earnings_year         = sum(earnings_q, na.rm = TRUE),
    offended_year         = as.integer(sum(offended_q)      > 0),
    .groups = "drop"
  ) %>%
  mutate(
    graduated_tip = as.integer(treatment == 1 & year >= end_year),
    rel_year      = as.integer(year - end_year),
    rel_year_bin  = as.integer(case_when(
      rel_year <= -4 ~ -4L,
      rel_year >=  4 ~  4L,
      TRUE           ~ as.integer(rel_year)
    ))
  )

cat("Annual panel:", nrow(da), "rows |", n_distinct(da$tip_id), "individuals\n")
cat("rel_year range:", min(da$rel_year), "to", max(da$rel_year), "\n\n")

da_window <- da %>% filter(rel_year <= 2)
cat("Windowed N:", nrow(da_window), "| Individuals:", n_distinct(da_window$tip_id), "\n\n")

# Helpers
tidy_cl <- function(model, label) {
  tidy(model, conf.int = TRUE) %>%
    filter(term %in% c("graduated_tip", "graduated_tip:cohort_2022")) %>%
    select(term, estimate, conf.low, conf.high, p.value) %>%
    mutate(model = label, OR = round(exp(estimate), 3),
           OR_lo = round(exp(conf.low), 3), OR_hi = round(exp(conf.high), 3),
           .before = term)
}

tidy_mp <- function(model, label) {
  tidy(model, conf.int = TRUE) %>%
    filter(term %in% c("treatment", "treatment:cohort_2022")) %>%
    select(term, estimate, conf.low, conf.high, p.value) %>%
    mutate(model = label, OR = round(exp(estimate), 3),
           OR_lo = round(exp(conf.low), 3), OR_hi = round(exp(conf.high), 3),
           .before = term)
}

# LPM — Individual + Year FEs

lpm_employ   <- feols(employed_full_year ~ graduated_tip + graduated_tip:cohort_2022 | tip_id + year,
                      cluster = ~tip_id, data = da_window)
lpm_partial  <- feols(employed_partial_year ~ graduated_tip + graduated_tip:cohort_2022 | tip_id + year,
                      cluster = ~tip_id, data = da_window)
lpm_offense  <- feols(offended_year ~ graduated_tip + graduated_tip:cohort_2022 | tip_id + year,
                      cluster = ~tip_id, data = da_window)
lpm_earn     <- feols(log(earnings_year + 1) ~ graduated_tip + graduated_tip:cohort_2022 | tip_id + year,
                      cluster = ~tip_id, data = da_window)

cat("\nEmployed (full-time):\n");  print(summary(lpm_employ))
cat("\nEmployed (partial):\n");    print(summary(lpm_partial))
cat("\nOffense:\n");               print(summary(lpm_offense))
cat("\nLog earnings:\n");          print(summary(lpm_earn))

# Spec 1 — Panel clogit: strata(tip_id) + year dummies, windowed ≤2yr

cl_employ <- clogit(
  employed_full_year ~ graduated_tip + graduated_tip:cohort_2022 +
    factor(year) + strata(tip_id),
  data = da_window, method = "efron"
)

cl_partial <- clogit(
  employed_partial_year ~ graduated_tip + graduated_tip:cohort_2022 +
    factor(year) + strata(tip_id),
  data = da_window, method = "efron"
)

cl_offense <- clogit(
  offended_year ~ graduated_tip + graduated_tip:cohort_2022 +
    factor(year) + strata(tip_id),
  data = da_window, method = "efron"
)

cat("\nResults (OR scale):\n")
print(bind_rows(tidy_cl(cl_employ,  "Employed (full-time)"),
                tidy_cl(cl_partial, "Employed (partial)"),
                tidy_cl(cl_offense, "Offense")), n = Inf)

cat("\nEmployed (full-time):\n"); print(summary(cl_employ))
cat("\nEmployed (partial):\n");   print(summary(cl_partial))
cat("\nOffense:\n");              print(summary(cl_offense))

cat("\nJoint test: net effect for 2022+ graduates\n")
for (nm in c("cl_employ", "cl_partial", "cl_offense")) {
  cat(nm, ":\n")
  print(linearHypothesis(get(nm), "graduated_tip + graduated_tip:cohort_2022 = 0"))
}

# Spec 2 — Matched-pair clogit: strata(subclass), person-level post-exit
cat("\n=== Spec 2: Matched-pair clogit strata(subclass), person-level ===\n")

person_level <- da %>%
  filter(rel_year >= 0, rel_year <= 1) %>%
  group_by(tip_id, treatment, cohort_2022) %>%
  summarise(
    employed_post         = as.integer(sum(employed_full_year)    > 0),
    employed_partial_post = as.integer(sum(employed_partial_year) > 0),
    any_offense_2yr       = as.integer(sum(offended_year)         > 0),
    .groups = "drop"
  ) %>%
  left_join(matched_ids %>% select(tip_id, subclass), by = "tip_id") %>%
  mutate(cohort_2022 = as.integer(cohort_2022), treatment = as.integer(treatment))

cat("Person-level N:", nrow(person_level), "\n")
cat("Offense events:", sum(person_level$any_offense_2yr), "\n\n")

cat("Discordant pairs (offense, pre-2022):\n")
person_level %>% filter(cohort_2022 == 0) %>%
  group_by(subclass) %>%
  summarise(s = sum(any_offense_2yr), n = n(), .groups = "drop") %>%
  filter(n == 2) %>% count(discordant = s == 1) %>% print()

cat("Discordant pairs (offense, 2022+):\n")
person_level %>% filter(cohort_2022 == 1) %>%
  group_by(subclass) %>%
  summarise(s = sum(any_offense_2yr), n = n(), .groups = "drop") %>%
  filter(n == 2) %>% count(discordant = s == 1) %>% print()

mp_employ  <- clogit(employed_post         ~ treatment + treatment:cohort_2022 + strata(subclass),
                     data = person_level, method = "efron")
mp_partial <- clogit(employed_partial_post ~ treatment + treatment:cohort_2022 + strata(subclass),
                     data = person_level, method = "efron")
mp_offense <- clogit(any_offense_2yr       ~ treatment + treatment:cohort_2022 + strata(subclass),
                     data = person_level, method = "efron")

cat("\nResults (OR scale):\n")
print(bind_rows(tidy_mp(mp_employ,  "Employed (full-time, 2yr)"),
                tidy_mp(mp_partial, "Employed (partial, 2yr)"),
                tidy_mp(mp_offense, "Any offense (2yr)")), n = Inf)

cat("\nEmployed (full-time):\n"); print(summary(mp_employ))
cat("\nEmployed (partial):\n");   print(summary(mp_partial))
cat("\nAny offense:\n");          print(summary(mp_offense))

cat("\nJoint test: net effect for 2022+ graduates\n")
for (nm in c("mp_employ", "mp_partial", "mp_offense")) {
  cat(nm, ":\n")
  print(linearHypothesis(get(nm), "treatment + treatment:cohort_2022 = 0"))
}

# Comparison table: LPM vs. Panel clogit
cat("\n=== LPM vs. Panel clogit (annual, windowed ≤2yr) ===\n")
cat(sprintf("%-36s  %8s %6s  %8s %6s\n", "Outcome", "LPM coef", "p", "Clogit OR", "p"))
cat(strrep("-", 72), "\n")

get_lpm <- function(model, term) {
  ct <- coef(summary(model))
  if (!term %in% rownames(ct)) return(c(NA, NA))
  c(ct[term, "Estimate"], ct[term, "Pr(>|t|)"])
}
get_cl <- function(model, term) {
  sm <- summary(model)$coefficients
  if (!term %in% rownames(sm)) return(c(NA, NA))
  c(exp(sm[term, "coef"]), sm[term, "Pr(>|z|)"])
}

for (row in list(
  list("Employed full — graduated_tip",      lpm_employ,  cl_employ,  "graduated_tip"),
  list("Employed full — × 2022+ cohort",     lpm_employ,  cl_employ,  "graduated_tip:cohort_2022"),
  list("Employed partial — graduated_tip",   lpm_partial, cl_partial, "graduated_tip"),
  list("Employed partial — × 2022+ cohort",  lpm_partial, cl_partial, "graduated_tip:cohort_2022"),
  list("Offense — graduated_tip",            lpm_offense, cl_offense, "graduated_tip"),
  list("Offense — × 2022+ cohort",           lpm_offense, cl_offense, "graduated_tip:cohort_2022")
)) {
  lpm <- get_lpm(row[[2]], row[[4]])
  cl  <- get_cl(row[[3]],  row[[4]])
  cat(sprintf("%-36s  %8.3f %6.3f  %8.3f %6.3f\n",
              row[[1]], lpm[1], lpm[2], cl[1], cl[2]))
}

# Event studies — annual
cat("\n=== Event studies (annual) ===\n")

da_pre  <- da %>% filter(cohort_2022 == 0)
da_post <- da %>% filter(cohort_2022 == 1)

es_employ   <- feols(employed_full_year ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                     cluster = ~tip_id, data = da)
es_partial  <- feols(employed_partial_year ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                     cluster = ~tip_id, data = da)
es_earn     <- feols(log(earnings_year + 1) ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                     cluster = ~tip_id, data = da)
es_offense  <- feols(offended_year ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                     cluster = ~tip_id, data = da)

es_employ_pre    <- feols(employed_full_year ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                          cluster = ~tip_id, data = da_pre)
es_partial_pre   <- feols(employed_partial_year ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                          cluster = ~tip_id, data = da_pre)
es_earn_pre      <- feols(log(earnings_year + 1) ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                          cluster = ~tip_id, data = da_pre)
es_offense_pre   <- feols(offended_year ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                          cluster = ~tip_id, data = da_pre)

es_employ_post   <- feols(employed_full_year ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                          cluster = ~tip_id, data = da_post)
es_partial_post  <- feols(employed_partial_year ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                          cluster = ~tip_id, data = da_post)
es_earn_post     <- feols(log(earnings_year + 1) ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                          cluster = ~tip_id, data = da_post)
es_offense_post  <- feols(offended_year ~ i(rel_year_bin, treatment, ref = -1) | tip_id + year,
                          cluster = ~tip_id, data = da_post)

# Combined plot — 4 outcomes
png("event_study_annual_combined.png", width = 2400, height = 600, res = 150)
par(mfrow = c(1, 4))
iplot(es_employ,  main = "Employment (full-time)", xlab = "Years relative to exit")
iplot(es_partial, main = "Employment (partial)",   xlab = "Years relative to exit")
iplot(es_earn,    main = "Log earnings",            xlab = "Years relative to exit")
iplot(es_offense, main = "Offense",                 xlab = "Years relative to exit")
par(mfrow = c(1, 1))
dev.off()
cat("Saved: event_study_annual_combined.png\n")

# Cohort split plot — 4 outcomes × 2 cohorts
png("event_study_annual_cohort_split.png", width = 3200, height = 1200, res = 150)
par(mfrow = c(2, 4))
iplot(es_employ_pre,   main = "Employment full — Pre-2022",    xlab = "Years relative to exit")
iplot(es_partial_pre,  main = "Employment partial — Pre-2022", xlab = "Years relative to exit")
iplot(es_earn_pre,     main = "Log earnings — Pre-2022",       xlab = "Years relative to exit")
iplot(es_offense_pre,  main = "Offense — Pre-2022",            xlab = "Years relative to exit")
iplot(es_employ_post,  main = "Employment full — 2022+",       xlab = "Years relative to exit")
iplot(es_partial_post, main = "Employment partial — 2022+",    xlab = "Years relative to exit")
iplot(es_earn_post,    main = "Log earnings — 2022+",          xlab = "Years relative to exit")
iplot(es_offense_post, main = "Offense — 2022+",               xlab = "Years relative to exit")
par(mfrow = c(1, 1))
dev.off()

# Pre-trend tests
pretrend_test <- function(model, label) {
  cn  <- names(coef(model))
  pre <- grep("rel_year_bin::-[2-9]|rel_year_bin::-[1][0-9]", cn, value = TRUE)
  if (length(pre) == 0) { cat(label, ": no pre-period terms\n"); return(invisible(NULL)) }
  test   <- linearHypothesis(model, paste(pre, "= 0"))
  pval   <- test$`Pr(>Chisq)`[2]
  result <- if (pval > 0.05) "PASS" else if (pval > 0.01) "MARGINAL" else "FAIL"
  cat(sprintf("%-55s  df=%d  p=%.3f  %s\n", label, length(pre), pval, result))
  invisible(test)
}

cat("\n=== Pre-trend tests (H0: all pre-period coefs = 0) ===\n")
pretrend_test(es_employ,        "Employment full (combined)")
pretrend_test(es_partial,       "Employment partial (combined)")
pretrend_test(es_earn,          "Log earnings (combined)")
pretrend_test(es_offense,       "Offense (combined)")
pretrend_test(es_employ_pre,    "Employment full (pre-2022)")
pretrend_test(es_partial_pre,   "Employment partial (pre-2022)")
pretrend_test(es_earn_pre,      "Log earnings (pre-2022)")
pretrend_test(es_offense_pre,   "Offense (pre-2022)")
pretrend_test(es_employ_post,   "Employment full (2022+)")
pretrend_test(es_partial_post,  "Employment partial (2022+)")
pretrend_test(es_earn_post,     "Log earnings (2022+)")
pretrend_test(es_offense_post,  "Offense (2022+)")
