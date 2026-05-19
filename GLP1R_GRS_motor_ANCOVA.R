# =============================================================================
# GLP1R Expression GRS × Treatment Interaction: Motor UPDRS Change
# Exenatide-PD3 Pharmacogenomics
#
# Analysis:  Baseline-adjusted ANCOVA with GxT interaction testing
# Outcome:   MDS-UPDRS Part III OFF-medication change (week 96 − baseline)
# Sample:    EUR-restricted (GP2 genetic ancestry inference, release 11)
# Models:    1) Continuous GRS (Z-scored)
#            2) Binary GRS (median split)
#            3) Stratified betas by treatment arm
# SEs:       OLS and HC3 heteroscedasticity-consistent robust SEs reported
#
# Author:    Rowan Gurney
# Date:      2026-03-05
# =============================================================================


# =============================================================================
# 00 — Packages
# =============================================================================

pkgs <- c("readxl", "readr", "dplyr", "stringr", "lmtest",
          "sandwich", "ggplot2")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}


# =============================================================================
# 01 — Paths  (edit to match your environment)
# =============================================================================

MASTER_XLSX  <- "[INSERT PATH]"
EUR_SAMPLES  <- "[INSERT PATH]"


# =============================================================================
# 02 — Load data and restrict to EUR analytic sample
# =============================================================================

# Load EUR sample list from GP2 .samples file
eur_list <- read_tsv(
  EUR_SAMPLES,
  col_names     = c("flag", "gp2_id"),
  show_col_types = FALSE
) %>%
  mutate(gp2_id = gsub("_s1$", "", gp2_id))

# Load master spreadsheet and apply EUR restriction
df <- read_xlsx(MASTER_XLSX, sheet = "Sheet1") %>%
  filter(gp2_id %in% eur_list$gp2_id, !is.na(patid)) %>%
  mutate(
    treat           = factor(str_trim(treat), levels = c("Placebo", "Exenatide")),
    GLP1R_GRS_SUM   = as.numeric(GLP1R_GRS_SUM),
    GLP1R_GRS_Z     = as.numeric(scale(GLP1R_GRS_SUM)),
    GLP1R_GRS_group = factor(
      ifelse(GLP1R_GRS_SUM >= median(GLP1R_GRS_SUM, na.rm = TRUE),
             "High", "Low"),
      levels = c("Low", "High")
    ),
    age             = as.numeric(age),
    sex             = factor(sex),
    diag_duration   = as.numeric(diag_duration),
    baseline_mds3   = as.numeric(mds3_off1),
    mds3_off_change = as.numeric(mds3_off10) - as.numeric(mds3_off1)
  )

cat("EUR sample:", n_distinct(df$gp2_id), "participants\n")

# Complete case analytic dataset
df_complete <- df %>%
  filter(complete.cases(
    mds3_off_change, GLP1R_GRS_Z, treat,
    baseline_mds3, age, sex, diag_duration
  ))

cat("Analytic N:", nrow(df_complete), "\n")
cat("Placebo:", sum(df_complete$treat == "Placebo"),
    "| Exenatide:", sum(df_complete$treat == "Exenatide"), "\n")


# =============================================================================
# 03 — Model: Continuous GRS × Treatment (primary)
# =============================================================================

cat("\n=== Model 1: Continuous GRS (per SD) × Treatment ===\n")

mod_cont <- lm(
  mds3_off_change ~ GLP1R_GRS_Z * treat +
    baseline_mds3 + age + sex + diag_duration,
  data = df_complete
)

summary(mod_cont)
confint(mod_cont)

cat("\n--- HC3 robust SEs ---\n")
coeftest(mod_cont, vcov = vcovHC(mod_cont, type = "HC3"))
coefci(mod_cont,   vcov = vcovHC(mod_cont, type = "HC3"))
