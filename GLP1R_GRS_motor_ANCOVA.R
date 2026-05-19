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

MASTER_XLSX  <- "D:/PhD/EXPD3_master.xlsx"
EUR_SAMPLES  <- "D:/PhD/Exenatide_PD3_genetics/nba/raw_genotypes/EUR/EXENATIDEPD3_EUR_release11.samples"


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
# 03 — Model 1: Continuous GRS × Treatment (primary)
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


# =============================================================================
# 04 — Model 2: Binary GRS (median split) × Treatment
# =============================================================================

cat("\n=== Model 2: Binary GRS (High vs Low) × Treatment ===\n")

mod_bin <- lm(
  mds3_off_change ~ GLP1R_GRS_group * treat +
    baseline_mds3 + age + sex + diag_duration,
  data = df_complete
)

summary(mod_bin)
confint(mod_bin)

cat("\n--- HC3 robust SEs ---\n")
coeftest(mod_bin, vcov = vcovHC(mod_bin, type = "HC3"))
coefci(mod_bin,   vcov = vcovHC(mod_bin, type = "HC3"))


# =============================================================================
# 05 — Stratified betas by treatment arm
# =============================================================================

cat("\n=== Model 3: Stratified betas by treatment arm ===\n")

run_stratified <- function(data, model_formula, prs_term, label) {
  exe_mod <- lm(model_formula, data = data %>% filter(treat == "Exenatide"))
  pla_mod <- lm(model_formula, data = data %>% filter(treat == "Placebo"))

  extract <- function(mod, arm) {
    ct <- summary(mod)$coefficients
    ci <- confint(mod)
    data.frame(
      Arm      = arm,
      Beta     = ct[prs_term, "Estimate"],
      SE       = ct[prs_term, "Std. Error"],
      CI_lower = ci[prs_term, 1],
      CI_upper = ci[prs_term, 2],
      t        = ct[prs_term, "t value"],
      p        = ct[prs_term, "Pr(>|t|)"],
      N        = nrow(mod$model)
    )
  }

  bind_rows(
    extract(exe_mod, "Exenatide"),
    extract(pla_mod, "Placebo")
  ) %>% mutate(Model = label)
}

strat_cont <- run_stratified(
  df_complete,
  mds3_off_change ~ GLP1R_GRS_Z + baseline_mds3 + age + sex + diag_duration,
  "GLP1R_GRS_Z",
  "Continuous (per SD)"
)

strat_bin <- run_stratified(
  df_complete,
  mds3_off_change ~ GLP1R_GRS_group + baseline_mds3 + age + sex + diag_duration,
  "GLP1R_GRS_groupHigh",
  "Binary (High vs Low)"
)

cat("Stratified: Continuous GRS\n")
print(strat_cont)

cat("\nStratified: Binary GRS\n")
print(strat_bin)


# =============================================================================
# 06 — Summary table
# =============================================================================

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
    ifelse(p < 0.001, formatC(p, format = "f", digits = 6),
      ifelse(p < 0.01, formatC(p, format = "f", digits = 4),
        formatC(p, format = "f", digits = 3))))
}

ct_cont      <- summary(mod_cont)$coefficients
ct_bin       <- summary(mod_bin)$coefficients
int_cont_row <- grep(":", rownames(ct_cont), value = TRUE)[1]
int_bin_row  <- grep(":", rownames(ct_bin),  value = TRUE)[1]

summary_df <- data.frame(
  Model = c("Continuous (GRS × Treatment)", "Binary (Group × Treatment)"),
  Beta  = c(ct_cont[int_cont_row, "Estimate"], ct_bin[int_bin_row, "Estimate"]),
  SE    = c(ct_cont[int_cont_row, "Std. Error"], ct_bin[int_bin_row, "Std. Error"]),
  t     = c(ct_cont[int_cont_row, "t value"], ct_bin[int_bin_row, "t value"]),
  p     = c(ct_cont[int_cont_row, "Pr(>|t|)"], ct_bin[int_bin_row, "Pr(>|t|)"])
) %>%
  mutate(p_fmt = fmt_p(p))

cat("\n=== GLP1R GRS × Treatment: MDS-UPDRS III change (ANCOVA summary) ===\n\n")
print(summary_df)


# =============================================================================
# 07 — Visualisation
# =============================================================================

# Scatter: continuous GRS vs UPDRS change, regression lines by arm
p1 <- ggplot(df_complete,
             aes(x = GLP1R_GRS_Z, y = mds3_off_change, colour = treat)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
  scale_colour_manual(values = c("Placebo" = "#4575b4",
                                  "Exenatide" = "#d73027")) +
  labs(
    title    = "GLP1R GRS × Treatment: MDS-UPDRS III Off Change",
    subtitle = "ANCOVA-adjusted for baseline UPDRS, age, sex, diagnosis duration",
    x        = "GLP1R GRS (z-scored)",
    y        = "MDS-UPDRS III Off Change (week 96 \u2212 baseline)",
    colour   = "Treatment"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

print(p1)

# Boxplot: UPDRS change by GRS group and treatment arm
p2 <- ggplot(
  df_complete %>% filter(!is.na(GLP1R_GRS_group)),
  aes(x    = interaction(GLP1R_GRS_group, treat, sep = "\n"),
      y    = mds3_off_change,
      fill = treat)
) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
  scale_fill_manual(values = c("Placebo" = "#4575b4",
                                "Exenatide" = "#d73027")) +
  labs(
    title = "Motor Change by GLP1R GRS Group and Treatment",
    x     = "",
    y     = "MDS-UPDRS III Off Change (week 96 \u2212 baseline)",
    fill  = "Treatment"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

print(p2)
