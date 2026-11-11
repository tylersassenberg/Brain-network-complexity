#lme_analyses.R
#performs group-level linear mixed-effects models across recurrence quantification analysis metrics,
#multiplex recurrence network metrics across modes and networks, and
#mode central frequency values before running experiment-wide FDR correction

#load packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(gt)

#folder containing outputs from run_MVMD.m and data_sorting.R
data_dir <- "output"

#helper function to safely relevel factors
safe_relevel <- function(vec, ref_level) {
  vec <- factor(vec)
  if (ref_level %in% levels(vec)) {
    return(relevel(vec, ref = ref_level))
  }
  return(vec)
}

#read and tidy data
modes <- 2:5
raw_list <- map(modes, ~ read.csv(file.path(data_dir, sprintf("MVMD_casnet_data_mode%d.csv", .x)))[, 2:290])

all_data <- map_dfr(raw_list, ~ mutate(.x, ID = row_number()))

tidy_data <- all_data %>%
  pivot_longer(
    cols = -ID,
    names_to = c("Network", "Metric", "Mode"),
    names_pattern = "^([^_]+)_(.*)_mode(\\d+)$"
  ) %>%
  mutate(
    Mode = factor(sprintf("%02d", as.numeric(Mode))),
    Network = factor(Network)
  )

raw_intermode <- read.csv(file.path(data_dir, "MVMD_intermode_MRN.csv"))

intermode_tidy_data <- raw_intermode %>%
  mutate(ID = row_number()) %>%
  pivot_longer(
    cols = starts_with("mode"),
    names_to = c("Mode_Pair", "Metric_Type", "Network"),
    names_pattern = "^(mode\\d+_mode\\d+)_(MI)_(.*)$"
  ) %>%
  mutate(
    Metric = paste(Metric_Type, Mode_Pair, sep = "_"),
    Network = factor(Network),
    Mode_Pair = factor(Mode_Pair)
  )

options(contrasts = c("contr.treatment", "contr.poly"))

#rqa whole cortex metrics models
rqa_metrics <- c("DET", "MEAN_dl", "ENT_dl", "LAM_vl", "TT_vl", "ENT_vl",
                 "GlobalClustering", "NetworkTransitivity", "GlobalEfficiency")

cortex_data <- tidy_data %>%
  filter(Metric %in% rqa_metrics) %>%
  group_by(ID, Metric, Mode) %>%
  summarize(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(Mode = safe_relevel(Mode, "02"))

cortex_results <- cortex_data %>%
  nest(data = -Metric) %>%
  mutate(
    model       = map(data, ~ lmer(value ~ Mode + (1|ID), data = .x)),
    anova_table = map(model, anova),
    tidy_anova  = map(anova_table, ~ suppressWarnings(tidy(.x))),
    tidy_coefs  = map(model, ~ tidy(.x, effects = "fixed"))
  )

#rqa network metrics models
rqa_data <- tidy_data %>%
  filter(Metric %in% rqa_metrics) %>%
  drop_na(value, Network, Mode) %>%
  mutate(
    Mode = safe_relevel(Mode, "02"),
    Network = safe_relevel(Network, "VisCent")
  )

rqa_results <- rqa_data %>%
  nest(data = -Metric) %>%
  mutate(
    data        = map(data, droplevels),
    model       = map(data, ~ lmer(value ~ Network + Mode + (1 | ID), data = .x)),
    anova_table = map(model, anova),
    tidy_anova  = map(anova_table, ~ suppressWarnings(tidy(.x))),
    tidy_coefs  = map(model, ~ tidy(.x, effects = "fixed"))
  )

#interlayer mutual information models
mode_data <- map_dfr(2:5, function(mode) {
  read.csv(file.path(data_dir, sprintf("MVMD_casnet_data_mode%d.csv", mode))) %>%
    select(2:290) %>%
    mutate(
      ID = row_number(),
      Mode = factor(sprintf("%02d", mode)),
      value = rowMeans(pick(154:289), na.rm = TRUE)
    ) %>%
    select(ID, Mode, value)
}) %>%
  mutate(Mode = relevel(Mode, ref = "02"))

mode_mod <- lmer(value ~ Mode + (1 | ID), data = mode_data)

net_cols <- seq(2, 98, by = 6)
net_names <- c(
  "VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB",
  "SalVentAttnA", "SalVentAttnB", "LimbicB", "LimbicA", "ContA", "ContB",
  "ContC", "DefaultA", "DefaultB", "DefaultC", "TempPar"
)

intermode <- read.csv(file.path(data_dir, "MVMD_intermode_MRN.csv"))

net_data <- map2_dfr(net_cols, net_names, function(col_start, net_name) {
  intermode %>%
    transmute(
      ID = row_number(),
      Network = net_name,
      value = rowMeans(pick(col_start:(col_start + 5)), na.rm = TRUE)
    )
}) %>%
  mutate(Network = relevel(factor(Network), ref = "VisCent"))

net_mod <- lmer(value ~ Network + (1 | ID), data = net_data)

#relative energy model
raw_energy <- read.csv(file.path(data_dir, "MVMD_network_relative_energy.csv"))

energy_tidy_data <- raw_energy %>%
  mutate(ID = row_number()) %>%
  pivot_longer(
    cols = -c(ID, subject_id),
    names_to = c("Network", "Mode"),
    names_pattern = "^([^_]+)_RE_mode(\\d+)$"
  ) %>%
  mutate(
    Mode = sprintf("%02d", as.numeric(Mode)),
    Network = safe_relevel(Network, "VisCent")
  ) %>%
  filter(Mode %in% c("02", "03", "04", "05")) %>%
  mutate(Mode = safe_relevel(Mode, "02"))

energy_mod <- lmer(value ~ Network + Mode + (1 | ID), data = energy_tidy_data)


#aggregate results and apply fdr correction
# Combine all ANOVA results
all_anova_results <- bind_rows(
  cortex_results %>% mutate(Analysis = "Whole Cortex") %>% unnest(tidy_anova),
  rqa_results    %>% mutate(Analysis = "RQA Network") %>% unnest(tidy_anova),
  tibble(Analysis = "Mode IMI",       Metric = "IMI_Mode",       tidy(anova(mode_mod))),
  tibble(Analysis = "Network IMI",    Metric = "IMI_Network",    tidy(anova(net_mod))),
  tibble(Analysis = "Relative Energy", Metric = "Relative_Energy", tidy(anova(energy_mod)))
) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr"))

# Combine all Fixed Effects results
all_fixed_effects <- bind_rows(
  cortex_results %>% mutate(Analysis = "Whole Cortex") %>% unnest(tidy_coefs),
  rqa_results    %>% mutate(Analysis = "RQA Network") %>% unnest(tidy_coefs),
  tibble(Analysis = "Mode IMI",       Metric = "IMI_Mode",       tidy(mode_mod, effects = "fixed")),
  tibble(Analysis = "Network IMI",    Metric = "IMI_Network",    tidy(net_mod, effects = "fixed")),
  tibble(Analysis = "Relative Energy", Metric = "Relative_Energy", tidy(energy_mod, effects = "fixed"))
) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr"))

#export tables
#omnibus ANOVA summary
anova_gt_table <- all_anova_results %>%
  select(Analysis, Metric, term, NumDF, DenDF, statistic, p.value, p.adj) %>%
  gt(groupname_col = "Analysis") %>%
  tab_header(
    title = "Omnibus F-Tests across Network and Mode Analyses",
  ) %>%
  cols_label(
    Metric = "Metric",
    term = "Effect",
    NumDF = "Num DF",
    DenDF = "Den DF",
    statistic = "F Statistic",
    p.value = "Raw p",
    p.adj = "FDR p"
  ) %>%
  fmt_number(columns = c(statistic, DenDF), decimals = 2) %>%
  fmt_scientific(columns = c(p.value, p.adj), decimals = 3) %>%
  tab_options(container.height = px(600), table.font.size = px(11))

#display and save
gtsave(anova_gt_table, filename = file.path(data_dir, "MVMD_Casnet_Omnibus_Anova_Results.html"))

#fixed effects summary (filtered to significant adjusted effects for readability)
coefs_gt_table <- all_fixed_effects %>%
  select(Analysis, Metric, term, estimate, std.error, statistic, p.value, p.adj) %>%
  gt(groupname_col = "Analysis") %>%
  tab_header(
    title = "Fixed Effects Parameter Estimates",
    subtitle = "Contrast differences against reference levels (Mode: 02, Network: VisCent)"
  ) %>%
  cols_label(
    Metric = "Metric",
    term = "Contrast Term",
    estimate = "Estimate (β)",
    std.error = "SE",
    statistic = "t Statistic",
    p.value = "Raw p",
    p.adj = "FDR p"
  ) %>%
  fmt_number(columns = c(estimate, std.error, statistic), decimals = 3) %>%
  fmt_scientific(columns = c(p.value, p.adj), decimals = 3) %>%
  tab_options(container.height = px(600), table.font.size = px(11))

gtsave(coefs_gt_table, filename = file.path(data_dir, "MVMD_Casnet_Fixed_Effects_Results.html"))