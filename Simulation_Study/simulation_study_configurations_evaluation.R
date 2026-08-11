# Evaluation of simulation study results ------------------------------------

library(tidyr)
library(dplyr)
library(ggplot2)
library(xtable)     # for LaTeX table export

# Set working directory to GitHub repo root ---------------------------------
setwd()
simulation_path <- "Simulation_Study/"

# Load simulation results (cases 1–4, various df and beta values) ----------
res_df <- read.csv(paste0(simulation_path, "output/results_systematic_simulation.csv"))
head(res_df)
res_df$dfh<-as.numeric(res_df$dfh)
# ---------------------------------------------------------------------------
# Causality-Test evaluation
# ---------------------------------------------------------------------------

# clear res_df where scenario is E or F and confounder is not heavier tailed 
res_df_1<-res_df[!(res_df$scenario %in% c("E", "F") &
                 pmin(res_df$df1, res_df$df2) == 2), ]
summary(res_df_1)
# Causality-Test without confounder ------------------------------------------------
causality_summary <- res_df_1 %>%
  mutate(decision = ifelse(p_val_causality <= 0.05, 1, 0)) %>%
  group_by(df1, df2, scenario) %>%
  summarise(
    H1_pct = 100 * mean(decision == 1),
    H0_pct = 100 * mean(decision == 0),
    .groups = "drop"
  )

print(causality_summary,n=44)
# Causality-Test with confounder --------------------------------------------------
causality_summary_conf <- res_df_1 %>%
  mutate(decision = ifelse(p_val_causality_conf <= 0.05, 1, 0)) %>%
  group_by(df1, df2, scenario) %>%
  summarise(
    H1_pct = 100 * mean(decision == 1),
    H0_pct = 100 * mean(decision == 0),
    .groups = "drop"
  )

print(causality_summary_conf,n=44)

# Significance Level evaluation on A and B
alpha_levels <- c(0.01, 0.025, 0.05, 0.10)

causality_alpha_summary <- lapply(alpha_levels, function(a) {
  res_df_1 %>%
    mutate(decision = ifelse(p_val_causality <= a, 1, 0),
           alpha = a) %>%
    group_by(df1, df2, scenario, alpha) %>%
    summarise(
      H1_pct = 100 * mean(decision == 1),
      H0_pct = 100 * mean(decision == 0),
      .groups = "drop"
    )
}) %>%
  bind_rows() %>%
  filter(scenario %in% c("A", "B"))

causality_alpha_wide <- causality_alpha_summary %>%
  select(df1, df2, scenario, alpha, H1_pct) %>%  # z.B. nur H1_pct
  pivot_wider(
    names_from = alpha,
    values_from = H1_pct,
    names_prefix = "alpha_"
  )

causality_alpha_wide

# ---------------------------------------------------------------------------
# Confounder test: basic type I/II evaluation
# ---------------------------------------------------------------------------

# Restrict to tail index combinations where confounder test is relevant

crit <- qnorm(0.95, mean = 0, sd = sqrt(1/12))

## 1) Without confounder (dfh NA)
conf_noH <- res_df %>%
  filter(!is.na(df1),
         !is.na(df2),
         is.na(dfh),
         df1 > df2) %>%
  mutate(decision = ifelse(Tc > crit, 1, 0)) %>%
  group_by(df1, df2, scenario) %>%
  summarise(
    reject_pct = 100 * mean(decision == 1),
    nonreject_pct = 100 * mean(decision == 0),
    .groups = "drop"
  )

## 2) With Confounder (dfh non NA)
conf_withH <- res_df_1 %>%
  filter(!is.na(df1),
         !is.na(df2),
         !is.na(dfh),
         df1 > df2) %>%
  mutate(decision = ifelse(Tc > crit, 1, 0)) %>%
  group_by(df1, df2, scenario) %>%
  summarise(
    reject_pct = 100 * mean(decision == 1),
    nonreject_pct = 100 * mean(decision == 0),
    .groups = "drop"
  )

## 3) Combine
conf_summary <- bind_rows(conf_noH, conf_withH)

conf_summary
# ---------------------------------------------------------------------------
# create scatter plots for confounder-test outcomes
# ---------------------------------------------------------------------------
res_df_1$conftest<-ifelse(res_df_1$Tc > qnorm(0.95, 0, sqrt(1 / 12)),1,0)
conf_plot_dat <- res_df_1 %>%
  filter(
    scenario %in% c("C", "D", "E", "F"),
    df1 %in% c(4),
    df2 %in% c(3)
  ) %>%
  mutate(
    conftest      = ifelse(Tc > qnorm(0.95, 0, sqrt(1 / 12)), 1, 0),
    conftest_label = ifelse(conftest == 0, "H[0]", "H[1]"),
    # Szenariolabels mit Klammern
    scenario_lab   = case_when(
      scenario == "C" ~ "(C)",
      scenario == "D" ~ "(D)",
      scenario == "E" ~ "(E)",
      scenario == "F" ~ "(F)"
    )
  )

p1<-ggplot(conf_plot_dat, aes(x = beta_h1, y = beta_h2, color = conftest_label)) +
  geom_point(size = 4) +
  scale_color_manual(
    values = c("H[0]" = "#1f78b4", "H[1]" = "orange"),
    labels = c("H[0]" = expression(H[0]), "H[1]" = expression(H[1])),
    name   = "Test Result"
  ) +
  facet_wrap(
    ~ factor(scenario_lab, levels = c("(C)", "(D)", "(E)", "(F)")),
    ncol = 2
  ) +
  labs(
    x = expression(beta["H 1"]),
    y = expression(beta["H 2"])
  ) +
  theme_minimal() +
  theme(text = element_text(size = 40))

ggsave(
  paste0("pics/conf_test_plot_k39.pdf"),
  plot   = p1,
  device = "pdf",
  width  = 14,
  height = 10,
  path   = simulation_path
)


conf_plot_dat_2 <- res_df_1 %>%
  filter(
    scenario %in% c("B"),
    df1 %in% c(3, 4),
    df2 %in% c(2, 3)
  ) %>%
  mutate(
    conftest       = ifelse(Tc > qnorm(0.95, 0, sqrt(1 / 12)), 1, 0),
    conftest_label = ifelse(conftest == 0, "H[0]", "H[1]"),
    df1_fac        = factor(df1, levels = c(3, 4),
                            labels = c(expression(alpha[1] == 3),
                                       expression(alpha[1] == 4))),
    df2_fac        = factor(df2, levels = c(2, 3),
                            labels = c(expression(alpha[2] == 2),
                                       expression(alpha[2] == 3)))
  )

p2 <- ggplot(conf_plot_dat_2,
             aes(x = beta1, y = beta_h2, color = conftest_label)) +
  geom_point(size = 4) +
  scale_color_manual(
    values = c("H[0]" = "#1f78b4", "H[1]" = "orange"),
    labels = c("H[0]" = expression(H[0]), "H[1]" = expression(H[1])),
    name   = "Test Result"
  ) +
  facet_grid(df1_fac ~ df2_fac, labeller = label_parsed) +
  labs(
    x = expression(beta[12]),
    y = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 40),
    axis.text.y  = element_blank(),  # keine Zahlen
    axis.ticks.y = element_blank()   # keine Ticks
  )
p2
ggsave(
  paste0("pics/conf_test_plot_scen_B_k39.pdf"),
  plot   = p2,
  device = "pdf",
  width  = 14,
  height = 10,
  path   = simulation_path
)

################# Tail Test  ################################################

res_df_1$indtest<-ifelse(res_df_1$Ti>55.44,1,0)

# Hoga test evaluation: counts by (case, df1, df2, indtest) -----------------
summary_df <- res_df_1 %>%
  group_by(scenario, df1, df2, indtest) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(scenario, df1, df2, indtest)

# Wide format: columns "0" and "1" with counts ------------------------------
summary_wide <- summary_df %>%
  pivot_wider(
    names_from   = indtest,
    values_from  = count,
    values_fill  = 0
  )

# Compute percentages of H0/H1 per (df1,df2,case) ---------------------------
summary_pct <- summary_wide %>%
  group_by(df1, df2, scenario) %>%
  summarise(
    total  = `0` + `1`,
    H0_pct = round(`0` / (`0` + `1`) * 100, 2),
    H1_pct = round(`1` / (`0` + `1`) * 100, 2),
    .groups = "drop"
  ) %>%
  select(df1, df2, scenario, H0_pct, H1_pct) %>%
  pivot_wider(
    names_from  = scenario,
    values_from = c(H0_pct, H1_pct),
    names_sep   = "_Scenario"
  )

# Reorder columns for LaTeX output ------------------------------------------
summary_pct <- summary_pct[, c(
  "df1", "df2",
  "H0_pct_ScenarioA", "H1_pct_ScenarioA",
  "H0_pct_ScenarioB", "H1_pct_ScenarioB",
  "H0_pct_ScenarioC", "H1_pct_ScenarioC",
  "H0_pct_ScenarioD", "H1_pct_ScenarioD",
  "H0_pct_ScenarioE", "H1_pct_ScenarioE",
  "H0_pct_ScenarioF", "H1_pct_ScenarioF"
)]

# LaTeX table for Hoga test -------------------------------------------------
print(
  xtable(summary_pct),
  include.rownames = FALSE,
  sanitize.text.function = identity
)

# Indicator for tail test based on normality of Hill estimator --------------
res_df_1$indtest2 <- as.numeric(res_df_1$p_val_hill < 0.05)

# Counts by (scenario, df1, df2, indtest2) --------------------------------------
summary_df <- res_df_1 %>%
  group_by(scenario, df1, df2, indtest2) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(scenario, df1, df2, indtest2)

summary_wide <- summary_df %>%
  pivot_wider(
    names_from   = indtest2,
    values_from  = count,
    values_fill  = 0
  )

print(summary_wide, n = 108)

summary_pct <- summary_wide %>%
  group_by(df1, df2, scenario) %>%
  summarise(
    total  = `0` + `1`,
    H0_pct = round(`0` / (`0` + `1`) * 100, 2),
    H1_pct = round(`1` / (`0` + `1`) * 100, 2),
    .groups = "drop"
  ) %>%
  select(df1, df2, scenario, H0_pct, H1_pct) %>%
  pivot_wider(
    names_from  = scenario,
    values_from = c(H0_pct, H1_pct),
    names_sep   = "_Scenario"
  )

summary_pct <- summary_pct[, c(
  "df1", "df2",
  "H0_pct_ScenarioA", "H1_pct_ScenarioA",
  "H0_pct_ScenarioB", "H1_pct_ScenarioB",
  "H0_pct_ScenarioC", "H1_pct_ScenarioC",
  "H0_pct_ScenarioD", "H1_pct_ScenarioD",
  "H0_pct_ScenarioE", "H1_pct_ScenarioE",
  "H0_pct_ScenarioF", "H1_pct_ScenarioF"
)]

# LaTeX table for standard Hill test ----------------------------------------
print(
  xtable(summary_pct),
  include.rownames = FALSE,
  sanitize.text.function = identity
)

