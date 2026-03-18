# Evaluation of simulation study results ------------------------------------

library(tidyr)
library(dplyr)
library(ggplot2)
library(xtable)     # for LaTeX table export

# Set working directory to GitHub repo root ---------------------------------
setwd()
simulation_path <- "Simulation_Study/"

# Load simulation results (cases 1–4, various df and beta values) ----------
res_df <- read.csv(paste0(simulation_path, "output/results_cases_k39.csv"))
head(res_df)

# ---------------------------------------------------------------------------
# Pre-Test evaluation 
# ---------------------------------------------------------------------------

# Pretest without confounder ------------------------------------------------
pre_test <- res_df[, c("case", "pretest")]
t1 <- table(pre_test)
round(prop.table(t1, margin = 1) * 100, 2)  # % of rejections per case

# Pretest with confounder --------------------------------------------------
pre_test <- res_df[, c("case", "pretestH")]
t1 <- table(pre_test)
round(prop.table(t1, margin = 1) * 100, 2)

# ---------------------------------------------------------------------------
# Confounder test: basic type I/II evaluation
# ---------------------------------------------------------------------------

# Restrict to tail index combinations where confounder test is relevant
conf_test <- res_df[res_df$df1 > res_df$df2 & res_df$dfh <= res_df$df2,
                    c("case", "conftest")]
t1 <- table(conf_test)
round(prop.table(t1, margin = 1) * 100, 2)


# ---------------------------------------------------------------------------
# Prepare subsets by scenario / case (1–4) for alpha_1>alpha_2
# ---------------------------------------------------------------------------

case1 <- res_df[res_df$case == 1 &
                  res_df$df1 %in% c(3, 4) &
                  res_df$df2 %in% c(2, 3) &
                  res_df$dfh <= res_df$df2, ]

case2 <- res_df[res_df$case == 2 &
                  res_df$df1 %in% c(3, 4) &
                  res_df$df2 %in% c(2, 3) &
                  res_df$dfh <= res_df$df2, ]

case3 <- res_df[res_df$case == 3 &
                  res_df$df1 %in% c(3, 4) &
                  res_df$df2 %in% c(2, 3) &
                  res_df$dfh <= res_df$df2, ]

case4 <- res_df[res_df$case == 4 &
                  res_df$df1 %in% c(3, 4) &
                  res_df$df2 %in% c(2, 3) &
                  res_df$dfh <= res_df$df2, ]

# Put all case-specific data frames into a named list -----------------------
cases <- list(case1 = case1, case2 = case2, case3 = case3, case4 = case4)

# Labels for dfh, df1, df2 used in the plots --------------------------------
dfh_labels <- c(
  "2" = expression(alpha[h] == 2),
  "3" = expression(alpha[h] == 3),
  "4" = expression(alpha[h] == 4)
)
df1_labels <- c(
  "2" = expression(alpha[1] == 2),
  "3" = expression(alpha[1] == 3),
  "4" = expression(alpha[1] == 4)
)
df2_labels <- c(
  "2" = expression(alpha[2] == 2),
  "3" = expression(alpha[2] == 3),
  "4" = expression(alpha[2] == 4)
)

# ---------------------------------------------------------------------------
# Function to create scatter plots for confounder-test outcomes
# ---------------------------------------------------------------------------

plot_conf_test <- function(data, case_name) {
  # Assign H0/H1 as labels based on conftest
  data <- data %>%
    mutate(
      conftest_label = ifelse(conftest == 0, "H[0]", "H[1]"),
      dfh_factor     = factor(dfh, levels = c(2, 3, 4),
                              labels = c("2", "3", "4"))
    )
  
  # Facet labels for dfh (parsed math expressions) --------------------------
  dfh_labels <- as_labeller(
    c(
      "2" = "alpha[h] == 2",
      "3" = "alpha[h] == 3",
      "4" = "alpha[h] == 4"
    ),
    default = label_parsed
  )
  
  # Plot 1: β_{H->X1} vs β_{H->X2}, faceted by dfh --------------------------
  p1 <- ggplot(data, aes(x = beta_h1, y = beta_h2, color = conftest_label)) +
    geom_point(size = 4) +
    scale_color_manual(
      values = c("H[0]" = "#1f78b4", "H[1]" = "orange"),
      labels = c("H[0]" = expression(H[0]), "H[1]" = expression(H[1])),
      name   = "Test Result"
    ) +
    facet_wrap(~ dfh_factor, labeller = dfh_labels) +
    labs(
      x = expression(beta["H 1"]),
      y = expression(beta["H 2"])
    ) +
    theme_minimal() +
    theme(text = element_text(size = 40))
  
  # Save plot 1 -------------------------------------------------------------
  ggsave(
    paste0("figures/tailh_", case_name, "_k39.pdf"),
    plot   = p1,
    device = "pdf",
    width  = 14,
    height = 5,
    path   = simulation_path
  )
  
  # Plot 2: same points but faceted by (df1, df2) grid ----------------------
  p2 <- ggplot(data, aes(x = beta_h1, y = beta_h2, color = conftest_label)) +
    geom_point(size = 4) +
    scale_color_manual(
      values = c("H[0]" = "#1f78b4", "H[1]" = "orange"),
      labels = c("H[0]" = expression(H[0]), "H[1]" = expression(H[1])),
      name   = "Test Result"
    ) +
    facet_grid(
      factor(df1, levels = c(2, 3, 4), labels = df1_labels) ~
        factor(df2, levels = c(2, 3, 4), labels = df2_labels),
      labeller = label_parsed
    ) +
    labs(
      x = expression(beta["H 1" ]),
      y = expression(beta["H 2"])
    ) +
    theme_minimal() +
    theme(text = element_text(size = 40))
  
  # Save plot 2 -------------------------------------------------------------
  ggsave(
    paste0("figures/tail12_", case_name, "_k39.pdf"),
    plot   = p2,
    device = "pdf",
    width  = 14,
    height = 10,
    path   = simulation_path
  )
}

# Loop over all four cases and create plots ---------------------------------
for (case_name in names(cases)) {
  plot_conf_test(cases[[case_name]], case_name)
}

################# Tail Test  ################################################



# Hoga test evaluation: counts by (case, df1, df2, indtest) -----------------
summary_df <- res_df %>%
  group_by(case, df1, df2, indtest) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(case, df1, df2, indtest)

# Wide format: columns "0" and "1" with counts ------------------------------
summary_wide <- summary_df %>%
  pivot_wider(
    names_from   = indtest,
    values_from  = count,
    values_fill  = 0
  )

# Compute percentages of H0/H1 per (df1,df2,case) ---------------------------
summary_pct <- summary_wide %>%
  group_by(df1, df2, case) %>%
  summarise(
    total  = `0` + `1`,
    H0_pct = round(`0` / (`0` + `1`) * 100, 2),
    H1_pct = round(`1` / (`0` + `1`) * 100, 2),
    .groups = "drop"
  ) %>%
  select(df1, df2, case, H0_pct, H1_pct) %>%
  pivot_wider(
    names_from  = case,
    values_from = c(H0_pct, H1_pct),
    names_sep   = "_Case"
  )

# Reorder columns for LaTeX output ------------------------------------------
summary_pct <- summary_pct[, c(
  "df1", "df2",
  "H0_pct_Case1", "H1_pct_Case1",
  "H0_pct_Case2", "H1_pct_Case2",
  "H0_pct_Case3", "H1_pct_Case3",
  "H0_pct_Case4", "H1_pct_Case4"
)]

# LaTeX table for Hoga test -------------------------------------------------
print(
  xtable(summary_pct),
  include.rownames = FALSE,
  sanitize.text.function = identity
)

# Indicator for tail test based on normality of Hill estimator --------------
res_df$indtest2 <- as.numeric(res_df$p_val < 0.05)

# Counts by (case, df1, df2, indtest2) --------------------------------------
summary_df <- res_df %>%
  group_by(case, df1, df2, indtest2) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(case, df1, df2, indtest2)

summary_wide <- summary_df %>%
  pivot_wider(
    names_from   = indtest2,
    values_from  = count,
    values_fill  = 0
  )

print(summary_wide, n = 108)

summary_pct <- summary_wide %>%
  group_by(df1, df2, case) %>%
  summarise(
    total  = `0` + `1`,
    H0_pct = round(`0` / (`0` + `1`) * 100, 2),
    H1_pct = round(`1` / (`0` + `1`) * 100, 2),
    .groups = "drop"
  ) %>%
  select(df1, df2, case, H0_pct, H1_pct) %>%
  pivot_wider(
    names_from  = case,
    values_from = c(H0_pct, H1_pct),
    names_sep   = "_Case"
  )

summary_pct <- summary_pct[, c(
  "df1", "df2",
  "H0_pct_Case1", "H1_pct_Case1",
  "H0_pct_Case2", "H1_pct_Case2",
  "H0_pct_Case3", "H1_pct_Case3",
  "H0_pct_Case4", "H1_pct_Case4"
)]

# LaTeX table for standard Hill test ----------------------------------------
print(
  xtable(summary_pct),
  include.rownames = FALSE,
  sanitize.text.function = identity

)
