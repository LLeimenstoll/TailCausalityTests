# Confounder Simulation Study -----------------------------------------------

library(dplyr)   # data manipulation (mutate, group_by, summarise)
library(ggplot2) # plotting
library(tidyr)   # reshaping 

# Set working directory to GitHub repo root ---------------------------------
setwd()
simulation_path <- "Simulation_Study/"

# Basic simulation setup ----------------------------------------------------
n      <- 10000                          # sample size per replication
vs     <- seq(0.15, 0.5, 0.05)           # exponents for k = n^v
ks     <- floor(n^vs)                    # tail sample sizes
rounds <- 10000                          # number of replications

# Combinations of tail indices (degrees of freedom of t): (df1, df2) --------
dfs <- list(c(3, 2), c(4, 2), c(4, 3))

results <- data.frame()  # will store df1, df2, k, beta, beta_h1, beta_h2, conftest
r <- 1                   # row index in results

# Main simulation loop ------------------------------------------------------
for (df in dfs) {
  print(df)
  
  for (i in 1:rounds) {
    print(i)
    
    # Potential confounder H (commented out here) -----------
    # H <- rt(n, 2)
    
    df1 <- df[1]
    df2 <- df[2]
    
    # Structural coefficient for X2 <- X1 -----------------------------------
    beta_1 <- runif(1, 0.1, 0.9)
    
    # Coefficients for confounder H  ------------------------
    # beta_h1 <- runif(1, 0.1, 0.9)
    # beta_h2 <- runif(1, 0.1, 0.9)
    
    # Heavy-tailed errors for X1, X2 ----------------------------------------
    e1 <- rt(n, df1)
    e2 <- rt(n, df2)
    
    # With confounder  -------------------------------------
    # X1 <- e1 + beta_h1 * H
    # X2 <- e2 + beta_1 * X1 + beta_h2 * H
    
    # Without confounder -------------------------------------
    X1 <- e1
    X2 <- e2 + beta_1 * X1
    
    # Loop over different k choices -----------------------------------------
    for (k in ks) {
      conftest <- 0  # indicator: 1 if confounder test rejects null
      
      # Nonparametric CTCs in both directions -------------------------------
      c12 <- causal_tail_coeff_basic(X1, X2, k = k, both_tails = FALSE)
      c21 <- causal_tail_coeff_basic(X2, X1, k = k, both_tails = FALSE)
      
      # Confounder test statistic: min(c12, c21) vs 0.5 ---------------------
      Tc <- sqrt(k) * (min(c12, c21) - 0.5)
      if (Tc > qnorm(0.95, 0, sqrt(1 / 12))) {
        conftest <- 1
      }
      
      # NOTE: beta_h1, beta_h2 are not defined in this (no confounder) block;
      #       they remain NA in this "without confounder" scenario.
      results[r, 1:7] <- c(df1, df2, k, beta_1, beta_h1, beta_h2, conftest)
      r <- r + 1
    }
  }
}

results

# Name the columns ----------------------------------------------------------
colnames(results) <- c("df1", "df2", "k", "beta", "beta_h1", "beta_h2", "conftest")

# Save results for the "no confounder" case ---------------------------------

write.csv(results,
          file = paste0(simulation_path, "output/results_confounder_test.csv"),
          row.names = FALSE)

# Load results for plotting -------------------------------------------------
results1 <- read.csv(file = paste0(simulation_path, "output/results_confounder_test.csv"))
results2 <- read.csv(file = paste0(simulation_path, "output/results_confounder_test_with_confs.csv"))

# ---------------------------------------------------------------------------
# Case B: No confounder. Plot Type I error of confounder test vs k.
# ---------------------------------------------------------------------------

# Create grouping variable based on (df1, df2) --------------------------------
results3 <- results1 %>%
  mutate(group = paste0("df1_", df1, "_df2_", df2))

summary_df <- results3 %>%
  group_by(group, k) %>%
  summarise(
    pct_1   = mean(conftest == 1),   # empirical Type I error (false positives)
    .groups = "drop"
  )

# Custom colors and legend labels for tail index combinations ----------------
cols_cust <- c("#1f78b4", "orange", "#B78BCF", "#a6cee3", "#fdbf6f", "#D4BDE4")
labels_cust <- c(
  expression(alpha[1] == 3 ~ "," ~ alpha[2] == 2),
  expression(alpha[1] == 4 ~ "," ~ alpha[2] == 2),
  expression(alpha[1] == 4 ~ "," ~ alpha[2] == 3)
)

p1 <- ggplot(summary_df,
             aes(x = k,
                 y = pct_1,
                 color = group)) +
  geom_line(linewidth = 2) +
  geom_point(size = 4) +
  scale_color_manual(values = cols_cust, labels = labels_cust) +
  labs(
    x     = "k",
    y     = "Type-I Error",
    color = "Tail Indices"
  ) +
  theme_minimal() +
  theme(text = element_text(size = 40))
p1

ggsave("pics/confounder_test_sim_caseB.pdf", plot = p1, device = "pdf",
       width = 14, height = 8, path = simulation_path)

# ---------------------------------------------------------------------------
# Case F: With confounder. Plot Type II error vs k.
# ---------------------------------------------------------------------------

results3 <- results2 %>%
  mutate(group = paste0("df1_", df1, "_df2_", df2))

summary_df <- results3 %>%
  group_by(group, k) %>%
  summarise(
    pct_1   = mean(conftest == 0),   # empirical Type II error (false negatives)
    .groups = "drop"
  )

cols_cust <- c("#1f78b4", "orange", "#B78BCF", "#a6cee3", "#fdbf6f", "#D4BDE4")

p2 <- ggplot(summary_df,
             aes(x = k,
                 y = pct_1,
                 color = group)) +
  geom_line(linewidth = 2) +
  geom_point(size = 4) +
  scale_color_manual(values = cols_cust, labels = labels_cust) +
  labs(
    x     = "k",
    y     = "Type-II Error",
    color = "Tail Indices"
  ) +
  theme_minimal() +
  theme(text = element_text(size = 40))
p2

ggsave("pics/confounder_test_sim_caseF.pdf", plot = p2, device = "pdf",
       width = 14, height = 8, path = simulation_path)
