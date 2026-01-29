# alpha-level behavior of pretests: simulation study ------------------------

library(dplyr)   # data manipulation (group_by, summarise)
library(ggplot2) # plotting
library(tidyr)   # pivot_longer

# Set working directory to GitHub repo root ---------------------------------
setwd()

simulation_path <- "Simulation_Study/"

# Load custom functions: CTC, tail estimator, etc. --------------------------
source("functions.R")

# Simulation setup ----------------------------------------------------------
n      <- 10000                    # sample size per simulation
k      <- floor(n^0.4)             # tail sample size for CTC
rounds <- 1000                     # number of Monte Carlo replications
alphas <- c(0.01, 0.025, 0.05, 0.1) # nominal significance levels to test

results <- data.frame()  # will hold df1, df2, k, pretest1, pretest2, pretest3, alpha
r <- 1                   # row index in results

# Main simulation loop ------------------------------------------------------
for (i in 1:rounds) {
  print(i)
  
  # Degrees of freedom for heavy-tailed errors (random per replication) ----
  df1 <- sample(c(2, 3, 4), 1)
  df2 <- sample(c(2, 3, 4), 1)
  
  e1 <- rt(n, df1)      # error series 1
  e2 <- rt(n, df2)      # error series 2
  
  # Random linear dependence X2 = e2 + beta * X1 ----------------------------
  beta <- runif(1, 0.1, 0.9)   # strength of dependence, for case 1 set to 0
  X1   <- e1
  X2   <- e2 + beta * X1

  # Nonparametric causal tail coefficients in both directions ---------------
  c12 <- causal_tail_coeff_basic(X1, X2, k = k, both_tails = FALSE)
  c21 <- causal_tail_coeff_basic(X2, X1, k = k, both_tails = FALSE)
  
  # Determine "direction" by larger CTC: c12 > c21 ⇒ X1 → X2 ----------------
  if (c12 > c21) {
    dir <- 1
  } else {
    dir <- 2
    tmp <- X1
    X1  <- X2
    X2  <- tmp
  }
  
  # Permutation-based CTC Pre-Test -------------------------------------------
  p_val <- CTC_causality_permutation_test(X1, X2, R = 2000, k = k)$Pmc
  
  # Indicator variables for Pre-Tests at different alpha values --------------
  pretest1 <- 0  # permutation test (reject / not)
  pretest2 <- 0  # additional test based on normal approximation (not needed)
  pretest3 <- 0  # combined test (not needed)
  
  # Loop over desired nominal alpha levels ----------------------------------
  for (alpha in alphas) {
    
    # Pretest 1: permutation p-value < alpha? -------------------------------
    if (p_val < alpha) {
      pretest1 <- 1
    }
    
    # Pretest 2: normal approximation of c12, c21 vs 0.5 --------------------
    tc1 <- sqrt(k) * (c12 - 0.5)
    tc2 <- sqrt(k) * (c21 - 0.5)
    # Under null, Var(ctc) = 1/12 
    q1 <- qnorm(1 - alpha / 2, 0, sqrt(1 / 12))
    q2 <- qnorm(alpha / 2,      0, sqrt(1 / 12))
    
    pretest2 <- (tc1 < q2 | tc1 > q1) | (tc2 < q2 | tc2 > q1)
    
    # Pretest 3: intersection of pretest1 and pretest2 (unused) -------------
    pretest3 <- pretest1 * pretest2
    
    # Store one row per (simulation i, alpha) --------------------------------
    results[r, 1:7] <- c(df1, df2, k, pretest1, pretest2, pretest3, alpha)
    r <- r + 1
  }
}

results

# Name columns --------------------------------------------------------------
colnames(results) <- c("df1", "df2", "k", "pretest1", "pretest2", "pretest3", "alpha")

# Save results (case 3) -----------------------------------------------------
simulation_path <- "C:/Users/be6427/Documents/Paper_2_Markets and weather extremes/Paper_final/Simulation_Study/"
write.csv(results,
          file = paste0("output/results_pretest_alphas_case3.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Plot Type I error behavior (case 1 results)
# ---------------------------------------------------------------------------

# Read results from case 1 (null case) --------------------------------------
results <- read.csv(paste0(simulation_path, "output/results_pretest_alphas_case1.csv"))

# Convert pretest columns to long format ------------------------------------
long_df <- results %>%
  pivot_longer(
    cols      = starts_with("pretest"),
    names_to  = "test",
    values_to = "value"
  )

# Compute empirical rejection probability per alpha and pretest ------------
summary_df <- long_df %>%
  group_by(alpha, test) %>%
  summarise(
    pct_0   = mean(value == 1),   # fraction of rejections (Type I error)
    .groups = "drop"
  )

# Plot Type I error vs nominal alpha ----------------------------------------
ggplot(summary_df,
       aes(x = alpha,
           y = pct_0,
           color = test)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    x     = "alpha",
    y     = "Type I Error",
    color = "Test"
    # title = "Rejection rate per pretest and alpha (case 1)"
  ) +
  theme_minimal(base_size = 14)

# ---------------------------------------------------------------------------
# Plot Type II error behavior (case 3 results)
# ---------------------------------------------------------------------------

results <- read.csv(paste0(simulation_path, "output/results_pretest_alphas_case3.csv"))

long_df <- results %>%
  pivot_longer(
    cols      = starts_with("pretest"),
    names_to  = "test",
    values_to = "value"
  )

# Here value == 0 = no rejection (Type II error in non-null case) ----------
summary_df <- long_df %>%
  group_by(alpha, test) %>%
  summarise(
    pct_0   = mean(value == 0),   # fraction of failures to reject
    .groups = "drop"
  )

ggplot(summary_df,
       aes(x = alpha,
           y = pct_0,
           color = test)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    x     = "alpha",
    y     = "Type II Error",
    color = "Test"
    # title = "Non-rejection rate per pretest and alpha (case 3)"
  ) +
  theme_minimal(base_size = 14)
