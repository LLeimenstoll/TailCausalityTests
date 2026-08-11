# Simulation study to compare Pre-Test with Lingam

library(pcalg)    # for lingam() causal discovery
library(ggplot2)  # for plotting

# Set working directory to GitHub repo root ---------------------------------
setwd()
simulation_path <- "Simulation_Study/"

# Load custom functions: CTC, tail estimators, permutation tests, etc. ------
source("functions.R")

# Simulation setup ----------------------------------------------------------
n_vals   <- seq(3000, 10000, 1000)     # sample sizes
a_vals   <- c(2, 3, 4)                 # df for X1 (controls tail heaviness)
b_vals   <- c(2, 3, 4)                 # df for X2
rho_vals <- c(0, 0.1, 0.3, 0.5)        # strength of X1 -> X2 dependence
k        <- 25                         # tail sample size for CTC
rounds   <- 100                        # Monte Carlo replications

# Store results: two rows per (n,a,b,rho) for directions c1, c2 -------------
ctcs_1 <- matrix(
  nrow = 2 * length(n_vals) * length(a_vals) * length(b_vals) * length(rho_vals),
  ncol = 8
)
j <- 1  # row pointer

# Nested loops over tail indices, dependence, and sample size ---------------
for (b in b_vals) {
  for (a in a_vals) {
    for (rho in rho_vals) {
      for (n in n_vals) {
        
        c1s <- c()   # CTC X1 -> X2
        c2s <- c()   # CTC X2 -> X1
        lin1 <- c()  # LiNGAM detects X1 -> X2
        lin2 <- c()  # LiNGAM detects X2 -> X1
        num_caus <- 0  # number of times CTC pretest detects causality
        
        # Threshold for "extreme" X1 (used in construction of Y) ------------
        q <- qt(0.95, a)
        
        # Monte Carlo for given (a,b,rho,n) ---------------------------------
        for (r in 1:rounds) {
          # Heavy-tailed X1
          x <- rt(n, a)
          # Nonlinear, partly tail-driven relation: only if x > q,
          # Y contains rho * x, otherwise just noise
          y <- ifelse(x > q, 1, 0) * rho * x + rt(n, b)
          
          mat <- cbind(x, y)
          
          # Nonparametric CTC in both directions ----------------------------
          c1s[r] <- causal_tail_coeff_basic(x, y, both_tails = FALSE, k = k)
          c2s[r] <- causal_tail_coeff_basic(y, x, both_tails = FALSE, k = k)
          
          # LiNGAM on (X1, X2): check inferred direction --------------------
          L <- pcalg::lingam(mat)
          lin1[r] <- as.numeric(isTRUE(L$Bpruned[2, 1] > 0))  # 1 if X1 -> X2
          lin2[r] <- as.numeric(isTRUE(L$Bpruned[1, 2] > 0))  # 1 if X2 -> X1
          
          # Permutation CTC pretest: does CTC detect asymmetry? -------------
          test <- CTC_causality_permutation_test(x, y, k = k, R = 200)
          if (test$Pmc < 0.05) {
            num_caus <- num_caus + 1
          }
        }
        
        # Store averages across repetitions for X1 -> X2 (c1) ---------------
        ctcs_1[j, 1] <- mean(c1s)              # mean CTC(X1 -> X2)
        ctcs_1[j, 2] <- n
        ctcs_1[j, 3] <- a
        ctcs_1[j, 4] <- b
        ctcs_1[j, 5] <- rho
        ctcs_1[j, 6] <- "c1"                   # direction label
        ctcs_1[j, 7] <- mean(lin1)             # LiNGAM correctness (X1 -> X2)
        ctcs_1[j, 8] <- num_caus / rounds      # fraction of CTC pretest rejections
        
        # Store averages for X2 -> X1 (c2) ----------------------------------
        ctcs_1[j + 1, 1] <- mean(c2s)          # mean CTC(X2 -> X1)
        ctcs_1[j + 1, 2] <- n
        ctcs_1[j + 1, 3] <- a
        ctcs_1[j + 1, 4] <- b
        ctcs_1[j + 1, 5] <- rho
        ctcs_1[j + 1, 6] <- "c2"
        ctcs_1[j + 1, 7] <- mean(lin2)         # LiNGAM correctness (X2 -> X1)
        ctcs_1[j + 1, 8] <- 1 - num_caus / rounds  # fraction of "no causality" in this direction
        
        j <- j + 2
        print(j)
      }
    }
  }
}

ctcs_1 <- as.data.frame(ctcs_1)
ctcs   <- ctcs_1

# Name columns: CTC estimate, sample size, df(a,b), rho, direction, LiNGAM, pretest
colnames(ctcs_1) <- c("ctc", "n", "a", "b", "rho", "dir", "linDir", "test")

# Save & reload results -----------------------------------------------------
write.csv(ctcs_1, paste0(simulation_path, "output/results_lingam_tests.csv"))
ctcs_1 <- read.csv(paste0(simulation_path, "output/results_lingam_tests.csv"))

# ---------------------------------------------------------------------------
# Prepare data for plotting
# ---------------------------------------------------------------------------

# Combined grouping for direction × rho, e.g. "c1_0.1"
ctcs_1$group <- interaction(ctcs_1$dir, ctcs_1$rho, sep = "_")

ctcs_1$group <- factor(
  ctcs_1$group,
  levels = c("c1_0.1", "c1_0.3", "c1_0.5", "c1_0",
             "c2_0.1", "c2_0.3", "c2_0.5", "c2_0")
)

# Ensure numeric types ------------------------------------------------------
ctcs_1$ctc    <- as.numeric(ctcs_1$ctc)
ctcs_1$test   <- as.numeric(ctcs_1$test)
ctcs_1$n      <- as.numeric(ctcs_1$n)
ctcs_1$a      <- as.numeric(ctcs_1$a)
ctcs_1$b      <- as.numeric(ctcs_1$b)
ctcs_1$linDir <- as.numeric(ctcs_1$linDir)

# Integer versions of df for facet labels -----------------------------------
ctcs_1$a1i <- round(ctcs_1$a)  # df for X1
ctcs_1$a2i <- round(ctcs_1$b)  # df for X2

# For plotting, focus on direction c1 (X1 -> X2) ----------------------------
ctcs_1 <- ctcs_1[ctcs_1$dir == "c1", ]

# ---------------------------------------------------------------------------
# Plot 1: power of CTC pretest (fraction of causal discoveries) ------------
# ---------------------------------------------------------------------------

cols_cust <- c("#1f78b4", "orange", "#B78BCF", "black")
pd <- position_dodge(0.00)

p1 <- ggplot(ctcs_1, aes(x = n, y = test, color = group)) +
  geom_line(alpha = .7, size = 0.8, position = pd) +
  geom_point(size = 2, position = pd, aes(fill = group)) +
  facet_grid(a2i ~ a1i, labeller = label_both) +  # rows: b, cols: a
  scale_color_manual(
    name   = expression("Causal Direction"),
    values = cols_cust,
    labels = c(
      expression("X"["1"] %->% "X"["2"] ~ "," ~ beta["1" %->% "2"] == 0.1),
      expression("X"["1"] %->% "X"["2"] ~ "," ~ beta["1" %->% "2"] == 0.3),
      expression("X"["1"] %->% "X"["2"] ~ "," ~ beta["1" %->% "2"] == 0.5),
      expression("X"["1"] %->% "X"["2"] ~ "," ~ beta["1" %->% "2"] == 0)
    )
  ) +
  theme_minimal() +
  labs(
    x     = "Sample Size",
    y     = "Percentage of Causal Discoveries",
    color = "Group (c, rho)"
  ) +
  ylim(c(0, 1)) +
  theme(
    panel.spacing = unit(1, "lines"),              # space between facets
    panel.border  = element_rect(color = "black", fill = NA),
    text          = element_text(size = 20)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  guides(fill = "none")

p1

# ---------------------------------------------------------------------------
# Plot 2: LiNGAM power (fraction of correct X1 -> X2 detections) ------------
# ---------------------------------------------------------------------------

p2 <- ggplot(ctcs_1, aes(x = n, y = linDir, color = group)) +
  geom_line(alpha = .7, size = 0.8, position = pd) +
  geom_point(size = 2, position = pd, aes(fill = group)) +
  facet_grid(a2i ~ a1i, labeller = label_both) +
  scale_color_manual(
    name   = expression("Causal Direction"),
    values = cols_cust,
    labels = c(
      expression("X"["1"] %->% "X"["2"] ~ "," ~ beta["1" %->% "2"] == 0.1),
      expression("X"["1"] %->% "X"["2"] ~ "," ~ beta["1" %->% "2"] == 0.3),
      expression("X"["1"] %->% "X"["2"] ~ "," ~ beta["1" %->% "2"] == 0.5),
      expression("X"["1"] %->% "X"["2"] ~ "," ~ beta["1" %->% "2"] == 0)
    )
  ) +
  theme_minimal() +
  labs(
    x     = "Sample Size",
    y     = "Percentage of Causal Discoveries",
    color = "Group (c, rho)"
  ) +
  ylim(c(0, 1)) +
  theme(
    panel.spacing = unit(1, "lines"),
    panel.border  = element_rect(color = "black", fill = NA),
    text          = element_text(size = 20)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  guides(fill = "none")

p2

# Save plots ----------------------------------------------------------------
ggsave("pics/pretest_pic.png", plot = p1, device = "png",
       width = 14, height = 8, path = simulation_path)
ggsave("pics/lingam_pic.png", plot = p2, device = "png",
       width = 14, height = 8, path = simulation_path)