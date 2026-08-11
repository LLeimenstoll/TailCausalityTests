# Example: Causal Tail Coefficient (CTC) simulation -------------------------

library(ggplot2)    # plotting
library(grid)       # low-level grid graphics (used for theme tweaks)

# Set working directory to GitHub repo root ---------------------------------
setwd()
simulation_path <- "Simulation_Study/"

# Load custom functions: CTC, tail estimators, etc. -------------------------
source("functions.R")

########### Simulations: convergence of CTC in a simple model ###############

# Simulation grid parameters ------------------------------------------------
n_vals  <- seq(3000, 10000, 1000)   # sample sizes
a_vals  <- c(2, 3, 4)               # df for X1 (tail index ~ 1/a)
b_vals  <- c(2, 3, 4)               # df for X2 (tail index ~ 1/b)
rho_vals <- c(0.1, 0.3, 0.5)        # linear dependence coefficient
k       <- 25                       # fixed k for CTC estimation

# Matrix to store mean CTCs over replications:
# rows: (c1, c2, delta) × all (n, a, b, rho) combinations
# columns: ctc, n, a, b, rho, dir (dir in {c1,c2,delta})
ctcs <- matrix(
  nrow = 3 * length(n_vals) * length(a_vals) * length(b_vals) * length(rho_vals),
  ncol = 6
)
j <- 1  # row pointer

# Nested loops over tail indices (a,b), dependence rho, and sample size n ---
for (b in b_vals) {
  for (a in a_vals) {
    for (rho in rho_vals) {
      for (n in n_vals) {
        
        c1s <- c()  # CTC X1 -> X2
        c2s <- c()  # CTC X2 -> X1
        ds  <- c()  # difference c1 - c2
        
        # Monte Carlo replications for given (n,a,b,rho) --------------------
        for (r in 1:50) {
          h <- 0                          # rt(n, 2) with confounder term 
          x <- rt(n, a) + 0.5 * h        
          y <- rt(n, b) + rho * x + 0.5 * h
          # y depends linearly on x plus possible confounder
          
          c1s[r] <- causal_tail_coeff_basic(x, y, both_tails = FALSE, k = k)
          c2s[r] <- causal_tail_coeff_basic(y, x, both_tails = FALSE, k = k)
          ds[r]  <- c1s[r] - c2s[r]
        }
        
        # Store means over replications for CTC in both directions and delta-
        ctcs[j,   1] <- mean(c1s)
        ctcs[j,   2] <- n
        ctcs[j,   3] <- a
        ctcs[j,   4] <- b
        ctcs[j,   5] <- rho
        ctcs[j,   6] <- "c1"
        
        ctcs[j+1, 1] <- mean(c2s)
        ctcs[j+1, 2] <- n
        ctcs[j+1, 3] <- a
        ctcs[j+1, 4] <- b
        ctcs[j+1, 5] <- rho
        ctcs[j+1, 6] <- "c2"
        
        ctcs[j+2, 1] <- mean(ds)
        ctcs[j+2, 2] <- n
        ctcs[j+2, 3] <- a
        ctcs[j+2, 4] <- b
        ctcs[j+2, 5] <- rho
        ctcs[j+2, 6] <- "delta"
        
        j <- j + 3
        print(j)
      }
    }
  }
}

ctcs

### Prepare data for plotting and saving ####################################

ctcs <- as.data.frame(ctcs)
colnames(ctcs) <- c("ctc", "n", "a", "b", "rho", "dir")

# Save simulation results (here for one setting) ------
saveRDS(ctcs, file = paste0(simulation_path, "output/deltas_t", ".rds"))

# Read precomputed results (here: t case) -------------------
ctcs <- readRDS(paste0(simulation_path, "output/deltas_t.rds"))

# Create combined factor for direction × rho to get 9 groups ----------------
ctcs$group <- interaction(ctcs$dir, ctcs$rho, sep = "_")

ctcs$group <- factor(
  ctcs$group,
  levels = c("c1_0.1", "c1_0.3", "c1_0.5",
             "c2_0.1", "c2_0.3", "c2_0.5",
             "delta_0.1", "delta_0.3", "delta_0.5")
)

# Convert columns to numeric where appropriate ------------------------------
ctcs$ctc <- as.numeric(ctcs$ctc)
ctcs$n   <- as.numeric(ctcs$n)
ctcs$a   <- as.numeric(ctcs$a)
ctcs$b   <- as.numeric(ctcs$b)

# Integer versions of a,b for facet labels ----------------------------------
ctcs$a1i <- round(ctcs$a)   # df for X1
ctcs$a2i <- round(ctcs$b)   # df for X2

# Colors for 9 groups (c1,c2,delta × 3 rho values) -------------------------
farben <- c("#1f78b4", "orange", "#B78BCF",
            "#a6cee3", "#fdbf6f", "#D4BDE4",
            "#1f78b4", "orange", "#B78BCF")

# Shapes: circles for c1,c2; crosses for delta ------------------------------
shapes <- c(21, 21, 21, 21, 21, 21, 4, 4, 4)

# Custom legend labels in math notation -------------------------------------
labels_cust <- c(
  expression(hat(Gamma)[X[1] %->% X[2]] ~ "," ~ beta["1" %->% "2"] == 0.1),
  expression(hat(Gamma)[X[1] %->% X[2]] ~ "," ~ beta["1" %->% "2"] == 0.3),
  expression(hat(Gamma)[X[1] %->% X[2]] ~ "," ~ beta["1" %->% "2"] == 0.5),
  expression(hat(Gamma)[X[2] %->% X[1]] ~ "," ~ beta["1" %->% "2"] == 0.1),
  expression(hat(Gamma)[X[2] %->% X[1]] ~ "," ~ beta["1" %->% "2"] == 0.3),
  expression(hat(Gamma)[X[2] %->% X[1]] ~ "," ~ beta["1" %->% "2"] == 0.5),
  expression(hat(Delta)[X[1] %->% X[2]] ~ "," ~ beta["1" %->% "2"] == 0.1),
  expression(hat(Delta)[X[1] %->% X[2]] ~ "," ~ beta["1" %->% "2"] == 0.3),
  expression(hat(Delta)[X[1] %->% X[2]] ~ "," ~ beta["1" %->% "2"] == 0.5)
)

pd <- position_dodge(0.00)

# Plot: CTC vs sample size, faceted by tail indices a (X1) and b (X2) -------
p1 <- ggplot(ctcs, aes(x = n, y = ctc, color = group)) +
  geom_line(alpha = .7, size = 0.8, position = pd) +
  geom_point(size = 2,
             position = pd, aes(color = group, fill = group, shape = group)) +
  facet_grid(a2i ~ a1i, labeller = label_both) +  # rows: df(X2), cols: df(X1)
  scale_shape_manual(values = shapes,
                     name   = expression("Causal Direction"),
                     labels = labels_cust) +
  scale_color_manual(name   = expression("Causal Direction"),
                     values = farben,
                     labels = labels_cust) +
  scale_fill_manual(name   = expression("Causal Direction"),
                    values = farben,
                    labels = labels_cust) +
  theme_minimal() +
  labs(
    x     = "Sample Size",
    y     = "Causal Tail Coefficient",
    color = "Group"
  ) +
  ylim(c(0, 1)) +
  theme(
    panel.spacing = unit(1, "lines"),           # space between facets
    panel.border  = element_rect(color = "black", fill = NA),
    text          = element_text(size = 20)
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1,   linetype = "dashed", color = "black")

p1

# Save plot as PNG ----------------------------------------------------------
ggsave("pics/simulations_conv_t_conf.png", plot = p1, device = "png",
       width = 14, height = 8, path = simulation_path)