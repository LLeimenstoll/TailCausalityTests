# Finance

# Load required libraries --------------------------------------------------
library(tidyr)        # data reshaping (pivot_longer, etc.)
library(ggplot2)      # plotting
library(gridExtra)    # arranging multiple ggplots
library(pcalg)        # causal discovery algorithms
library(quantmod)     # financial data (getSymbols, etc.)
library(causalXtreme) # causal tail coefficients, etc.
library( fGarch )     # fit ARMA-GARCH models

# Set directories -----------------------------------------------------------
setwd()
application_path <- "Application/"

# Source custom functions (user-defined) -----------------------------------
source("functions.R")   # contains Pareto/Hill, CTC functions, etc.

# Get data ------------------------------------------------------------------
options("getSymbols.warning4.0" = FALSE)
options("getSymbols.yahoo.warning" = FALSE)

startdate <- '2010-07-18'
enddate   <- "2024-12-04"

# Read Bitcoin price series (Coin Metrics) ---------------------------------
BTC <- read.csv(paste0(application_path, "data/coin-metrics.csv"), sep = ";")
BTC$Time <- as.Date(BTC$Time)
BTC <- xts(BTC[, 2], order.by = BTC[, 1])  # close price as xts
colnames(BTC) <- c("BTC.Close")

# S&P 500 and VIX from Yahoo Finance ---------------------------------------
getSymbols("^GSPC", from = startdate, to = enddate,
           warnings = FALSE, auto.assign = TRUE)
getSymbols("^VIX",  from = startdate, to = enddate,
           warnings = FALSE, auto.assign = TRUE)
VIX <- na.omit(VIX)

# MSCI Europe from local CSV -----------------------------------------------
msci_eu <- read.csv(paste0(application_path, "data/MSCI_eu.csv"), sep = ";")

# Compute log returns -------------------------------------------------------
# BTC: log-differences of prices
BTC$BTC.Return <- 0
for (i in 2:nrow(BTC)) {
  BTC$BTC.Return[i] <- log(as.numeric(BTC$BTC.Close[i]) /
                             as.numeric(BTC$BTC.Close[i - 1]))
}

# S&P 500: log-return from open to close
GSPC$GSPC.Return <- 0
for (i in 1:nrow(GSPC)) {
  GSPC$GSPC.Return[i] <- log(GSPC$GSPC.Close[i] / GSPC$GSPC.Open[i])
}

# VIX: log-return from open to close
VIX$VIX.Return <- 0
for (i in 1:nrow(VIX)) {
  VIX$VIX.Return[i] <- log(VIX$VIX.Close[i] / VIX$VIX.Open[i])
}

# MSCI Europe: log-differences of index level
msci_eu$return <- 0
for (i in 2:nrow(msci_eu)) {
  msci_eu$return[i] <- log(as.numeric(msci_eu$MSCI.Europe.Index[i]) /
                             as.numeric(msci_eu$MSCI.Europe.Index[i - 1]))
}

# Fit GARCH models and extract residuals -----------------------------------

fitA1 <- garchFit(data ~ arma(1, 1) + garch(1, 1), as.numeric(GSPC$GSPC.Return),
                  trace=F, cond.dist = "sstd") # S&P 500, based on best fit
                                                      
fitA2 <- garchFit(data ~ arma(2, 0) + garch(1, 1), as.numeric(BTC$BTC.Return),
                  trace = F, cond.dist = "sstd")  # BTC, based on best fit

GSPC_res <- residuals(fitA1, standardize = FALSE)
BTC_res  <- residuals(fitA2, standardize = FALSE)

# Attach residuals to datasets
GSPC$GSPC.res <- GSPC_res
BTC$BTC.res   <- BTC_res

# Merge S&P and BTC on common dates ----------------------------------------
data <- merge(GSPC, BTC, by = 0, all = FALSE)
df   <- as.data.frame(data)

# Scatter plot of residuals (S&P vs BTC) -----------------------------------
plot(as.numeric(df$GSPC.res), as.numeric(df$BTC.res),
     xlab = "S&P 500", ylab = "Bitcoin", title(main = "Returns"))

# QQ-plot comparing returns ----------------------------------
n  <- nrow(df)
xq <- quantile(df$GSPC.Return, probs = seq(0, 1, length.out = n))
yq <- quantile(df$BTC.Return,  probs = seq(0, 1, length.out = n))
qs <- cbind.data.frame(xq, yq)

qq_plot <- ggplot(qs, aes(x = xq, y = yq)) +
  geom_point(size = 4, col = "#1f78b4") +
  geom_abline(slope = 1, intercept = 0,
              color = "black", linetype = "dashed", linewidth = 1) +
  theme_bw() +
  theme(
    text = element_text(size = 35),
    legend.position = c(0.01, 0.98),
    legend.justification = c(0, 1)
  ) +
  xlab("S&P 500 Quantiles") +
  ylab("BTC Quantiles");qq_plot

ggsave("figures/finance_qq.pdf", plot = qq_plot, device = "pdf",
       width = 10, height = 10, path = application_path)

############# Equal tail index test (BTC vs S&P) ############################
# Data frame for collecting test results -----------------------------------
test_results <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(test_results) <- c("tail", "k_sp", "k_btc", "test_val", "sp_est", "btc_est")

# Split residuals into positive and negative parts (absolute for left tail) -
sp_pos  <- df$GSPC.res[df$GSPC.res > 0]
btc_pos <- df$BTC.res[df$BTC.res > 0]
sp_neg  <- -df$GSPC.res[df$GSPC.res < 0]
btc_neg <- -df$BTC.res[df$BTC.res < 0]

# Shifts for Pareto-type tail modeling -------------------------------------
# (positive tail: shift to upper endpoint; negative: likewise)
shift_sp_pos  <- max(sp_pos)
shift_sp_neg  <- max(sp_neg)
shift_btc_pos <- quantile(btc_pos, 0.95)
shift_btc_neg <- quantile(btc_neg, 0.95)

# Hill/Pareto estimators and plots -----------------------------------------
p1 <- pareto_hill(sp_pos,  shift = shift_sp_pos)
p2 <- pareto_hill(sp_neg,  shift = shift_sp_neg)
p3 <- pareto_hill(btc_pos, shift = shift_btc_pos)
p4 <- pareto_hill(btc_neg, shift = shift_btc_neg)
ggsave("figures/hill_plots_sp_neg.pdf",  plot = p2$plot, device = "pdf",
       width = 14, height = 10, path = application_path)
ggsave("figures/hill_plots_btc_pos.pdf", plot = p3$plot, device = "pdf",
       width = 14, height = 10, path = application_path)
ggsave("figures/hill_plots_btc_neg.pdf", plot = p4$plot, device = "pdf",
       width = 14, height = 10, path = application_path)

# Optimal k selection for Hill estimator ------------------------------------
k_sp_pos  <- find_optimal_k(sp_pos  + shift_sp_pos)$k_star
k_sp_neg  <- find_optimal_k(sp_neg  + shift_sp_neg)$k_star
k_btc_pos <- find_optimal_k(btc_pos + shift_btc_pos)$k_star
k_btc_neg <- find_optimal_k(btc_neg + shift_btc_neg)$k_star

# Two-sample tail index equality tests -------------------------------------
t1 <- tail_test(sp_pos  + shift_sp_pos,
                btc_pos + shift_btc_pos,
                min(k_sp_pos, k_btc_pos), 0.2)

t2 <- tail_test(sp_neg  + shift_sp_neg,
                btc_neg + shift_btc_neg,
                min(k_sp_neg, k_btc_neg), 0.2)

# Hoga tail index estimators at p = 1 --------------------------------------
h_sp_pos  <- hoga_estimator_vec(c(1), sp_pos  + shift_sp_pos,
                                min(k_sp_pos, k_btc_pos))
h_btc_pos <- hoga_estimator_vec(c(1), btc_pos + shift_btc_pos,
                                min(k_sp_pos, k_btc_pos))
h_sp_neg  <- hoga_estimator_vec(c(1), sp_neg  + shift_sp_neg,
                                min(k_sp_neg, k_btc_neg))
h_btc_neg <- hoga_estimator_vec(c(1), btc_neg + shift_btc_neg,
                                min(k_sp_neg, k_btc_neg))


# Collect test results in summary table ------------------------------------
test_results[1, 1] <- "right"         # right tail
test_results[1, 2] <- k_sp_pos
test_results[1, 3] <- k_btc_pos
test_results[1, 4] <- t1[1]           # test statistic
test_results[1, 5] <- h_sp_pos        # S&P tail index estimate
test_results[1, 6] <- h_btc_pos       # BTC tail index estimate

test_results[2, 1] <- "left"          # left tail
test_results[2, 2] <- k_sp_neg
test_results[2, 3] <- k_btc_neg
test_results[2, 4] <- t2[1]
test_results[2, 5] <- h_sp_neg
test_results[2, 6] <- h_btc_neg

test_results

############################# Causal Tail Coefficient (CTC) ################

# k according to n^0.4 ------------------------------------------------------
k1 <- floor(nrow(df)^0.4)

# CTC for upper tail (GSPC -> BTC and reverse) -----------------------------
ctc1a <- causal_tail_coeff_basic(as.numeric(df$GSPC.res),
                                 as.numeric(df$BTC.res),
                                 both_tails = FALSE, k = k1)
ctc1b <- causal_tail_coeff_basic(as.numeric(df$BTC.res),
                                 as.numeric(df$GSPC.res),
                                 both_tails = FALSE, k = k1)

# CTC for lower tail (min=TRUE; negative tail) ------------------------------
ctc2a <- causal_tail_coeff_basic(as.numeric(df$GSPC.res),
                                 as.numeric(df$BTC.res),
                                 both_tails = FALSE, k = k1, min = TRUE)
ctc2b <- causal_tail_coeff_basic(as.numeric(df$BTC.res),
                                 as.numeric(df$GSPC.res),
                                 both_tails = FALSE, k = k1, min = TRUE)

ctc1a
ctc1b
ctc1a - ctc1b  # directional difference

# Pre-Test and Confounder-Test for upper-tail causality --------------------------------
set.seed(123)
perm1 <- CTC_causality_permutation_test(as.numeric(df$GSPC.res),
                                        as.numeric(df$BTC.res),
                                        k = k1, R = 2000)
quantile(perm1$diff12p, 0.975)        # 97.5% quantile
tconf <- sqrt(k1) * (ctc1b - 0.5); tconf  # confounder-Test statistic

# Pre-Test and Confounder-Test for lower-tail causality --------------------------------
ctc2a
ctc2b
ctc2a - ctc2b

set.seed(123)
perm1 <- CTC_causality_permutation_test(-as.numeric(df$GSPC.res),
                                        -as.numeric(df$BTC.res),
                                        k = k1, R = 2000)
quantile(perm1$diff12p, 0.975)
tconf <- sqrt(k1) * (ctc2b - 0.5); tconf

# LiNGAM causal discovery on residuals --------------------------------------
mat <- df[, c("GSPC.res", "BTC.res")] %>%
  as.matrix()

pcalg::lingam(mat)                       
causal_discovery(mat, method = c("direct_lingam"))# direct LiNGAM model

####################  Convergence study (no confounder) #####################

min    <- TRUE                         # TRUE to consider lower tail, FALSE upper
n_vals <- seq(500, nrow(df), 100)      # subsample sizes
ks     <- floor(c(nrow(df)^0.4, nrow(df)^0.5))  # two k choices

# Matrix to store mean CTCs and deltas --------------------------------------
ctcs <- matrix(nrow = length(n_vals) * 3 * length(ks), ncol = 3)
i    <- 1

for (k1 in ks) {
  for (n in n_vals) {
    print(n)
    ctc_1 <- 0
    ctc_2 <- 0
    delta <- 0
    # Monte Carlo repetitions over random subsamples ------------------------
    for (r in 1:100) {
      s <- sample(nrow(df), n)
      x <- as.numeric(df[s, "GSPC.res"])
      y <- as.numeric(df[s, "BTC.res"])
      ctc_1[r] <- causal_tail_coeff_basic(x, y,
                                          k = k1, both_tails = FALSE,
                                          min = min)
      ctc_2[r] <- causal_tail_coeff_basic(y, x,
                                          k = k1, both_tails = FALSE,
                                          min = min)
      delta[r] <- ctc_1[r] - ctc_2[r]
    }
    # Average over repetitions ----------------------------------------------
    ctcs[i, 2]     <- mean(ctc_1)
    ctcs[i + 1, 2] <- mean(ctc_2)
    ctcs[i + 2, 2] <- mean(delta)
    
    ctcs[i, 1]     <- paste0("S&P500 \u2192 BTC (k=", k1, ")")
    ctcs[i + 1, 1] <- paste0("BTC \u2192 S&P500 (k=", k1, ")")
    ctcs[i + 2, 1] <- paste0("delta (k=", k1, ")")
    
    ctcs[i, 3]     <- paste0(n)
    ctcs[i + 1, 3] <- paste0(n)
    ctcs[i + 2, 3] <- paste0(n)
    
    i <- i + 3
  }
}

ctcs_1 <- ctcs

# Prepare data frame for plotting -------------------------------------------
ctcs <- data.frame(ctcs)
colnames(ctcs) <- c("group", "CTC", "n")
ctcs$n   <- as.numeric(ctcs$n)
ctcs$CTC <- as.numeric(ctcs$CTC)

ctcs$group <- factor(
  ctcs$group,
  levels = c(
    paste0("S&P500 \u2192 BTC (k=", ks[1], ")"),
    paste0("BTC \u2192 S&P500 (k=", ks[1], ")"),
    paste0("S&P500 \u2192 BTC (k=", ks[2], ")"),
    paste0("BTC \u2192 S&P500 (k=", ks[2], ")"),
    paste0("delta (k=", ks[1], ")"),
    paste0("delta (k=", ks[2], ")")
  )
)

pd    <- position_dodge(0.00)
cols  <- c("#1f78b4", "#a6cee3", "orange", "#fdbf6f",
           "#1f78b4", "orange")
shapes <- c(21, 21, 21, 21, 4, 4)

# Custom legend labels with math notation ----------------------------------
labels_cust <- c(
  bquote(hat(Gamma)["S&P 500" %->% "BTC"] * "," ~ k == .(ks[1])),
  bquote(hat(Gamma)["BTC" %->% "S&P 500"] * "," ~ k == .(ks[1])),
  bquote(hat(Gamma)["S&P 500" %->% "BTC"] * "," ~ k == .(ks[2])),
  bquote(hat(Gamma)["BTC" %->% "S&P500"] * "," ~ k == .(ks[2])),
  bquote(hat(Delta)["S&P 500" %->% "BTC"] * "," ~ k == .(ks[1])),
  bquote(hat(Delta)["S&P 500" %->% "BTC"] * "," ~ k == .(ks[2]))
)

# Convergence plot for CTC (no confounder) ---------------------------------
conv_plot <- ggplot(ctcs, aes(n, CTC, color = group)) +
  geom_line(alpha = .7, size = 1.3, position = pd) +
  geom_point(size = 3.5,
             position = pd, aes(color = group, fill = group, shape = group)) +
  theme_bw() +
  scale_shape_manual(values = shapes,
                     name   = expression("Causal Direction"),
                     labels = labels_cust,
                     guide  = guide_legend(ncol = 3)) +
  scale_color_manual(name   = expression("Causal Direction"),
                     values = cols,
                     labels = labels_cust,
                     guide  = guide_legend(ncol = 3)) +
  scale_fill_manual(name   = expression("Causal Direction"),
                    values = cols,
                    labels = labels_cust,
                    guide  = guide_legend(ncol = 3)) +
  labs(color = 'Causal Direction') +
  xlab("Sample Size") +
  ylab("Causal Tail Coefficient") +
  ylim(c(0, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1,   linetype = "dashed", color = "black") +
  theme(text = element_text(size = 26)) +
  theme(legend.position = c(0.48, 0.86))

conv_plot

ggsave("figures/finance_conv_plot_left.pdf", plot = conv_plot, device = "pdf",
       width = 14, height = 10, path = application_path)

################################################################################
######## Include confounders: VIX and MSCI Europe #############################

# Fit GARCH for MSCI EU and VIX --------------------------------------------

fitA3 <- garchFit(data ~ arma(1, 2) + garch(1, 1), as.numeric(msci_eu$return),
                  trace=F, cond.dist = "sstd") # MSCI, based on best fit

fitA4 <- garchFit(data ~ arma(0, 2) + garch(1, 1), as.numeric(VIX$VIX.Return),
                  trace=F, cond.dist = "sstd") # VIX, based on best fit

msci_eu_res <- residuals(fitA3, standardize = FALSE)
VIX_res     <- residuals(fitA4, standardize = FALSE)

msci_eu$IEUR.res <- msci_eu_res
# Align VIX residuals with price series (shifted by one) ----
VIX$VIX.res      <- c(0, VIX_res[-length(VIX_res)])

# Prepare merged data sets with confounders ---------------------------------
rownames(msci_eu) <- msci_eu$Date
data_ieur <- merge(df, msci_eu, by = 0, all = FALSE)
data_vix  <- merge(df, VIX,     by = 0, all = FALSE)

# Conditional CTC with VIX as confounder (upper tail)------------------------
k1   <- floor(nrow(data_vix)^0.4)
ctc3a <- LGPD_causal_tail_coeff(data_vix$GSPC.res,
                                data_vix$BTC.res,
                                H = data_vix$VIX.res, k = k1)
ctc3b <- LGPD_causal_tail_coeff(data_vix$BTC.res,
                                data_vix$GSPC.res,
                                H = data_vix$VIX.res, k = k1)
# Conditional CTC with VIX (lower tail)--------------------------------------
ctc4a <- LGPD_causal_tail_coeff(-data_vix$GSPC.res,
                                -data_vix$BTC.res,
                                H = -data_vix$VIX.res, k = k1)
ctc4b <- LGPD_causal_tail_coeff(-data_vix$BTC.res,
                                -data_vix$GSPC.res,
                                H = -data_vix$VIX.res, k = k1)

# Conditional CTC with MSCI Europe as confounder (upper tail)---------------------------
k2   <- floor(nrow(data_ieur)^0.4)
ctc5a <- LGPD_causal_tail_coeff(data_ieur$GSPC.res,
                                data_ieur$BTC.res,
                                H = data_ieur$IEUR.res, k = k2)
ctc5b <- LGPD_causal_tail_coeff(data_ieur$BTC.res,
                                data_ieur$GSPC.res,
                                H = data_ieur$IEUR.res, k = k2)

# Conditional CTC with MSCI Europe (lower tail) -----------------------
ctc6a <- LGPD_causal_tail_coeff(-data_ieur$GSPC.res,
                                -data_ieur$BTC.res,
                                H = -data_ieur$IEUR.res,
                                k = k2, threshold_q = 0.94)
ctc6b <- LGPD_causal_tail_coeff(-data_ieur$BTC.res,
                                -data_ieur$GSPC.res,
                                H = -data_ieur$IEUR.res,
                                k = k2, threshold_q = 0.94)

# Pre-Test and Confounder-Test with confounders ---------------------------------------
ctc4a
ctc4b
ctc4a - ctc4b

set.seed(123)
perm1 <- CTC_causality_permutation_test(-as.numeric(df$GSPC.res),
                                        -as.numeric(df$BTC.res),
                                        H = -data_vix$VIX.res,
                                        k = k1, R = 2000)
quantile(perm1$diff12p, 0.975)
tconf <- sqrt(k1) * (ctc4b - 0.5); tconf

ctc6a
ctc6b
ctc6a - ctc6b

set.seed(123)
perm1 <- CTC_causality_permutation_test(-as.numeric(df$GSPC.res),
                                        -as.numeric(df$BTC.res),
                                        H = -data_ieur$IEUR.res,
                                        k = k2, threshold_q = 0.94, R = 2000)
quantile(perm1$diff12p, 0.975)
tconf <- sqrt(k2) * (ctc6b - 0.5); tconf

########### Delta (permutation) plots #######################################

# Custom legend labels for conditional CTC with MSCI EU -------------------
labels_cust <- c(
  bquote(hat(Gamma)["S&P 500" %->% "BTC|MSCI EU"] * "," ~ k == .(ks[1])),
  bquote(hat(Gamma)["BTC" %->% "S&P 500|MSCI EU"] * "," ~ k == .(ks[1])),
  bquote(hat(Gamma)["S&P 500" %->% "BTC|MSCI EU"] * "," ~ k == .(ks[2])),
  bquote(hat(Gamma)["BTC" %->% "S&P500|MSCI EU"] * "," ~ k == .(ks[2])),
  bquote(hat(Delta)["S&P 500" %->% "BTC|MSCI EU"] * "," ~ k == .(ks[1])),
  bquote(hat(Delta)["S&P 500" %->% "BTC|MSCI EU"] * "," ~ k == .(ks[2]))
)

# Example: permutation with MSCI EU (commented variant retained) -----------
k1 <- floor(nrow(data_ieur)^0.4)
k2 <- floor(nrow(data_ieur)^0.5)
h  <- -as.numeric(data_ieur[, "IEUR.res"])

set.seed(123)
perm_res  <- CTC_causality_permutation_test(-data_ieur$GSPC.res,
                                            -data_ieur$BTC.res,
                                            k = k1, R = 2000,
                                            H = h, threshold_q = 0.94)
set.seed(123)
perm_res2 <- CTC_causality_permutation_test(-data_ieur$GSPC.res,
                                            -data_ieur$BTC.res,
                                            k = k2, R = 2000,
                                            H = h, threshold_q = 0.94)
# Here: unconditional permutation for right tail ---------------------------
k1 <- floor(nrow(df)^0.4)
k2 <- floor(nrow(df)^0.5)

set.seed(123)
perm_res  <- CTC_causality_permutation_test(df$GSPC.res,
                                            df$BTC.res,
                                            k = k1, R = 2000)
set.seed(123)
perm_res2 <- CTC_causality_permutation_test(df$GSPC.res,
                                            df$BTC.res,
                                            k = k2, R = 2000)

perm <- cbind.data.frame(perm_res$diff12p, perm_res2$diff12p)
colnames(perm) <- c("k1", "k2")

perm_long <- perm |>
  pivot_longer(
    cols      = everything(),
    names_to  = "group",
    values_to = "diff"
  )

# Histogram for k1 ---------------------------------------------------------
bs1 <- ggplot(perm_long, aes(x = diff, fill = group)) +
  geom_histogram(
    data = subset(perm_long, group == "k1"),
    aes(y = after_stat(density), fill = "k1"),
    binwidth = 0.02, color = "black", alpha = 1
  ) +
  geom_vline(xintercept = perm_res$diff12, color = "black",
             linetype = "dashed", linewidth = 1) +
  annotate("text",
           x = perm_res$diff12, y = 10,
           label = expression(hat(Delta)["   S&P 500" %->% "BTC"]),
           vjust = 1.2, hjust = 1.05, color = "black", size = 10) +
  scale_fill_manual(
    name   = NULL,
    values = c("k1" = "#1f78b4"),
    labels = c("k1" = "k=26")
  ) +
  labs(x = expression("Bootstrapped " ~ hat(Delta)), y = "Density") +
  theme_bw() +
  theme(
    text = element_text(size = 26),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  scale_x_continuous(limits = c(-0.3, 0.3)) +
  scale_y_continuous(limits = c(0, 10))

bs1

# Histogram for k2 ---------------------------------------------------------
bs2 <- ggplot(perm_long, aes(x = diff, fill = group)) +
  geom_histogram(
    data = subset(perm_long, group == "k2"),
    aes(y = after_stat(density), fill = "k2"),
    binwidth = 0.02, color = "black", alpha = 1
  ) +
  geom_vline(xintercept = perm_res2$diff12, color = "black",
             linetype = "dashed", linewidth = 1) +
  annotate("text",
           x = perm_res2$diff12, y = 10,
           label = expression(hat(Delta)["   S&P 500" %->% "BTC"]),
           vjust = 1.2, hjust = 1.05, color = "black", size = 10) +
  scale_fill_manual(
    name   = NULL,
    values = c("k2" = "orange"),
    labels = c("k2" = "k=60")
  ) +
  labs(x = expression("Bootstrapped " ~ hat(Delta)), y = "Density") +
  theme_bw() +
  theme(
    text = element_text(size = 26),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  scale_x_continuous(limits = c(-0.3, 0.3)) +
  scale_y_continuous(limits = c(0, 10))

bs2

bs_plot <- grid.arrange(bs1, bs2, nrow = 2, ncol = 1); bs_plot

ggsave("figures/finance_conf_delta_plot_right.pdf", plot = bs_plot, device = "pdf",
       width = 14, height = 10, path = application_path)

####################  Convergence with confounder ###########################
n_vals<-seq(1500,nrow(data_ieur),100)
ctcs <- matrix(nrow = length(n_vals) * 3 * length(ks), ncol = 3)
i    <- 1

for (k1 in ks) {
  for (n in n_vals) {
    print(n)
    ctc_1 <- 0
    ctc_2 <- 0
    for (r in 1:100) {
      s <- sample(nrow(data_ieur), n)
      x <- -as.numeric(data_ieur[s, "GSPC.res"])
      y <- -as.numeric(data_ieur[s, "BTC.res"])
      h <- -as.numeric(data_ieur[s, "IEUR.res"])
      ctc_1[r] <- LGPD_causal_tail_coeff(x, y, H = h, k = k1)
      ctc_2[r] <- LGPD_causal_tail_coeff(y, x, H = h, k = k1)
      delta[r] <- ctc_1[r] - ctc_2[r]
    }
    ctcs[i, 2]     <- mean(ctc_1)
    ctcs[i + 1, 2] <- mean(ctc_2)
    ctcs[i + 2, 2] <- mean(delta)
    
    ctcs[i, 1]     <- paste0("S&P500 \u2192 BC (k=", k1, ")")
    ctcs[i + 1, 1] <- paste0("BC \u2192 S&P500 (k=", k1, ")")
    ctcs[i + 2, 1] <- paste0("delta (k=", k1, ")")
    
    ctcs[i, 3]     <- paste0(n)
    ctcs[i + 1, 3] <- paste0(n)
    ctcs[i + 2, 3] <- paste0(n)
    
    i <- i + 3
  }
}

ctcs <- data.frame(ctcs)
colnames(ctcs) <- c("group", "CTC", "n")
ctcs$n   <- as.numeric(ctcs$n)
ctcs$CTC <- as.numeric(ctcs$CTC)

ctcs$group <- factor(
  ctcs$group,
  levels = c(
    "S&P500 \u2192 BC (k=26)", "BC \u2192 S&P500 (k=26)",
    "S&P500 \u2192 BC (k=60)", "BC \u2192 S&P500 (k=60)",
    "delta (k=26)","delta (k=60)")
)


# Reuse pd, shapes, cols, labels_cust from before --------------------------
conv_plot <- ggplot(ctcs, aes(n, CTC, color = group)) +
  geom_line(alpha = .7, size = 1.3, position = pd) +
  geom_point(size = 3.5,
             position = pd, aes(color = group, fill = group, shape = group)) +
  theme_bw() +
  scale_shape_manual(values = shapes,
                     name   = expression("Causal Direction"),
                     labels = labels_cust, guide = guide_legend(ncol = 3)) +
  scale_color_manual(name   = expression("Causal Direction"),
                     values = cols,
                     labels = labels_cust,
                     guide  = guide_legend(ncol = 3)) +
  scale_fill_manual(name   = expression("Causal Direction"),
                    values = cols,
                    labels = labels_cust,
                    guide  = guide_legend(ncol = 3)) +
  labs(color = 'Causal Direction') +
  xlab("Sample Size") +
  ylab("Causal Tail Coefficient") +
  ylim(c(0, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1,   linetype = "dashed", color = "black") +
  theme(text = element_text(size = 26)) +
  theme(legend.position = c(0.48, 0.86))

conv_plot

ggsave("figures/finance_conf_conv_plot_left.pdf", plot = conv_plot, device = "pdf",

       width = 14, height = 10, path = application_path)
