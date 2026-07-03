# Match Swiss train delay data with weather (precipitation) data ------------

# Load required packages ----------------------------------------------------
library(dplyr)      # data manipulation (group_by, summarise, etc.)
library(tidyr)      # data reshaping (used later with pivot_longer)
library(ggplot2)    # plotting
library(pcalg)      # causal discovery (LiNGAM, etc.)
library(causalXtreme)
# Set working directory to GitHub repo root ---------------------------------
setwd()

# Define application path (relative to working directory) -------------------
application_path <- "Application/"

# Load custom functions (tail causality, Hill, etc.) ------------------------
source("functions.R")

# Read combined train and weather data --------------------------------------
data_all <- read.csv(paste0(application_path, "data/data_combined_train_weather.csv"))

data_all_2 <- data_all

# Extract date and calendar information from 'key' --------------------------
data_all_2$date    <- as.Date(substr(data_all_2$key, 1, 10), format = "%d.%m.%Y")
data_all_2$weekday <- weekdays(data_all_2$date)
data_all_2$month   <- substr(data_all_2$key, 4, 5)

head(data_all_2)

################ Analyse Data ###############################################

# Simple scatter plots: rain vs delay at different aggregation levels -------
plot(data_all_2$rain,    data_all_2$delay)
plot(data_all_2$rain_6,  data_all_2$delay)
plot(data_all_2$rain_12, data_all_2$delay)

# Mean and quantile of delay by month and weekday ---------------------------
ggplot(data_all_2, aes(x = month,   y = delay)) +
  stat_summary(fun = "mean", geom = "bar", fill = "#1874CD", col = "black")

ggplot(data_all_2, aes(x = month,   y = delay)) +
  stat_summary(fun = "quantile", geom = "bar", fill = "#1874CD", col = "black")

ggplot(data_all_2, aes(x = weekday, y = delay)) +
  stat_summary(fun = "mean", geom = "bar", fill = "#1874CD", col = "black")

ggplot(data_all_2, aes(x = weekday, y = delay)) +
  stat_summary(fun = "quantile", geom = "bar", fill = "#1874CD", col = "black")

# Mean daily precipitation by month (12h window) ----------------------------
ggplot(data_all_2, aes(x = month, y = rain_12)) +
  stat_summary(fun = "mean", geom = "bar", fill = "#1874CD", col = "black")

# Inspect “special days” with unusually high delays -------------------------
quantile(data_all_2$delay, 0.95)          # 95th percentile of delays
high_delays <- data_all_2[data_all_2$delay > 4, ]

ggplot(high_delays) +
  geom_bar(aes(date))                     # count of high-delay days by date

high_delays$d1 <- as.factor(high_delays$date)
summary(high_delays$d1)

# Visualize the frequency of observations by month (YYYY-MM) ----------------
ggplot(data_all_2) +
  geom_bar(aes(substring(date, 1, 7)))

high_delays$d1 <- as.factor(high_delays$date)

########### Data Cleaning ###################################################

# Swiss public holidays (to be excluded) ------------------------------------
public_holidays <- c(
  "2021-01-01","2021-04-02","2021-04-05","2021-05-01","2021-05-13",
  "2021-05-24","2021-08-01","2021-12-25","2021-12-26",
  "2022-01-01","2022-04-15","2022-04-18","2022-05-01","2022-05-26",
  "2022-06-06","2022-08-01","2022-12-25","2022-12-26",
  "2023-01-01","2023-04-07","2023-04-10","2023-05-01","2023-05-18",
  "2023-05-29","2023-08-01","2023-12-25","2023-12-26",
  "2024-01-01","2024-03-29","2024-04-01","2024-05-01","2024-05-09",
  "2024-05-20"
)
public_holidays <- as.Date(public_holidays, format = "%Y-%m-%d")

days_to_exclude <- c(public_holidays)

# Remove public holidays ----------------------------------------------------
data_all_3 <- subset(data_all_2, !(data_all_2$date %in% days_to_exclude))

# Exclude Saturdays and Sundays (keep Monday–Friday only) -------------------
data_all_3 <- subset(data_all_3,
                     !(data_all_3$weekday %in% c("Samstag", "Sonntag")))

# Keep only one travel direction, e.g. Zurich → Bern ("ZB") -----------------
data_all_3 <- subset(data_all_3, data_all_3$direction == "ZB")

summary(data_all_3)

# Remove rows with any missing values ---------------------------------------
data_all_3 <- na.omit(data_all_3)

############# Plot Histograms ###############################################

# Histogram of delays -------------------------------------------------------
ggplot(data_all_3, aes(x = delay)) +
  geom_histogram(aes(y = ..density..),
                 binwidth = 0.5,
                 color = "black", fill = "#1874CD") +
  theme_bw() +
  xlab("Delay in minutes for trains from Zurich to Bern")

# Histogram of hourly precipitation (Bern) ----------------------------------
ggplot(data_all_3, aes(x = rain)) +
  geom_histogram(aes(y = ..density..),
                 binwidth = 0.2,
                 color = "black", fill = "#1874CD") +
  theme_bw() +
  ylim(c(0, 0.15)) +
  xlab("Total hourly precipiation in mm in Bern")

# Scatter plot: precipitation vs delay --------------------------------------
ggplot(data_all_3, aes(rain, delay)) +
  geom_point(color = "black", fill = "#1874CD",
             stroke = 0.1, size = 2, shape = 21) +
  theme_bw() +
  xlab("Hourly preciptiation in mm") +
  ylab("Delay in minutes")

############# Monthly averages of delays and rain ###########################
monthly <- data_all_3 %>%
  group_by(month) %>%
  summarise(
    mean_delay  = mean(delay),
    mean_rain   = mean(rain_zh),
    mean_rain_6 = mean(rain_6_zh),
    mean_rain_12 = mean(rain_12)
  )

ggplot(monthly, aes(x = month, y = mean_rain)) +
  stat_summary(fun = "mean", geom = "bar", fill = "#1874CD", col = "black")

# Focus on “rainy” months (here: May–August) --------------------------------
rainy_months <- c("05", "06", "07", "08")
data_all_4   <- data_all_3[data_all_3$month %in% rainy_months, ]


##########################################
# Pre-Test and Confounder-Test --------------------------------------------

method1 <- "first"
k1      <- floor(nrow(data_all_4)^0.4)

set.seed(123)
t1 <- CTC_causality_permutation_test(data_all_4$rain_zh,
                                     data_all_4$delay,
                                     k = k1, method = method1, R = 2000)
set.seed(123)
t2 <- CTC_causality_permutation_test(data_all_4$rain_3_zh,
                                     data_all_4$delay,
                                     k = k1, method = method1, R = 2000)

t1_crit <- quantile(t1$diff12p, 0.975); t1_crit
t2_crit <- quantile(t2$diff12p, 0.975); t2_crit

t1_crit < t1$diff12; t1$diff12     # check significance for 1h rain
t2_crit < t2$diff12; t2$diff12     # check significance for 3h rain


# Confounder Test Statistic--------------------------------------------------
sqrt(k1) * (t1$ctc21 - 0.5)
sqrt(k1) * (t2$ctc21 - 0.5)

################# Equal Tail Test (delay vs rain) ###########################
test_results <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(test_results) <- c("tail", "k_delay", "k_rain",
                            "test_val", "delay_est", "rain_est")

delay_pos  <- data_all_4$delay[data_all_4$delay > 0]
rain_pos   <- data_all_4$rain_zh[data_all_4$rain_zh > 0]
rain3_pos  <- data_all_4$rain_3_zh[data_all_4$rain_3_zh > 0]

shift_delay  <- 0
shift_rain   <- quantile(rain_pos,  0.95)
shift_rain3  <- quantile(rain3_pos, 0.95)

# Hill/Pareto plots for each variable --------------------------------------
p1 <- pareto_hill(delay_pos, shift = shift_delay)
p2 <- pareto_hill(rain_pos,  shift = shift_rain)
p3 <- pareto_hill(rain3_pos, shift = shift_rain3)

ggsave(paste0("figures/hill_plots_train_delay.pdf"), plot = p1$plot,
       device = "pdf", width = 14, height = 10, path = application_path)
ggsave(paste0("figures/hill_plots_train_rain.pdf"),  plot = p2$plot,
       device = "pdf", width = 14, height = 10, path = application_path)
ggsave(paste0("figures/hill_plots_train_rain3.pdf"), plot = p3$plot,
       device = "pdf", width = 14, height = 10, path = application_path)

# Optimal k for delay and rain tails ----------------------------------------
k_delay <- find_optimal_k(delay_pos + shift_delay)$k_star;  k_delay
k_rain  <- find_optimal_k(rain_pos  + shift_rain)$k_star;   k_rain
k_rain3 <- find_optimal_k(rain3_pos + shift_rain3)$k_star;  k_rain3

# Two-sample equal tail index tests -----------------------------------------
t1 <- tail_test(delay_pos + shift_delay,
                rain_pos  + shift_rain,
                min(k_delay, k_rain), 0.2);  t1
t2 <- tail_test(delay_pos + shift_delay,
                rain3_pos + shift_rain3,
                min(k_delay, k_rain3), 0.2); t2



# Store results -------------------------------------------------------------
test_results[1, 1] <- "rain"
test_results[1, 2] <- k_delay
test_results[1, 3] <- k_rain
test_results[1, 4:6] <- t1

test_results[2, 1] <- "rain3"
test_results[2, 2] <- k_delay
test_results[2, 3] <- k_rain3
test_results[2, 4:6] <- t2

test_results

#########################################
# QQ-Plot: precipitation (3h) vs delays ------------------------------------
n  <- nrow(data_all_4)
xq <- quantile(data_all_4$rain_3_zh,
               probs = seq(0, 1, length.out = n))
yq <- quantile(data_all_4$delay,
               probs = seq(0, 1, length.out = n))
qs <- cbind.data.frame(xq, yq)

qq_plot <- ggplot(qs, aes(x = xq, y = yq)) +
  geom_point(size = 4, col = "#1f78b4") +
  geom_abline(slope = 1, intercept = 0,
              color = "black", linetype = "dashed", size = 1) +
  theme_bw() +
  theme(
    text = element_text(size = 35),
    legend.position = c(0.01, 0.98),
    legend.justification = c(0, 1)
  ) +
  xlab("Precipitation Quantiles") +
  ylab("Delay Quantiles");qq_plot

ggsave(paste0("figures/train_qq.pdf"), plot = qq_plot, device = "pdf",
       width = 10, height = 10, path = application_path)

##################### LiNGAM and causal discovery ###########################
# LiNGAM (pcalg) on 'mat' (defined earlier) ---------------------------------
pcalg::lingam(mat[,1:3])

causal_discovery(mat[, 1:3], method = c("direct_lingam"))


################### Convergence study for CTC (3h rain vs delay) ###########
n_vals <- seq(500, nrow(data_all_4), 100)
ks     <- floor(c(nrow(data_all_4)^0.4, nrow(data_all_4)^0.5))
method <- "first"

ctcs <- matrix(nrow = length(n_vals) * 3 * length(ks), ncol = 3)
i    <- 1

for (k1 in ks) {
  for (n in n_vals) {
    print(n)
    ctc_1 <- 0
    ctc_2 <- 0
    delta <- 0
    for (r in 1:100) {
      s <- sort(sample(nrow(data_all_4), n))
      x <- as.numeric(data_all_4[s, "rain_3_zh"])
      y <- as.numeric(data_all_4[s, "delay"])
      ctc_1[r] <- causal_tail_coeff_basic(x, y,
                                          k = k1, both_tails = FALSE,
                                          method = method)
      ctc_2[r] <- causal_tail_coeff_basic(y, x,
                                          k = k1, both_tails = FALSE,
                                          method = method)
      delta[r] <- ctc_1[r] - ctc_2[r]
    }
    
    ctcs[i, 2]     <- mean(ctc_1)
    ctcs[i + 1, 2] <- mean(ctc_2)
    ctcs[i + 2, 2] <- mean(delta)
    
    ctcs[i, 1]     <- paste0("precipitation \u2192 delay (k=", k1, ")")
    ctcs[i + 1, 1] <- paste0("delay \u2192 precipitation (k=", k1, ")")
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
    paste0("precipitation \u2192 delay (k=", ks[1], ")"),
    paste0("delay \u2192 precipitation (k=", ks[1], ")"),
    paste0("precipitation \u2192 delay (k=", ks[2], ")"),
    paste0("delay \u2192 precipitation (k=", ks[2], ")"),
    paste0("delta (k=", ks[1], ")"),
    paste0("delta (k=", ks[2], ")")
  )
)

pd     <- position_dodge(0.00)
cols   <- c("#1f78b4", "#a6cee3", "orange", "#fdbf6f",
            "#1f78b4", "orange")
shapes <- c(21, 21, 21, 21, 4, 4)

# Legend labels for 3h precipitation vs delay -------------------------------
labels_cust_3 <- c(
  bquote(hat(Gamma)["precipitation 3h" %->% "delay"] * "," ~ k == .(ks[1])),
  bquote(hat(Gamma)["delay" %->% "precipitation 3h"] * "," ~ k == .(ks[1])),
  bquote(hat(Gamma)["precipitation 3h" %->% "delay"] * "," ~ k == .(ks[2])),
  bquote(hat(Gamma)["delay" %->% "precipitation 3h"] * "," ~ k == .(ks[2])),
  bquote(hat(Delta)["precipitation 3h" %->% "delay"] * "," ~ k == .(ks[1])),
  bquote(hat(Delta)["precipitation 3h" %->% "delay"] * "," ~ k == .(ks[2]))
)
# Legend labels for precipitation vs delay -------------------------------
labels_cust <- c(
  bquote(hat(Gamma)["precipitation" %->% "delay"] * "," ~ k == .(ks[1])),
  bquote(hat(Gamma)["delay" %->% "precipitation"] * "," ~ k == .(ks[1])),
  bquote(hat(Gamma)["precipitation" %->% "delay"] * "," ~ k == .(ks[2])),
  bquote(hat(Gamma)["delay" %->% "precipitation"] * "," ~ k == .(ks[2])),
  bquote(hat(Delta)["precipitation" %->% "delay"] * "," ~ k == .(ks[1])),
  bquote(hat(Delta)["precipitation" %->% "delay"] * "," ~ k == .(ks[2]))
)

conv_plot <- ggplot(ctcs, aes(n, CTC, color = group)) +
  geom_line(alpha = .7, size = 1.3, position = pd) +
  geom_point(size = 3.5,
             position = pd,
             aes(color = group, fill = group, shape = group)) +
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
  theme(legend.position = c(0.45, 0.85))

conv_plot

ggsave(paste0("figures/train_", method, "_3h_conv_plot.pdf"),
       plot = conv_plot, device = "pdf",
       width = 14, height = 10, path = application_path)

#### Delta (permutation) plots ##############################################
k1 <- floor(nrow(data_all_4)^0.4)
k2 <- floor(nrow(data_all_4)^0.5)

set.seed(123)
perm_res  <- CTC_causality_permutation_test(data_all_4$rain_3_zh,
                                            data_all_4$delay,
                                            k = k1, method = "first", R = 2000)
set.seed(123)
perm_res2 <- CTC_causality_permutation_test(data_all_4$rain_3_zh,
                                            data_all_4$delay,
                                            k = k2, method = "first", R = 2000)

perm <- cbind.data.frame(perm_res$diff12p, perm_res2$diff12p)
colnames(perm) <- c("k1", "k2")

perm_long <- perm |>
  pivot_longer(
    cols      = everything(),
    names_to  = "group",
    values_to = "diff"
  )

# Hist for k1 ---------------------------------------------------------------
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
           label = expression(hat(Delta)[precipitation~"3h" %->% delay]),
           vjust = 1.2, hjust = 1.05, color = "black", size = 10) +
  scale_fill_manual(
    name   = NULL,
    values = c("k1" = "#1f78b4"),
    labels = c("k1" = "k=27")
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

# Hist for k2 ---------------------------------------------------------------
bs2 <- ggplot(perm_long, aes(x = diff, fill = group)) +
  geom_histogram(
    data = subset(perm_long, group == "k2"),
    aes(y = after_stat(density), fill = "k2"),
    binwidth = 0.02, color = "black", alpha = 1
  ) +
  geom_vline(xintercept = perm_res2$diff12, color = "black",
             linetype = "dashed", linewidth = 1) +
  annotate("text",
           x = perm_res$diff12, y = 10,
           label = expression(hat(Delta)[precipitation~"3h" %->% delay]),
           vjust = 1.2, hjust = 1.05, color = "black", size = 10) +
  scale_fill_manual(
    name   = NULL,
    values = c("k2" = "orange"),
    labels = c("k2" = "k=63")
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

bs_plot <- grid.arrange(bs1, bs2, nrow = 2, ncol = 1);bs_plot

ggsave(paste0("figures/train_3h_", method, "_delta_plot.pdf"),
       plot = bs_plot, device = "pdf",

       width = 14, height = 10, path = application_path)
