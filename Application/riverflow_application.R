# riverflows and precipitation

# Load required libraries ---------------------------------------------------
library(dplyr)     # data manipulation (group_by, summarise, etc.)
library(ggplot2)   # plotting
library(pcalg)     # causal discovery (LiNGAM, etc.)
library(knitr)     # table formatting (kable)

# Set working directory to GitHub repo root ---------------------------------
setwd()

# Define application path (relative to working directory) -------------------
application_path <- "Application/"

# Source custom functions (tail causality, Hill estimator, etc.) -----------
source("functions.R")

# Read river flow and precipitation data -----------------------------------
Donau <- read.csv2(paste0(application_path, "data/passau.csv"),
                   header = TRUE, sep = ";")  # Danube at Passau
Passau_prec <- read.table(paste0(application_path, "data/passau_precipitation.txt"),
                          sep = ";", header = TRUE)

Main <- read.csv2(paste0(application_path, "data/wurzburg.csv"),
                  header = TRUE, sep = ";")   # Main at Wuerzburg
Wurzburg_prec <- read.table(paste0(application_path, "data/wurzburg_precipitation.txt"),
                            sep = ";", header = TRUE)

SF <- read.csv2(paste0(application_path, "data/schweinfurth.csv"),
                header = TRUE, sep = ";")     # Main at Schweinfurt (upstream)

# ----------------------------------------------------------------------------
# Date conversion for river series
# ----------------------------------------------------------------------------
rivers <- list(Donau, Main, SF)

rivers <- lapply(rivers, function(df) {
  # Convert German date format "dd.mm.yyyy" to Date
  df$date <- as.Date(format(df[, "Datum"], format = "%d.%m.%Y"),
                     format = "%d.%m.%Y")
  df$Datum <- NULL  # drop original character date column
  return(df)
})

head(rivers[[1]])
tail(rivers[[2]])

# ----------------------------------------------------------------------------
# Date conversion and basic filtering for precipitation series
# ----------------------------------------------------------------------------
precipitation <- list(Passau_prec, Wurzburg_prec)

precipitation <- lapply(precipitation, function(df) {
  # Convert numeric date "YYYYMMDD" to Date
  df$dates_converted <- as.Date(format(df[,"MESS_DATUM"], format = "%Y%m%d"),
                                format = "%Y%m%d")
  df[,"date"] <- as.Date(df[,"dates_converted"], format = "%d.%m.%Y")
  df$MESS_DATUM     <- NULL
  df$dates_converted <- NULL

  # Keep only precipitation and date
  df <- df[, c("RS", "date")]    # RS: daily precipitation
  # Extract month as character ("01", "02", ...)
  df$month <- format(as.Date(df$date, format = "%d/%m/%Y"), "%m")

  # Average precipitation per month (for exploratory purposes) --------------
  grouped_data <- df %>%
    group_by(month) %>%
    summarise(mean_value = mean(RS))
  plot(grouped_data$month, grouped_data$mean_value)

  # Keep only late spring and summer months: May–August ---------------------
  df <- df[df$month %in% c("05", "06", "07", "08"), ]
  return(df)
})

head(precipitation[[1]])

# Quick sanity checks on date ranges ---------------------------------------
for (p in precipitation) {
  print(summary(p$date))
}
for (p in rivers) {
  print(summary(p$date))
}

# ----------------------------------------------------------------------------
# Merge precipitation with river flows 
# ----------------------------------------------------------------------------

# Passau: precipitation + Danube discharge
df1 <- merge(precipitation[[1]],
             rivers[[1]][, c("date", "Mittelwert")],  # Mittelwert = discharge
             by = "date")

# Wuerzburg: precipitation + Main discharge (Wuerzburg + Schweinfurt)
df2 <- merge(precipitation[[2]],
             rivers[[2]][, c("date", "Mittelwert")],
             by = "date")
df2 <- merge(df2,
             rivers[[3]][, c("date", "Mittelwert")],
             by = "date")

summary(df1$date)
summary(df2$date)

############ QQ plot: precipitation vs river discharge #######################
n  <- nrow(df2)
xq <- quantile(df2$RS,         probs = seq(0, 1, length.out = n)) # precip
yq <- quantile(df2$Mittelwert, probs = seq(0, 1, length.out = n)) # discharge
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
  ylab("River Discharge Quantiles")

ggsave(paste0("figures/river_qq.pdf"), plot = qq_plot, device = "pdf",
       width = 10, height = 10, path = application_path)

###################### Build analysis data frames ###########################
# dfs[[1]] = Passau/Donau, dfs[[2]] = Wuerzburg/Main (+Schweinfurth)
dfs <- list(Donau = df1, Main = df2)

dfs <- lapply(dfs, function(df) {
  # Standardize column names:
  #  date, precipitation (RS), month, river (=discharge at station)
  colnames(df) <- c("date", "precipitation", "month", "river")

  # Create lagged precipitation variables (1–8 days)
  df$precipitation_1 <- df$precipitation
  df$precipitation_2 <- df$precipitation
  df$precipitation_3 <- df$precipitation
  df$precipitation_4 <- df$precipitation
  df$precipitation_5 <- df$precipitation
  df$precipitation_6 <- df$precipitation
  df$precipitation_7 <- df$precipitation
  df$precipitation_8 <- df$precipitation


  # Create lags: precipitation_t-1,...,t-8 and change in discharge ----------
  for (i in 9:nrow(df)) {
    df$precipitation_1[i] <- df$precipitation[(i - 1)]
    df$precipitation_2[i] <- df$precipitation[(i - 2)]
    df$precipitation_3[i] <- df$precipitation[(i - 3)]
    df$precipitation_4[i] <- df$precipitation[(i - 4)]
    df$precipitation_5[i] <- df$precipitation[(i - 5)]
    df$precipitation_6[i] <- df$precipitation[(i - 6)]
    df$precipitation_7[i] <- df$precipitation[(i - 7)]
    df$precipitation_8[i] <- df$precipitation[(i - 8)]
  }
  return(df)
})

head(dfs[[1]])

# For Main data, add second river column (Schweinfurth) renamed as river2 ---
colnames(dfs[[2]]) <- c("date", "precipitation", "month",
                        "river", "river2",
                        "precipitation_1", "precipitation_2",
                        "precipitation_3", "precipitation_4",
                        "precipitation_5", "precipitation_6",
                        "precipitation_7", "precipitation_8")

# CTC without confounder for Main -------------------------
df<-dfs[[2]]
k1 <- floor(nrow(df)^0.4)

c1<-causal_tail_coeff_basic(df$precipitation_1, df$river, k = k1, both_tails = FALSE)
c2<-causal_tail_coeff_basic(df$river, df$precipitation_1, k = k1, both_tails = FALSE)
set.seed(123)
t1 <- CTC_causality_permutation_test(df$precipitation_1, df$river,
                                     k = k1, method = method, R = 2000)
round(quantile(t1$diff12p, 0.975), 4)
t1$diff12
t_conf <- sqrt(k1) * (c2 - 0.5);t_conf
################# Tail-causality tests for precipitation lags #####################
method <- "first"               
i      <- 1
var_1  <- "river"                # response variable

table_res <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(table_res) <- c("river", "ctc1", "ctc2", "delta", "deltacrit",
                         "t_conf", "lingam")

# Precipitation lag variables to test --------------------------------------
lags <- c("precipitation_1", "precipitation_2", "precipitation_3",
          "precipitation_4", "precipitation_5", "precipitation_6",
          "precipitation_7", "precipitation_8")

for (df in dfs) {
  k1 <- floor(nrow(df)^0.4)      # tail sample size

  for (j in 1:length(lags)) {
    var <- lags[j]

    if (i == 2) {
      # For Main: include confounder river2 in LGPD causal tail coefficient -
      c1 <- LGPD_causal_tail_coeff(df[, var], df[, var_1],
                                   H = df$river2, k = k1, method = method)
      c2 <- LGPD_causal_tail_coeff(df[, var_1], df[, var],
                                   H = df$river2, k = k1, method = method)
      set.seed(123)
      t1 <- CTC_causality_permutation_test(df[, var], df[, var_1],
                                           H = df$river2, k = k1,
                                           method = method, R = 2000)
    } else {
      # For Donau: no confounder; use basic causal tail coefficient ----------
      c1 <- causal_tail_coeff_basic(df[, var], df[, var_1],
                                    k = k1, both_tails = FALSE,
                                    method = method)
      c2 <- causal_tail_coeff_basic(df[, var_1], df[, var],
                                    k = k1, both_tails = FALSE,
                                    method = method)
      set.seed(123)
      t1 <- CTC_causality_permutation_test(df[, var], df[, var_1],
                                           k = k1, method = method, R = 2000)
    }

    # Test statistic for Confounder-Test
    t_conf <- sqrt(k1) * (c2 - 0.5)

    # LiNGAM 
    mat  <- cbind(df[, var], df[, var_1])
    ling <- pcalg::lingam(mat)

    # Row index in result table
    index <- (i - 1) * length(lags) + j

    table_res[index, 1] <- names(dfs)[i]
    table_res[index, 2] <- round(c1, 4)
    table_res[index, 3] <- round(c2, 4)
    table_res[index, 4] <- round(t1$diff12, 4)
    table_res[index, 5] <- round(quantile(t1$diff12p, 0.975), 4)
    table_res[index, 6] <- round(t_conf, 4)
    table_res[index, 7] <- round(ling$Bpruned[2, 1], 4)
    print(index)
  }
  i <- i + 1
}

table_res

# LaTeX table for river results (ctc1,...,lingam) ----------------------
kable(table_res[table_res$river == "Main", 2:7],
      format = "latex", booktabs = TRUE)
kable(table_res[table_res$river == "Donau", 2:7],
      format = "latex", booktabs = TRUE)

############# Equal tail index test for rain vs discharge ###################
i <- 1
test_results <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(test_results) <- c("river", "k_rain", "k_river",
                            "test_val", "rain_est", "river_est")

for (df in dfs) {
  df_tmp     <- df
  rain_pos   <- df_tmp$precipitation[df_tmp$precipitation > 0]
  discharge_pos <- df_tmp$river

  # Shifts for tail estimation (rain to upper endpoint, discharge no shift)
  shift_rain      <- max(rain_pos)
  shift_discharge <- 0

  # Hill/Pareto plots for both variables -----------------------------------
  p1 <- pareto_hill(rain_pos, shift = shift_rain);       plot(p1$plot)
  p2 <- pareto_hill(discharge_pos, shift = shift_discharge); plot(p2$plot)

  ggsave(paste0("figures/hill_plots_rain_", names(dfs)[i], ".pdf"),
         plot  = p1$plot, device = "pdf",
         width = 14, height = 10, path = application_path)
  ggsave(paste0("figures/hill_plots_river_", names(dfs)[i], ".pdf"),
         plot  = p2$plot, device = "pdf",
         width = 14, height = 10, path = application_path)

  # Optimal k for both variables -------------------------------------------
  k_rain      <- find_optimal_k(rain_pos + shift_rain)$k_star
  k_discharge <- find_optimal_k(discharge_pos + shift_discharge)$k_star

  # Two-sample tail index equality test ------------------------------------
  t1 <- tail_test(rain_pos + shift_rain,
                  discharge_pos + shift_discharge,
                  min(k_rain, k_discharge), 0.2)


  # Fill summary table -----------------------------------------------------
  test_results[i, 1] <- names(dfs)[i]
  test_results[i, 2] <- k_rain
  test_results[i, 3] <- k_discharge
  test_results[i, 4:6] <- t1

  i <- i + 1
}

test_results

############ Convergence study (with confounder for Main) ###################
data   <- dfs[2]$Main  # for Main; use dfs[1]$Donau for Danube
n_vals <- seq(500, nrow(data), 100)
ks     <- floor(c(nrow(data)^0.4, nrow(data)^0.5))
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
      s <- sort(sample(nrow(data), n))
      x <- as.numeric(data[s, "precipitation_1"])
      y <- as.numeric(data[s, "river"])
      h <- as.numeric(data[s, "river2"])
      thres <- min(0.95, (1 - (k1 + 1) / n))  # adjust threshold to n,k1

      # Conditional CTC with confounder h (river2)
      ctc_1[r] <- LGPD_causal_tail_coeff(x, y, H = h, k = k1,
                                         method = method,
                                         threshold_q = thres)
      ctc_2[r] <- LGPD_causal_tail_coeff(y, x, H = h, k = k1,
                                         method = method,
                                         threshold_q = thres)
      # For unconditional versions, see commented lines:
      # ctc_1[r] <- causal_tail_coeff_basic(x, y, k = k1,
      #                                     both_tails = FALSE, method = method)
      # ctc_2[r] <- causal_tail_coeff_basic(y, x, k = k1,
      #                                     both_tails = FALSE, method = method)
      delta[r] <- ctc_1[r] - ctc_2[r]
    }

    # Store means over repetitions -----------------------------------------
    ctcs[i, 2]     <- mean(ctc_1)
    ctcs[i + 1, 2] <- mean(ctc_2)
    ctcs[i + 2, 2] <- mean(delta)

    ctcs[i, 1]     <- paste0("precipitation \u2192 river discharge (k=", k1, ")")
    ctcs[i + 1, 1] <- paste0("river discharge \u2192 precipitation (k=", k1, ")")
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
    paste0("precipitation \u2192 river discharge (k=", ks[1], ")"),
    paste0("river discharge \u2192 precipitation (k=", ks[1], ")"),
    paste0("precipitation \u2192 river discharge (k=", ks[2], ")"),
    paste0("river discharge \u2192 precipitation (k=", ks[2], ")"),
    paste0("delta (k=", ks[1], ")"),
    paste0("delta (k=", ks[2], ")")
  )
)

pd     <- position_dodge(0.00)
cols   <- c("#1f78b4", "#a6cee3", "orange", "#fdbf6f",
            "#1f78b4", "orange")
shapes <- c(21, 21, 21, 21, 4, 4)

# Legend labels for Danube and Main (not both used simultaneously) ---------
labels_cust <- c(
  bquote(hat(Gamma)["precipitation" %->% "Danube"] * "," ~ k == .(ks[1])),
  bquote(hat(Gamma)["Danube" %->% "precipitation"] * "," ~ k == .(ks[1])),
  bquote(hat(Gamma)["precipitation" %->% "Danube"] * "," ~ k == .(ks[2])),
  bquote(hat(Gamma)["Danube" %->% "precipitation"] * "," ~ k == .(ks[2])),
  bquote(hat(Delta)["precipitation" %->% "Danube"] * "," ~ k == .(ks[1])),
  bquote(hat(Delta)["precipitation" %->% "Danube"] * "," ~ k == .(ks[2]))
)

labels_cust2 <- c(
  bquote(hat(Gamma)["precipitation" %->% "Main"] * "," ~ k == .(ks[1])),
  bquote(hat(Gamma)["Main" %->% "precipitation"] * "," ~ k == .(ks[1])),
  bquote(hat(Gamma)["precipitation" %->% "Main"] * "," ~ k == .(ks[2])),
  bquote(hat(Gamma)["Main" %->% "precipitation"] * "," ~ k == .(ks[2])),
  bquote(hat(Delta)["precipitation" %->% "Main"] * "," ~ k == .(ks[1])),
  bquote(hat(Delta)["precipitation" %->% "Main"] * "," ~ k == .(ks[2]))
)

# Convergence plot ----------------------------------------------------------
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
  theme(legend.position = c(0.45, 0.86))

conv_plot

ggsave(paste0("figures/river_conv_plot_main_", method, ".pdf"),
       plot = conv_plot, device = "pdf",
       width = 14, height = 10, path = application_path)

############################ Delta Histogram ################################
data <- dfs[1]$Donau     # use Donau here; change to dfs[2]$Main for Main
var1 <- "precipitation_1"
var2 <- "river"

k1 <- floor(nrow(data)^0.4)
k2 <- floor(nrow(data)^0.5)
#h  <- as.numeric(data[, "river2"])

set.seed(123)
perm_res  <- CTC_causality_permutation_test(data[, var1], data[, var2],
                                            k = k1, method = "first", R = 2000) # ,H=h
set.seed(123)
perm_res2 <- CTC_causality_permutation_test(data[, var1], data[, var2],
                                            k = k2, method = "first", R = 2000) # ,H=h

perm <- cbind.data.frame(perm_res$diff12p, perm_res2$diff12p)
colnames(perm) <- c("k1", "k2")

perm_long <- perm |>
  pivot_longer(
    cols      = everything(),
    names_to  = "group",
    values_to = "diff"
  )

# Histogram for k1 ----------------------------------------------------------
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
           label = expression(hat(Delta)[precipitation %->% Danube]),
           vjust = 1.2, hjust = 1.05, color = "black", size = 10) +
  scale_fill_manual(
    name   = NULL,
    values = c("k1" = "#1f78b4"),
    labels = c("k1" = "k=24")
  ) +
  labs(x = expression("Bootstrapped " ~ hat(Delta)), y = "Density") +
  theme_bw() +
  theme(
    text = element_text(size = 26),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  scale_x_continuous(limits = c(-0.35, 0.35)) +
  scale_y_continuous(limits = c(0, 10))
bs1

# Histogram for k2 ----------------------------------------------------------
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
           label = expression(hat(Delta)[precipitation %->% Danube]),
           vjust = 1.2, hjust = 1.05, color = "black", size = 10) +
  scale_fill_manual(
    name   = NULL,
    values = c("k2" = "orange"),
    labels = c("k2" = "k=54")
  ) +
  labs(x = expression("Bootstrapped " ~ hat(Delta)), y = "Density") +
  theme_bw() +
  theme(
    text = element_text(size = 26),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  scale_x_continuous(limits = c(-0.35, 0.35)) +
  scale_y_continuous(limits = c(0, 10))
bs2

bs_plot <- grid.arrange(bs1, bs2, nrow = 2, ncol = 1); bs_plot

ggsave(paste0("figures/river_bs_plot_danube_delta_", method, ".pdf"),
       plot = bs_plot, device = "pdf",
       width = 14, height = 10, path = application_path)

# Precipitation histogram (marginal distribution) ---------------------------
hist <- ggplot(data, aes(x = precipitation)) +
  geom_histogram(aes(y = ..density..),
                 bins = 30, fill = "#1f78b4",
                 color = "black", alpha = 0.6) +
  labs(x = "Precipitation", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 26))

ggsave(paste0("figures/precipiation_hist.pdf"), plot = hist, device = "pdf",
       width = 8, height = 10, path = application_path)

