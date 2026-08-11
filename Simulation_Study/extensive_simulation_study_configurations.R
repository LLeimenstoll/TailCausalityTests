#Simulation Study for test evaluation for different configurations A-F
setwd()
simulation_path <- "Simulation_Study/"

source("functions.R")
source("testing_strategy.R")

# -----------------------------
# Settings
# -----------------------------
n <- 10000
k1 <- floor(n^0.4)
rep_per_setting <- 1000

df_values <- c(2, 3, 4)

# Combinations of df
df_grid <- expand.grid(df1 = df_values, df2 = df_values)

# Configurations A-F
scenario_grid <- data.frame(
  scenario = c("A", "B", "C", "D", "E", "F"),
  conf_type = c("none", "none", "light", "light", "heavy", "heavy"),
  causal = c(0, 1, 0, 1, 0, 1),
  stringsAsFactors = FALSE
)

results_list <- list()
row_counter <- 1

# -----------------------------
# Simulation
# -----------------------------
for (i in 1:nrow(df_grid)) {

  df1 <- df_grid$df1[i]
  df2 <- df_grid$df2[i]

  for (j in 1:nrow(scenario_grid)) {

    scenario <- scenario_grid$scenario[j]
    conf_type <- scenario_grid$conf_type[j]
    causal <- scenario_grid$causal[j]

    # Set dfh depending on df1 and df2 and configuration
    if (conf_type == "none") {
      dfh <- NA
      conf_indicator <- 0
    } else if (conf_type == "light") {
      dfh <- max(df1, df2) + 1
      conf_indicator <- 1
    } else if (conf_type == "heavy") {
      dfh <- max(min(df1, df2) - 1,2)
      conf_indicator <- 1

      # Check
      if (dfh <= 1) {
        warning(paste("Skipping invalid heavy-tail setting: df1 =", df1,
                      "df2 =", df2, "gives dfh =", dfh))
        next
      }
    }

    for (r in 1:rep_per_setting) {
      cat("\r iteration:",r, " out of ", rep_per_setting)
      # Structural coefficients
      beta_h1 <- runif(1, 0.1, 0.9)
      beta_h2 <- runif(1, 0.1, 0.9)
      beta_1  <- runif(1, 0.1, 0.9)

      # Error terms
      e1 <- rt(n, df1)
      e2 <- rt(n, df2)

      # Confounder
      if (conf_indicator == 1) {
        H <- rt(n, dfh)
      } else {
        H <- rep(0, n)
      }

      # Generate data
      X1 <- e1 + conf_indicator * beta_h1 * H
      X2 <- e2 + causal * beta_1 * X1 + conf_indicator * beta_h2 * H

      # Run tests
      res_main <- testing_strategy(X1, X2, k = k1)

      # Run tests with known confounder
      if (conf_indicator == 1) {
        res_conf <- testing_strategy(X1, X2, H = H, k = k1)
      } else {
        res_conf <- rep(NA, 5)
      }

      # Report tail relation
      tail_relation <- ifelse(df1 > df2, "df1>df2",
                              ifelse(df1 < df2, "df1<df2", "df1=df2"))

      # Save
      results_list[[row_counter]] <- data.frame(
        replication = r,
        scenario = scenario,
        conf_type = conf_type,
        causal = causal,
        tail_relation = tail_relation,
        df1 = df1,
        df2 = df2,
        dfh = dfh,
        beta1 = beta_1,
        beta_h1 = beta_h1,
        beta_h2 = beta_h2,

        dir = res_main[1],
        Tc = res_main[2],
        c12 = res_main[3],
        c21 = res_main[4],
        Ti = res_main[5],
        hill1 = res_main[6],
        hill2 = res_main[7],
        p_val_hill = res_main[8],
        p_val_causality = res_main[9],

        dirH = res_conf[1],
        c12H = res_conf[2],
        c21H = res_conf[3],
        p_val_causality_conf = res_conf[4],
        Tc_conf = res_conf[5]

      )

      row_counter <- row_counter + 1
    }

    cat("Finished setting:",
        "scenario =", scenario,
        ", df1 =", df1,
        ", df2 =", df2,
        ", dfh =", dfh, "\n")
  }
}

# Combine
df_results <- do.call(rbind, results_list)

# Export
write.csv(
  df_results,
  file = paste0("results_systematic_simulation.csv"),
  row.names = FALSE
)

head(df_results)
