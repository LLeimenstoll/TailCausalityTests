# functions.R: helper functions for tail indices, causal tail coefficients and tests

# Load required packages ----------------------------------------------------
library(qrmtools)   # Hill_estimator and risk tools
library(maxLik)     # maximum likelihood (used in constrained GPD fits)
library(gridExtra)  # arranging multiple ggplots
library(ggplot2)    # plotting for diagnostic plots (Pareto–Hill)
library(evd)        # extreme value distributions (fpot, pgpd)
library(ismev)      # extreme value modeling (gpd.fit)

#---------------------------------------------------------------------------
# Hill-based test for equality of tail indices (normal approximation)
#---------------------------------------------------------------------------

hill_test <- function(x1, x2, k) {
  # x1, x2: samples from two distributions
  # k:     tail sample size used in Hill/Hoga estimator

  # Hoga tail index estimators (gamma_hat) at t = 1
  g1 <- hoga_estimator_vec(1, x1, k)
  g2 <- hoga_estimator_vec(1, x2, k)

  # Asymptotic variances (gamma^2 / k)
  var1 <- g1^2 / k
  var2 <- g2^2 / k

  # Wald-type statistic for difference in tail indices
  z <- (g1 - g2) / sqrt(var1 + var2)

  # Two-sided p-value under normal approximation
  pval <- 2 * (1 - pnorm(abs(z)))

  list(z = z, p_value = pval)
}

#---------------------------------------------------------------------------
# CTC_causality_permutation_test:
# Pre-test for causal tail coefficients (Pasche et al., 2023), adjusted p-value
#---------------------------------------------------------------------------

CTC_causality_permutation_test <- function(
    X1, X2, H = NULL,
    k = floor(NROW(X1)^0.4),
    R = 10000,
    constrained_fit = FALSE,
    threshold_q = 0.95,
    both_tails = FALSE,
    parametric_F1 = TRUE,
    method = "average"
) {
  n <- NROW(X1)          # sample size

  # If k not specified, default to n^0.4
  if (is.null(k)) {
    k <- floor(n^0.4)
  }

  # Basic checks ------------------------------------------------------------
  if (threshold_q >= 1 | threshold_q <= 0) {
    stop("threshold_q must be between 0 and 1.")
  }
  if (k <= 1 | k > floor((1 - threshold_q) * n)) {
    stop("k must be greater than 1 and smaller than the number of threshold exceedences.")
  }
  if (both_tails) {
    stop("Both tails CTC not possible in CTC_causality_permutation_test.")
  }
  if (!parametric_F1) {
    stop("Non-Parametric F1 not possible for the LGPD CTC in CTC_causality_permutation_test.")
  }

  # Ranks of X1 and X2 ------------------------------------------------------
  r1 <- rank(X1, ties.method = method)
  r2 <- rank(X2, ties.method = method)

  # Transform to (approximate) uniforms F1(X1), F2(X2) ----------------------
  if (is.null(H)) {
    # Non-parametric CTC: just rank-based empirical CDF
    FX1 <- r1 / n
    FX2 <- r2 / n
  } else {
    # LGPD CTC with covariate H
    Hcs <- scale(H)                     # center and scale confounder
    u2  <- quantile(X2, threshold_q)    # thresholds for POT
    u1  <- quantile(X1, threshold_q)

    if (constrained_fit) {
      # Constrained GPD fits using maxLik (monotonicity / positivity)

      # X2 | H
      fit2 <- gpd_constraints_fit_maxLik(
        X2, threshold = u2, ydat = as.matrix(Hcs), sigl = 1, show = FALSE,
        constraints = list(
          ineqA = matrix(
            c(1, min(Hcs), 0,
              1, max(Hcs), 0),
            nrow = 2, ncol = 3, byrow = TRUE
          ),
          ineqB = matrix(0, nrow = 2, ncol = 1, byrow = TRUE)
        )
      )
      FX2 <- pgpd(X2, loc = u2,
                  scale = fit2$mle[1] + fit2$mle[2] * Hcs,
                  shape = fit2$mle[3])
      FX2 <- (FX2 * (1 - threshold_q) + threshold_q) * (X2 >= u2) +
        r2 / n * (X2 < u2)

      # X1 | H
      fit1 <- gpd_constraints_fit_maxLik(
        X1, threshold = u1, ydat = as.matrix(Hcs), sigl = 1, show = FALSE,
        constraints = list(
          ineqA = matrix(
            c(1, min(Hcs), 0,
              1, max(Hcs), 0),
            nrow = 2, ncol = 3, byrow = TRUE
          ),
          ineqB = matrix(0, nrow = 2, ncol = 1, byrow = TRUE)
        )
      )
      FX1 <- pgpd(X1, loc = u1,
                  scale = fit1$mle[1] + fit1$mle[2] * Hcs,
                  shape = fit1$mle[3])
      FX1 <- (FX1 * (1 - threshold_q) + threshold_q) * (X1 >= u1) +
        r1 / n * (X1 < u1)
    } else {
      # Post-fit corrected LGPD CTC (unconstrained gpd.fit)

      # X2 | H
      fit2 <- gpd.fit(X2, threshold = u2,
                      ydat = as.matrix(Hcs), sigl = 1, show = FALSE)
      FX2 <- pgpd(
        X2, loc = u2,
        scale = pmax(fit2$mle[1] + fit2$mle[2] * Hcs, 0.0001)[1:length(Hcs)],
        shape = fit2$mle[3]
      )
      FX2 <- (FX2 * (1 - threshold_q) + threshold_q) * (X2 >= u2) +
        r2 / n * (X2 < u2)

      # X1 | H
      fit1 <- gpd.fit(X1, threshold = u1,
                      ydat = as.matrix(Hcs), sigl = 1, show = FALSE)
      FX1 <- pgpd(
        X1, loc = u1,
        scale = pmax(fit1$mle[1] + fit1$mle[2] * Hcs, 0.0001)[1:length(Hcs)],
        shape = fit1$mle[3]
      )
      FX1 <- (FX1 * (1 - threshold_q) + threshold_q) * (X1 >= u1) +
        r1 / n * (X1 < u1)
    }
  }

  # Causal tail coefficients ctc12: X1 -> X2; ctc21: X2 -> X1 ---------------
  ctc12 <- mean(FX2[FX1 > (1 - k / n)])
  ctc21 <- mean(FX1[FX2 > (1 - k / n)])
  diff12 <- ctc12 - ctc21

  # Permutation step: generate null distribution under no causal asymmetry --
  ctc12p <- rep(as.double(NA), R)
  ctc21p <- rep(as.double(NA), R)

  for (i in 1:R) {
    FX1p <- rep(as.double(NA), n)
    FX2p <- rep(as.double(NA), n)

    # Per observation, randomly choose which variable gets FX1 or FX2
    permbool <- sample(x = c(TRUE, FALSE), size = n, replace = TRUE)

    FX1p[!permbool] <- FX1[!permbool]
    FX1p[ permbool] <- FX2[ permbool]
    FX2p[!permbool] <- FX2[!permbool]
    FX2p[ permbool] <- FX1[ permbool]

    ctc12p[i] <- mean(FX2p[FX1p > (1 - k / n)])
    ctc21p[i] <- mean(FX1p[FX2p > (1 - k / n)])
  }

  diff12p <- ctc12p - ctc21p

  # Two-sided Monte Carlo p-value -------------------------------------------
  Pmc <- (1 + sum(abs(diff12p) >= abs(diff12))) / (R + 1)

  return(list(
    Pmc    = Pmc,
    ctc12  = ctc12,
    ctc21  = ctc21,
    diff12 = diff12,
    ctc12p = ctc12p,
    ctc21p = ctc21p,
    diff12p = diff12p
  ))
}

#---------------------------------------------------------------------------
# Hoga (2017) tail index estimator
#---------------------------------------------------------------------------

hoga_estimator_vec <- function(t_vec, x_vals, k) {
  # t_vec: vector of fractions in (0,1]; tail index estimated for each t
  # x_vals: sample (assumed heavy-tailed)
  # k: base tail sample size at t = 1

  n <- length(x_vals)
  results <- numeric(length(t_vec))

  for (i in seq_along(t_vec)) {
    t  <- t_vec[i]
    nt <- floor(n * t)
    kt <- floor(k * t)

    # Use the first nt observations; sort them in ascending order
    x_nt        <- x_vals[1:nt]
    x_sorted_nt <- sort(x_nt)

    # Threshold = (nt - kt)-th order statistic
    x_kt <- x_sorted_nt[nt - kt]

    # Hill-style log differences for the top kt observations
    logs <- log(x_sorted_nt[(nt - kt + 1):nt] / x_kt)
    est  <- mean(logs)

    results[i] <- est
  }
  return(results)
}

#---------------------------------------------------------------------------
# Integrand for Hoga's index test
#---------------------------------------------------------------------------

integrator <- function(t, x, y, k) {
  # t: integration variable in [t0,1]
  # x, y: two samples
  # k: tail parameter
  val <- t^2 * (
    (hoga_estimator_vec(t, x, k) - hoga_estimator_vec(t, y, k)) -
      (hoga_estimator_vec(1, x, k) - hoga_estimator_vec(1, y, k))
  )^2
  return(val)
}

#---------------------------------------------------------------------------
# Index test for equality of tail indices (Hoga)
#---------------------------------------------------------------------------

tail_test <- function(x, y, k, t0) {
  # x, y: samples
  # k: tail parameter
  # t0: lower limit of integration (e.g. 0.2)

  hill_x <- hoga_estimator_vec(1, x, k)
  hill_y <- hoga_estimator_vec(1, y, k)

  num <- (hill_x - hill_y)^2

  denum <- integrate(
    integrator, lower = t0, upper = 1,
    x = x, y = y, k = k, subdivisions = 1000
  )$value

  res <- num / denum
  return(c(res, hill_x, hill_y))  # test statistic and two tail index estimates
}

#---------------------------------------------------------------------------
# Pareto–Hill diagnostic plots to choose shift parameter
#---------------------------------------------------------------------------

pareto_hill <- function(x, shift = 0) {
  # x: sample of interest
  # shift: additive constant to improve Pareto tail behavior (for x + shift)

  n <- length(x)
  j <- seq(1, n, 1)

  # Pareto quantiles log((n+1)/j)
  pareto_quantiles <- log((n + 1) / j)

  # Sort data in descending order
  xis <- sort(x, decreasing = TRUE)
  par_df <- data.frame(pq = pareto_quantiles, x = log(xis))

  # Plot of log(X_(j)) vs Pareto quantiles ----------------------------------
  y_max <- par_df$x[n]
  y_min <- par_df$x[1]

  p1 <- ggplot(par_df, aes(x = pq, y = x)) +
    geom_point() +
    theme_bw() +
    theme(text = element_text(size = 35)) +
    xlab(expression(log((n + 1) / j))) +
    ylab("") +
    ylim(c(y_max, y_min))

  # Hill estimator over k = 2,...,n ----------------------------------------
  hill <- c()
  for (i in 2:n) {
    hill[i - 1] <- hoga_estimator_vec(1, xis, i)
  }

  # Shape parameter from GPD fit over varying k ----------------------------
  shapes <- c()
  for (i in 2:n) {
    q     <- 1 - i / n
    thres <- quantile(x, q, na.rm = TRUE)

    fit <- tryCatch(
      {
        evd::fpot(x, thres, start = list(scale = 5, shape = 1))$estimate[2]
      },
      error = function(e) {
        NA
      }
    )
    shapes[i - 1] <- fit
  }

  df <- data.frame(
    k     = seq(2, n, 1),
    hill  = hill,
    shape = shapes
  )

  y_max2 <- max(
    quantile(df$hill,  0.99, na.rm = TRUE),
    quantile(df$shape, 0.99, na.rm = TRUE)
  )

  # Plot Hill vs. k and GPD shape vs. k ------------------------------------
  p2 <- ggplot(df, aes(k, hill)) +
    geom_point() +
    geom_point(aes(k, shape), colour = "#1f78b4") +
    ylim(c(-0.1, y_max2)) +
    theme_bw() +
    ylab("") +
    theme(text = element_text(size = 35))

  # Repeat same diagnostics for shifted data x + shift ----------------------
  x_sh  <- x + shift
  xis_sh <- sort(x_sh, decreasing = TRUE)
  par_df_sh <- data.frame(pq = pareto_quantiles, x = log(xis_sh))

  y_max_sh <- par_df_sh$x[n]
  y_min_sh <- par_df_sh$x[1]

  p4 <- ggplot(par_df_sh, aes(x = pq, y = x)) +
    geom_point() +
    theme_bw() +
    theme(text = element_text(size = 35)) +
    xlab(expression(log((n + 1) / j))) +
    ylab("") +
    ylim(c(y_max_sh, y_min_sh))

  hill_sh <- c()
  for (i in 2:n) {
    hill_sh[i - 1] <- hoga_estimator_vec(1, x_sh, i)
  }

  df_sh <- data.frame(
    k     = seq(2, n, 1),
    hill  = hill_sh,
    shape = shapes
  )

  y_max_sh2 <- max(
    quantile(df_sh$hill,  0.99, na.rm = TRUE),
    quantile(df_sh$shape, 0.99, na.rm = TRUE)
  )

  p5 <- ggplot(df_sh, aes(k, hill)) +
    geom_point() +
    geom_point(aes(k, shape), colour = "#1f78b4") +
    ylim(c(-0.1, y_max_sh2)) +
    theme_bw() +
    ylab("") +
    theme(text = element_text(size = 35))

  # Arrange plots in a 2×2 grid --------------------------------------------
  pgrid <- grid.arrange(p1, p2, p4, p5, nrow = 2, ncol = 2)

  return(list(df = df, df_shift = df_sh, plot = pgrid))
}

#---------------------------------------------------------------------------
# Find “optimal” k for Hoga/Hill index test (sup-norm heuristic)
#---------------------------------------------------------------------------

find_optimal_k <- function(X) {
  # X: sample (assumed heavy-tailed, sorted internally)

  X_sorted <- sort(X, decreasing = TRUE)
  n <- length(X_sorted)

  # Search range for k
  k_min <- max(floor(0.05 * n), 50)
  k_max <- floor(n^0.8)

  # Objective: sup_j |X_(j+1) - model(j,k)| with model based on hoga estimator
  objective <- function(k) {
    gamma_hat <- hoga_estimator_vec(1, X_sorted, k)
    x_k       <- X_sorted[k + 1]
    diffs <- sapply(1:k_max, function(j) {
      var_hat <- (k / j)^gamma_hat * x_k
      abs(X_sorted[j + 1] - var_hat)
    })
    return(max(diffs))
  }

  k_vals   <- k_min:k_max
  obj_vals <- sapply(k_vals, objective)

  k_opt <- k_vals[which.min(obj_vals)]

  return(list(
    k_star         = k_opt,
    objective_value = min(obj_vals),
    k_grid         = k_vals,
    objective_grid = obj_vals
  ))
}

#---------------------------------------------------------------------------
# GPD_causal_tail_coeff (Pasche et al., 2023): unconditional CTC
#---------------------------------------------------------------------------

GPD_causal_tail_coeff <- function(
    X1, X2, k = NULL,
    threshold_q = 0.9,
    parametric_F1 = TRUE,
    both_tails = FALSE
) {
  # X1, X2: numeric vectors
  # k: tail sample size (if NULL, n^0.4)
  # threshold_q: quantile for POT threshold
  # parametric_F1: TRUE -> parametric for both X1 and X2; FALSE -> nonparametric X1
  # both_tails: not implemented

  n <- NROW(X1)
  if (is.null(k)) {
    k <- floor(n^0.4)
  }

  # Checks ------------------------------------------------------------------
  if (threshold_q >= 1 | threshold_q <= 0) {
    stop("threshold_q must be between 0 and 1.")
  }
  if (k <= 1 | k > floor((1 - threshold_q) * n)) {
    stop("k must be greater than 1 and smaller than the number of threshold exceedences.")
  }
  if (both_tails) {
    stop("both_tails not yet possible.")
  }

  r1 <- rank(X1, ties.method = "first")
  r2 <- rank(X2, ties.method = "first")

  # Fit GPD to X2 and compute F2(X2) ---------------------------------------
  u2   <- quantile(X2, threshold_q)
  fit2 <- fpot(X2, threshold = u2)
  FX2  <- pgpd(X2, loc = u2,
               scale = fit2$estimate[1],
               shape = fit2$estimate[2])
  FX2  <- (FX2 * (1 - threshold_q) + threshold_q) * (X2 >= u2) +
    r2 / n * (X2 < u2)

  if (parametric_F1) {
    # Parametric X1 as well
    u1   <- quantile(X1, threshold_q)
    fit1 <- fpot(X1, threshold = u1)
    FX1  <- pgpd(X1, loc = u1,
                 scale = fit1$estimate[1],
                 shape = fit1$estimate[2])
    FX1  <- (FX1 * (1 - threshold_q) + threshold_q) * (X1 >= u1) +
      r1 / n * (X1 < u1)

    # CTC (not necessarily exactly k exceedances)
    ctc <- mean(FX2[FX1 > (1 - k / n)])
  } else {
    # Nonparametric F1: use ranks for X1
    ctc <- 1 / k * sum(FX2[r1 > n - k])
  }

  return(ctc)
}

#---------------------------------------------------------------------------
# LGPD_causal_tail_coeff (Pasche et al., 2023): conditional on H
#---------------------------------------------------------------------------
LGPD_causal_tail_coeff <- function(X1, X2, H, k = NULL, threshold_q = 0.95,
                                   parametric_F1=TRUE, return_MLEs=FALSE,
                                   both_tails = FALSE, opposite=FALSE,
                                   method="average"){
  #Post-fit corrected H-conditional LGPD causal tail coefficient estimator

  n <- NROW(X1) # number of observations
  if(is.null(k)){
    k<-floor(n ^ 0.4)
  }
  #Checks
  if (threshold_q>=1 | threshold_q<=0) {
    stop("threshold_q must be between 0 and 1.")
  }
  if (k <= 1 | k>floor((1-threshold_q)*n) ) {
    stop("k must be greater than 1 and smaller than the number of threshold exceedences.")
  }
  if (both_tails) {
    stop("both_tails not yet possible.")
  }
  Hcs <- scale(H) #Centering and scaling H
  r1 <- rank(X1, ties.method = method) # ranks of X1
  r2 <- rank(X2, ties.method = method) # ranks of X2
  if(opposite){
    r1<-max(r1)-r1+1
  }
  #Fit POT GPD model and compute the fitted cdf in the data points for X2
  u2 <- quantile(X2,threshold_q)
  fit2 <- gpd.fit(X2,threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE)
  FX2 <- pgpd(X2, loc=u2, scale=pmax(fit2$mle[1]+fit2$mle[2]*Hcs,0.0001)[1:length(Hcs)], shape=fit2$mle[3])
  FX2 <- (FX2*(1-threshold_q)+threshold_q)*(X2>=u2) + r2/n * (X2<u2)
  #Compute the parametric causal tail coefficient
  if(parametric_F1){
    #Fit POT GPD model and compute the fitted cdf in the data points for X1
    u1 <- quantile(X1,threshold_q)
    fit1 <- gpd.fit(X1,threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE)
    FX1 <- pgpd(X1, loc=u1, scale=pmax(fit1$mle[1]+fit1$mle[2]*Hcs,0.0001)[1:length(Hcs)], shape=fit1$mle[3])
    FX1 <- (FX1*(1-threshold_q)+threshold_q)*(X1>=u1) + r1/n * (X1<u1)
    #Compute the parametric causal tail coefficient
    ctc <- mean(FX2[FX1 > (1-k/n)]) #not always exactly k excesses for parametric approach
    if(return_MLEs){
      return(list(ctc=ctc, F2_scale0=fit2$mle, F2_scale1=fit2$mle, F2_shape=fit2$mle,
                  se_F2_scale0=fit2$se, se_F2_scale1=fit2$se, se_F2_shape=fit2$se,
                  F1_scale0=fit1$mle, F1_scale1=fit1$mle, F1_shape=fit1$mle,
                  se_F1_scale0=fit1$se, se_F1_scale1=fit1$se, se_F1_shape=fit1$se))
    }
  }else{
    ctc <- 1/k * sum(FX2[r1 > n - k])
    if(return_MLEs){
      return(list(ctc=ctc, F2_scale0=fit2$mle, F2_scale1=fit2$mle, F2_shape=fit2$mle,
                  se_F2_scale0=fit2$se, se_F2_scale1=fit2$se, se_F2_shape=fit2$se,
                  F1_scale0=NULL, F1_scale1=NULL, F1_shape=NULL,
                  se_F1_scale0=NULL, se_F1_scale1=NULL, se_F1_shape=NULL))
    }
  }
  return(ctc)
}
#---------------------------------------------------------------------------
# Nonparametric CTC estimator (Gnecco et al., 2021), add of min for lower tail
#---------------------------------------------------------------------------

causal_tail_coeff_basic <- function(v1, v2, k = floor(n ^ 0.4), to_rank = TRUE,
                                    both_tails = TRUE, opposite = FALSE, min=FALSE, method="average"){
  # number of observations
  n <- NROW(v1)

  # check k
  if (k <= 1 | k >= n) {
    stop("k must be greater than 1 and smaller than n.")
  }

  # rank variables?
  if (to_rank){
    if(method=="average"){
      r1 <- rank(v1, ties.method = "average")
      r2 <- rank(v2, ties.method ="average")
    } else{
      r1 <- rank(v1, ties.method = "first")
      r2 <- rank(v2, ties.method ="first")
    }


  } else{
    r1 <- v1
    r2 <- v2
  }
  if(opposite){
    r2 <- max(r2)-r2+1
  }
  if(min){
    r1<-max(r1)-r1+1
    r2 <- max(r2)-r2+1
  }
  # compute causal tail coefficient
  if (both_tails){
    k <- (k %/% 2) * 2
    1 / (k * n) * sum(2 * abs(r2[r1 > n - k / 2 | r1 <= k / 2] - (n + 1) / 2))
  } else{
    #print(length(r2[r1 > n - k]))
    k1<-length(r2[r1 > n - k])
    1 / (k1 * n) * sum(r2[r1 > n - k])

  }
}
