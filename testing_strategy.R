# Strategy to test for causal tail direction and confounding ----------------
testing_strategy <- function(X1, X2, H = NULL, k = floor(NROW(X1)^0.4)) {
  # X1, X2 : numeric vectors (two variables whose tail-causality we test)
  # H      : optional confounder. If NULL → unconditional; else H-conditional.
  # k      : tail sample size (default ~ n^0.4)
  #
  # Returns (unconditional case):
  #   c(dir, pretest, conftest, indtest, c12, c21, Ti[2], Ti[3], T2)
  #   where:
  #     dir      = 1 if c12 > c21 (X1 → X2), 2 otherwise (after swap)
  #     pretest  = 1 if Pre-Test is significant at 5%
  #     conftest = 1 if Confounder-Test rejects
  #     indtest  = 1 if Hoga index test rejects equal tails
  #     c12, c21 = estimated CTCs
  #     Ti[2:3]  = Hoga tail index estimates for X1 and X2
  #     T2       = p-value of Hill-based z-test
  #
  # Conditional (H ≠ NULL) case returns:
  #   c(dir, pretest, pretest2, conftest, c12, c21)
  
  pretest  <- 0   # 0/1: significance of permutation CTC pre-test
  conftest <- 0   # 0/1: confounding present (both tails > 0.5)
  indtest  <- 0   # 0/1: Hoga index test rejects equal tail index
  pretest_conf <- c()
  
  # -------------------------------------------------------------------------
  # Case 1: No confounder (unconditional CTC)
  # -------------------------------------------------------------------------
  if (is.null(H)) {
    
    # Nonparametric causal tail coefficients in both directions -------------
    c12 <- causal_tail_coeff_basic(X1, X2, k = k, both_tails = FALSE)
    c21 <- causal_tail_coeff_basic(X2, X1, k = k, both_tails = FALSE)
    
    # Decide "direction" by comparing CTCs: larger ⇒ more causal influence --
    if (c12 > c21) {
      dir <- 1            # X1 → X2
    } else {
      dir <- 2            # X2 → X1, so we swap to standardize orientation
      tmp <- X1
      X1  <- X2
      X2  <- tmp
    }
    
    # Pre-Test for causal asymmetry -----------------------------
    p_val <- CTC_causality_permutation_test(X1, X2, k = k, R = 2000)$Pmc
    if (p_val < 0.05) {
      pretest <- 1
    }
    
    # Normal-approximation test of CTC against 0.5 --------------------------
    tc1 <- sqrt(k) * (c12 - 0.5)
    tc2 <- sqrt(k) * (c21 - 0.5)
    alpha <- 0.05
    # Under null, asymptotic variance ~ 1/12
    q1 <- qnorm(1 - alpha / 2, 0, sqrt(1 / 12))
    q2 <- qnorm(alpha / 2,      0, sqrt(1 / 12))
    
    # pretest2 = 1 if either CTC significantly differs from 0.5 -------------
    pretest2 <- (tc1 < q2 | tc1 > q1) | (tc2 < q2 | tc2 > q1)
    
    # Confounding test: is min(c12, c21) significantly above 0.5? -----------
    Tc <- sqrt(k) * (min(c12, c21) - 0.5)
    if (Tc > qnorm(0.95, 0, sqrt(1 / 12))) {
      conftest <- 1
    }
    
    # Tail index equality test (Hoga index-test) ----------------------------
    # k_opt could be chosen via find_optimal_k; here fixed at 500.
    # k_opt_1 <- find_optimal_k(X1)$k_star
    # k_opt_2 <- find_optimal_k(X2)$k_star
    k_opt <- 500  # min(k_opt_1, k_opt_2)
    
    Ti <- tail_test(X1, X2, k_opt, 0.2)
    # Ti[1] = test statistic; threshold 55.44 from Hoga
    if (Ti[1] > 55.44) {
      indtest <- 1
    }
    
    # Hill-based test p-value for tail index equality -----------------------
    T2 <- hill_test(X1, X2, 500)$p_value
    
    return(c(
      dir, pretest, conftest, indtest,
      c12, c21, Ti[2], Ti[3], T2
    ))
    
  } else {
    # -----------------------------------------------------------------------
    # Case 2: With confounder H (conditional LGPD CTC)
    # -----------------------------------------------------------------------
    
    # Conditional CTC with permutation test ---------------------------------
    pretest_conf <- CTC_causality_permutation_test(
      X1, X2, H = H, k = k, R = 2000
    )
    c12 <- pretest_conf$ctc12
    c21 <- pretest_conf$ctc21
    
    # Direction via comparison of conditional CTCs --------------------------
    if (c12 > c21) {
      dir <- 1
    } else {
      dir <- 2
      tmp <- X1
      X1  <- X2
      X2  <- tmp
    }
    
    # Permutation-based pretest for asymmetry -------------------------------
    if (pretest_conf$Pmc < 0.05) {
      pretest <- 1
    }
    
    # Confounding-Test in the conditional setting  --------------
    Tc <- sqrt(k) * (min(c12, c21) - 0.5)
    if (Tc > qnorm(0.95, 0, sqrt(1 / 12))) {
      conftest <- 1
    }
    
    return(c(dir, pretest, conftest, c12, c21))
  }
}