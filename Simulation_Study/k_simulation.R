# Simulation Study to find the best k

# ------------------------------------------------------------------------------
# This file contains a modified version of code originally written by:
#   Nicola Gnecco, Nicolai Meinshausen, Jonas Peters, Sebastian Engelke (2022),
#   available at: https://github.com/nicolagnecco/causalXtreme
#
# Original license: GNU General Public License v3.0 (GPL-3.0)
#
# Modifications by: Lisa Leimenstoll, 2025
#   - add second tail index and pareto distribution in simulate_data function
#   - add pareto distribution in simulate_noise
#   - in simulation_0 change arguments, add second tail index
#   - add second tail index and pareto distribution in my_args
#   - change of plot labels and adjust for two tail indices
#   - add additional simualtion study for percentage of wrong causal inference 
#     between two variables


library(causalXtreme)
library(tidyverse)
library(doParallel)
library(doRNG)
library(rngtools)
library(tictoc)
library(EnvStats)
library(latex2exp)
###################################################

############## simulation functions ############################################
simulate_data <- function(n, p, prob_connect,
                          distr = c("student_t", "gaussian", "log_normal", "pareto"),
                          tail_index = 2.5,tail_index_2= 1.5, has_confounder = FALSE,
                          is_nonlinear = FALSE, has_uniform_margins = FALSE){
  
  if (p <= 1 | n <= 1){
    stop("n and p must be larger than 1!")
  }
  
  distr <- match.arg(distr)
  
  # Simulate random DAG
  dag <- random_dag(p = p, prob_connect = prob_connect,caus_order=1:p)
  
  # Add confounders to DAG
  if (has_confounder){
    ll <- add_random_confounders(dag, prob_confound = 2 / (3 * (p - 1)))
    dag <- ll$dag_confounders
    pos_confounders <- ll$pos_confounders
    p <- p + length(pos_confounders)
  } else{
    dag <- dag
    pos_confounders <- integer(0)
  }
  
  # Create random adjacency matrix
  adj_mat <- random_coeff(dag)
  
  # Simulate the noise variables
  noise <- cbind(simulate_noise(n, p/2, distr, tail_index),simulate_noise(n, p/2, distr, tail_index_2))
  
  # Simulate data
  if (is_nonlinear){
    dataset <- nonlinear_scm(adj_mat, noise)
  } else {
    dataset <- t(solve(diag(p) - t(adj_mat), t(noise)))
  }
  
  # Marginally transform each variable?
  if (has_uniform_margins){
    dataset <- apply(dataset, 2, uniform_margin)
  }
  
  # Remove confounders
  if (has_confounder){
    if (length(pos_confounders) > 0){
      dataset <- dataset[, -pos_confounders]
    }
  }
  
  # Return list
  ll <- list()
  ll$dataset <- dataset
  ll$dag <- dag
  ll$pos_confounders <- pos_confounders
  return(ll)
}

pick_elements <- function(vec, prob){
  r <- stats::rbinom(n = length(vec), size = 1, prob = prob)
  vec[r == 1]
}


#' Inverse mirror uniform
#'
#' Produces the quantile of a mirrored uniform distribution
#' associated to the probability \code{prob}. The mirrored uniform
#' distribution has support [-max, -min] \eqn{\cup} [min, max].
#'
#' @param prob Numeric --- between 0 and 1. The probability associated
#' to the desired quantile.
#' @param min Numeric --- above 0. The lower bound of the positive half of
#' the support.
#' @param max Numeric --- above 0. The upper bound of the positive half of
#' the support. Note that \code{max} must be strictly greater than
#' \code{min}.
#'
#' @return Numeric --- between \code{-max} and \code{max}. The quantile
#' associated to the probability \code{prob}.
#' @noRd
inverse_mirror_uniform <- function(prob, min, max){
  
  # check if max < min
  if (!(max > min)){
    stop(paste("The maximum value, max, must be strictly greater",
               "than the minimum value, min."))
  }
  
  # check if min > 0 and max > 0
  if (!(min > 0 & max > 0)){
    stop(paste("Both the minimum and the maximum values, min and max,",
               "must be positive."))
  }
  
  if (prob == 0){
    return(-max)
  } else if (prob > 0 & prob < 1 / 2) {
    return(2 * prob * (max - min) - max)
  } else if (prob >= 1 / 2 & prob < 1){
    return( (2 * prob - 1) * (max - min) + min)
  } else if (prob == 1){
    return(max)
  } else {
    stop("The probability, prob, must be between 0 and 1.")
  }
}


#' Sample from uniform family
#'
#' Sample n elements from a uniform distribution with lower and upper
#' bound equal to \code{min} and \code{max}, respectively.
#' If \code{mirror == TRUE}, the elements are sampled from a mirrored
#' uniform distribution, see \code{\link{inverse_mirror_uniform}}.
#'
#' @param n Integer. The number of elements to sample.
#' @param min Numeric. The lower bound of the distribution support.
#' If \code{mirror == TRUE}, this represents the lower bound of
#' the positive half of the support.
#' @param max Numeric. The upper bound of the distribution support.
#' If \code{mirror == TRUE}, this represents the upper bound of
#' the positive half of the support.
#' Note that \code{max} must be strictly greater than \code{min}.
#' @param mirror Logical. Should the elements be sampled from a
#' mirrored uniform distribution? If \code{mirror == TRUE},
#' \code{min} and \code{max} must be strictly greater than 0.
#'
#' @return Numeric vector. A vector with the sampled elements.
#' @noRd
sample_uniform <- function(n, min, max, mirror = FALSE){
  
  # check if max < min
  if (!(max > min)){
    stop(paste("The maximum value, max, must be strictly greater",
               "than the minimum value, min."))
  }
  
  if (mirror == TRUE){
    
    purrr::map_dbl(stats::runif(n),
                   inverse_mirror_uniform, min = min, max = max)
    
  } else {
    
    stats::runif(n = n, min = min, max = max)
  }
}


#' Simulate random DAG
#'
#' Simulates a directed acyclic graph (DAG) and returns its adjacency matrix.
#' Copyright (c) 2013 Jonas Peters \email{peters@@math.ku.dk}.
#' All rights reserved.
#'
#' @param p Integer --- greater than 0. Number of nodes.
#' @param prob_connect Numeric --- between 0 and 1. The probability that an edge
#' \eqn{i {\rightarrow} j} is added to the DAG.
#' @param caus_order Numeric vector. The causal order of the DAG.
#' If the argument is not provided it is generated randomly.
#'
#' @return Square binary matrix. A matrix representing the random DAG.
#' @keywords internal
random_dag <- function(p, prob_connect,
                       caus_order = sample(p, p, replace = FALSE)){
  
  # check inputs
  if (p < 1){
    stop("The number of nodes p must be greater than 0.")
  }
  
  if (!missing(caus_order)){
    if (length(caus_order) != p){
      stop(paste("The number of nodes p does not match with the length of",
                 "the causal order."))
    }
  }
  
  # get random dag
  dag <- matrix(0, nrow = p, ncol = p)
  
  if (p > 1){
    for (i in 1:(p - 1)){
      elms <- pick_elements(caus_order[(i + 1):p], prob_connect)
      dag[caus_order[i], elms] <- 1
    }
  }
  dag
}


#' Sample random coefficients
#'
#' Sample random coefficients from uniform distribution
#' for the given DAG \code{dag}.
#' Copyright (c) 2013 Jonas Peters \email{peters@@math.ku.dk}.
#' All rights reserved.
#'
#' @inheritParams compute_caus_order
#' @param lb Numeric. The lower bound of the coefficient that can be
#' sampled.
#' @param ub Numeric. The lower bound of the coefficient that can be
#' sampled. It must be stricly greater than \code{lb}.
#' @param two_intervals Logical. Should the coefficient be sampled
#' from two symmetric uniform distributions?
#' If \code{two_intervals == TRUE}, \code{lb} and \code{ub} must be
#' positive.
#' @return Square numeric matrix. The adjacency matrix of the underlying
#' DAG \code{dag}.
#' @noRd
random_coeff <- function(dag, lb = 0.1, ub = 0.9, two_intervals = TRUE){
  
  # check if dag is a (non-weighted) adjacency matrix
  if (!all(dag %in% c(0, 1))){
    stop("The entries of dag must be either 0 or 1.")
  }
  
  adj_mat <- matrix(0, nrow = NROW(dag), ncol = NCOL(dag))
  num_coeff <- sum(dag)
  adj_mat[dag == 1] <- sample_uniform(num_coeff, lb, ub, two_intervals)
  adj_mat
}



#' Simulate noise observations
#'
#' Sample \code{n} observations for \code{p} independent noise variables
#' from a certain distribution \code{distr}.
#' @param n Positive integer. The number of observations, must be larger
#' than 1.
#' @param p Positive integer. The number of variables, must be larger
#' than 1.
#' @param distr Character. The distribution of the noise. It is one of:
#' \itemize{
#' \item \code{student_t}, in this case the user has to specify the
#' \code{tail_index}, i.e., the degrees of freedom,
#' \item \code{gaussian},
#' \item \code{log_normal}.
#' }
#' @param tail_index Positive numeric. The tail index, i.e., degrees
#' of freedom, of the noise.
#' @return Numeric matrix. Dataset matrix with \code{n}
#' rows (observations) and \code{p} columns (variables).
#' @keywords internal
simulate_noise <- function(n, p, distr = c("student_t", "gaussian",
                                           "log_normal", "pareto"), tail_index){
  
  if (p < 1 | n <= 1){
    stop("n and p must be larger than 1!")
  }
  
  distr <- match.arg(distr)
  
  switch(distr,
         "student_t" = {
           
           noise <- array(stats::rt(n * p, df = tail_index), dim = c(n, p))
           
         },
         "gaussian" = {
           
           noise <- array(stats::rnorm(n * p), dim = c(n, p))
         },
         "log_normal" = {
           
           noise <- array(stats::rlnorm(n * p), dim = c(n, p))
           
         },
         "pareto" = {
           
           noise <- array(EnvStats::rpareto(n * p, location =1, shape=tail_index), dim = c(n, p))
           
         })
  
  return(noise)
}


simulation_settings <- function(){
  ## void -> tibble
  ## returns a tibble with simulation settings
  
  nexp <- 50
  my_args <- list(
    experiment = 1:nexp,
    n = c(5e2, 1e3, 1e4),
    p = c(4, 7, 10, 15, 20, 30, 50),
    distr = c('student_t'),
    tail_index = c(1.5, 2.5, 3.5),
    has_confounder = c(F, T),
    is_nonlinear = c(F, T),
    has_uniform_margins = c(F, T))
  
  my_args <- expand.grid(my_args, stringsAsFactors = F) %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(prob_connect = min(5/(p - 1), 1/2)) %>%
    filter(has_confounder + is_nonlinear + has_uniform_margins <= 1) %>%
    rowid_to_column("id") %>%
    select(id, everything())
  
  return(my_args)
}

method_settings <- function(){
  ## void -> tibble
  ## returns a tibble with all the used methods
  method_tbl <- tibble(method = c("ease", "ica_lingam",
                                  "direct_lingam", "pc_rank", "random"))
  
  return(method_tbl)
}

set_simulations <- function(experiment_ids=1:NROW(simulation_settings()),
                            method_nms=method_settings()$method,
                            seed){
  ## numeric_vector character_vector integer -> list
  ## prepares the settings for the simulations
  ## INPUTS:
  ## - experiment_ids: a numeric vector that specifies the rows of the tibble generated
  ## by calling the function simulation_settings().
  ## - method_ids: a character vector that specifies the method to use. This is
  ## a subset of methods from the output of method_settings().
  ## - seed: an integer seed that determines the outcome of the simulation.
  ##
  ## RETURNS:
  ## The function returns a list made of:
  ## - simulation_arguments: a tibble generated by the function simulation_settings(),
  ## where only the rows in experiments_ids are kept.
  ## - method_arguments: a character vectors generated by the function method_settings(),
  ## where only the methods in method_nms are kept.
  ## NOTE:
  ## The simulations are fully repeatable. That is, even after running the full
  ## simulation, you can repeat it over a subset of arguments/methods.
  ## This is possible because this function assigns to each simulation run
  ## and method a unique random seed. These random seeds are generated with
  ## L'Ecuyer RNG method and are independent of each other.
  
  # set simulation options
  simulation_arguments <- simulation_settings()
  m <- NROW(simulation_arguments)
  if (any(!(experiment_ids %in% 1:m))){
    stop(paste("Argument experiment_ids must contain integers between 1 and ",
               m, ".", sep = ""))
  }
  
  experiment_ids <- sort(unique(experiment_ids))
  
  # set methods
  method_arguments <- method_settings()
  k <- NROW(method_arguments)
  if (any(!(method_nms %in% method_arguments$method))){
    stop(paste("Argument method_nms must contains a subset of these methods: ",
               paste(method_arguments, collapse = ", "), ".", sep = ""))
  }
  
  # create independent RNG streams with L'Ecuyer method
  rng <- RNGseq(m + m * k, seed = seed)
  
  # select RNG streams for simulations
  sims_rng_stream <- rng[1:m]
  simulation_arguments$rng <- sims_rng_stream
  
  # select RNG streams for methods
  meth_rng_stream <- list()
  
  for (j in 1:k){
    # cat(j*m + experiment_ids, '\n')
    meth_rng_stream[[j]] <- rng[j * m + experiment_ids]
  }
  
  method_arguments$rng <- meth_rng_stream
  
  # filter simulations
  simulation_arguments <- simulation_arguments %>%
    filter(id %in% experiment_ids)
  
  # filter methods
  method_arguments <- method_arguments %>%
    filter(method %in% method_nms)
  
  # return list
  list(simulation_arguments=simulation_arguments,
       method_arguments=method_arguments)
  
}

get_method_args <- function(method = c("ease", "lingam", "ica_lingam",
                                       "direct_lingam",
                                       "pc", "pc_rank", "random"), n){
  ## character integer -> list
  ## produces a list with the arguments used by the method in the simulation
  method <- match.arg(method)
  switch(method,
         ease = {
           list(k = floor(n ** 0.4))
         },
         ica_lingam = {
           list(contrast_fun = "logcosh")
         },
         direct_lingam = {
           list()
         },
         pc_rank = {
           list(alpha = 5e-4)
         },
         random = {
           list()
         },
         lingam = {
           list(contrast_fun = "logcosh")
         })
}

causal_discovery_wrapper <- function(dat, method, set_args){
  time1 <- Sys.time()
  ll <- causal_discovery(dat, method, set_args)
  ll$time <- Sys.time() - time1
  
  return(ll)
}

wrapper_sim <- function(i, j, sims_args, method_args, trimdata = FALSE,
                        meth_args = FALSE){
  
  # Checks
  if (meth_args){
    if (!("set_args" %in% names(method_args))){
      stop("When meth_args = TRUE, the tibble `method_args` must contain the
           column `set_args`.")
    }
  } else {
    if ("set_args" %in% names(method_args)){
      stop("When meth_args = FALSE, the tibble `method_args` must *not* contain the
           column `set_args`.")
    }
  }
  
  # Set up index variables
  m <- nrow(sims_args)
  k <- nrow(method_args)
  cat("Simulation", i, "out of", m, "--- Inner iteration", j, "out of", k, "\n")
  
  # Simulate data
  current_exp <- sims_args[i, ]
  args_simulate <- current_exp %>% select(-experiment, -id, -rng)
  
  rng_sims <- current_exp$rng[[1]]
  rngtools::setRNG(rng_sims)
  X <- do.call(simulate_data, args_simulate)
  
  if (trimdata){
    X$dataset <- trim_data(X$dataset)
  }
  
  n <- NROW(X$dataset)
  
  
  # Run method
  method <- method_args$method[j]
  
  if (meth_args){
    set_args <- method_args$set_args[[j]]
  } else {
    set_args <- get_method_args(method, n)
  }
  
  rng_method <- method_args$rng[[j]][[i]]
  rngtools::setRNG(rng_method)
  
  out <- causal_discovery_wrapper(dat = X$dataset,
                                  method = method,
                                  set_args = set_args)
  
  # Evaluate algorithms
  algo_result <- causal_metrics(simulated_data = X,
                                estimated_graphs = out[1:2])
  
  # Return tibble
  res <- current_exp %>% select(id) %>%
    bind_cols(tibble(method = method)) %>%
    bind_cols(tibble(sid = algo_result$sid,
                     shd = algo_result$shd,
                     time = out$time))
  
  if (meth_args){
    # create tibble with argument name-value pair
    if (!length(set_args)){
      args_tbl <- tibble(arg_name = "", arg_value = NA)
    } else {
      args_tbl <- tibble(arg_name = names(set_args)[1],
                         arg_value = set_args[[1]])
    }
    
    res <- res %>%
      bind_cols(args_tbl)
  }
  
  # return value
  return(res)
}

trim_data <- function(data){
  ## numeric_matrix -> numeric_matrix
  ## keeps observations in the tails
  
  n <- nrow(data)
  p <- ncol(data)
  data_trimmed <- data %>%
    as_tibble() %>%
    mutate(id = 1:n()) %>%
    gather(key = "variab",
           value = "value",
           -id) %>%
    group_by(variab) %>%
    mutate(uq = quantile(value, .9),
           lq = quantile(value, .1),
           bulk = value > lq & value < uq,
           tail = !bulk) %>%
    select(id, variab, tail) %>%
    ungroup() %>%
    spread(variab, tail)
  
  is_in_tail <- (apply(data_trimmed[, -1], 1, sum) >= sqrt(p))
  is_in_bulk <- !(is_in_tail)
  
  c(sum(is_in_bulk)/n)
  
  tail_id <- which(is_in_tail)
  bulk_id <- which(is_in_bulk)
  bulk_id <- sample(x = bulk_id, size = length(bulk_id) / 2)
  
  ids <- sample(c(tail_id))
  
  data[ids, ]
}
#################################### simulation_0 #############################


simulation_0 <- function(log_file,
                         result_file,
                         is_demo = FALSE){
  
  
  ## ROBUSTNESS OF K ####
  set.seed(1436) # Gutenberg's press
  
  # simulation settings
  nexp <- 80
  my_args <- list(
    experiment = 1:nexp,
    n = c(5e2, 1e3, 1e4),
    p = c(4,8, 10,16, 20,30, 50),
    distr = c("student_t"),#c("pareto"),
    tail_index = c(2,3,4),
    tail_index_2 = c(2,3,4),
    has_confounder = c(F),
    is_nonlinear = c(F),
    has_uniform_margins = c(F))
  
  n <- c("500" = 5e2, "1000" = 1e3, "10000" = 1e4)
  root <- seq(.2, .69, by = 5e-2)
  n_root <- as_tibble(cbind(sapply(n, function(n) floor(n ** root)), root)) %>%
    gather('500', '1000', '10000', key = 'n', value = 'k') %>%
    mutate(n = as.numeric(n))
  
  
  my_args <- expand.grid(my_args, stringsAsFactors = F) %>%
    as_tibble() %>%
    left_join(n_root, by = 'n') %>%
    rowwise() %>%
    mutate(prob_connect = min(5/(p - 1), 1/2)) %>%
    filter(has_confounder + is_nonlinear + has_uniform_margins <= 1) %>%
    rowid_to_column("id")
  
  niter <- NROW(my_args)
  
  if (is_demo){
    niter <- 100
  }
  
  # create clusters
  cores <-  if(Sys.info()["user"] == "elvis"){
    detectCores()
  } else {
    detectCores() - 1
  }
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl)
  
  # update each worker's environment
  clusterEvalQ(cl, {
    #library(causalXtreme)
    library(doParallel)
    library(tidyverse)
    source("Simulation_Study/simulation_functions.R")
  }
  )
  
  
  # Loop through all simulations
  tic()
  sink(file = log_file)
  cat("**** Simulation 0 **** \n")
  
  ll <- foreach(i = 1:niter, .combine = rbind) %dorng% {
    cat("Experiment number", i, "out of", niter, "\n",
        file = log_file, append = TRUE)
    
    # Generate data
    current_exp <- my_args[i, ]
    args_simulate <- current_exp %>% select(-experiment, -id, -root, -k)
    X <- do.call(simulate_data, args_simulate)
    
    # Methods arguments
    args_methods <- tibble(method = c("ease"),
                           set_args = list(
                             list(k = current_exp$k)))
    
    # Run algorithms
    out <- pmap(args_methods, causal_discovery_wrapper, dat = X$dataset)
    time_elapsed <- map(out, "time") %>% unlist()
    algo_results <- map(out, magrittr::extract, c("est_g", "est_cpdag"))
    
    # Evaluate algorithms
    algo_evaluations <- algo_results %>%
      map_dfr(causal_metrics, simulated_data = X)
    
    # Collect results
    res <- algo_evaluations %>%
      mutate(time = time_elapsed) %>%
      mutate(method = c("ease")) %>%
      mutate(id = current_exp$id) %>%
      left_join(current_exp, by = "id") %>%
      group_by_at(vars(-method, -time, -sid, -shd)) %>%
      nest()
    
    
  }
  sink()
  toc()
  stopCluster(cl)
  closeAllConnections()
  
  # save results
  saveRDS(ll, file = result_file)
  
}


##################### Create data #############################################

setwd()

file.name <- paste("Simulation_Study/output/k_sim_SID", ".rds", sep = "")

simulation_0(log_file = "Simulation_Study/output/sims_0.txt", 
             result_file = file.name,
             is_demo = FALSE)

###################### Create plot ###########################################
# produce chart
sim0_file = "Simulation_Study/output/k_sim_SID.rds"


# Set plotting theme
theme_set(theme_bw() +
            theme(plot.background = element_blank(),
                  legend.background = element_blank()))

### CONSTANTS ####
SIMULATION_K <- sim0_file

SID_KVARYING <- "Simulation_Study/pics/k_sim_SID.pdf"
tolPalette <- c(tolBlue = "#4477AA",
                tolRed = "#EE6677",
                tolGreen = "#228833",
                tolYellow = "#CCBB44",
                tolCyan = "#66CCEE",
                tolPurple = "#AA3377",
                tolGrey = "#BBBBBB",
                tolLightg ="lightgreen",
                tolVio = "violet")

tolPalette <- c(
  # Group 1: Blues/Cyans (cool family)
  "#1f78b4",
  "#66c2e0",
  "#a6cee3",
  "#ff8c42",
  "orange",
  "#fdbf6f",
  "#8E4585",
  "#B78BCF",
  "#D4BDE4"
)


## Robustness of k across different settings wrt SID ----
# Import data
dat <- read_rds(SIMULATION_K) %>% unnest(cols = c(data)) %>%
  filter(has_confounder == F, is_nonlinear == F,
         has_uniform_margins == F)

res <- dat %>%
  group_by(tail_index, root, tail_index_2) %>%
  summarise(mean_sid = mean(sid),
            N = n(),
            se = sd(sid)/sqrt(N)) %>%
  ungroup() %>%
  mutate(tail_diff = as.factor(paste0(tail_index,"_",tail_index_2)))%>%
  mutate(tail_index = as.factor(tail_index)) %>%
  mutate(tail_index_2 = as.factor(tail_index_2))

res$tail_diff<-paste0(substr(res$tail_diff,1,1),",",substr(res$tail_diff,3,3))
pd <- position_dodge(0.00)
desired_order<-c("2,2" ,"3,2", "4,2", "2,3", "3,3", "4,3", "2,4" ,"3,4" ,"4,4")
res$tail_diff <- factor(res$tail_diff, levels = desired_order)


g <- ggplot(res) +
  geom_line(aes(x = root, y = mean_sid, color = tail_diff),alpha = .7, 
            size = 1.2, position = pd) +
  geom_point(aes(x = root, y = mean_sid, color = tail_diff, fill=tail_diff),
             size = 2.5, 
             position = pd) +  
  scale_color_manual(values = unname(tolPalette)) +
  scale_fill_manual(values = unname(tolPalette)) +
  xlab(TeX("Fractional exponent $\\nu$ of $k=n^\\nu")) +
  ylab("Structural Intervention Distance") +
  labs(color = TeX("Tail index $\\alpha_1$,$\\alpha_2$")) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)
  ) + guides(fill=FALSE); g

ggsave(SID_KVARYING, g, width = 7.5, height = 5, units = c("in"))


######## Percentage of wrong causal inference between two variables ###########

# Grid of exponents ν for k = n^ν -------------------------------------------
nus   <- c(seq(.2, .69, by = 5e-2))  # ν ∈ [0.2, 0.69] in steps of 0.05
ns    <- c(500, 1000, 10000)         # sample sizes
times <- 10000                       # Monte Carlo repetitions per (n,ν,a,b)

# Tail indices (via df of t-distribution) for X and Y ------------------------
as <- c(2, 3, 4)  # df for X
bs <- c(2, 3, 4)  # df for Y

# Result matrix: one row per (ν,a,b), storing count of wrong directions ------
ffs <- matrix(nrow = length(nus) * length(as) * length(bs), ncol = 4)
r   <- 1  # row index in ffs

for (nu in nus) {
  print(nu)
  for (a in as) {
    for (b in bs) {
      f <- 0  # counter for wrong causal assignments over all ns and runs
      for (n in ns) {
        print(n)
        for (t in 1:times) {
          k1   <- n^nu                  # tail threshold k = n^ν
          beta <- runif(1, 0.1, 0.9)    # strength of linear dependence
          
          # Generate heavy-tailed X, Y with linear dependence Y = beta*X + noise
          x <- rt(n, a)
          y <- beta * x + rt(n, b)
          
          # Estimate causal tail coefficients in both directions -------------
          c1 <- causal_tail_coeff_basic(x, y, k = k1, both_tails = FALSE) # X -> Y
          c2 <- causal_tail_coeff_basic(y, x, k = k1, both_tails = FALSE) # Y -> X
          
          # If CTC suggests Y -> X (c2 > c1), this is a wrong direction here --
          if (c2 > c1) {
            f <- f + 1
          }
        }
      }
      
      # Store results for this (ν, a, b) -------------------------------------
      ffs[r, 1] <- nu
      ffs[r, 2] <- a
      ffs[r, 3] <- b
      ffs[r, 4] <- f     # total number of wrong decisions
      
      r <- r + 1
    }
  }
}

# Convert to data frame and normalise by number of runs ----------------------
ff <- as.data.frame(ffs)
colnames(ff) <- c("root", "tail_index", "tail_index_2", "error")

# Combine tail indices into a single label for plotting
ff$tail_diff <- paste0(ff$tail_index, ",", ff$tail_index_2)

# mean_sid: proportion of wrong causal inference over all n and repetitions
ff$error <- ff$error / (length(ns) * times)

write.csv(ff, file = "Simulation_Study/output/k_sim_wrong_percentage.csv",
          row.names = FALSE)

############# Create plot ##################################################

wrong_dir_data <- read.csv("Simulation_Study/output/k_sim_wrong_percentage.csv")

g <- ggplot(wrong_dir_data) +
  geom_line(
    aes(x = root, y = error, color = tail_diff),
    alpha = .7, size = 1.2, position = pd
  ) +
  geom_point(
    aes(x = root, y = error, color = tail_diff, fill = tail_diff),
    size = 2.5, position = pd
  ) +
  # Custom color palette for tail-index combinations (tolPalette must exist)
  scale_color_manual(values = unname(tolPalette)) +
  scale_fill_manual(values = unname(tolPalette)) +
  xlab(TeX("Fractional exponent $\\nu$ of $k=n^\\nu$")) +
  ylab("Percentage of wrong causal inference") +
  labs(color = TeX("Tail index $\\alpha_1,\\alpha_2$")) +
  theme(
    axis.text  = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  guides(fill = FALSE)

g

ggsave("Simulation_Study/pics/k_sim_wrong_percentage.pdf", g,
       width = 7.5, height = 5, units = "in")
