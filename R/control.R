
# Control functions --------

## Rejection --------

#' abc_control.rejection
#' @description Create list of control parameters for rejection ABC sampler
#' @export
#' @param n numeric number of proposals in output
#' @param epsilon numeric proposals corresponding to distances under this number are kept
#' @param delta numeric control number of new proposals to be resampled if distances are greater than epsilon
abc_control.rejection <- function(
  n = 1000,
  epsilon = 0.05,
  delta = 1e-3){

  stopifnot(n > 0 & n %% 1 == 0)
  stopifnot(epsilon >= 0 %% is.numeric(epsilon))
  stopifnot(delta >= 0 %% is.numeric(delta))

  return(as.list(environment(), all=TRUE))
}

abc_output.rejection <- function(
  include_dist = FALSE){

  return(as.list(environment(), all=TRUE))
}

## RABC -----------

#' abc_control.RABC
#' @description Create list of control parameters for Replenishment ABC sampler
#' @export
#' @param n numeric number of proposals in output
#' @param a numeric fraction of proposals to drop at each iteration
#' @param c1 controls number of MCMC steps to make at each iteration
#' @param R starting value of number of MCMC steps to make at each iteration
#' @param eps_final numeric algorithm will finish when proposals have distances below this number
#' @param pacc_final numeric algorithm will finish when acceptance proportion is below this number
#' @param d_eps_final numeric algorithm will stop when percent difference in epsilon between iterations is below this number
#' @param num_drop numeric number of proposals to drop at each iteration
#' @param num_keep numeric number of proposals to keep at each iteration
#' @param n_param numeric for internal use, don't change
#' @param cov_func function method for computing covariance
#' @param prior_eval evaluator of prior function
abc_control.RABC <- function(
  n = 1000,
  a = 0.5,
  c1 = 0.01,
  R = 10,
  eps_final = 0,
  pacc_final = 0.02,
  d_eps_final = 0,
  num_drop = ceiling(n * a),
  num_keep = n - num_drop,
  n_param = NA,
  cov_func = ifelse(n_param == 1, stats::var, stats::cov),
  prior_eval = function(x,y){return(1)}){

  stopifnot(n > 2 & n %% 1 == 0)
  stopifnot(a > 0 & a < 1)
  stopifnot(is.numeric(c1))
  stopifnot(R %% 1 == 0)
  stopifnot(eps_final >= 0 %% is.numeric(eps_final))
  stopifnot(d_eps_final >= 0 %% is.numeric(d_eps_final))
  stopifnot(num_drop %% 1 == 0)
  stopifnot(num_keep >= 2 && num_keep %% 1 == 0)
  stopifnot(num_keep == n - num_drop)

  return(as.list(environment(), all=TRUE))
}

abc_output.RABC <- function(
  include_dist = FALSE,
  print_output = TRUE){

  return(as.list(environment(), all=TRUE))
}

evidence_control <- function(
  n_param = NA,
  cov_func = ifelse(n_param == 1, stats::var, stats::cov)
  ){
  return(as.list(environment(), all = TRUE))
}



