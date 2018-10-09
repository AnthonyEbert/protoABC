
# Control functions --------

## Rejection --------

abc_control.rejection <- function(
  n = 1000,
  epsilon = 0.05,
  delta = 1e-3){

  return(as.list(environment(), all=TRUE))
}

abc_output.rejection <- function(
  include_dist = FALSE){

  return(as.list(environment(), all=TRUE))
}

## RABC -----------

abc_control.RABC <- function(
  n = 1000,
  a = 0.5,
  c1 = 0.01,
  R = 10,
  eps_final = 0,
  pacc_final = 0.02,
  d_eps_final = 0,
  num_drop = n * a,
  num_keep = n - num_drop,
  n_param = NA,
  cov_func = ifelse(n_param == 1, stats::var, stats::cov),
  prior_eval = function(x,y){return(1)}){

  return(as.list(environment(), all=TRUE))
}

abc_output.RABC <- function(
  include_dist = FALSE,
  print_output = TRUE){

  return(as.list(environment(), all=TRUE))
}
