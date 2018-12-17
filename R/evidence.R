
#' @export
evidence <- function(prior_eval, kernel, posterior_sample, kern_args, control = list()){

  cl <- NULL

  parallel <- ifelse(is.null(cl), 1,
                     ifelse("cluster" %in% class(cl), 2,
                            ifelse(cl == "mclapply", 3,
                                   ifelse(cl == "test", 4, 5)))
  )

  lfunc <- make_lfunc(parallel, cl)

  control$n_param <- dim(posterior_sample)[2]
  control <- do.call("evidence_control", control)

  rw_cov <- control$cov_func(posterior_sample)

  rmulti <- function(mean, sigma){
    param = mvtnorm::rmvnorm(mean = mean, sigma, n = 1)
    value = mvtnorm::dmvnorm(param, mean = mean, sigma)
    output <- c(param, value)
  }

  param <- apply(
    posterior_sample,
    1,
    rmulti,
    sigma = rw_cov
  )

  param_s <- matrix(t(param), ncol = control$n_param + 1)

  param <- param_s[,-I(control$n_param + 1), drop = FALSE]

  kern_values <- lfunc(as.matrix(param), 1, kernel, kern_args)

  output <- apply(param, 1, prior_eval) * kern_values / mean(param_s[,control$n_param + 1])


  return(output)
}
