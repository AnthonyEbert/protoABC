
<!-- README.md is generated from README.Rmd. Please edit that file -->
protoABC
========

The goal of protoABC is to provide an interface to ABC inference which is as flexible as possible.

Installation
------------

You can install protoABC from github with:

``` r
# install.packages("devtools")
devtools::install_github("AnthonyEbert/protoABC", auth_token = "0f7acf8a9c7faa1c678ce5fd8afb195badbca24b")
```

Example
-------

``` r

library(protoABC)

# Infer mean parameter for normal distribution from some observed data

distance_args <- rnorm(1000, 3, 1.5)

summary(distance_args)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  -1.878   1.941   2.929   2.957   3.957   7.998

prior <- function(n){data.frame(mu = rnorm(n, 5))}

distance <- function(theta, distance_args){
  sim <- rnorm(1000, theta)
  output <- abs(mean(sim) - mean(distance_args))
  return(output)
}

abc_post_1 <- abc_start(
  prior,
  distance,
  distance_args,
  method = "rejection",
  control = list(epsilon = 0.1)
)

summary(abc_post_1)
#>        mu       
#>  Min.   :2.792  
#>  1st Qu.:2.915  
#>  Median :2.969  
#>  Mean   :2.965  
#>  3rd Qu.:3.017  
#>  Max.   :3.134
```

``` r
# Infer mean and standard deviation for normal distribution given some observed data

prior <- function(n){
  data.frame(mu = runif(n, 2, 4), sd = rgamma(n, 1, 1))
}

prior_eval <- function(theta){
  prior_value <- dunif(theta["mu"], 2, 4) * dgamma(theta["sd"], 1, 1)
  
  return(prior_value)
}


distance <- function(theta, distance_args){
  sim <- rnorm(1000, theta["mu"], theta["sd"])
  output <- sqrt( (mean(sim) - mean(distance_args))^2 + (sd(sim) - sd(distance_args))^2)
  return(output)
}

abc_post_2 <- abc_start(
  prior,
  distance,
  distance_args,
  method = "RABC",
  control = list(n = 1000, prior_eval = prior_eval, pacc_final = 0.1), 
  output_control = list(print_output = FALSE)
)
#> Warning in rnorm(1000, theta["mu"], theta["sd"]): NAs produced

#> Warning in rnorm(1000, theta["mu"], theta["sd"]): NAs produced

summary(abc_post_2)
#>        mu              sd       
#>  Min.   :2.805   Min.   :1.437  
#>  1st Qu.:2.923   1st Qu.:1.521  
#>  Median :2.958   Median :1.545  
#>  Mean   :2.959   Mean   :1.546  
#>  3rd Qu.:2.993   3rd Qu.:1.569  
#>  Max.   :3.121   Max.   :1.666
```
