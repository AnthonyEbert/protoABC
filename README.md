
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

# Infer mean parameter for normal distribution given that observed sample average is 3

prior <- function(n){data.frame(mu = rnorm(n, 5))}

distance <- function(theta){
  data <- rnorm(1000, theta)
  output <- abs(mean(data) - 3)
  return(output)
}

abc_post_1 <- abc_start(
  prior,
  distance,
  method = "rejection",
  control = list(epsilon = 0.1)
)

summary(abc_post_1)
#>        mu       
#>  Min.   :2.823  
#>  1st Qu.:2.961  
#>  Median :3.013  
#>  Mean   :3.010  
#>  3rd Qu.:3.063  
#>  Max.   :3.159

# Infer mean and standard deviation for normal distribution given that observed sample average is 3 and observed standard deviation estimate is 1

prior <- function(n){
  data.frame(mu = runif(n, 2, 4), sdp = rgamma(n, 1, 1))
}

distance <- function(theta){
  data <- rnorm(1000, theta[1], theta[2])
  output <- sqrt( (mean(data) - 3)^2 + (sd(data) - 1)^2)
  return(output)
}

abc_post_2 <- abc_start(
  prior,
  distance,
  method = "rejection",
  control = list(epsilon = 1)
)

summary(abc_post_2)
#>        mu             sdp          
#>  Min.   :2.021   Min.   :0.005847  
#>  1st Qu.:2.624   1st Qu.:0.382364  
#>  Median :3.008   Median :0.694683  
#>  Mean   :3.020   Mean   :0.767758  
#>  3rd Qu.:3.418   3rd Qu.:1.083556  
#>  Max.   :3.976   Max.   :1.991111
```
