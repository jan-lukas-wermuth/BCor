
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package BCor

## Overview

<!-- badges: start -->
<!-- badges: end -->

This package accompanies the paper *Measuring Dependence between Events*
by Marc-Oliver Pohle, Timo Dimitriadis and Jan-Lukas Wermuth.

It provides functions which compute estimators for correlation
coefficients that are specific to $2 \times 2$ contingency tables.
Additionally, the functions *YuleQ*, *Cole* and *Phi* yield the lower
and upper bounds of the respective confidence intervals. Furthermore, we
offer the corresponding inverses, i.e. functions that compute the
underlying $2 \times 2$ contingency table given the value of the
dependence measure and the marginal frequencies.

## Installation

To install the package from GitHub and load it, run the following `R`
code:

``` r
install.packages("devtools")
library(devtools)
install_github("jan-lukas-wermuth/BCor")
library(BCor)
```

## Example

We provide a simple example illustrating the use of correlation
coefficients for $2 \times 2$ contingency tables as well as their
inverses. Confidence intervals are only included for Phi, Cole and
Yule’s Q.

``` r
# matrix with x[1,1] = 10, x[2,1] = 20, x[1,2] = 30, x[2,2] = 5
x <- matrix(c(10, 20, 30, 5), ncol = 2) 

BCor::Phi(x)
#> # A tibble: 1 × 3
#>      Phi CI_lower CI_upper
#>    <dbl>    <dbl>    <dbl>
#> 1 -0.537   -0.711   -0.300
BCor::Cole(x)
#> # A tibble: 1 × 3
#>        C CI_lower CI_upper
#>    <dbl>    <dbl>    <dbl>
#> 1 -0.629   -0.843    -0.24
BCor::YuleQ(x)
#> # A tibble: 1 × 3
#>        Q CI_lower CI_upper
#>    <dbl>    <dbl>    <dbl>
#> 1 -0.846   -0.952   -0.559
BCor::YuleQ(x, g = 0.5)
#> # A tibble: 1 × 1
#>        Q
#>    <dbl>
#> 1 -0.552
BCor::Tetrachoric(x)
#> [1] -0.7559687
BCor::Oddsratio(x)
#> [1] 0.08333333

C <- 0.5
Q <- 0.5

# a 4 x 1 vector of marginals c(R1, R2, C1, C2)
m <- c(10, 40, 20, 30) 

BCor::Cole.inv(C, m)
#>      [,1] [,2]
#> [1,]    7    3
#> [2,]   13   27
BCor::YuleQ.inv(Q, m)
#>           [,1]      [,2]
#> [1,]  6.139991  3.860009
#> [2,] 13.860009 26.139991
```
