
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package BCor

## Overview

<!-- badges: start -->
<!-- badges: end -->

This package accompanies the paper *Measuring Dependence between Events*
by Marc-Oliver Pohle, Timo Dimitriadis and Jan-Lukas Wermuth.

It provides functions which compute estimators for correlation
coefficients that are specific to $2 \times 2$ contingency tables and
fulfill the characteristic of attainability, i.e. that attain the value
1 for perfectly positively dependent (comonotonic) Bernoulli
distributions and -1 for perfectly negatively dependent
(countermonotonic) Bernoulli distributions. Additionally, we offer the
corresponding inverses, i.e. functions that compute the underlying
$2 \times 2$ contingency table given the value of the dependence measure
and the marginal frequencies.

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

We provide a simple example illustrating the use of attainable
correlations for $2 \times 2$ contingency tables as well as their
inverses.

``` r
x <- matrix(c(10, 20, 30, 5), ncol = 2)

Cole(x)
#> # A tibble: 1 × 3
#>        C CI_lower CI_upper
#>    <dbl>    <dbl>    <dbl>
#> 1 -0.629    -0.81   -0.208
YuleQ(x)
#> # A tibble: 1 × 3
#>        Q CI_lower CI_upper
#>    <dbl>    <dbl>    <dbl>
#> 1 -0.846   -0.952   -0.559
Tetrachoric(x)
#> [1] -0.7559687

C <- 0.5
Q <- 0.5
m <- c(10, 40, 20, 30)

Cole.inv(C, m)
#>      [,1] [,2]
#> [1,]    7    3
#> [2,]   13   27
YuleQ.inv(Q, m)
#>           [,1]      [,2]
#> [1,]  6.139991  3.860009
#> [2,] 13.860009 26.139991
```
