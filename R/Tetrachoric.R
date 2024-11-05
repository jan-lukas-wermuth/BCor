#' Tetrachoric Correlation
#'
#' `Tetrachoric()` computes the tetrachoric correlation.
#'
#' @param X a numeric vector, matrix or data frame with values of 0 or 1. Alternatively, a 2 x 2 matrix / contingency table can be supplied. In that case, the upper left corner shall contain the joint success probability (or frequency). If a 3-dimensional vector of probabilities is supplied, the order (p, q, r) shall be respected. p denotes the success probability of the row variable, q the success probability of the column variable and r the joint success probability.
#' @param Y NULL (default) or a n x 1 numeric vector, matrix or data frame with values of 0 or 1 and compatible dimensions to X.
#'
#' @return The value of the tetrachoric correlation.
#' @export
#'
#' @references
#' - \insertRef{pohle2024measuringdependenceevents}{BCor}
#' - \insertRef{pearson1900}{BCor}
#'
#' @examples
#' # Insert a contingency table with frequencies in form of a matrix.
#' x <- matrix(c(10, 20, 30, 5), ncol = 2)
#' Tetrachoric(x)
#'
#' # Insert a contingency table with relative frequencies in form of a matrix.
#' x <- matrix(c(0.2, 0.1, 0.4, 0.3), ncol = 2)
#' Tetrachoric(x)
#'
#' # Insert two vectors of observations.
#' x <- c(0,1,1,1,1,0,1,0,0,1)
#' y <- c(1,0,1,1,0,0,0,0,1,0)
#' Tetrachoric(x, y)
#'
#' # Insert two marginal success probabilities (p for the row variable and q for the column variable)
#' # and a joint success probability r.
#' p <- 0.6
#' q <- 0.3
#' r <- 0.2
#' Tetrachoric(c(p, q, r))
Tetrachoric <- function (X, Y = NULL) {
  if (!is.null(Y))
    X <- table(X, Y)
  if(is.null(Y) && length(X) == 3){
    p <- X[1]
    q <- X[2]
    r <- X[3]
    X <- matrix(NA, ncol = 2, nrow = 2)
    X[1,1] <- r
    X[1,2] <- p - r
    X[2,1] <- q - r
    X[2,2] <- 1 - p - q + r
  }
  stopifnot(prod(dim(X)) == 4 || length(X) == 4)
  a <- X[1, 1]
  b <- X[1, 2]
  c <- X[2, 1]
  d <- X[2, 2]

  # Check whether we have co- or countermonotonicity and report boundary value
  if (any(dplyr::near(X / (a + b + c + d), matrix(c(0, 0, 0, 0), ncol = 2)))) {
    if (dplyr::near(a / (a + b + c + d), 0) | dplyr::near(d / (a + b + c + d), 0)) {
      return(-1)
    }
    else if (dplyr::near(b / (a + b + c + d), 0) |
             dplyr::near(c / (a + b + c + d), 0)) {
      return(1)
    }
  }
  # Compute tetrachoric correlation implicitly via the standard definition
  else {
    f <-
      function(rho) {
        a / (a + b + c + d) - mvtnorm::pmvnorm(
          lower = c(stats::qnorm(1 - (a + b) / (a + b + c + d)), stats::qnorm(1 - (a + c) / (a +
                                                                                               b + c + d))),
          upper = c(Inf, Inf),
          mean = c(0, 0),
          corr = matrix(c(1, rho, rho, 1), ncol = 2)
        )
      }
    result <- stats::uniroot(f, c(-1, 1), tol = .Machine$double.eps ^ 0.5)
    return(result$root)
  }
}
