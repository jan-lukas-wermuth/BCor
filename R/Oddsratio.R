#' Odds Ratio
#'
#' `Oddsratio()` computes the odds ratio.
#'
#' @param X a numeric vector, matrix or data frame with values of 0 or 1. Alternatively, a 2 x 2 matrix / contingency table can be supplied. In that case, the upper left corner shall contain the joint success probability (or frequency). If a 3-dimensional vector of probabilities is supplied, the order (p, q, r) shall be respected. p denotes the success probability of the row variable, q the success probability of the column variable and r the joint success probability.
#' @param Y NULL (default) or a n x 1 numeric vector, matrix or data frame with values of 0 or 1 and compatible dimensions to X.
#'
#' @return The value of the odds ratio.
#' @export
#'
#' @examples
#' # Insert a contingency table with frequencies in form of a matrix.
#' x <- matrix(c(10, 20, 30, 5), ncol = 2)
#' Oddsratio(x)
#'
#' # Insert a contingency table with relative frequencies in form of a matrix.
#' x <- matrix(c(0.2, 0.1, 0.4, 0.3), ncol = 2)
#' Oddsratio(x)
#'
#' # Insert two vectors of observations.
#' x <- c(0,1,1,1,1,0,1,0,0,1)
#' y <- c(1,0,1,1,0,0,0,0,1,0)
#' Oddsratio(x, y)
#'
#' # Insert two marginal success probabilities (p for the row variable and q for the column variable)
#' # and a joint success probability r.
#' p <- 0.6
#' q <- 0.3
#' r <- 0.2
#' Oddsratio(c(p, q, r))
Oddsratio <- function (X, Y = NULL) {
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
  result <- (a * d) / (b * c)
  return(result)
}
