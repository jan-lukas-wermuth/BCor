#' Odds Ratio
#'
#' `Oddsratio()` computes the oddsratio
#'
#' @param X a n x 1 numeric vector, matrix or data frame. If a contingency table is supplied, the upper left corner shall contain the joint success probability (or frequency). If a 3-dimensional vector of probabilities is supplied, the order (p, q, r) shall be respected.
#' @param Y NULL (default) or a n x 1 numeric vector, matrix or data frame and compatible dimensions to X.
#'
#' @return The value of the oddsratio
#' @export
#'
#' @examples
#' X <- matrix(c(10, 20, 30, 5), ncol = 2)
#' Oddsratio(X)
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
