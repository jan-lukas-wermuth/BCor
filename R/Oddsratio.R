#' Odds Ratio
#'
#' `Oddsratio()` computes the oddsratio
#'
#' @param X a n x 1 numeric vector, matrix or data frame.
#' @param Y a n x 1 numeric vector, matrix or data frame.
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
  stopifnot(prod(dim(X)) == 4 || length(X) == 4)
  a <- X[1, 1]
  b <- X[1, 2]
  c <- X[2, 1]
  d <- X[2, 2]
  result <- (a * d) / (b * c)
  return(result)
}
