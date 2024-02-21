#' Tetrachoric Correlation
#'
#' `Tetrachoric()` computes the tetrachoric correlation.
#'
#' @param X a n x 1 numeric vector, matrix or data frame.
#' @param Y a n x 1 numeric vector, matrix or data frame.
#'
#' @return The value of the tetrachoric correlation.
#' @export
#'
#' @examples
#' X <- matrix(c(10, 20, 30, 5), ncol = 2)
#' Tetrachoric(X)
Tetrachoric <- function (X, Y = NULL) {
  if (!is.null(Y))
    X <- table(X, Y)
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
