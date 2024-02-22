#' Phi
#'
#' `Phi()` computes Phi and corresponding confidence intervals.
#'
#' @param X a n x 1 numeric vector, matrix or data frame.
#' @param Y a n x 1 numeric vector, matrix or data frame.
#' @param alpha confidence level for the returned confidence interval. NA yields Cole's correlation without confidence intervals.
#' @param covar data assumptions: iid ("iid"), heteroskedasticity ("HC") or heteroskedasticity and autocorrelation ("HAC").
#' @param n sample size. Only necessary if a contingency table of probabilities is provided.
#'
#' @return The value of Phi together with the specified confidence interval.
#' @export
#'
#' @examples
#' x <- matrix(c(10, 20, 30, 5), ncol = 2)
#' Phi(x)
Phi <- function (X, Y = NULL, alpha = 0.95, covar = "iid") {
  if (alpha == NA){
    if (!is.null(Y))
      X <- table(X, Y)
    stopifnot(prod(dim(X)) == 4 || length(X) == 4)
    a <- X[1, 1]
    b <- X[1, 2]
    c <- X[2, 1]
    d <- X[2, 2]
    Phi <- (a * d - b * c)/(sqrt((a + b) * (a + c) * (b + d) * (c + d)))
    res <- dplyr::tribble(~Phi
                          #--
                          Phi)
    return(res)
  }
  if (is.null(Y)){
    stopifnot(prod(dim(X)) == 4 || length(X) == 4)
    colnames(X) <- c(1,0)
    rownames(X) <- c(1,0)
    if (near(sum(X), 1)){
      X <- X * n
    }
    x <- as.numeric(epitools::expand.table(X)[,1]) - 1
    y <- as.numeric(epitools::expand.table(X)[,2]) - 1
  } else {
    x <- X
    y <- Y
  }

  if (length(x) != length(y)) {
    stop("X and Y must have the same length")
  }
  n <- length(x)

  if ((alpha < 0.5) | (alpha >= 1)) {
    warning("alpha must be between 0.5 and 1. The results are based on alpha=0.95 instead.")
    alpha <- 0.95
  }

  ### Compute the estimators
  df <- dplyr::tibble(X = x, Y = y) %>% dplyr::mutate(XY = x * y)
  df_means <- df %>% dplyr::summarize_all(mean)
  p_est <- df_means$X
  q_est <- df_means$Y
  r_est <- df_means$XY
  cov <- r_est - p_est * q_est

  Phi <- stats::cor(x, y)
  PhiZ <- atanh(Phi)

  ### Check whether the estimated probabilities are on the boundaries for inference
  boundary_eps <- 10 ^ (-8)
  boundary_conditions <- (any(c(p_est, q_est, r_est) < boundary_eps | c(p_est, q_est, r_est) > 1 - boundary_eps) | (r_est > min(p_est, q_est) - boundary_eps) | (r_est < p_est + q_est - 1 + boundary_eps))
  boundary_case <- F
  if (boundary_conditions) {
    boundary_case <- T
    warning("One of the probabilities is on a boundary: Do not report inference!")
  }

  # If we are on the boundary, do not report inference!
  if (boundary_case) {
    res <- dplyr::tribble(~Phi
                          #--
                          Phi)
    return(res)
  } else {
    # Estimate the covariance matrix Omega, either by HC, HAC or the sample cov
    if (covar == "iid") {
      Omega <- stats::var(df)
    } else if (covar == "HAC") {
      Omega <- n * sandwich::vcovHAC(stats::lm(as.matrix(df) ~ 1))
    } else if (covar == "HC") {
      Omega <- n * sandwich::vcovHC(stats::lm(as.matrix(df) ~ 1))
    } else stop("Please insert a valid option for computing the covariance matrix!")

  ### Inference
  L_est <- L(p_est, q_est, r_est)
  PhiZ_Var <- as.numeric(t(L_est) %*% Omega %*% L_est * 1 / (1 - Phi^2)^2) / n
  PhiZ_CI <- c(tanh(PhiZ + stats::qnorm((1 - alpha) / 2) * sqrt(PhiZ_Var)),
               tanh(PhiZ - stats::qnorm((1 - alpha) / 2) * sqrt(PhiZ_Var)))

  res <- dplyr::tribble(~Phi, ~CI_lower, ~CI_upper,
                 #--|--|--
                 Phi, PhiZ_CI[1], PhiZ_CI[2])
  }
  return(res)
}


#' @keywords internal
L <- function(p, q, r) {
  denom <- sqrt(p * (1 - p) * q * (1 - q))
  mult <- (p * q - r) / (2 * p * (1 - p) * q * (1 - q))
  L1 <- (-q + mult * (1 - p) * q * (1 - q) - p * q * (1 - q)) / denom
  L2 <- (-p + mult * p * (1 - p) * (1 - q) - p * (1 - p) * q) / denom
  L3 <- 1 / denom
  return(c(L1, L2, L3))
}

