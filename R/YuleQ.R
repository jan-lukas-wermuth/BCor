#' Yule's Q
#'
#' `YuleQ()` computes Yule's Q and corresponding confidence intervals.
#'
#' @param X a n x 1 numeric vector, matrix or data frame with values of 0 or 1.
#' @param Y a n x 1 numeric vector, matrix or data frame with values of 0 or 1.
#' @param g a scalar between 0 and 1.
#' @param alpha confidence level for the returned confidence interval. FALSE yields Yule's Q without confidence intervals.
#' @param Fisher Indicator whether confidence intervals should be computed by using the Fisher Transformation. Default is TRUE.
#' @param covar data assumptions: iid ("iid"), heteroskedasticity ("HC") or heteroskedasticity and autocorrelation ("HAC").
#' @param n sample size. Only necessary if a contingency table of probabilities is provided and confidence intervals are desired.
#'
#' @return The value of Yule's Q or any of its related measures. g = 0.5 yields Yule's Y and g = 0.75 yields Digby's H.
#' @export
#'
#' @examples
#' x <- matrix(c(10, 20, 30, 5), ncol = 2)
#' YuleQ(x)
YuleQ <- function (X, Y = NULL, g = 1, alpha = 0.95, Fisher = TRUE, covar = "iid", n) {
  if (isFALSE(alpha)){
    if (!is.null(Y)){
      X <- table(X, Y)
    }
    stopifnot(prod(dim(X)) == 4 || length(X) == 4)
    a <- as.numeric(X[1, 1])
    b <- as.numeric(X[1, 2])
    c <- as.numeric(X[2, 1])
    d <- as.numeric(X[2, 2])
    Q <- (a * d - b * c)/(a * d + b * c)
    res <- dplyr::tribble(~Q,
                          #--
                          Q)
    return(res)
  }
  if (is.null(Y)){
    stopifnot(prod(dim(X)) == 4 || length(X) == 4)
    colnames(X) <- c(1,0)
    rownames(X) <- c(1,0)
    if (dplyr::near(sum(X), 1)){
      X <- X * n
    }
    x <- as.numeric(epitools::expand.table(X)[,1]) - 1
    y <- as.numeric(epitools::expand.table(X)[,2]) - 1
  } else {
    x <- X
    y <- Y
  }

  if (isFALSE(all(unique(x) %in% c(0,1))) || isFALSE(all(unique(y) %in% c(0,1)))){
    stop("Please insert data with values of 0 or 1!")
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

  Q <- cov / (r_est * (1 - p_est - q_est + r_est) + (q_est - r_est) * (p_est - r_est))
  QZ <- atanh(Q)

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
    res <- dplyr::tribble(~Q,
                          #--
                          Q)
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

  ### Inference for Yule's Q with Fisher transformation
  if (isTRUE(Fisher)){
    H_est <- H(p_est, q_est, r_est)
    QZ_Var <- as.numeric(t(H_est) %*% Omega %*% H_est) / n
    QZ_CI <- c(tanh(QZ + stats::qnorm((1 - alpha) / 2) * sqrt(QZ_Var)),
              tanh(QZ - stats::qnorm((1 - alpha) / 2) * sqrt(QZ_Var)))
    res <- dplyr::tribble(~Q, ~CI_lower, ~CI_upper,
                  #--|--|--
                  Q, QZ_CI[1], QZ_CI[2])
  }
  ### Inference for Yule's Q without Fisher transformation
  else {
    G_est <- G(p_est, q_est, r_est)
    Q_Var <- as.numeric(t(G_est) %*% Omega %*% G_est) / n
    # Confidence Interval for Q
    Q_CI <- c(Q + stats::qnorm((1 - alpha) / 2) * sqrt(Q_Var),
              Q - stats::qnorm((1 - alpha) / 2) * sqrt(Q_Var))
    res <- dplyr::tribble(~Q, ~CI_lower, ~CI_upper,
                          #--|--|--
                          Q, Q_CI[1], Q_CI[2])
  }
  }
  return(res)
}

# Helper functions for the Jacobians used for the Delta-Method
# Yule's Q with Fisher transformation
#' @keywords internal
H <- function(p, q, r) {
  denom1 <- -2 * (p - r) * (q - r)
  denom2 <- -2 * (1 - p - q + r)
  H1 <- (q - r) / denom1 + 1 / denom2
  H2 <- (p - r) / denom1 + 1 / denom2
  H3 <- (2 * r - p - q) / denom1 - (1 - p - q + 2 * r) / (r * denom2)
  return(c(H1, H2, H3))
}

# Yule's Q
#' @keywords internal
G <- function(p, q, r) {
  denom <- (p * (q - 2 * r) + r * (1 - 2 * q + 2 * r)) ^ 2
  G1 <- 2 * (q - 1) * r * (q - r) / denom
  G2 <- 2 * (p - 1) * r * (p - r) / denom
  G3 <- -(2 * (q * p ^ 2 + p * q * (-1 + q - 2 * r) + r ^ 2)) / denom
  return(c(G1, G2, G3))
}

