#' Phi Coefficient
#'
#' `Phi()` computes Phi and corresponding confidence intervals.
#'
#' @param X a n x 1 numeric vector, matrix or data frame with values of 0 or 1. If a contingency table is supplied, the upper left corner shall contain the joint success probability (or frequency). If a 3-dimensional vector of probabilities is supplied, the order (p, q, r) shall be respected.
#' @param Y NULL (default) or a n x 1 numeric vector, matrix or data frame with values of 0 or 1 and compatible dimensions to X.
#' @param alpha confidence level for the returned confidence interval. FALSE yields Phi without confidence intervals.
#' @param Fisher Indicator whether confidence intervals should be computed by using the Fisher Transformation. Default is TRUE.
#' @param covar data assumptions: iid ("iid"), heteroskedasticity ("HC") or heteroskedasticity and autocorrelation ("HAC").
#' @param n sample size. Only necessary if probabilities are provided and confidence intervals are desired.
#'
#' @return The value of Phi together with the specified confidence interval.
#' @export
#'
#' @references
#' - \insertRef{pohle2024measuringdependenceevents}{BCor}
#' - \insertRef{pearson1900}{BCor}
#' - \insertRef{boas1909}{BCor}
#' - \insertRef{yule1912}{BCor}
#'
#' @examples
#' # Insert a contingency table with frequencies in form of a matrix.
#' x <- matrix(c(10, 20, 30, 5), ncol = 2)
#' Phi(x)
#'
#' # Insert a contingency table with relative frequencies in form of a matrix.
#' # In that case, disable confidence intervals via alpha = FALSE or supply a sample size n.
#' x <- matrix(c(0.2, 0.1, 0.4, 0.3), ncol = 2)
#' Phi(x, alpha = FALSE)
#'
#' # Insert two vectors of observations.
#' x <- c(0,1,1,1,1,0,1,0,0,1)
#' y <- c(1,0,1,1,0,0,0,0,1,0)
#' Phi(x, y)
#'
#' # Insert two marginal success probabilities (p for the row variable and q for the column variable)
#' # and a joint success probability r. In that case, disable confidence intervals via alpha = FALSE
#' # or supply a sample size n.
#' p <- 0.6
#' q <- 0.3
#' r <- 0.2
#' Phi(c(p, q, r), alpha = FALSE)
Phi <- function (X, Y = NULL, alpha = 0.95, Fisher = TRUE, covar = "iid", n = 10000) {
  if (isFALSE(alpha)){
    if (!is.null(Y)){
      X <- table(X, Y)
    }
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
    a <- as.numeric(X[1, 1])
    b <- as.numeric(X[1, 2])
    c <- as.numeric(X[2, 1])
    d <- as.numeric(X[2, 2])
    Phi <- (a * d - b * c)/sqrt((a + b) * (a + c) * (b + d) * (c + d))
    res <- dplyr::tribble(~Phi,
                          #--
                          Phi)
    return(res)
  }
  if (is.null(Y)){
    if (length(X) == 3){
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
    colnames(X) <- c(1,0)
    rownames(X) <- c(1,0)
    if (dplyr::near(sum(X), 1)){
      X <- X * n
    }
    x <- as.numeric(as.character(epitools::expand.table(X)[,1]))
    y <- as.numeric(as.character(epitools::expand.table(X)[,2]))
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
    res <- dplyr::tribble(~Phi,
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

  ### Inference for Phi with Fisher Transformation
  L_est <- L(p_est, q_est, r_est)
  if (isTRUE(Fisher)){
    PhiZ_Var <- as.numeric(t(L_est) %*% Omega %*% L_est / (1 - Phi^2)^2) / n
    PhiZ_CI <- c(tanh(PhiZ + stats::qnorm((1 - alpha) / 2) * sqrt(PhiZ_Var)),
                 tanh(PhiZ - stats::qnorm((1 - alpha) / 2) * sqrt(PhiZ_Var)))
    res <- dplyr::tribble(~Phi, ~CI_lower, ~CI_upper,
                          #--|--|--
                          Phi, PhiZ_CI[1], PhiZ_CI[2])
  } else { # Inference for Phi without Fisher Transformation
    Phi_Var <- as.numeric(t(L_est) %*% Omega %*% L_est) / n
    Phi_CI <- c(Phi + stats::qnorm((1 - alpha) / 2) * sqrt(Phi_Var),
                Phi - stats::qnorm((1 - alpha) / 2) * sqrt(Phi_Var))
    res <- dplyr::tribble(~Phi, ~CI_lower, ~CI_upper,
                          #--|--|--
                          Phi, Phi_CI[1], Phi_CI[2])
  }
  }
  return(res)
}

# Helper functions for the Jacobians used for the Delta-Method
#' @keywords internal
L <- function(p, q, r) {
  denom <- sqrt(p * (1 - p) * q * (1 - q))
  mult <- (p * q - r) / (2 * p * (1 - p) * q * (1 - q))
  L1 <- (-q + mult * (1 - p) * q * (1 - q) - p * q * (1 - q)) / denom
  L2 <- (-p + mult * p * (1 - p) * (1 - q) - p * (1 - p) * q) / denom
  L3 <- 1 / denom
  return(c(L1, L2, L3))
}

