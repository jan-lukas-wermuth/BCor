#' Cole's Correlation
#'
#' `Cole()` computes Cole's correlation and corresponding confidence intervals.
#'
#' @param X a n x 1 numeric vector, matrix or data frame with values of 0 or 1. If a contingency table is supplied, the upper left corner shall contain the joint success probability (or frequency). If a 3-dimensional vector of probabilities is supplied, the order (p, q, r) shall be respected.
#' @param Y NULL (default) or a n x 1 numeric vector, matrix or data frame with values of 0 or 1 and compatible dimensions to X.
#' @param alpha confidence level for the returned confidence interval. FALSE yields Cole's C without confidence intervals.
#' @param Fisher Indicator whether confidence intervals should be computed by using the Fisher Transformation. Default is TRUE.
#' @param covar data assumptions: iid ("iid"), heteroskedasticity ("HC") or heteroskedasticity and autocorrelation ("HAC").
#' @param m_rep number of Monte Carlo replications used for approximating the limiting distribution of Cole's C.
#' @param c_seq sequence of C_0 to be tested for computing the confidence interval. Default is the sequence from -0.999 to 0.999 with distance 0.001.
#' @param n sample size. Only necessary if probabilities are provided and confidence intervals are desired.
#' @param details_inference logical. If TRUE, the output includes all 3 p-values that arise during the testing procedure that yields confidence intervals.
#'
#' @return The value of Cole's correlation.
#' @export
#'
#' @examples
#' x <- matrix(c(10, 20, 30, 5), ncol = 2)
#' Cole(x)
Cole <- function (X, Y = NULL, alpha = 0.95, Fisher = TRUE, covar = "iid", m_rep = 10000, c_seq = NA, n = 10000, details_inference = FALSE) {
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
    if (a*d-b*c >= 0){
      C <- (a*d-b*c) / min((a + b) * (b + d), (a + c) * (c + d))
    }
    else C <- (a*d-b*c) / min((a + b) * (a + c), (c + d) * (b + d))
    res <- dplyr::tribble(~C,
                          #--
                          C)
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

  m_plus <- min(p_est, q_est) - p_est * q_est
  m_minus <- p_est * q_est - max(0, p_est + q_est - 1)
  C <- 1 * (cov >= 0) * cov / m_plus + 1 * (cov < 0) * cov / m_minus
  CZ <- atanh(C)

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
    res <- dplyr::tribble(~C,
                          #--
                          C)
    return(res)
  } else {
    # Estimate the covariance matrix Omega, either by HC, HAC or the sample cov
    if (covar == "iid") {
      Omega <- stats::var(df)
    } else if (covar == "HAC") {
      Omega <- n * sandwich::vcovHAC(stats::lm(as.matrix(df) ~ 1))
    } else if (covar == "HC"){
      Omega <- n * sandwich::vcovHC(stats::lm(as.matrix(df) ~ 1))
    } else stop("Please insert a valid option for computing the covariance matrix!")

  ### Inference
  # Set the sequence of possible values c_seq here for the CI through inverted tests
  if (any(is.na(c_seq))) {
    c_seq <- seq(-1, 1, length.out = 2001) %>% utils::head(-1) %>% utils::tail(-1)
  }

  # "Pre-test" for sigma <= 0 and sigma >= 0
  var_sigma <- c(-q_est,-p_est, 1) %*% Omega %*% c(-q_est,-p_est, 1) / n
  T_stat_sigma0 <- cov / sqrt(var_sigma)
  p_val_ctrue_pos <- stats::pnorm(T_stat_sigma0)
  p_val_ctrue_neg <- 1 - stats::pnorm(T_stat_sigma0)
  p_val_ctrue_0 <- 2 * stats::pnorm(-abs(T_stat_sigma0))

  p_vals_pretest <- c(rep(p_val_ctrue_neg, (length(c_seq) - 1) / 2), p_val_ctrue_0, rep(p_val_ctrue_pos, (length(c_seq) - 1) / 2))

  # "Pre-test" for p=q and p=1-q
  pv_eq_pre_pos <- 2 * stats::pnorm(-abs(p_est - q_est), mean = 0, sd = sqrt((Omega[1, 1] + Omega[2, 2] - 2 * Omega[1, 2]) / n))
  pv_eq_pre_neg <- 2 * stats::pnorm(-abs(p_est + q_est), mean = 1, sd = sqrt((Omega[1, 1] + Omega[2, 2] + 2 * Omega[1, 2]) / n))

  pv_eq_pre <- c(rep(pv_eq_pre_neg, (length(c_seq) - 1) / 2), NA, rep(pv_eq_pre_pos, (length(c_seq) - 1) / 2))

  # The case {cc=0}
  Delta <- c(-q_est,-p_est, 1)

  # The cases {cc>0 and p!=q} and {cc<0 and p!=1-q}
  Lambda_pos <- c(1 / m_plus,-cov / m_plus ^ 2)
  Lambda_neg <- c(1 / m_minus,-cov / m_minus ^ 2)
  Jh_est_pos <- Jh_pos(p_est, q_est, r_est)
  Jh_est_neg <- Jh_neg(p_est, q_est, r_est)
  C_Var_pos <- as.numeric(t(Lambda_pos) %*% t(Jh_est_pos) %*% Omega %*% Jh_est_pos %*% Lambda_pos) / n
  C_Var_neg <- as.numeric(t(Lambda_neg) %*% t(Jh_est_neg) %*% Omega %*% Jh_est_neg %*% Lambda_neg) / n

  # The cases {cc>0 and p=q} and {cc<0 and p=1-q}
  Jf <- rbind(c(1, 0, 0),
              c(0, 1, 0),
              c(q_est, p_est, 0),
              c(-q_est, -p_est, 1))
  V_Var <- Jf %*% Omega %*% t(Jf)

  level_seq <- seq(0, 1, length.out = 1001) %>% utils::head(-1) %>% utils::tail(-1)

  # Draw from the (estimated) asymptotic distributions...
  mV <- MASS::mvrnorm(m_rep, mu = rep(0, 4), Sigma = V_Var)

  # ... in the positive case
  sigma_mpos_vec <- rbind(mV[, 4], (mV[, 1] - mV[, 3]) - 1 * (mV[, 1] > mV[, 2]) * (mV[, 1] - mV[, 2]))
  Cdiff_pos_samples <- as.numeric(Lambda_pos %*% sigma_mpos_vec) # Samples from sqrt(n) (hat C_n - cc)

  # ... in the negative case
  sigma_mneg_vec <- rbind(mV[, 4], mV[, 3] - 1 * (mV[, 1] > -mV[, 2]) * (mV[, 1] + mV[, 2]))
  Cdiff_neg_samples <- as.numeric(Lambda_neg %*% sigma_mneg_vec)   # Samples from sqrt(n) (hat C_n - cc)

  eps <- 0.000001 # Small number such that we avoid cc < 0 in the later if-clauses

  # Initialize P-values
  pvals <- matrix(rep(NA, 2 * length(c_seq)), ncol = 2)

  # Fisher-transformed confidence intervals
  if (isTRUE(Fisher)){
    CZ_Var_pos <- C_Var_pos / (1 - C^2)^2
    CZ_Var_neg <- C_Var_neg / (1 - C^2)^2

    # Draw from the (estimated) asymptotic distributions in the positive case
    CZdiff_pos_samples <- Cdiff_pos_samples / (1 - C^2) # Fisher-transformed Samples from sqrt(n) (hat C_n - cc)

    # Draw from the (estimated) asymptotic distributions in the negative case
    CZdiff_neg_samples <- Cdiff_neg_samples / (1 - C^2) # Fisher-transformed Samples from sqrt(n) (hat C_n - cc)

    ### Construct a (possibly conservative) confidence interval for C through inverted tests
    pvals <- matrix(rep(NA, 2 * length(c_seq)), ncol = 2)

    for (c_index in 1:length(c_seq)) {
      ccZ <- atanh(c_seq[c_index])

      # The case {cc < 0}
      if (ccZ < -eps) {
        # The case {p != 1-q}
        TZ_stat <- (CZ - ccZ) / sqrt(CZ_Var_neg)
        pvals[c_index, 1] <- 2 * stats::pnorm(-abs(TZ_stat))

        # The case {p == 1-q}
        # Compute the quantile level corresponding to cc, and the respective p-value, this is simply a resampling p-value
        ccZ_quantile_level <- mean(CZ - CZdiff_neg_samples / sqrt(n) <= ccZ)
        pvals[c_index, 2] <- 2 * min(ccZ_quantile_level, 1 - ccZ_quantile_level)

        # The case {cc > 0}
      } else if (ccZ > eps) {
        # The case {p != q}
        TZ_stat <- (CZ - ccZ) / sqrt(CZ_Var_pos)
        pvals[c_index, 1] <- 2 * stats::pnorm(-abs(TZ_stat))

        # The case {p == q}
        # Compute the quantile level corresponding to cc, and the respective p-value
        ccZ_quantile_level <- mean(CZ - CZdiff_pos_samples / sqrt(n) <= ccZ)
        pvals[c_index, 2] <- 2 * min(ccZ_quantile_level, 1 - ccZ_quantile_level)

        # The case {cc == 0}
      } else {
        # Compare with the sum of two truncated normals:
        # This test draws on samples from the asymptotic distribution
        Z1 <- stats::rnorm(m_rep)
        CZdiff_zero_samples <- sqrt(as.numeric(t(Delta) %*% Omega %*% Delta) / (1 - C^2)^2) * (Z1 * 1 * (Z1 >= 0) / m_plus + Z1 * 1 * (Z1 < 0) / m_minus)
        pvals[c_index, 2] <- mean(abs(CZdiff_zero_samples) > sqrt(n) * abs(CZ - ccZ)) / 2 # Divide by 2 because in the merging of p-values we (wrongly) multiply by 2. No Bonferroni necessary!
      }
    }
  }
  # Non-Fisher-transformed confidence intervals
  else {
    ### Construct a (possibly conservative) confidence interval for C through inverted tests
    for (c_index in 1:length(c_seq)) {
      cc <- c_seq[c_index]

      # The case {cc < 0}
      if (cc < -eps) {
        # The case {p != 1-q}
        T_stat <- (C - cc) / sqrt(C_Var_neg)
        pvals[c_index, 1] <- 2 * stats::pnorm(-abs(T_stat))

        # The case {p == 1-q}
        # Compute the quantile level corresponding to cc, and the respective p-value, this is simply a resampling p-value
        cc_quantile_level <- mean(C - Cdiff_neg_samples / sqrt(n) <= cc)
        pvals[c_index, 2] <- 2 * min(cc_quantile_level, 1 - cc_quantile_level)

        # The case {cc > 0}
      } else if (cc > eps) {
        # The case {p != q}
        T_stat <- (C - cc) / sqrt(C_Var_pos)
        pvals[c_index, 1] <- 2 * stats::pnorm(-abs(T_stat))

        # The case {p == q}
        # Compute the quantile level corresponding to cc, and the respective p-value
        cc_quantile_level <- mean(C - Cdiff_pos_samples / sqrt(n) <= cc)
        pvals[c_index, 2] <- 2 * min(cc_quantile_level, 1 - cc_quantile_level)

        # The case {cc == 0}
      } else {
        # Compare with the sum of two truncated normals:
        # This test draws on samples from the asymptotic distribution
        Z1 <- stats::rnorm(m_rep)
        Cdiff_zero_samples <- sqrt(as.numeric(t(Delta) %*% Omega %*% Delta)) * (Z1 * 1 * (Z1 >= 0) / m_plus + Z1 * 1 * (Z1 < 0) / m_minus)
        pvals[c_index, 2] <- mean(abs(Cdiff_zero_samples) > sqrt(n) * abs(C - cc)) / 2 # Divide by 2 because in the merging of p-values we (wrongly) multiply by 2. No Bonferroni necessary!
      }
    }
  }
  # Merge the p-values
  pval_min_pretest <- pmin(1, 2 * apply(cbind(p_vals_pretest, apply(cbind(2 * apply(cbind(pvals[,2], pv_eq_pre), 1, min, na.rm = T), pvals[,1]), 1, max, na.rm = T)), 1, min, na.rm=T))

  # Define the lower and upper limits of the confidence intervals for C
  lower_vals_C <- min(c_seq[which(pval_min_pretest >= 1 - alpha)])
  upper_vals_C <- max(c_seq[which(pval_min_pretest >= 1 - alpha)])

  # If no rejections occur, select -1 and 1 instead of Inf and -Inf
  lower_vals_C <- ifelse(lower_vals_C == Inf,-1, lower_vals_C)
  upper_vals_C <- ifelse(upper_vals_C == -Inf, 1, upper_vals_C)

  # Define the confidence interval
  C_CI <- c(lower_vals_C, upper_vals_C)

  res <- dplyr::tribble(~C, ~CI_lower, ~CI_upper,
                        #--|--|--
                        C, C_CI[1], C_CI[2])
  }
  if (details_inference){
    details <- dplyr::tibble(c_seq = c_seq, pv_uneq = pvals[,1], pv_eq = pvals[,2], pv_pretest = p_vals_pretest, pv_eq_pretest = pv_eq_pre)
    return(list(res, details))
  } else return(res)
}


#' @keywords internal
Jh_pos <- function(p, q, r) {
  Jh <- cbind(c(-q,-p, 1),
              c(1 * (p < q) - q, 1 * (p > q) - p, 0))
  return(Jh)
}

#' @keywords internal
Jh_neg <- function(p, q, r) {
  Jh <- cbind(c(-q,-p, 1),
              c(q - 1 * (p + q > 1), p - 1 * (p + q > 1), 0))
  return(Jh)
}

