#' Inverse Cole's C for 2 x 2 contingency tables
#'
#' `Cole.inv()` computes a 2 x 2 contingency table based on the marginal frequencies or probabilities and Cole's C.
#'
#' @param C Cole's C
#' @param m a 4 x 1 vector of marginals c(R1, R2, C1, C2)
#'
#' @return a matrix corresponding to the inputs of a 2x2 contingency table
#' @export
#'
#' @references
#' - \insertRef{pohle2024measuringdependenceevents}{BCor}
#' - \insertRef{Cole1949}{BCor}
#'
#' @examples
#' C <- 0.5
#' m <- c(10, 40, 20, 30)
#' Cole.inv(C, m)
Cole.inv <- function (C, m)
{
  # check whether C and m are numeric and whether m is a reasonable vector
  if (!(is.numeric(C) && is.numeric(m) && sum(m[1:2]) == sum(m[3:4]))){
    stop("`C` and `m` must be numeric and the marginal frequencies have to sum to the same number", call. = FALSE)
  }
  n <- sum(m) / 2
  R1 <- m[1] / n
  R2 <- m[2] / n
  C1 <- m[3] / n
  C2 <- m[4] / n

  R <- matrix(ncol = 2, nrow = 2)
  if (C >= 0){
    R[1, 1] <- C * min(R1 * C2, R2 * C1) + R1 * C1
    R[1, 2] <- R1 - R[1, 1]
    R[2, 1] <- C1 - R[1, 1]
    R[2, 2] <- R2 - C1 + R[1, 1]
  }
  else {
    R[1, 1] <- C * min(R1 * C1, R2 * C2) + R1 * C1
    R[1, 2] <- R1 - R[1, 1]
    R[2, 1] <- C1 - R[1, 1]
    R[2, 2] <- R2 - C1 + R[1, 1]
  }
  return(round(R * n, digits = 10))
}

#' Inverse Yule's Q for 2 x 2 contingency tables
#'
#' `YuleQ.inv()` computes a 2 x 2 contingency table based on the marginal frequencies or probabilities and Yule's Q.
#'
#' @param Q Yule's Q
#' @param m a 4 x 1 vector of marginals c(R1, R2, C1, C2)
#'
#' @return a matrix corresponding to the inputs of a 2x2 contingency table
#' @export
#'
#' @references
#' - \insertRef{pohle2024measuringdependenceevents}{BCor}
#' - \insertRef{yule1900}{BCor}
#'
#' @examples
#' Q <- 0.5
#' m <- c(10, 40, 20, 30)
#' YuleQ.inv(Q, m)
YuleQ.inv <- function (Q, m)
{
  # check whether Q and m are numeric and whether m is a reasonable vector
  if (!(is.numeric(Q) && is.numeric(m) && sum(m[1:2]) == sum(m[3:4]))){
    stop("`Q` and `m` must be numeric and the marginal frequencies have to sum to the same number", call. = FALSE)
  }
  n <- sum(m) / 2
  R1 <- m[1] / n
  R2 <- m[2] / n
  C1 <- m[3] / n
  C2 <- m[4] / n
  # For Q = 0, we have
  R_1 <- matrix(ncol = 2, nrow = 2)
  if (Q == 0){
    R_1[1, 1] <- R1 * C1
    R_1[1, 2] <- R1 - R_1[1, 1]
    R_1[2, 1] <- C1 - R_1[1, 1]
    R_1[2, 2] <- R2 - C1 + R_1[1, 1]
    R_1 <- round(R_1, digits = 10)
    return(R_1 * n)
  }
  # For Q \ne 0, we have a quadratic equation with two solutions.
  else{
    # Solution 1:
    R_2 <- matrix(ncol = 2, nrow = 2)
    R_2[1, 1] <- - (Q * (1 - 2 * R1 - 2 * C1) - 1) / (4 * Q) + sqrt(((Q * (1 - 2 * R1 - 2 * C1) - 1) ^ 2) / (16 * Q ^ 2) - (Q + 1) / (2 * Q) * R1 * C1)
    R_2[1, 2] <- R1 - R_2[1, 1]
    R_2[2, 1] <- C1 - R_2[1, 1]
    R_2[2, 2] <- R2 - C1 + R_2[1, 1]

    # Round in order to avoid small negative probabilities that only result due to computational inaccuracy
    R_2 <- round(R_2, digits = 10)

    # Solution 2:
    R_3 <- matrix(ncol = 2, nrow = 2)
    R_3[1, 1] <- - (Q * (1 - 2 * R1 - 2 * C1) - 1) / (4 * Q) - sqrt(((Q * (1 - 2 * R1 - 2 * C1) - 1) ^ 2) / (16 * Q ^ 2) - (Q + 1) / (2 * Q) * R1 * C1)
    R_3[1, 2] <- R1 - R_3[1, 1]
    R_3[2, 1] <- C1 - R_3[1, 1]
    R_3[2, 2] <- R2 - C1 + R_3[1, 1]

    # Round in order to avoid small negative probabilities that only result due to computational inaccuracy
    R_3 <- round(R_3, digits = 10)

    # Only one of the two solutions should deliver exclusively positive joint probabilities
    if (all(R_2 >= matrix(c(0, 0, 0, 0), ncol = 2))){
      return(R_2 * n)
    }
    else if (all(R_3 >= matrix(c(0, 0, 0, 0), ncol = 2))){
      return(R_3 * n)
    }
    else {print("We have a problem. Probably there is some approximation error.")}
  }
}










