% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2x2_inverse.R
\name{Cole.inv}
\alias{Cole.inv}
\title{Cole.inv}
\usage{
Cole.inv(C, m)
}
\arguments{
\item{C}{Cole's C}

\item{m}{a 4 x 1 vector of marginals c(R1, R2, C1, C2)}
}
\value{
a matrix corresponding to the inputs of a 2x2 contingency table
}
\description{
Cole.inv
}
\examples{
C <- 0.5
m <- c(10, 40, 20, 30)
Cole.inv(C, m)
}
