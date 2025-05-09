% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/updates.R
\name{update_lambdas}
\alias{update_lambdas}
\title{Update Lambda Values}
\usage{
update_lambdas(K, x, f, delta, alpha, lambda)
}
\arguments{
\item{K}{The precision matrix.}

\item{x}{The data matrix.}

\item{f}{The current factor values.}

\item{delta}{The regularization parameter.}

\item{alpha}{The current values of alpha.}

\item{lambda}{The current values of lambda.}
}
\value{
A matrix of updated lambda values.
}
\description{
This function updates the lambda values based on the current values of alpha and f.
}
\examples{
set.seed(123)
x <- matrix(rnorm(20), 5, 4)
f <- rnorm(5)
alpha <- rnorm(4)
lambda <- rnorm(4)
K <- diag(4)
update_lambdas(K, x, f, 1, alpha, lambda)
}
