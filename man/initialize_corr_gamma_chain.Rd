% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gibbs_sampler_helpers.R
\name{initialize_corr_gamma_chain}
\alias{initialize_corr_gamma_chain}
\title{Initialize correlation matrices chain}
\usage{
initialize_corr_gamma_chain(steps, p, L, start_corr_gamma, has_start)
}
\arguments{
\item{steps}{Number of MCMC steps}

\item{p}{Number of variables}

\item{L}{Number of groups}

\item{start_corr_gamma}{Starting values for correlation matrices}

\item{has_start}{Logical indicating if starting values are provided}
}
\value{
Initialized correlation matrices array
}
\description{
Initializes the chain for group correlation matrices.
}
\keyword{internal}
