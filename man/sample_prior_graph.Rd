% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priors.R
\name{sample_prior_graph}
\alias{sample_prior_graph}
\title{Propose new graph and accept with Metropolis ratio}
\usage{
sample_prior_graph(graph, calculate_prior_ratio, prior_params)
}
\arguments{
\item{graph}{A graph represented as an adjacency matrix.}

\item{calculate_prior_ratio}{A function to calculate the prior ratio.}

\item{prior_params}{A list with parameters needed by the prior function.}
}
\value{
The updated graph.
}
\description{
Proposes a new graph by toggling a random edge and applies the Metropolis-Hastings rule.
}
\examples{
g <- diag(0, 3)
GSFM:::sample_prior_graph(g, bern_prior_ratio, list(p = 3, phi = 0.5))
}
