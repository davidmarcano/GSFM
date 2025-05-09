% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/updates.R
\name{update_mus_list}
\alias{update_mus_list}
\title{Update Group Means in Lattice Model}
\usage{
update_mus_list(Kw, mus, f, n)
}
\arguments{
\item{Kw}{Current spatial precision matrix (L x L)}

\item{mus}{Current mean parameters (vector of length L)}

\item{f}{List of factor scores (length L, each element vector of size n_l)}

\item{n}{Vector of sample sizes for each group (length L)}
}
\value{
Updated mean parameters (vector of length L)
}
\description{
Updates the mean parameters (mus) for each group in the lattice model using Gibbs sampling.
Handles both single-group and multi-group cases.
}
