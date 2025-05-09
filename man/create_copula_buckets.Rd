% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gibbs_sampler.R
\name{create_copula_buckets}
\alias{create_copula_buckets}
\title{Create Copula Buckets for Non-Continuous Data}
\usage{
create_copula_buckets(xdf)
}
\arguments{
\item{xdf}{Data frame or matrix of observations for a single group}
}
\value{
List containing:
\itemize{
\item lessthan - Array \link{n x p x num_breakpoints} indicating values less than each breakpoint
\item greaterthan - Array \link{n x p x num_breakpoints} indicating values greater than each breakpoint
\item obs_buckets - Matrix \link{n x p} indicating bucket membership for each observation
}
}
\description{
Generates data structures needed for copula modeling of non-continuous variables
by creating breakpoints and bucket indicators.
}
