library(stats)

#' Compute the number of edges in an undirected graph
#'
#' This function computes the total number of edges in an undirected graph represented as an adjacency matrix.
#' It counts the number of 1s in the upper triangle of the matrix (excluding the diagonal), which represents the edges.
#'
#' @param g An undirected graph represented as a square adjacency matrix.
#'
#' @return The number of edges in the graph.
#'
#' @examples
#' g <- matrix(c(0, 1, 0,
#'               1, 0, 1,
#'               0, 1, 0), nrow = 3, byrow = TRUE)
#' GSFM:::size(g)
size <- function(g) {
  sum(g[upper.tri(g)])
}

#' Uniform prior ratio
#'
#' Computes the prior ratio between two graphs under a uniform prior. Always returns 1.
#'
#' @param g1 A graph represented as an adjacency matrix.
#' @param g2 A graph represented as an adjacency matrix.
#' @param prior_params A list of prior parameters (not used in this case).
#'
#' @return A constant value of 1.
#' @export
#'
#' @examples
#' g1 <- diag(0, 3)
#' g2 <- matrix(c(0, 1, 0,
#'                1, 0, 0,
#'                0, 0, 0), 3, 3)
#' unif_prior_ratio(g1, g2, list())
unif_prior_ratio <- function(g1, g2, prior_params) {
  1
}

#' Bernoulli prior ratio
#'
#' Computes the prior ratio between two graphs under a Bernoulli prior.
#'
#' @param g1 A graph represented as an adjacency matrix.
#' @param g2 A graph represented as an adjacency matrix.
#' @param prior_params A list containing the parameter \code{phi}.
#'
#' @return A numeric value representing the prior ratio.
#' @export
#'
#' @examples
#' g1 <- diag(0, 3)
#' g2 <- matrix(c(0, 1, 0,
#'                1, 0, 0,
#'                0, 0, 0), 3, 3)
#' bern_prior_ratio(g1, g2, list(phi = 0.5))
bern_prior_ratio <- function(g1, g2, prior_params) {
  phi <- prior_params$phi
  diff <- size(g1) - size(g2)
  phi^diff * (1 - phi)^(-diff)
}

#' Beta prior ratio
#'
#' Computes the prior ratio between two graphs under a Beta prior.
#'
#' @param g1 A graph represented as an adjacency matrix.
#' @param g2 A graph represented as an adjacency matrix.
#' @param prior_params A list with parameters \code{a}, \code{b}, and \code{m}.
#'
#' @return A numeric value representing the prior ratio.
#' @export
#'
#' @examples
#' g1 <- diag(0, 3)
#' g2 <- matrix(c(0, 1, 0,
#'                1, 0, 0,
#'                0, 0, 0), 3, 3)
#' beta_prior_ratio(g1, g2, list(a = 1, b = 1, m = 3))
beta_prior_ratio <- function(g1, g2, prior_params) {
  a <- prior_params$a
  b <- prior_params$b
  m <- prior_params$m
  beta(a + size(g1), b + m - size(g1)) / beta(a + size(g2), b + m - size(g2))
}

#' Hierarchical prior ratio
#'
#' Computes the prior ratio under a hierarchical prior using binomial coefficients.
#'
#' @param g1 A graph represented as an adjacency matrix.
#' @param g2 A graph represented as an adjacency matrix.
#' @param prior_params A list with parameter \code{m}.
#'
#' @return A numeric value representing the prior ratio.
#' @export
#'
#' @examples
#' g1 <- diag(0, 3)
#' g2 <- matrix(c(0, 1, 0,
#'                1, 0, 0,
#'                0, 0, 0), 3, 3)
#' hierarchy_prior_ratio(g1, g2, list(m = 3))
hierarchy_prior_ratio <- function(g1, g2, prior_params) {
  m <- prior_params$m
  choose(m, size(g2)) / choose(m, size(g1))
}

#' Propose new graph and accept with Metropolis ratio
#'
#' Proposes a new graph by toggling a random edge and applies the Metropolis-Hastings rule.
#'
#' @param graph A graph represented as an adjacency matrix.
#' @param calculate_prior_ratio A function to calculate the prior ratio.
#' @param prior_params A list with parameters needed by the prior function.
#'
#' @return The updated graph.
#'
#' @examples
#' g <- diag(0, 3)
#' GSFM:::sample_prior_graph(g, bern_prior_ratio, list(p = 3, phi = 0.5))
sample_prior_graph <- function(graph, calculate_prior_ratio, prior_params) {
  p <- prior_params$p
  candidate_edge <- sample(1:p, 2)
  i <- candidate_edge[1]
  j <- candidate_edge[2]

  candidate_graph <- graph
  candidate_graph[i, j] <- 1 - graph[i, j]
  candidate_graph[j, i] <- candidate_graph[i, j]

  ratio <- calculate_prior_ratio(candidate_graph, graph, prior_params)

  if (stats::rbinom(1, 1, min(1, ratio))) {
    graph <- candidate_graph
  }

  graph
}

#' Generic MCMC sampler
#'
#' Runs a Metropolis-Hastings sampler for graphs.
#'
#' @param graph A graph represented as an adjacency matrix.
#' @param prior_params A list of prior parameters including \code{p} and \code{iter}.
#' @param prior_ratio_fn A prior ratio function.
#'
#' @return A list with the final graph and edge probabilities.
#'
#' @examples
#' g <- diag(0, 3)
#' GSFM:::run_mc_sampler(g, list(p = 3, iter = 10, phi = 0.5), bern_prior_ratio)
run_mc_sampler <- function(graph, prior_params, prior_ratio_fn) {
  p <- prior_params$p
  iter <- prior_params$iter
  prior_params$m <- choose(p, 2)

  graphs_sum <- graph
  for (i in 1:iter) {
    graph <- sample_prior_graph(graph, prior_ratio_fn, prior_params)
    graphs_sum <- graphs_sum + graph
  }

  list(graph = graph, edge_probs = graphs_sum / iter)
}

#' MCMC sampler with uniform prior
#'
#' @inheritParams run_mc_sampler
#' @return A list with the final graph and edge probabilities.
#' @export
#'
#' @examples
#' g <- diag(0, 3)
#' mc_unif_prior(g, list(p = 3, iter = 10))
mc_unif_prior <- function(graph, prior_params) run_mc_sampler(graph, prior_params, unif_prior_ratio)

#' MCMC sampler with Bernoulli prior
#'
#' @inheritParams run_mc_sampler
#' @return A list with the final graph and edge probabilities.
#' @export
#'
#' @examples
#' g <- diag(0, 3)
#' mc_bernoulli_prior(g, list(p = 3, iter = 10, phi = 0.5))
mc_bernoulli_prior <- function(graph, prior_params) run_mc_sampler(graph, prior_params, bern_prior_ratio)

#' MCMC sampler with Beta prior
#'
#' @inheritParams run_mc_sampler
#' @return A list with the final graph and edge probabilities.
#' @export
#'
#' @examples
#' g <- diag(0, 3)
#' mc_betas_prior(g, list(p = 3, iter = 10, a = 1, b = 1, m = 3))
mc_betas_prior <- function(graph, prior_params) run_mc_sampler(graph, prior_params, beta_prior_ratio)

#' MCMC sampler with hierarchical prior
#'
#' @inheritParams run_mc_sampler
#' @return A list with the final graph and edge probabilities.
#' @export
#'
#' @examples
#' g <- diag(0, 3)
#' mc_hierarchical_prior(g, list(p = 3, iter = 10, m = 3))
mc_hierarchical_prior <- function(graph, prior_params) run_mc_sampler(graph, prior_params, hierarchy_prior_ratio)

#' Sample graph from uniform prior
#'
#' Samples a random graph from the uniform prior.
#'
#' @param graph Not used.
#' @param prior_params A list with parameter \code{p}.
#'
#' @return A sampled graph as an adjacency matrix.
#'
#' @examples
#' GSFM:::sample_unif_prior(NULL, list(p = 3))
sample_unif_prior <- function(graph, prior_params) {
  p <- prior_params$p
  upper <- stats::rbinom((p - 1) * p / 2, 1, 0.5)
  g <- matrix(0, p, p)
  g[upper.tri(g)] <- upper
  g + t(g)
}

#' Sample graph from Bernoulli prior
#'
#' Samples a random graph under a Bernoulli prior.
#'
#' @param graph Not used.
#' @param prior_params A list with \code{p} and \code{phi}.
#'
#' @return A sampled graph as an adjacency matrix.
#'
#' @examples
#' GSFM:::sample_bernoulli_prior(NULL, list(p = 3, phi = 0.5))
sample_bernoulli_prior <- function(graph, prior_params) {
  p <- prior_params$p
  phi <- prior_params$phi
  upper <- stats::rbinom((p - 1) * p / 2, 1, phi)
  g <- matrix(0, p, p)
  g[upper.tri(g)] <- upper
  g + t(g)
}
