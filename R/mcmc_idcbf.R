library(BDgraph)

#' @title Multi-group Identifiable Double Conditional Bayes Factor Sampler
#' @description This function runs a MCMC multi-group version of the Identifiable Direct Double Conditional Bayes Factor (IDCBF) sampler
#' for Gaussian Graphical Models. It is based on the DCBF algorithm from Hinne et. al.
#'
#' @param sample_cov_list A list of length \code{L} (number of groups), where each element is the sample covariance matrix for a group.
#' @param ns A vector of sample sizes for each group. If not provided, the function assumes it is embedded in \code{prior_params}.
#' @param method Graphical model type. Currently only \code{"ggm"} is supported.
#' @param algorithm Graph estimation algorithm. Currently only \code{"ggmdcbf"} is supported.
#' @param iter Total number of MCMC iterations.
#' @param burnin Number of initial iterations to discard (burn-in period).
#' @param prior A function to sample an initial graph from the prior.
#' @param prior_params Parameters passed to the \code{prior} function.
#' @param prior_ratio A function returning the prior ratio between two graphs.
#' @param df_prior Degrees of freedom for the G-Wishart prior (must be >= 3).
#' @param identifiable Logical. If \code{TRUE}, ensures that sampled graphs are identifiable.
#' @param fast Logical. If \code{TRUE}, expects \code{num_edges_sample} to control sampling efficiency.
#' @param num_edges_sample Optional. Number of edges to sample per iteration if \code{fast = TRUE}.
#' @param nu Optional hyperparameter passed to DCBF inner sampler.
#' @param theta Optional hyperparameter passed to DCBF inner sampler.
#'
#' @return A list of length \code{iter - burnin}, where each element is itself a list of length \code{L} containing the sampled graph and precision matrix for each group.
#' @export
mcmc_idcbf <- function(sample_cov_list,
                               ns = NULL,
                               method = "ggm",
                               algorithm = "ggmdcbf",
                               iter = 1000,
                               burnin = iter / 2,
                               prior = mc_unif_prior,
                               prior_params = NULL,
                               prior_ratio = NULL,
                               df_prior = 3,
                               identifiable = FALSE,
                               fast = FALSE,
                               num_edges_sample = NULL,
                               nu = NULL,
                               theta = NULL) {

  if (method != "ggm") stop("Currently the only 'method' supported is 'ggm'.")
  if (algorithm != "ggmdcbf") stop("Currently the only 'algorithm' supported is 'ggmdcbf'.")
  if (df_prior < 3) stop("'df_prior' must be >= 3.")
  if (iter < burnin) stop("Number of iterations must exceed number of burn-in iterations.")
  if (is.null(prior_params)) stop("Missing required prior parameters.")
  if (is.null(prior_ratio)) stop("Please specify a prior ratio function.")
  if (fast && is.null(num_edges_sample)) stop("Specify the number of edges to sample when using 'fast = TRUE'.")

  burnin <- floor(burnin)
  L <- length(sample_cov_list)

  example_cov <- sample_cov_list[[1]]
  p <- ncol(example_cov)

  delta <- df_prior
  D <- diag(p)

  Gs <- vector("list", L)
  Ks <- vector("list", L)
  alphas <- vector("list", L)
  deltan <- numeric(L)

  for (l in 1:L) {
    graph <- matrix(0, p, p)
    G <- prior(graph, prior_params)$graph

    if (identifiable) {
      while (!is_identifiable(G)) {
        G <- prior(graph, prior_params)$graph
      }
    }

    Gs[[l]] <- G
    Ks[[l]] <- BDgraph::rgwish(1, G, delta, D)
    deltan[l] <- delta + ns[l]
  }

  results <- vector("list", iter)

  if (method == "ggm" && algorithm == "ggmdcbf") {
    for (k in 1:iter) {
      if (k %% 100 == 0) cat("Iteration:", k, "\n")

      results_groups <- vector("list", L)
      for (l in 1:L) {
        results_groups[[l]] <- idcbf_step(
          Gs = Gs,
          D = D,
          S = sample_cov_list[[l]],
          K = Ks[[l]],
          p = p,
          delta = delta,
          deltan = deltan[l],
          nu = nu[l],
          theta = theta,
          l = l,
          prior = prior,
          prior_params = prior_params,
          prior_ratio = prior_ratio,
          identifiable = identifiable
        )
        Gs[[l]] <- results_groups[[l]]$G
        Ks[[l]] <- results_groups[[l]]$K
        alphas[[l]] <- results_groups[[l]]$alphas
      }

      results[[k]] <- results_groups
    }
  }

  return(results[-seq_len(burnin)])
}
