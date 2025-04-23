# This file implements a Gibbs sampler for a lattice-based model that handles
# both continuous and non-continuous (copula/probit) data cases. It's separated from
# lattice.R to allow for variable-sized matrices per group using lists.

# Load required libraries and source files
library(MASS)
library(BDgraph)
library(invgamma)
library(pheatmap)
library(igraph)
library(matrixcalc)
library(parallel)

lbw <- FALSE  # Global flag for low birth weight mode

#' Update Group Parameters in Lattice Model
#'
#' Performs all parameter updates for a single group in the lattice model, including
#' graph structure, precision matrix, factor scores, loadings, and intercepts.
#'
#' @param l Group index to update
#' @param alphal Current alpha parameters (intercepts) for all groups
#' @param lambdal Current lambda parameters (loadings) for all groups
#' @param gl Current graph structures for all groups
#' @param kl Current precision matrices for all groups
#' @param nu Current nu parameters
#' @param theta Current theta matrix
#' @param mul Current mu parameters (factor means) for all groups
#' @param delta_group Current delta parameters for all groups
#' @param fl Current factor scores for all groups
#' @param x List of observed data matrices for all groups
#' @param n Vector of sample sizes for each group
#' @param p Number of variables
#' @param Dl Prior precision matrix for graph learning
#' @param deltal Degrees of freedom for G-Wishart prior
#' @param identifiable Logical indicating if identifiability constraints should be enforced
#' @param non_continuous Logical indicating non-continuous data
#' @param prior Graph prior function
#' @param prior_params Parameters for graph prior
#' @param prior_ratio Prior ratio function for graph proposals
#' @param alpha Hyperparameter for theta/gamma updates
#' @param beta Hyperparameter for theta/gamma updates
#' @param w Hyperparameter for theta/gamma updates
#' @param a Hyperparameter for nu updates
#' @param b Hyperparameter for nu updates
#' @param c Hyperparameter for delta updates
#' @param d Hyperparameter for delta updates
#' @param n0 Hyperparameter for alpha updates
#' @param s Current MCMC iteration
#' @param copy_x Original observed data (for non-continuous case)
#' @param zs Current latent variables (for non-continuous case)
#' @param xpred Current predicted values
#' @param cumulativeMargProb Cumulative marginal probabilities (for copula)
#' @param copula Logical indicating copula model
#' @param copula_lessthan_buckets Copula less-than indicators
#' @param copula_greaterthan_buckets Copula greater-than indicators
#' @param obs_buckets Observation buckets for copula
#' @return List containing updated parameters for the group:
#' \itemize{
#'   \item lambda_chain - Updated loadings
#'   \item f_temp - Updated factor scores
#'   \item graph_chain - Updated graph structure
#'   \item K_chain - Updated precision matrix
#'   \item alpha_chain - Updated intercepts
#'   \item tau_chain - Updated delta parameter
#'   \item acceptprob_chain - Acceptance probabilities
#'   \item latent_zs - Updated latent variables (if non-continuous)
#'   \item corr_gamma_chain - Updated correlation matrix
#'   \item xpred - Updated predicted values
#' }
update_group_parameters <- function(l, alphal, lambdal, gl, kl, nu, theta,
                                    mul, delta_group, fl, x, n, p, Dl, deltal,
                                    identifiable, non_continuous, prior, prior_params,
                                    prior_ratio, alpha, beta, w, a, b, c, d, n0, s,
                                    copy_x = NULL, zs = NULL, xpred = NULL,
                                    cumulativeMargProb = NULL, copula = FALSE,
                                    copula_lessthan_buckets = NULL,
                                    copula_greaterthan_buckets = NULL,
                                    obs_buckets = NULL) {

  # Calculate shifted data
  shifted <- calculate_shifted_x(x[[l]], alphal[,l], lambdal[,l], fl[[l]])
  S <- (n[l]-1)*cov(shifted)

  # Update graph and precision matrix using IDCBF
  dcbf <- idcbf_step(gl, Dl, S, kl[,,l], p, deltal, deltal + n[l],
                         nu, theta, l, prior, prior_params, prior_ratio, identifiable)

  Gl <- dcbf$G
  Kl <- dcbf$K
  acceptprob_chain <- dcbf$alphas

  # Handle non-continuous case
  corr_gamma <- cov2cor(solve(Kl))
  if(non_continuous) {
    Kl <- solve(corr_gamma)
  }

  # Update factors, alphas and lambdas
  fl_update <- update_factors_spatial(Kl, x[[l]], lambdal[,l], alphal[,l], mul[l])
  alphas <- update_alphas(Kl, x[[l]], n[l], n0, fl_update, lambdal[,l])

  lambdas <- update_lambdas(Kl, x[[l]], fl_update, delta_group[l], alphas, lambdal[,l])

  # Ensure lambda convergence
  l_count <- 0
  while(lambdas[1] == Inf) {
    lambdas <- update_lambdas(Kl, x[[l]], fl_update, delta_group[l], alphas, lambdal[,l])
    l_count <- l_count + 1
    if(l_count > 1000) stop("Gibbs sampler non-convergence")
  }

  deltas <- update_delta(c, d, lambdas, Kl)

  # Handle copula case
  latent_zs <- xpred <- NULL
  if(non_continuous) {
    latent_zs <- update_zs(Kl, Gl, copy_x[[l]], lambdas, alphas, fl_update,
                           x[[l]], copula, copula_lessthan_buckets[[l]],
                           copula_greaterthan_buckets[[l]], obs_buckets[[l]])

    xpred <- randomMVNSampler(1, p, corr_gamma, cumulativeMargProb[,,l],
                              alphas, lambdas, fl_update)
  }

  return(list(
    lambda_chain = lambdas,
    f_temp = fl_update,
    graph_chain = Gl,
    K_chain = Kl,
    alpha_chain = alphas,
    tau_chain = deltas,
    acceptprob_chain = acceptprob_chain,
    latent_zs = latent_zs,
    corr_gamma_chain = corr_gamma,
    xpred = xpred
  ))
}


#' Create Copula Buckets for Non-Continuous Data
#'
#' Generates data structures needed for copula modeling of non-continuous variables
#' by creating breakpoints and bucket indicators.
#'
#' @param xdf Data frame or matrix of observations for a single group
#' @return List containing:
#' \itemize{
#'   \item lessthan - Array [n x p x num_breakpoints] indicating values less than each breakpoint
#'   \item greaterthan - Array [n x p x num_breakpoints] indicating values greater than each breakpoint
#'   \item obs_buckets - Matrix [n x p] indicating bucket membership for each observation
#' }
create_copula_buckets <- function(xdf) {
  n <- nrow(xdf)
  p <- ncol(xdf)

  # Define breakpoints
  copula_breakpoints <- c(0, 25.5, 30.5, 35.5, 40.5, 45.5)
  binary_breakpoints <- c(0, 1)

  # Initialize result containers
  observation_buckets <- matrix(NA, n, p)
  copula_lessthan <- array(NA, dim = c(n, p, length(copula_breakpoints)))
  copula_greaterthan <- array(NA, dim = c(n, p, length(copula_breakpoints)))

  for (k in 1:p) {
    xvec <- xdf[, k]
    breakpoints <- if (k == 1 && !lbw) copula_breakpoints else binary_breakpoints
    num_breakpoints <- length(breakpoints)

    # Compute less-than and greater-than indicators
    for (i in 1:num_breakpoints) {
      right_endpoint <- ifelse(i == num_breakpoints, Inf, breakpoints[i + 1])
      copula_lessthan[, k, i] <- xvec < breakpoints[i]
      copula_greaterthan[, k, i] <- xvec > right_endpoint
    }

    # Assign observation buckets
    left_open <- !(k == 1 && !lbw)
    observation_buckets[, k] <- findInterval(xvec, breakpoints,
                                             rightmost.closed = FALSE,
                                             left.open = left_open) + 1
  }

  return(list(
    lessthan = copula_lessthan,
    greaterthan = copula_greaterthan,
    obs_buckets = observation_buckets
  ))
}

#' Gibbs Sampler for Lattice Model
#'
#' Main function implementing the Gibbs sampler for the lattice-based model that handles
#' both continuous and non-continuous (copula/probit) data cases.
#'
#' @param x List of observed data matrices (one per group)
#' @param steps Number of MCMC iterations
#' @param W Spatial adjacency matrix
#' @param prior Graph prior function
#' @param prior_params Parameters for graph prior
#' @param prior_ratio Prior ratio function for graph proposals
#' @param identifiable Logical indicating if identifiability constraints should be enforced (default TRUE)
#' @param start_lambda Starting values for lambda parameters (loadings)
#' @param start_alpha Starting values for alpha parameters (intercepts)
#' @param start_mus Starting values for mu parameters (factor means)
#' @param start_kw Starting value for Kw matrix (spatial precision)
#' @param start_f Starting values for factor scores
#' @param start_delta Starting values for delta parameters
#' @param start_gl Starting values for graph structures
#' @param start_kl Starting values for precision matrices
#' @param start_nu Starting value for nu parameter
#' @param start_theta Starting value for theta matrix
#' @param start_gamma Starting value for gamma matrix
#' @param start_corr_gamma Starting values for correlation matrices
#' @param start Logical indicating if starting values are provided (default FALSE)
#' @param non_continuous Logical indicating non-continuous data (default FALSE)
#' @param z Starting values for latent variables (required if non_continuous=TRUE)
#' @param copula Logical indicating copula model (default FALSE)
#' @return List containing MCMC chains for all parameters:
#' \itemize{
#'   \item lambda_chain - Loadings
#'   \item f_chain - Factor scores
#'   \item graph_chain - Graph structures
#'   \item K_chain - Precision matrices
#'   \item alpha_chain - Intercepts
#'   \item tau_chain - Delta parameters
#'   \item acceptprob_chain - Acceptance probabilities
#'   \item mu_chain - Factor means
#'   \item kw_chain - Spatial precision matrices
#'   \item theta_chain - Theta matrices
#'   \item gamma_chain - Gamma matrices
#'   \item nu_chain - Nu parameters
#'   \item delta_chain - Delta parameters
#'   \item z_chain - Latent variables (if non_continuous)
#'   \item corr_gamma_chain - Correlation matrices
#'   \item xpred - Predicted values
#' }
#' @export
gibbs_sampler_lattice <- function(x, steps, W, prior, prior_params, prior_ratio,
                                  identifiable = TRUE, start_lambda = NULL,
                                  start_alpha = NULL, start_mus = NULL,
                                  start_kw = NULL, start_f = NULL, start_delta = NULL,
                                  start_gl = NULL, start_kl = NULL, start_nu = NULL,
                                  start_theta = NULL, start_gamma = NULL,
                                  start_corr_gamma = NULL, has_start = FALSE,
                                  non_continuous = FALSE, z = NULL, copula = FALSE) {

  # Validate inputs
  if(non_continuous && is.null(z)) stop("Latent start required for non-continuous model")
  if(copula && !non_continuous) stop("Copula requires non-continuous variables")

  # Initialize dimensions
  L <- length(x)
  n <- sapply(x, nrow)
  p <- ncol(x[[1]])

  # Set hyperparameters
  params <- set_hyperparameters(p, L, W)

  # Initialize chains
  chains <- initialize_chains(steps, L, p, n, has_start, start_lambda, start_alpha,
                              start_mus, start_kw, start_f, start_delta, start_gl,
                              start_kl, start_nu, start_theta, start_gamma,
                              start_corr_gamma, non_continuous, z, params)


  # Initialize copula structures if needed
  copula_structs <- NULL
  if(copula) {
    copula_structs <- initialize_copula_structures(x, L)
  }

  # Initialize parallel processing
  cl <- setup_parallel_processing()

  # Main Gibbs sampling loop
  for(s in 2:steps) {
    # Update group parameters in parallel
    results <- update_groups_parallel(cl, L, chains, x, n, p, params,
                                      non_continuous, copula, copula_structs, s)

    # Update chains with results
    chains <- update_chains(chains, results, L, s)

    # Update shared parameters
    # if(!non_continuous) {
    #   chains <- update_shared_parameters(chains, params, s)
    # }

    # Update means and precision
    chains$mu_chain[s,] <- update_mus_list(chains$kw_chain[s-1,,],
                                           chains$mu_chain[s-1,],
                                           chains$f_temp, n)
    chains$kw_chain[s,,] <- update_kw(params$deltaw, params$rho, W,
                                      chains$mu_chain[s,])

    # Save progress periodically
    if(s %% 200 == 0) save_progress(chains, s)
  }

  stopCluster(cl)
  return(chains)
}

