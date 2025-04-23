library(BDgraph)

#' Initialize Latent Variables for Copula Model
#'
#' Initializes latent variables (zs) for a copula model with mixed variable types.
#' For ordinal variables, maintains rank ordering through truncation bounds.
#' For binary variables, uses simple truncation at 0.
#'
#' @param x Observed data matrix (n x p)
#' @param var_types Character vector specifying variable types:
#'        "ordinal" or "binary" for each column (length p)
#' @param init_sd Standard deviation for initialization normal distribution (default = 10)
#' @return Matrix of initialized latent variables (n x p)
#'
#' @details
#' For ordinal variables:
#' - Maintains observed rank ordering through truncation bounds
#' - Uses a large initialization variance (init_sd) to allow flexibility
#' For binary variables:
#' - Values initialized as N(0, init_sd) truncated below/above 0 depending on observed value
#'
#' @examples
#' # Binary data only
#' bin_data <- matrix(rbinom(100, 1, 0.5), ncol = 5)
#' zs <- init_zs_copula(bin_data, rep("binary", 5))
#'
#' # Mixed ordinal/binary
#' mixed_data <- cbind(sample(1:5, 100, replace = TRUE),
#'                    matrix(rbinom(400, 1, 0.5), ncol = 4))
#' zs <- init_zs_copula(mixed_data, c("ordinal", rep("binary", 4)))
#' @importFrom truncnorm rtruncnorm
init_zs_copula <- function(x, var_types, init_sd = 10) {
  # Input validation
  if (missing(var_types)) {
    stop("var_types argument must be specified")
  }
  if (ncol(x) != length(var_types)) {
    stop("Length of var_types must match number of columns in x")
  }
  if (!all(var_types %in% c("ordinal", "binary"))) {
    stop("var_types must contain only 'ordinal' or 'binary'")
  }
  if (any(is.na(x))) {
    stop("Input matrix x contains missing values")
  }

  n <- nrow(x)
  p <- ncol(x)
  zs <- matrix(NA, nrow = n, ncol = p)

  for (j in 1:p) {
    is_ordinal <- var_types[j] == "ordinal"

    for (i in 1:n) {
      if (is_ordinal) {
        # Handle ordinal variable with rank-preserving initialization
        current_value <- x[i, j]

        # Find bounds based on already initialized values
        lb <- if (any(x[1:i, j] < current_value, na.rm = TRUE)) {
          max(zs[x[1:i, j] < current_value, j], na.rm = TRUE)
        } else -Inf

        ub <- if (any(x[1:i, j] > current_value, na.rm = TRUE)) {
          min(zs[x[1:i, j] > current_value, j], na.rm = TRUE)
        } else Inf

        zs[i, j] <- truncnorm::rtruncnorm(1, a = lb, b = ub, mean = 0, sd = init_sd)

      } else {
        # Handle binary variable
        if (x[i, j] == 1) {
          zs[i, j] <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = 0, sd = init_sd)
        } else {
          zs[i, j] <- truncnorm::rtruncnorm(1, a = -Inf, b = 0, mean = 0, sd = init_sd)
        }
      }
    }
  }

  return(zs)
}

#' Generate Synthetic Data with Flexible Variable Types
#'
#' Creates simulated data for testing models with mixed variable types (continuous, ordinal, binary).
#'
#' @param L Number of groups (default: 3)
#' @param n Vector of sample sizes for each group (default: c(100, 100, 100))
#' @param var_types Character vector specifying variable types:
#'        "continuous", "ordinal", or "binary" (length p)
#' @param rho Spatial correlation parameter (default: 0.5)
#' @param W Optional adjacency matrix for spatial structure (default: NULL)
#' @param graph_density Density of true graph structure (default: 0.3)
#' @return List containing:
#' \itemize{
#'   \item data - List of observed data matrices
#'   \item latent - List of latent variables (for non-continuous cases) or NULL
#'   \item true_params - List of true parameter values used in generation
#' }
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom BDgraph rgwish
#' @importFrom stats rnorm rbinom
generate_synthetic_data <- function(L = 3, n = c(100, 100, 100),
                                    var_types = c("continuous", "binary", "binary"),
                                    rho = 0.5, W = NULL,
                                    graph_density = 0.3) {
  # Validate inputs
  p <- length(var_types)
  if (!all(var_types %in% c("continuous", "ordinal", "binary"))) {
    stop("var_types must contain only 'continuous', 'ordinal' or 'binary'")
  }
  if (length(n) != L) {
    stop("Length of n must equal L")
  }
  if (rho <= -1 || rho >= 1) {
    stop("rho must be between -1 and 1")
  }

  # Generate proximity weights matrix if not provided
  if (is.null(W)) {
    W <- matrix(0, L, L)
    if (L > 1) {
      for (i in 1:(L-1)) {
        W[i, i+1] <- W[i+1, i] <- 1  # Create chain structure
      }
    }
  } else if (nrow(W) != L || ncol(W) != L) {
    stop("W must be L x L matrix")
  }

  # Generate true parameters
  delta_mu <- 3
  I_L <- diag(L)
  true_kw <- if (L == 1) matrix(1) else
    BDgraph::rgwish(1, W, b = delta_mu, D = (delta_mu - 2) * solve(I_L - rho * W)) # this needs to be fixed
  true_mus <- MASS::mvrnorm(1, rep(0, L), solve(true_kw))

  # Initialize storage
  true_gl <- true_kl <- true_lambda <- true_alpha <- vector("list", L)
  x <- z <- vector("list", L)

  for (l in 1:L) {
    # Generate sparse graph structure
    true_gl[[l]] <- matrix(0, p, p)
    true_gl[[l]][upper.tri(true_gl[[l]])] <- rbinom(p*(p-1)/2, 1, graph_density)
    true_gl[[l]] <- true_gl[[l]] + t(true_gl[[l]])

    # Generate precision matrix
    true_kl[[l]] <- BDgraph::rgwish(1, true_gl[[l]], 3, diag(p))

    # Generate loadings and intercepts
    true_lambda[[l]] <- abs(rnorm(p, 1, 0.2))
    true_alpha[[l]] <- rnorm(p, 0, 0.5)

    # Generate factors
    f <- rnorm(n[l], true_mus[l], 1)

    # Generate latent continuous variables
    latent <- t(true_alpha[[l]] + true_lambda[[l]] %*% t(f) +
                  t(MASS::mvrnorm(n[l], rep(0, p), solve(true_kl[[l]]))))

    # Convert to specified variable types
    x[[l]] <- matrix(NA, n[l], p)
    z[[l]] <- if (any(var_types != "continuous")) latent else NULL

    for (j in 1:p) {
      if (var_types[j] == "continuous") {
        x[[l]][, j] <- latent[, j]
      } else if (var_types[j] == "ordinal") {
        # Create 5-level ordinal variable
        breaks <- qnorm(seq(0, 1, length.out = 6), mean(latent[, j]), sd(latent[, j]))
        x[[l]][, j] <- as.numeric(cut(latent[, j], breaks = breaks))
      } else { # binary
        x[[l]][, j] <- ifelse(latent[, j] > 0, 1, 0)
      }
    }
  }

  list(
    data = x,
    latent = if (any(var_types != "continuous")) z else NULL,
    true_params = list(
      gl = true_gl,
      kl = true_kl,
      mus = true_mus,
      kw = true_kw,
      lambda = true_lambda,
      alpha = true_alpha,
      rho = rho,
      W = W,
      var_types = var_types
    )
  )
}
#'
#' #' Test the IDCBF graph learning
#' test_idcbf <- function(x, true_gl, true_kl, p, L, n, steps = 100) {
#'   set.seed(123)
#'   results <- list()
#'
#'   sample_cov_list <- lapply(1:L, function(l){(n[l] - 1) * cov(x[[l]])})
#'
#'   prior <- GSFM::mc_hierarchical_prior
#'   prior_params <- list(iter = 1000, p = p)
#'   prior_params$m <- choose(p, 2)
#'   prior_ratio <- GSFM::hierarchy_prior_ratio
#'
#'   burnin <- round(steps/5)
#'   results <- mcmc_idcbf(sample_cov_list,
#'                         ns = n,
#'                         iter = steps,
#'                         burnin = burnin,
#'                         prior = prior,
#'                         prior_params = prior_params,
#'                         prior_ratio = prior_ratio,
#'                         df_prior = 3,
#'                         identifiable = TRUE,
#'                         nu = rep(1, L),
#'                         theta = 1)
#'   for(i in 1:(steps-burnin)) {
#'     for(j in 1:L) {
#'       results[[i]][[j]]$accuracy <- results[[i]][[j]]$G == true_gl[[j]]
#'     }
#'   }
#'   return(results)
#' }
#'
#' #' Test the full Gibbs sampler
#' test_gibbs_sampler <- function(x, L, n, p, W = NULL, steps = 500, continuous = TRUE,
#'                                z = NULL, copula = FALSE) {
#'   # Generate proximity weights matrix if not provided
#'   if (is.null(W)) {
#'     W <- matrix(0, L, L)
#'     if (L > 1) {
#'       for (i in 1:(L-1)) {
#'         W[i, i+1] <- W[i+1, i] <- 1  # Create chain structure
#'       }
#'     }
#'   } else if (nrow(W) != L || ncol(W) != L) {
#'     stop("W must be L x L matrix")
#'   }
#'
#'   # Run Gibbs sampler
#'   results <- gibbs_sampler_lattice(
#'     x = x,
#'     steps = steps,
#'     W = W,
#'     prior = prior_bernoulli,
#'     prior_params = list(prob = 0.1),
#'     prior_ratio = function(g1, g2) 1,
#'     identifiable = TRUE,
#'     non_continuous = !continuous,
#'     z = z,
#'     copula = copula,
#'     has_start = FALSE
#'   )
#'
#'   return(results)
#' }
#'
#' #' Test the spatial parameter estimation
#' test_spatial_parameters <- function(gibbs_results, true_params) {
#'   # Calculate posterior means
#'   post_mean_kw <- apply(gibbs_results$kw_chain, c(2,3), mean)
#'   post_mean_mu <- colMeans(gibbs_results$mu_chain)
#'
#'   # Compare with true values
#'   list(
#'     kw_correlation = cor(c(post_mean_kw), c(true_params$kw)),
#'     mu_correlation = cor(post_mean_mu, true_params$mus),
#'     kw_mae = mean(abs(post_mean_kw - true_params$kw)),
#'     mu_mae = mean(abs(post_mean_mu - true_params$mus))
#'   )
#' }
#'
#' #' Test the loading estimation
#' test_loadings <- function(gibbs_results, true_params) {
#'   L <- dim(gibbs_results$lambda_chain)[3]
#'   p <- dim(gibbs_results$lambda_chain)[2]
#'
#'   # Calculate posterior means
#'   post_mean_lambda <- apply(gibbs_results$lambda_chain, c(2,3), mean)
#'   post_mean_alpha <- apply(gibbs_results$alpha_chain, c(2,3), mean)
#'
#'   # Compare with true values
#'   lambda_cor <- alpha_cor <- numeric(L)
#'   for(l in 1:L) {
#'     lambda_cor[l] <- cor(post_mean_lambda[,l], true_params$lambda[[l]])
#'     alpha_cor[l] <- cor(post_mean_alpha[,l], true_params$alpha[[l]])
#'   }
#'
#'   list(
#'     lambda_correlation = mean(lambda_cor),
#'     alpha_correlation = mean(alpha_cor),
#'     lambda_mae = mean(abs(post_mean_lambda - do.call(cbind, true_params$lambda))),
#'     alpha_mae = mean(abs(post_mean_alpha - do.call(cbind, true_params$alpha)))
#'   )
#' }
#'
#'
#' #' Run comprehensive tests for the lattice model
#' run_comprehensive_tests <- function() {
#'   # Set seed for reproducibility
#'   set.seed(123)
#'
#'   # Test Case 1: Continuous data
#'   cat("Running tests for continuous data...\n")
#'   cont_data <- generate_synthetic_data(L = 3, n = c(100, 100, 100),
#'                                        var_types = rep("continuous", 5))
#'
#'   # Test IDCBF separately
#'   idcbf_results <- test_idcbf(
#'     cont_data$data,
#'     cont_data$true_params$gl,
#'     cont_data$true_params$kl,
#'     p = 5, L = 3, n = c(100, 100, 100)
#'   )
#'
#'   cat("IDCBF Test Results (Continuous):\n")
#'   print(sapply(idcbf_results, function(x) x$accuracy))
#'
#'   # Test full Gibbs sampler
#'   gibbs_results <- test_gibbs_sampler(
#'     x = cont_data$data,
#'     L = 3,
#'     n = c(100, 100, 100),
#'     p = 5,
#'     W = cont_data$true_params$W,
#'     steps = 500,
#'     continuous = TRUE
#'   )
#'
#'   # Test spatial parameters
#'   spatial_test <- test_spatial_parameters(gibbs_results, cont_data$true_params)
#'   cat("\nSpatial Parameter Test (Continuous):\n")
#'   print(spatial_test)
#'
#'   # Test loadings
#'   loading_test <- test_loadings(gibbs_results, cont_data$true_params)
#'   cat("\nLoading Test (Continuous):\n")
#'   print(loading_test)
#'
#'   # Test Case 2: Discrete data with copula
#'   cat("\nRunning tests for discrete data with copula...\n")
#'   disc_data <- generate_synthetic_data(L = 3, n = c(100, 100, 100), p = 5, continuous = FALSE)
#'
#'   gibbs_results_disc <- test_gibbs_sampler(
#'     x = disc_data$data,
#'     L = 3,
#'     n = c(100, 100, 100),
#'     p = 5,
#'     steps = 500,
#'     continuous = FALSE,
#'     z = disc_data$latent,
#'     copula = TRUE
#'   )
#'
#'   # Test loadings for discrete case
#'   loading_test_disc <- test_loadings(gibbs_results_disc, disc_data$true_params)
#'   cat("\nLoading Test (Discrete):\n")
#'   print(loading_test_disc)
#'
#'   # Return all results
#'   invisible(list(
#'     continuous = list(
#'       data = cont_data,
#'       dcbf_results = dcbf_results,
#'       gibbs_results = gibbs_results,
#'       tests = list(
#'         spatial = spatial_test,
#'         loadings = loading_test
#'       )
#'     ),
#'     discrete = list(
#'       data = disc_data,
#'       gibbs_results = gibbs_results_disc,
#'       tests = list(
#'         loadings = loading_test_disc
#'       )
#'     )
#'   ))
#' }
#'
#' # Run the comprehensive tests
#' test_results <- run_comprehensive_tests()
#'
#' #' Plot graph estimation results
#' plot_graph_results <- function(estimated, true, title = "") {
#'   par(mfrow = c(1, 2))
#'   pheatmap::pheatmap(estimated, cluster_rows = FALSE, cluster_cols = FALSE,
#'                      main = paste("Estimated", title))
#'   pheatmap::pheatmap(true, cluster_rows = FALSE, cluster_cols = FALSE,
#'                      main = paste("True", title))
#'   par(mfrow = c(1, 1))
#' }
#'
#' #' Plot trace plots for parameters
#' plot_traces <- function(chain, true_value = NULL, param_name = "") {
#'   plot(chain, type = "l", main = paste("Trace plot for", param_name),
#'        ylab = param_name, xlab = "Iteration")
#'   if(!is.null(true_value)) {
#'     abline(h = true_value, col = "red", lwd = 2)
#'   }
#' }
#'
#' #' Plot posterior distributions
#' plot_posteriors <- function(chain, true_value = NULL, param_name = "") {
#'   hist(chain, main = paste("Posterior distribution for", param_name),
#'        xlab = param_name, col = "lightblue")
#'   if(!is.null(true_value)) {
#'     abline(v = true_value, col = "red", lwd = 2)
#'   }
#' }
#'
#'
#' # Run tests (already done above)
#' # test_results <- run_comprehensive_tests()
#'
#' # Visualize results for continuous case
#' l <- 1  # Look at first group
#'
#' # Graph structure comparison
#' plot_graph_results(
#'   test_results$continuous$gibbs_results$g_chain[500,,,l],
#'   test_results$continuous$data$true_params$gl[[l]],
#'   "Graph Structure"
#' )
#'
#' # Precision matrix comparison
#' plot_graph_results(
#'   cov2cor(test_results$continuous$gibbs_results$k_chain[500,,,l]),
#'   cov2cor(test_results$continuous$data$true_params$kl[[l]]),
#'   "Precision Matrix (Correlation)"
#' )
#'
#' # Trace plots for spatial parameters
#' plot_traces(
#'   test_results$continuous$gibbs_results$kw_chain[,1,1],
#'   test_results$continuous$data$true_params$kw[1,1],
#'   "Spatial Precision (1,1)"
#' )
#'
#' # Posterior distribution for loadings
#' plot_posteriors(
#'   test_results$continuous$gibbs_results$lambda_chain[,1,l],
#'   test_results$continuous$data$true_params$lambda[[l]][1],
#'   "First Loading"
#' )
