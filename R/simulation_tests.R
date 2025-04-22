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

#' Generate Synthetic Data for Lattice Model
#'
#' Creates simulated data for testing lattice models, handling both continuous and discrete cases,
#' and properly supporting single-group (L=1) scenarios.
#'
#' @param L Number of groups (default: 3)
#' @param n Vector of sample sizes for each group (default: c(100, 100, 100))
#' @param p Number of variables (default: 5)
#' @param continuous Logical for continuous data (TRUE) or discrete (FALSE) (default: TRUE)
#' @param rho Spatial correlation parameter (default: 0.5)
#' @param W Optional adjacency matrix for spatial structure (default: NULL)
#' @return List containing:
#' \itemize{
#'   \item data - List of observed data matrices
#'   \item latent - List of latent variables (for discrete case) or NULL
#'   \item true_params - List of true parameter values used in generation
#' }
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom BDgraph rgwish
generate_synthetic_data <- function(L = 3, n = c(100, 100, 100), p = 5,
                                    continuous = TRUE, rho = 0.5, W = NULL) {

  # Validate inputs
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
      # Create simple chain structure for L > 1
      for (i in 1:(L-1)) {
        W[i, i+1] <- W[i+1, i] <- 1
      }
    }
  } else {
    if (nrow(W) != L || ncol(W) != L) {
      stop("W must be L x L matrix")
    }
  }

  # Generate true parameters
  true_params <- list()

  # 1. Generate spatial parameters (handle L=1 case specially)
  I_L <- diag(L)
  if (L == 1) {
    true_kw <- matrix(1, 1, 1)  # Single group case
  } else {
    true_kw <- (L - 1) * solve(I_L - rho * W)
  }
  true_mus <- MASS::mvrnorm(1, rep(0, L), solve(true_kw))

  # 2. Generate group-specific parameters
  true_gl <- true_kl <- true_lambda <- true_alpha <- vector("list", L)
  x <- z <- vector("list", L)

  for (l in 1:L) {
    # Generate graph structure (sparse upper triangular)
    true_gl[[l]] <- matrix(0, p, p)
    true_gl[[l]][upper.tri(true_gl[[l]])] <- rbinom(p*(p-1)/2, 1, 0.3)
    true_gl[[l]] <- true_gl[[l]] + t(true_gl[[l]])

    # Generate precision matrix
    true_kl[[l]] <- BDgraph::rgwish(1, true_gl[[l]], 3, diag(p))

    # Generate loadings and intercepts
    true_lambda[[l]] <- abs(rnorm(p, 1, 0.2))
    true_alpha[[l]] <- rnorm(p, 0, 0.5)

    # Generate factors
    f <- rnorm(n[l], true_mus[l], 1)

    # Generate data
    if (continuous) {
      # Continuous case
      epsilon <- MASS::mvrnorm(n[l], rep(0, p), solve(true_kl[[l]]))
      x[[l]] <- t(true_alpha[[l]] + true_lambda[[l]] * f + t(epsilon))
    } else {
      # Discrete case (probit model)
      z[[l]] <- matrix(NA, n[l], p)
      x[[l]] <- matrix(NA, n[l], p)

      latent <- t(true_alpha[[l]] + true_lambda[[l]] * f +
                    t(MASS::mvrnorm(n[l], rep(0, p), solve(true_kl[[l]])))

                  # First variable is ordinal, others binary
                  x[[l]][, 1] <- as.numeric(cut(latent[, 1],
                                                breaks = c(-Inf, -1, 0, 1, Inf)))
                  if (p > 1) {
                    x[[l]][, 2:p] <- ifelse(latent[, 2:p] > 0, 1, 0)
                  }
                  z[[l]] <- latent
    }
  }

  list(
    data = x,
    latent = if (!continuous) z else NULL,
    true_params = list(
      gl = true_gl,
      kl = true_kl,
      mus = true_mus,
      kw = true_kw,
      lambda = true_lambda,
      alpha = true_alpha,
      rho = rho,
      W = W
    )
  )
}

#' Test the DCBF graph learning
test_dcbf <- function(x, true_gl, true_kl, p, L, n, steps = 100) {
  results <- list()

  for(l in 1:L) {
    # Prepare data
    shifted <- scale(x[[l]], scale = FALSE)
    S <- (n[l]-1) * cov(shifted)

    # Run DCBF
    dcbf_result <- ggm_dcbf_multi(
      G = matrix(0, p, p),  # Start with empty graph
      D = diag(p),
      S = S,
      K = diag(p),
      p = p,
      delta = 3,
      deltan = 3 + n[l],
      nu = matrix(1, p, p),
      theta = matrix(0.5, L, L),
      l = l,
      prior = prior_bernoulli,  # Need to define this
      prior_params = list(prob = 0.1),
      prior_ratio = function(g1, g2) 1,
      identifiable = TRUE
    )

    # Compare with true graph
    results[[l]] <- list(
      estimated_graph = dcbf_result$G,
      true_graph = true_gl[[l]],
      precision = dcbf_result$K,
      true_precision = true_kl[[l]],
      accuracy = mean(dcbf_result$G == true_gl[[l]])
    )
  }

  return(results)
}

#' Test the full Gibbs sampler
test_gibbs_sampler <- function(x, L, n, p, steps = 500, continuous = TRUE,
                               z = NULL, copula = FALSE) {
  # Generate default W matrix if not provided
  W <- matrix(0, L, L)
  W[1,2] <- W[2,1] <- W[2,3] <- W[3,2] <- 1

  # Run Gibbs sampler
  results <- gibbs_sampler_lattice(
    x = x,
    steps = steps,
    W = W,
    prior = prior_bernoulli,
    prior_params = list(prob = 0.1),
    prior_ratio = function(g1, g2) 1,
    identifiable = TRUE,
    non_continuous = !continuous,
    z = z,
    copula = copula
  )

  return(results)
}

#' Test the spatial parameter estimation
test_spatial_parameters <- function(gibbs_results, true_params) {
  # Calculate posterior means
  post_mean_kw <- apply(gibbs_results$kw_chain, c(2,3), mean)
  post_mean_mu <- colMeans(gibbs_results$mu_chain)

  # Compare with true values
  list(
    kw_correlation = cor(c(post_mean_kw), c(true_params$kw)),
    mu_correlation = cor(post_mean_mu, true_params$mus),
    kw_mae = mean(abs(post_mean_kw - true_params$kw)),
    mu_mae = mean(abs(post_mean_mu - true_params$mus))
  )
}

#' Test the loading estimation
test_loadings <- function(gibbs_results, true_params) {
  L <- dim(gibbs_results$lambda_chain)[3]
  p <- dim(gibbs_results$lambda_chain)[2]

  # Calculate posterior means
  post_mean_lambda <- apply(gibbs_results$lambda_chain, c(2,3), mean)
  post_mean_alpha <- apply(gibbs_results$alpha_chain, c(2,3), mean)

  # Compare with true values
  lambda_cor <- alpha_cor <- numeric(L)
  for(l in 1:L) {
    lambda_cor[l] <- cor(post_mean_lambda[,l], true_params$lambda[[l]])
    alpha_cor[l] <- cor(post_mean_alpha[,l], true_params$alpha[[l]])
  }

  list(
    lambda_correlation = mean(lambda_cor),
    alpha_correlation = mean(alpha_cor),
    lambda_mae = mean(abs(post_mean_lambda - do.call(cbind, true_params$lambda))),
    alpha_mae = mean(abs(post_mean_alpha - do.call(cbind, true_params$alpha)))
  )
}


#' Run comprehensive tests for the lattice model
run_comprehensive_tests <- function() {
  # Set seed for reproducibility
  set.seed(123)

  # Test Case 1: Continuous data
  cat("Running tests for continuous data...\n")
  cont_data <- generate_synthetic_data(L = 3, n = c(100, 100, 100), p = 5, continuous = TRUE)

  # Test DCBF separately
  dcbf_results <- test_dcbf(
    cont_data$data,
    cont_data$true_params$gl,
    cont_data$true_params$kl,
    p = 5, L = 3, n = c(100, 100, 100)
  )

  cat("DCBF Test Results (Continuous):\n")
  print(sapply(dcbf_results, function(x) x$accuracy))

  # Test full Gibbs sampler
  gibbs_results <- test_gibbs_sampler(
    x = cont_data$data,
    L = 3,
    n = c(100, 100, 100),
    p = 5,
    steps = 500,
    continuous = TRUE
  )

  # Test spatial parameters
  spatial_test <- test_spatial_parameters(gibbs_results, cont_data$true_params)
  cat("\nSpatial Parameter Test (Continuous):\n")
  print(spatial_test)

  # Test loadings
  loading_test <- test_loadings(gibbs_results, cont_data$true_params)
  cat("\nLoading Test (Continuous):\n")
  print(loading_test)

  # Test Case 2: Discrete data with copula
  cat("\nRunning tests for discrete data with copula...\n")
  disc_data <- generate_synthetic_data(L = 3, n = c(100, 100, 100), p = 5, continuous = FALSE)

  gibbs_results_disc <- test_gibbs_sampler(
    x = disc_data$data,
    L = 3,
    n = c(100, 100, 100),
    p = 5,
    steps = 500,
    continuous = FALSE,
    z = disc_data$latent,
    copula = TRUE
  )

  # Test loadings for discrete case
  loading_test_disc <- test_loadings(gibbs_results_disc, disc_data$true_params)
  cat("\nLoading Test (Discrete):\n")
  print(loading_test_disc)

  # Return all results
  invisible(list(
    continuous = list(
      data = cont_data,
      dcbf_results = dcbf_results,
      gibbs_results = gibbs_results,
      tests = list(
        spatial = spatial_test,
        loadings = loading_test
      )
    ),
    discrete = list(
      data = disc_data,
      gibbs_results = gibbs_results_disc,
      tests = list(
        loadings = loading_test_disc
      )
    )
  ))
}

# Run the comprehensive tests
test_results <- run_comprehensive_tests()

#' Plot graph estimation results
plot_graph_results <- function(estimated, true, title = "") {
  par(mfrow = c(1, 2))
  pheatmap::pheatmap(estimated, cluster_rows = FALSE, cluster_cols = FALSE,
                     main = paste("Estimated", title))
  pheatmap::pheatmap(true, cluster_rows = FALSE, cluster_cols = FALSE,
                     main = paste("True", title))
  par(mfrow = c(1, 1))
}

#' Plot trace plots for parameters
plot_traces <- function(chain, true_value = NULL, param_name = "") {
  plot(chain, type = "l", main = paste("Trace plot for", param_name),
       ylab = param_name, xlab = "Iteration")
  if(!is.null(true_value)) {
    abline(h = true_value, col = "red", lwd = 2)
  }
}

#' Plot posterior distributions
plot_posteriors <- function(chain, true_value = NULL, param_name = "") {
  hist(chain, main = paste("Posterior distribution for", param_name),
       xlab = param_name, col = "lightblue")
  if(!is.null(true_value)) {
    abline(v = true_value, col = "red", lwd = 2)
  }
}


# Run tests (already done above)
# test_results <- run_comprehensive_tests()

# Visualize results for continuous case
l <- 1  # Look at first group

# Graph structure comparison
plot_graph_results(
  test_results$continuous$gibbs_results$g_chain[500,,,l],
  test_results$continuous$data$true_params$gl[[l]],
  "Graph Structure"
)

# Precision matrix comparison
plot_graph_results(
  cov2cor(test_results$continuous$gibbs_results$k_chain[500,,,l]),
  cov2cor(test_results$continuous$data$true_params$kl[[l]]),
  "Precision Matrix (Correlation)"
)

# Trace plots for spatial parameters
plot_traces(
  test_results$continuous$gibbs_results$kw_chain[,1,1],
  test_results$continuous$data$true_params$kw[1,1],
  "Spatial Precision (1,1)"
)

# Posterior distribution for loadings
plot_posteriors(
  test_results$continuous$gibbs_results$lambda_chain[,1,l],
  test_results$continuous$data$true_params$lambda[[l]][1],
  "First Loading"
)
