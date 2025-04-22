#' Generate synthetic data for lattice model testing
#'
#' @param L Number of groups
#' @param n Vector of sample sizes for each group
#' @param p Number of variables
#' @param continuous Logical, whether to generate continuous data
#' @param rho Spatial correlation parameter
#' @param W Adjacency matrix for spatial structure
#' @return List containing generated data and true parameters
generate_synthetic_data <- function(L = 3, n = c(100, 100, 100), p = 5,
                                    continuous = TRUE, rho = 0.5, W = NULL) {

  # Generate proximity weights matrix if not provided
  if(is.null(W)) {
    W <- matrix(0, L, L)
    W[1,2] <- W[2,1] <- 1
    W[2,3] <- W[3,2] <- 1
  }

  # Generate true parameters
  true_params <- list()

  # 1. Generate spatial parameters
  I_L <- diag(L)
  true_kw <- (L - 1) * solve(I_L - rho * W)
  true_mus <- MASS::mvrnorm(1, rep(0, L), solve(true_kw))

  # 2. Generate group-specific parameters
  true_gl <- true_kl <- true_lambda <- true_alpha <- vector("list", L)
  x <- z <- vector("list", L)

  for(l in 1:L) {
    # Generate graph structure (triangular)
    true_gl[[l]] <- matrix(0, p, p)
    true_gl[[l]][upper.tri(true_gl[[l]])] <- rbinom(p*(p-1)/2, 1, 0.3)
    true_gl[[l]] <- true_gl[[l]] + t(true_gl[[l]])

    # Generate precision matrix
    true_kl[[l]] <- rgwish(1, true_gl[[l]], 3, diag(p))

    # Generate loadings and intercepts
    true_lambda[[l]] <- abs(rnorm(p, 1, 0.2))
    true_alpha[[l]] <- rnorm(p, 0, 0.5)

    # Generate factors
    f <- rnorm(n[l], true_mus[l], 1)

    # Generate continuous or discrete data
    if(continuous) {
      # Continuous case
      epsilon <- MASS::mvrnorm(n[l], rep(0,p), solve(true_kl[[l]]))
      x[[l]] <- t(true_alpha[[l]] + true_lambda[[l]] * f + t(epsilon))
    } else {
      # Discrete case (probit model)
      z[[l]] <- matrix(NA, n[l], p)
      x[[l]] <- matrix(NA, n[l], p)

      latent <- t(true_alpha[[l]] + true_lambda[[l]] * f +
                    t(MASS::mvrnorm(n[l], rep(0,p), solve(true_kl[[l]]))))

      for(j in 1:p) {
        if(j == 1) {
          # First variable is ordinal (age-like)
          breaks <- c(-Inf, -1, 0, 1, Inf)
          x[[l]][,j] <- as.numeric(cut(latent[,j], breaks))
        } else {
          # Other variables are binary
          x[[l]][,j] <- ifelse(latent[,j] > 0, 1, 0)
        }
      }
      z[[l]] <- latent
    }
  }

  list(
    data = x,
    latent = if(!continuous) z else NULL,
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
