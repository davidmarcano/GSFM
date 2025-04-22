library(MASS)
library(stats)
library(invgamma)
library(tmvtnorm)
library(truncnorm)

#' Update Latent Factor Scores in Spatial Model
#'
#' Updates the latent factor scores (f) for observations in a spatial model using Gibbs sampling.
#' The function calculates conditional posterior distributions for each factor score and samples
#' from these distributions.
#'
#' @param K Current precision matrix (p x p)
#' @param x Data matrix (n x p) for the current group
#' @param lambda Current loading vector (length p)
#' @param alpha Current intercept vector (length p)
#' @param mu Current mean parameter for the factors
#'
#' @return A vector of updated factor scores (length n)
#'
#' @details
#' For each observation j, the function:
#' \enumerate{
#'   \item Calculates the conditional posterior mean for f[j]
#'   \item Calculates the conditional posterior variance
#'   \item Samples a new f[j] from the normal distribution
#'   \item Centers the final factor scores
#' }
#'
#' @examples
#' # Example usage:
#' # updated_f <- update_factors_spatial(K, x, lambda, alpha, mu)
#'
#' @export
update_factors_spatial <- function(K, x, lambda, alpha, mu) {
  # Precompute constant terms
  lambda_sq <- t(lambda) %*% lambda
  lambda_K_lambda <- (t(lambda) %*% solve(K) %*% lambda) / lambda_sq
  sigma_f <- sqrt(lambda_K_lambda / (lambda_sq + lambda_K_lambda))

  # Initialize updated factors
  n_obs <- nrow(x)
  updated_f <- numeric(n_obs)

  # Update each factor score
  for (j in 1:n_obs) {
    # Calculate conditional posterior mean
    numerator <- (t(lambda) %*% (x[j, ] - alpha)) + lambda_K_lambda * mu
    posterior_mean <- numerator / (lambda_sq + lambda_K_lambda)

    # Sample from conditional posterior
    updated_f[j] <- rnorm(1, mean = posterior_mean, sd = sigma_f)
  }

  # Center the factor scores
  updated_f <- updated_f - mean(updated_f)

  return(updated_f)
}


#' Update Group Means in Lattice Model
#'
#' Updates the mean parameters (mus) for each group in the lattice model using Gibbs sampling.
#' Handles both single-group and multi-group cases.
#'
#' @param Kw Current spatial precision matrix (L x L)
#' @param mus Current mean parameters (vector of length L)
#' @param f List of factor scores (length L, each element vector of size n_l)
#' @param n Vector of sample sizes for each group (length L)
#' @return Updated mean parameters (vector of length L)
update_mus_list <- function(Kw, mus, f, n) {
  L <- length(mus)
  new_mus <- numeric(L)

  if (L == 1) {
    # Special case for single group
    mu <- -sum(f[[1]])
    vars <- 1 / (n + Kw[1,1])
    new_mus <- rnorm(1, mean = mu, sd = sqrt(vars))
  } else {
    # General case for multiple groups
    for (l in 1:L) {
      # Calculate conditional mean and variance
      spatial_effect <- sum(Kw[l, -l] * mus[-l]) / Kw[l, l]
      data_effect <- sum(f[[l]])
      mu <- -spatial_effect + data_effect
      vars <- 1 / (n[l] + Kw[l, l])

      # Sample new mu
      new_mus[l] <- rnorm(1, mean = mu, sd = sqrt(vars))
    }
  }

  return(new_mus)
}

#' Update Spatial Precision Matrix
#'
#' Updates the spatial precision matrix Kw using a G-Wishart distribution.
#' Ensures numerical symmetry and handles single-group case.
#'
#' @param deltaw Degrees of freedom for G-Wishart prior
#' @param rho Spatial autocorrelation parameter
#' @param W Adjacency matrix for spatial structure (L x L)
#' @param mus Current mean parameters (vector of length L)
#' @return Updated spatial precision matrix (L x L)
update_kw <- function(deltaw, rho, W, mus) {
  L <- length(mus)
  I_L <- diag(1, L, L)

  # Calculate scale matrix D
  D <- (deltaw - 2) * solve(I_L - rho * W) + mus %*% t(mus)

  # Ensure numerical symmetry
  D[lower.tri(D)] <- t(D)[lower.tri(D)]

  # Handle single-group case
  if (L == 1) {
    W <- matrix(W, 1, 1)
  }

  # Sample from G-Wishart distribution
  rgwish(1, W, b = deltaw + 1, D = D)
}

#' Update the Factors
#'
#' This function updates the factors (f) based on the current values of lambda and alpha.
#'
#' @param K The precision matrix.
#' @param x The data matrix with rows representing observations.
#' @param lambda The current values of lambda.
#' @param alpha The current values of alpha.
#'
#' @return A numeric vector of updated factors (f).
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(20), 5, 4)
#' lambda <- rnorm(4)
#' alpha <- rnorm(4)
#' K <- diag(4)
#' update_factors(K, x, lambda, alpha)
update_factors <- function(K, x, lambda, alpha) {
  lamlam <- t(lambda) %*% lambda
  lamk <- (t(lambda) %*% solve(K) %*% lambda) / lamlam
  sigma <- sqrt(lamk / (lamlam + lamk))
  updated_f <- numeric(nrow(x))
  for (j in 1:nrow(x)) {
    mu <- (t(lambda) %*% (x[j, ] - alpha)) / (lamlam + lamk)
    updated_f[j] <- stats::rnorm(1, mean = mu, sd = sigma)
  }
  return(updated_f - mean(updated_f))
}

#' Update the Alphas
#'
#' This function updates the alphas (alpha) based on the current values of lambda and f.
#'
#' @param K The precision matrix.
#' @param x The data matrix.
#' @param n The number of observations.
#' @param n0 The prior parameter for alpha.
#' @param f The current factor values.
#' @param lambda The current values of lambda.
#'
#' @return A numeric vector of updated alphas.
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(20), 5, 4)
#' lambda <- rnorm(4)
#' f <- rnorm(5)
#' K <- diag(4)
#' update_alphas(K, x, 5, 1, f, lambda)
update_alphas <- function(K, x, n, n0, f, lambda) {
  p <- ncol(x)
  covariance <- solve(n * K + n0 * diag(rep(1, p)))
  x_hat <- x - t(lambda %*% t(f))
  new_alphas <- MASS::mvrnorm(1, mu = covariance %*% K %*% colSums(x_hat), Sigma = covariance)
  return(new_alphas)
}

#' Update Lambda Values
#'
#' This function updates the lambda values based on the current values of alpha and f.
#'
#' @param K The precision matrix.
#' @param x The data matrix.
#' @param f The current factor values.
#' @param delta The regularization parameter.
#' @param alpha The current values of alpha.
#' @param lambda The current values of lambda.
#'
#' @return A matrix of updated lambda values.
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(20), 5, 4)
#' f <- rnorm(5)
#' alpha <- rnorm(4)
#' lambda <- rnorm(4)
#' K <- diag(4)
#' update_lambdas(K, x, f, 1, alpha, lambda)
update_lambdas <- function(K, x, f, delta, alpha, lambda) {
  n <- nrow(x)
  p <- ncol(x)
  x_hat <- sweep(x, 1, alpha, "-")
  deltafsum <- 1 / delta + sum(f^2)
  mu <- 1 / deltafsum * colSums(f * x_hat)
  covariance <- solve(deltafsum * K)
  Hprecision <- deltafsum * K
  Hprecision[lower.tri(Hprecision)] <- t(Hprecision)[lower.tri(Hprecision)]
  samp <- MASS::mvrnorm(1, mu = mu, Sigma = covariance)
  resamples <- 100
  while ((samp[1] < 0) | is.nan(samp[1]) | is.infinite(samp[1])) {
    print("Truncated sampling failed")
    samp <- tmvtnorm::rtmvnorm(1, mean = mu, H = Hprecision, lower = c(0, rep(-Inf, p - 1)), algorithm = "gibbs")
    resamples <- resamples - 1
    if (resamples == 0) {
      print("Reached 100 resamples, breaking out")
      samp <- MASS::mvrnorm(1, mu = mu, Sigma = covariance)
      break
    }
  }
  return(matrix(samp, 1))
}

#' Update Delta
#'
#' This function updates the delta parameter based on the current values of lambda and the precision matrix K.
#'
#' @param c The shape parameter for the inverse gamma prior.
#' @param d The scale parameter for the inverse gamma prior.
#' @param lambda The current values of lambda.
#' @param K The precision matrix.
#'
#' @return A scalar representing the updated delta value.
#' @export
#'
#' @examples
#' set.seed(123)
#' lambda <- rnorm(4)
#' K <- diag(4)
#' update_delta(2, 1, lambda, K)
update_delta <- function(c, d, lambda, K) {
  p <- length(lambda)
  if (is.numeric(lambda)) {
    lambda <- matrix(lambda, 1, p)
  }
  B <- chol(K)
  lambda_hat <- B %*% t(lambda)
  shape <- c + 0.5 * p
  scale <- 0.5 * (c * d + t(lambda_hat) %*% lambda_hat)
  return(invgamma::rinvgamma(1, shape = shape, rate = scale))
}

#' Calculate the Shifted X Matrix
#'
#' This function computes the shifted X matrix by subtracting alpha and lambda * f from each row of the data.
#'
#' @param x The data matrix.
#' @param alpha The current values of alpha.
#' @param lambda The current values of lambda.
#' @param f The current factor values.
#'
#' @return The shifted data matrix.
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(20), 5, 4)
#' alpha <- rnorm(4)
#' lambda <- rnorm(4)
#' f <- rnorm(5)
#' calculate_shifted_x(x, alpha, lambda, f)
calculate_shifted_x <- function(x, alpha, lambda, f) {
  n <- nrow(x)
  x_hat <- sweep(x, 1, alpha + lambda * f, "-")
  return(x_hat)
}

#' Update the Z Values
#'
#' This function updates the latent variable z using truncated normal or copula-based sampling.
#'
#' @param K The precision matrix.
#' @param G The group indicator matrix.
#' @param x The data matrix.
#' @param lambda The current values of lambda.
#' @param alpha The current values of alpha.
#' @param f The current factor values.
#' @param z The current values of z.
#' @param copula A logical value indicating if copula sampling should be used.
#' @param lessthan_buckets A matrix of indices for values less than the bucket.
#' @param greaterthan_buckets A matrix of indices for values greater than the bucket.
#' @param buckets A matrix of bucket indices for each observation.
#'
#' @return A matrix of updated z values.
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 5; p <- 4
#' K <- diag(p)
#' G <- matrix(1, p, p)
#' x <- matrix(rbinom(n * p, 1, 0.5), n, p)
#' lambda <- rnorm(p)
#' alpha <- rnorm(p)
#' f <- rnorm(n)
#' z <- matrix(rnorm(n * p), n, p)
#' buckets <- matrix(sample(1:3, n * p, TRUE), n, p)
#' lessthan_buckets <- array(TRUE, dim = c(n, p, 3))
#' greaterthan_buckets <- array(TRUE, dim = c(n, p, 3))
#' update_zs(K, G, x, lambda, alpha, f, z, FALSE, lessthan_buckets, greaterthan_buckets, buckets)
update_zs <- function(K, G, x, lambda, alpha, f, z, copula = FALSE,
                      lessthan_buckets, greaterthan_buckets, buckets) {
  n <- nrow(z)
  p <- ncol(z)
  new_z <- matrix(0, n, p)

  if (copula) {
    for (i in 1:p) {
      zi <- new_z[, i]
      sigma <- 1 / K[i, i]
      unique_buckets <- sort(unique(buckets[, i]))
      num_buckets <- length(unique_buckets)
      lbs <- rep(-Inf, num_buckets)
      ubs <- rep(Inf, num_buckets)

      for (b in 1:num_buckets) {
        bucket <- unique_buckets[b]
        lessthan_current <- zi[lessthan_buckets[, i, bucket]]
        greaterthan_current <- zi[greaterthan_buckets[, i, bucket]]
        if (length(lessthan_current) > 0) lbs[b] <- max(lessthan_current)
        if (length(greaterthan_current) > 0) ubs[b] <- min(greaterthan_current)
      }

      for (j in 1:n) {
        current_bucket <- buckets[j, i]
        mu <- alpha[i] + lambda[i] * f[j] -
          sum((K[i, ] / K[i, i] * (z[j, ] - alpha - lambda * f[j]))[G[i, ] == 1])
        new_z[j, i] <- truncnorm::rtruncnorm(1, a = lbs[current_bucket], b = ubs[current_bucket], mean = mu, sd = sigma)
      }
    }
  } else {
    for (j in 1:n) {
      for (i in 1:p) {
        mu <- alpha[i] + lambda[i] * f[j] -
          sum((K[i, ] / K[i, i] * (z[j, ] - alpha - lambda * f[j]))[G[i, ] == 1])
        sigma <- sqrt(1 / K[i, i])
        lower_bound <- ifelse(x[j, i], 0, -Inf)
        upper_bound <- ifelse(x[j, i], Inf, 0)
        new_z[j, i] <- truncnorm::rtruncnorm(1, a = lower_bound, b = upper_bound, mean = mu, sd = sigma)
      }
    }
  }

  return(new_z)
}

#' Update Lambda Values with Spike-and-Slab Prior
#'
#' @inheritParams update_lambdas
#' @param pi_lambda The probability of lambda being non-zero.
#'
#' @return A matrix of updated lambda values.
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(20), 5, 4)
#' f <- rnorm(5)
#' alpha <- rnorm(4)
#' lambda <- rnorm(4)
#' K <- diag(4)
#' update_lambdas_spikeslab(K, x, f, 1, alpha, lambda, 0.8)
update_lambdas_spikeslab <- function(K, x, f, delta, alpha, lambda, pi_lambda) {
  n <- nrow(x)
  p <- ncol(x)
  x_hat <- sweep(x, 1, alpha, "-")
  deltafsum <- 1 / delta + sum(f^2)
  mu <- 1 / deltafsum * colSums(f * x_hat)
  covariance <- solve(deltafsum * K)
  Hprecision <- deltafsum * K
  Hprecision[lower.tri(Hprecision)] <- t(Hprecision)[lower.tri(Hprecision)]
  samp <- tmvtnorm::rtmvnorm(1, mean = mu, H = Hprecision, lower = c(0, rep(-Inf, p - 1)), algorithm = "gibbs")

  resamples <- 100
  while ((samp[1] < 0) | is.nan(samp[1]) | is.infinite(samp[1])) {
    print("Truncated sampling failed")
    samp <- tmvtnorm::rtmvnorm(1, mean = mu, H = Hprecision, lower = c(0, rep(-Inf, p - 1)), algorithm = "gibbs")
    resamples <- resamples - 1
    if (resamples == 0) {
      print("Reached 100 resamples, breaking out")
      samp <- MASS::mvrnorm(1, mu = mu, Sigma = covariance)
      break
    }
  }

  zeroes_lam <- stats::rbinom(p, 1, pi_lambda)
  samp[which(zeroes_lam == 0)] <- 0

  return(samp)
}

#' Update Probability of Lambda Being Non-Zero
#'
#' @param lambda The current values of lambda.
#'
#' @return A scalar representing the updated probability of lambda being non-zero.
#' @export
#'
#' @examples
#' set.seed(123)
#' lambda <- c(0, 1.2, 0, 0.8)
#' update_pi_prob(lambda)
update_pi_prob <- function(lambda) {
  shape1 <- 1 + sum(lambda != 0)
  shape2 <- 1 + length(lambda) - sum(lambda != 0)
  return(stats::rbeta(1, shape1, shape2))
}
