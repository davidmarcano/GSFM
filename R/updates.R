library(MASS)
library(stats)
library(invgamma)
library(tmvtnorm)
library(truncnorm)

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
