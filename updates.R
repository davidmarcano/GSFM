source("./priors.R")
source("./identifiability.R")

library(MASS)
library(invgamma)
library(parallel)
library(tmvtnorm)
library(truncnorm)
library(statmod)


# Function to update the factors (f) based on current values of lambda and alpha
update_factors <- function(K, x, lambda, alpha) {
  lamlam <- t(lambda) %*% lambda  # Calculate lambda' * lambda
  lamk <- (t(lambda) %*% solve(K) %*% lambda) / lamlam  # Calculate a factor for lambda
  sigma <- sqrt(lamk / (lamlam + lamk))  # Calculate the standard deviation for the normal distribution
  
  updated_f <- numeric(nrow(x))  # Initialize the updated factors vector
  
  for (j in 1:nrow(x)) {
    # Compute the mean for the normal distribution for each observation
    mu <- (t(lambda) %*% (x[j, ] - alpha)) / (lamlam + lamk)
    # Sample a new factor value from the normal distribution
    updated_f[j] <- rnorm(1, mean = mu, sd = sigma)
  }
  
  # Center the factors by subtracting their mean
  return(updated_f - mean(updated_f))
}

# Function to update the alphas (alpha) based on current values of lambda and f
update_alphas <- function(K, x, n, n0, f, lambda) {
  p <- ncol(x)  # Number of features (dimensions)
  covariance <- solve(n * K + n0 * diag(rep(1, p)))  # Compute the covariance matrix
  x_hat <- sweep(x, 1, lambda %*% t(f), "-")  # Center the data using lambda and f
  
  # Sample new alphas from a multivariate normal distribution
  new_alphas <- mvrnorm(1, mu = covariance %*% K %*% colSums(x_hat), Sigma = covariance)
  return(new_alphas)
}

# Function to update lambda values
update_lambdas <- function(K, x, f, delta, alpha, lambda) {
  n <- nrow(x)
  p <- ncol(x)
  x_hat <- sweep(x, 1, alpha, "-")  # Center the data by subtracting alpha
  
  deltafsum <- 1 / delta + sum(f^2)  # Calculate the scaling factor for the lambda update
  mu <- 1 / deltafsum * colSums(f * x_hat)  # Compute the mean for the multivariate normal
  covariance <- solve(deltafsum * K)  # Compute the covariance matrix
  Hprecision <- deltafsum * K  # Compute the precision matrix for truncated sampling
  
  Hprecision[lower.tri(Hprecision)] <- t(Hprecision)[lower.tri(Hprecision)]  # Ensure symmetry
  
  # Sample from the truncated multivariate normal distribution
  samp <- mvrnorm(1, mu = mu, Sigma = covariance)
  
  resamples <- 100  # Number of resampling attempts
  while ((samp[1] < 0) | is.nan(samp[1]) | is.infinite(samp[1])) {
    print("Truncated sampling failed")
    samp <- rtmvnorm(1, mean = mu, H = Hprecision, lower = c(0, rep(-Inf, p - 1)), algorithm = "gibbs")
    resamples <- resamples - 1
    if (resamples == 0) {
      print("Reached 100 resamples, breaking out")
      samp <- mvrnorm(1, mu = mu, Sigma = covariance)  # Fall back to regular normal sampling
      break
    }
  }
  
  return(matrix(samp, 1))
}

# Function to update delta
update_delta <- function(c, d, lambda, K) {
  p <- length(lambda)
  B <- chol(K)  # Cholesky decomposition of the precision matrix K
  lambda_hat <- B %*% t(lambda)  # Transform lambda
  shape <- c + 0.5 * p  # Shape parameter for the inverse gamma distribution
  scale <- 0.5 * (c * d + t(lambda_hat) %*% lambda_hat)  # Scale parameter for inverse gamma
  return(rinvgamma(1, shape = shape, rate = scale))  # Sample from inverse gamma distribution
}

# Function to calculate the shifted X matrix (X - alpha - lambda * f)
calculate_shifted_x <- function(x, alpha, lambda, f) {
  n <- nrow(x)
  x_hat <- sweep(x, 1, alpha + lambda * f, "-")  # Subtract alpha + lambda * f from each row
  return(x_hat)
}

# Function to update z values (latent variable)
update_zs <- function(K, G, x, lambda, alpha, f, z, copula = FALSE,
                      lessthan_buckets, greaterthan_buckets, buckets) {
  n <- nrow(z)
  p <- ncol(z)
  
  new_z <- matrix(0, n, p)  # Initialize the new z matrix
  
  if (copula) {
    # Handle copula case: truncating normal samples within buckets
    for (i in 1:p) {
      zi <- new_z[, i]
      sigma <- 1 / K[i, i]  # Standard deviation for the normal distribution
      
      # Define bounds for each bucket
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
      
      # Update z values for each observation using truncated normal distribution
      for (j in 1:n) {
        current_bucket <- buckets[j, i]
        mu <- alpha[i] + lambda[i] * f[j] - sum((K[i, ] / K[i, i] * (z[j, ] - alpha - lambda * f[j]))[G[i, ] == 1])
        new_z[j, i] <- rtruncnorm(1, a = lbs[current_bucket], b = ubs[current_bucket], mean = mu, sd = sigma)
      }
    }
  } else {
    # Standard case without copula: using regular truncation
    for (j in 1:n) {
      for (i in 1:p) {
        mu <- alpha[i] + lambda[i] * f[j] - sum((K[i, ] / K[i, i] * (z[j, ] - alpha - lambda * f[j]))[G[i, ] == 1])
        sigma <- sqrt(1 / K[i, i])
        
        # Define bounds for truncation based on the value of x
        lower_bound <- ifelse(x[j, i], 0, -Inf)
        upper_bound <- ifelse(x[j, i], Inf, 0)
        
        new_z[j, i] <- rtruncnorm(1, a = lower_bound, b = upper_bound, mean = mu, sd = sigma)
      }
    }
  }
  
  return(new_z)
}

# --- Functions to update spike-and-slab prior parameters --- #

# Function to update lambda values with a spike-and-slab prior
update_lambdas_spikeslab <- function(K, x, f, delta, alpha, lambda, pi_lambda) {
  n <- nrow(x)
  p <- ncol(x)
  x_hat <- sweep(x, 1, alpha, "-")
  
  deltafsum <- 1 / delta + sum(f^2)  # Scaling factor for update
  mu <- 1 / deltafsum * colSums(f * x_hat)  # Compute the mean
  covariance <- solve(deltafsum * K)  # Covariance matrix
  Hprecision <- deltafsum * K  # Precision matrix
  
  Hprecision[lower.tri(Hprecision)] <- t(Hprecision)[lower.tri(Hprecision)]  # Ensure symmetry
  
  samp <- rtmvnorm(1, mean = mu, H = Hprecision, lower = c(0, rep(-Inf, p - 1)), algorithm = "gibbs")
  
  resamples <- 100  # Number of resampling attempts
  while ((samp[1] < 0) | is.nan(samp[1]) | is.infinite(samp[1])) {
    print("Truncated sampling failed")
    samp <- rtmvnorm(1, mean = mu, H = Hprecision, lower = c(0, rep(-Inf, p - 1)), algorithm = "gibbs")
    resamples <- resamples - 1
    if (resamples == 0) {
      print("Reached 100 resamples, breaking out")
      samp <- mvrnorm(1, mu = mu, Sigma = covariance)  # Fall back to regular normal sampling
      break
    }
  }
  
  # Apply spike-and-slab prior by setting some lambdas to zero
  zeroes_lam <- rbinom(p, 1, pi_lambda)  # Binomial sampling for spike-and-slab
  samp[which(zeroes_lam == 0)] <- 0  # Set some lambdas to zero
  
  return(samp)
}

# Function to update the probability of lambda being non-zero (pi_lambda)
update_pi_prob <- function(lambda) {
  shape1 <- 1 + sum(lambda != 0)  # Shape parameter for the beta distribution (successes)
  shape2 <- 1 + length(lambda) - sum(lambda != 0)  # Shape parameter for the beta distribution (failures)
  return(rbeta(1, shape1, shape2))  # Sample from the beta distribution
}


# --- Functions to update GDP prior parameters --- #

# Calculate lambda_hat for GDP prior
calculate_lambda_hat_gdp <- function(psi_hat, f, x_hat_col) {
  s <- sum(f * x_hat_col)  # Sum over f and x_hat
  return(psi_hat * s)  # Return the product for lambda_hat
}

# Update lambdas for GDP prior
update_lambdas_gdp <- function(K, x, f, psi, alpha) {
  new_lambdas <- numeric(ncol(x))  # Initialize new lambda vector
  p <- ncol(x)
  x_hat <- sweep(x, 1, alpha, "-")  # Center the data by subtracting alpha
  
  # Compute covariance and mean for the update
  covar <- solve(sum(f^2) * K + solve(diag(psi)))
  mu <- as.numeric(covar %*% (K %*% colSums(f * x_hat)))
  
  # Sample new lambdas from a truncated normal distribution
  new_lambdas <- rtmvnorm(1, mean = mu, sigma = covar, lower = c(0, rep(-Inf, p - 1)), algorithm = "gibbs")
  
  resamples <- 100  # Number of resampling attempts
  while ((new_lambdas[1] < 0) | is.nan(new_lambdas[1]) | is.infinite(new_lambdas[1])) {
    print("Truncated sampling failed")
    new_lambdas <- rtmvnorm(1, mean = mu, sigma = covar, lower = c(0, rep(-Inf, p - 1)), algorithm = "gibbs")
    resamples <- resamples - 1
    if (resamples == 0) {
      print("Reached 100 resamples, breaking out")
      new_lambdas <- mvrnorm(1, mu = mu, Sigma = covar)  # Fall back to regular normal sampling
      break
    }
  }
  return(new_lambdas)
}
