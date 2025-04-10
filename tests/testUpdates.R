library(testthat)  # Load testthat for unit testing
source("../updates.R")

# Test 1: Test the update_factors function
test_that("update_factors returns a vector of correct length", {
  K <- matrix(1, 3, 3)  # A simple 3x3 matrix for K
  x <- matrix(1, 5, 3)  # A 5x3 matrix for x
  lambda <- c(1, 1, 1)  # A 3-element vector for lambda
  alpha <- c(1, 1, 1)  # A 3-element vector for alpha
  
  updated_factors <- update_factors(K, x, lambda, alpha)
  
  # Check if the result is a vector with length 5 (number of rows in x)
  expect_equal(length(updated_factors), 5)
})

# Test 2: Test the update_alphas function
test_that("update_alphas returns a vector of correct length", {
  K <- matrix(1, 3, 3)  # A 3x3 matrix for K
  x <- matrix(1, 5, 3)  # A 5x3 matrix for x
  n <- 5  # Number of samples
  n0 <- 1  # Prior parameter for alpha
  f <- c(1, 1, 1, 1, 1)  # A 5-element vector for f
  lambda <- c(1, 1, 1)  # A 3-element vector for lambda
  
  updated_alphas <- update_alphas(K, x, n, n0, f, lambda)
  
  # Check if the result is a vector with length 3 (number of columns in x)
  expect_equal(length(updated_alphas), 3)
})

# Test 3: Test the update_lambdas function
test_that("update_lambdas returns a vector of correct length", {
  K <- matrix(1, 3, 3)  # A 3x3 matrix for K
  x <- matrix(1, 5, 3)  # A 5x3 matrix for x
  f <- c(1, 1, 1, 1, 1)  # A 5-element vector for f
  delta <- 1  # A scalar for delta
  alpha <- c(1, 1, 1)  # A 3-element vector for alpha
  lambda <- c(1, 1, 1)  # A 3-element vector for lambda
  
  updated_lambdas <- update_lambdas(K, x, f, delta, alpha, lambda)
  
  # Check if the result is a vector with the same length as lambda
  expect_equal(length(updated_lambdas), 3)
})

# Test 4: Test the update_lambdas_spikeslab function
test_that("update_lambdas_spikeslab returns a vector of correct length", {
  K <- matrix(1, 3, 3)  # A 3x3 matrix for K
  x <- matrix(1, 5, 3)  # A 5x3 matrix for x
  f <- c(1, 1, 1, 1, 1)  # A 5-element vector for f
  delta <- 1  # A scalar for delta
  alpha <- c(1, 1, 1)  # A 3-element vector for alpha
  lambda <- c(1, 1, 1)  # A 3-element vector for lambda
  pi_lambda <- 0.5  # Probability of lambda being non-zero
  
  updated_lambdas <- update_lambdas_spikeslab(K, x, f, delta, alpha, lambda, pi_lambda)
  
  # Check if the result is a vector with length 3 (same as lambda)
  expect_equal(length(updated_lambdas), 3)
})

# Test 5: Test the update_pi_prob function
test_that("update_pi_prob returns a scalar between 0 and 1", {
  lambda <- c(1, 0, 1)  # A vector with some zero entries for lambda
  
  pi_prob <- update_pi_prob(lambda)
  
  # Check if the result is a scalar between 0 and 1
  expect_true(is.numeric(pi_prob))
  expect_true(pi_prob >= 0 && pi_prob <= 1)
})

# Test 6: Test the update_delta function
test_that("update_delta returns a positive scalar", {
  c <- 1  # Shape parameter for the inverse gamma distribution
  d <- 1  # Rate parameter for the inverse gamma distribution
  lambda <- c(1, 2, 3)  # A vector for lambda
  K <- matrix(1, 3, 3)  # A 3x3 matrix for K
  
  updated_delta <- update_delta(c, d, lambda, K)
  
  # Check if the result is a positive scalar
  expect_true(is.numeric(updated_delta))
  expect_true(updated_delta > 0)
})

# Test 7: Test the calculate_shifted_x function
test_that("calculate_shifted_x returns a matrix of correct dimensions", {
  x <- matrix(1, 5, 3)  # A 5x3 matrix for x
  alpha <- c(1, 1, 1)  # A 3-element vector for alpha
  lambda <- c(1, 1, 1)  # A 3-element vector for lambda
  f <- c(1, 1, 1, 1, 1)  # A 5-element vector for f
  
  shifted_x <- calculate_shifted_x(x, alpha, lambda, f)
  
  # Check if the result has the same number of rows as x and columns as x
  expect_equal(dim(shifted_x), c(5, 3))
})

# Test 8: Test the update_zs function
test_that("update_zs returns a matrix of correct dimensions", {
  K <- matrix(1, 3, 3)  # A 3x3 matrix for K
  G <- matrix(1, 3, 3)  # A 3x3 matrix for G
  x <- matrix(1, 5, 3)  # A 5x3 matrix for x
  lambda <- c(1, 1, 1)  # A 3-element vector for lambda
  alpha <- c(1, 1, 1)  # A 3-element vector for alpha
  f <- c(1, 1, 1, 1, 1)  # A 5-element vector for f
  z <- matrix(1, 5, 3)  # A 5x3 matrix for z
  
  # Simulate buckets for testing
  lessthan_buckets <- matrix(1, 5, 3)
  greaterthan_buckets <- matrix(1, 5, 3)
  buckets <- matrix(1, 5, 3)
  
  updated_z <- update_zs(K, G, x, lambda, alpha, f, z, copula = FALSE,
                         lessthan_buckets, greaterthan_buckets, buckets)
  
  # Check if the result is a matrix with the same dimensions as z
  expect_equal(dim(updated_z), c(5, 3))
})

# Test 9: Test the calculate_lambda_hat_gdp function
test_that("calculate_lambda_hat_gdp returns a scalar", {
  psi_hat <- 1  # A scalar for psi_hat
  f <- c(1, 1, 1, 1, 1)  # A 5-element vector for f
  x_hat_col <- c(1, 1, 1, 1, 1)  # A 5-element vector for x_hat_col
  
  lambda_hat <- calculate_lambda_hat_gdp(psi_hat, f, x_hat_col)
  
  # Check if the result is a scalar
  expect_true(is.numeric(lambda_hat))
})

# Test 10: Test the update_lambdas_gdp function
test_that("update_lambdas_gdp returns a vector of correct length", {
  K <- matrix(1, 3, 3)  # A 3x3 matrix for K
  x <- matrix(1, 5, 3)  # A 5x3 matrix for x
  f <- c(1, 1, 1, 1, 1)  # A 5-element vector for f
  psi <- c(1, 1, 1)  # A 3-element vector for psi
  alpha <- c(1, 1, 1)  # A 3-element vector for alpha
  
  updated_lambdas <- update_lambdas_gdp(K, x, f, psi, alpha)
  
  # Check if the result is a vector with length 3 (same as lambda)
  expect_equal(length(updated_lambdas), 3)
})
