library(testthat)  # Load testthat for unit testing
library(BDgraph)

p <- 5
egraph <- matrix(0, p, p)
D <- diag(1, p, p)

# Test 1: Test the update_factors function
test_that("update_factors returns a vector of correct length", {
  K <- BDgraph::rgwish(1, egraph, D = D)  # A simple pxp matrix for K
  x <- matrix(1, 5, p)  # A 5xp matrix for x
  lambda <- rep(1, times = p)  # A p-element vector for lambda
  alpha <- rep(1, p)  # A p-element vector for alpha

  updated_factors <- update_factors(K, x, lambda, alpha)

  # Check if the result is a vector with length 5 (number of rows in x)
  expect_equal(length(updated_factors), 5)
})

# Test 2: Test the update_alphas function
test_that("update_alphas returns a vector of correct length", {
  K <- BDgraph::rgwish(1, egraph, D = D)  # A pxp matrix for K
  x <- matrix(1, 5, p)  # A 5xp matrix for x
  n <- 5  # Number of samples
  n0 <- 1  # Prior parameter for alpha
  f <- c(1, 1, 1, 1, 1)  # A 5-element vector for f
  lambda <- rep(1, p)  # A p-element vector for lambda

  updated_alphas <- update_alphas(K, x, n, n0, f, lambda)

  # Check if the result is a vector with length p (number of columns in x)
  expect_equal(length(updated_alphas), p)
})

# Test p: Test the update_lambdas function
test_that("update_lambdas returns a vector of correct length", {
  K <- BDgraph::rgwish(1, egraph, D = D)  # A pxp matrix for K
  x <- matrix(1, 5, p)  # A 5xp matrix for x
  f <- c(1, 1, 1, 1, 1)  # A 5-element vector for f
  delta <- 1  # A scalar for delta
  alpha <- rep(1, p)  # A p-element vector for alpha
  lambda <- rep(1, p)  # A p-element vector for lambda

  updated_lambdas <- update_lambdas(K, x, f, delta, alpha, lambda)

  # Check if the result is a vector with the same length as lambda
  expect_equal(length(updated_lambdas), p)
})

# Test 4: Test the update_lambdas_spikeslab function
test_that("update_lambdas_spikeslab returns a vector of correct length", {
  K <- BDgraph::rgwish(1, egraph, D = D)  # A pxp matrix for K
  x <- matrix(1, 5, p)  # A 5xp matrix for x
  f <- c(1, 1, 1, 1, 1)  # A 5-element vector for f
  delta <- 1  # A scalar for delta
  alpha <- rep(1, p)  # A p-element vector for alpha
  lambda <- rep(1, p)  # A p-element vector for lambda
  pi_lambda <- 0.5  # Probability of lambda being non-zero

  updated_lambdas <- update_lambdas_spikeslab(K, x, f, delta, alpha, lambda, pi_lambda)

  # Check if the result is a vector with length p (same as lambda)
  expect_equal(length(updated_lambdas), p)
})

# Test 5: Test the update_pi_prob function
test_that("update_pi_prob returns a scalar between 0 and 1", {
  lambda <- rep(1, p)  # A vector with some zero entries for lambda

  pi_prob <- update_pi_prob(lambda)

  # Check if the result is a scalar between 0 and 1
  expect_true(is.numeric(pi_prob))
  expect_true(pi_prob >= 0 && pi_prob <= 1)
})

# Test 6: Test the update_delta function
test_that("update_delta returns a positive scalar", {
  c <- 1  # Shape parameter for the inverse gamma distribution
  d <- 1  # Rate parameter for the inverse gamma distribution
  lambda <- rep(1, p)  # A vector for lambda
  K <- BDgraph::rgwish(1, egraph, D = D)  # A pxp matrix for K

  updated_delta <- update_delta(c, d, lambda, K)

  # Check if the result is a positive scalar
  expect_true(is.numeric(updated_delta))
  expect_true(updated_delta > 0)
})

# Test 7: Test the calculate_shifted_x function
test_that("calculate_shifted_x returns a matrix of correct dimensions", {
  x <- matrix(1, 5, p)  # A 5xp matrix for x
  alpha <- rep(1, p)  # A p-element vector for alpha
  lambda <- rep(1, p)  # A p-element vector for lambda
  f <- c(1, 1, 1, 1, 1)  # A 5-element vector for f

  shifted_x <- calculate_shifted_x(x, alpha, lambda, f)

  # Check if the result has the same number of rows as x and columns as x
  expect_equal(dim(shifted_x), c(5, p))
})

# Test 8: Test the update_zs function
test_that("update_zs returns a matrix of correct dimensions", {
  K <- BDgraph::rgwish(1, egraph, D = D)  # A pxp matrix for K
  G <- matrix(1, p, p)  # A pxp matrix for G
  x <- matrix(1, 5, p)  # A 5xp matrix for x
  lambda <- rep(1, p)  # A p-element vector for lambda
  alpha <- rep(1, p)  # A p-element vector for alpha
  f <- c(1, 1, 1, 1, 1)  # A 5-element vector for f
  z <- matrix(1, 5, p)  # A 5xp matrix for z

  # Simulate buckets for testing
  lessthan_buckets <- matrix(1, 5, p)
  greaterthan_buckets <- matrix(1, 5, p)
  buckets <- matrix(1, 5, p)

  updated_z <- update_zs(K, G, x, lambda, alpha, f, z, copula = FALSE,
                         lessthan_buckets, greaterthan_buckets, buckets)

  # Check if the result is a matrix with the same dimensions as z
  expect_equal(dim(updated_z), c(5, p))
})

