# test-init_zs_copula.R
context("Copula Latent Variable Initialization - Improved Version")

test_that("Basic initialization works", {
  bin_data <- matrix(rbinom(100, 1, 0.5), ncol = 5)
  zs <- init_zs_copula(bin_data, rep("binary", 5))

  expect_equal(dim(zs), dim(bin_data))
  expect_true(is.numeric(zs))
  expect_false(any(is.na(zs)))
})

test_that("Mixed variable types are handled correctly", {
  mixed_data <- cbind(sample(1:5, 10), rbinom(10, 1, 0.5))
  zs <- init_zs_copula(mixed_data, c("ordinal", "binary"))

  # Check ordinal ranking
  expect_true(all(zs[mixed_data[,1] == 1, 1] < zs[mixed_data[,1] == 2, 1]))

  # Check binary constraints
  expect_true(all(zs[mixed_data[,2] == 0, 2] <= 0))
  expect_true(all(zs[mixed_data[,2] == 1, 2] >= 0))
})

test_that("Input validation works", {
  good_data <- matrix(rbinom(10, 1, 0.5), ncol = 2)

  # Missing var_types
  expect_error(init_zs_copula(good_data), "var_types argument must be specified")

  # Wrong length var_types
  expect_error(init_zs_copula(good_data, "binary"), "must match number of columns")

  # Invalid var_types
  expect_error(init_zs_copula(good_data, c("binary", "invalid")),
               "must contain only 'ordinal' or 'binary'")

  # NA values
  bad_data <- good_data
  bad_data[1,1] <- NA
  expect_error(init_zs_copula(bad_data, c("binary", "binary")),
               "contains missing values")
})

test_that("All ordinal initialization works", {
  ord_data <- cbind(sample(1:3, 10, replace = TRUE),
                    sample(1:4, 10, replace = TRUE))
  zs <- init_zs_copula(ord_data, c("ordinal", "ordinal"))

  # Check both columns maintain ordering
  expect_true(all(zs[ord_data[,1] == 1, 1] < zs[ord_data[,1] == 2, 1]))
  expect_true(all(zs[ord_data[,2] == 1, 2] < zs[ord_data[,2] == 2, 2]))
})

test_that("Edge cases are handled", {
  # Single row
  single_row <- matrix(c(1,0), nrow =
