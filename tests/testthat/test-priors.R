library(testthat)

# ---- Setup ----
set.seed(123)
p <- 5
iter <- 1000
initial_graph <- matrix(0, p, p)

# ---- Prior Parameter Lists ----
unif_params <- list(p = p, iter = iter)
bern_params <- list(p = p, iter = iter, phi = 0.3)
beta_params <- list(p = p, iter = iter, a = 1, b = 1)
hier_params <- list(p = p, iter = iter)

# ---- Helper Function ----
expect_edge_probs <- function(result, label) {
  test_that(paste(label, "- edge_probs is correct shape and numeric"), {
    expect_true(is.matrix(result$edge_probs))
    expect_equal(dim(result$edge_probs), c(p, p))
    expect_type(result$edge_probs, "double")
  })
}

# ---- Test Each Prior ----

test_that("Uniform prior MCMC runs and returns expected output", {
  result_unif <- mc_unif_prior(initial_graph, unif_params)
  expect_edge_probs(result_unif, "Uniform")
})

test_that("Bernoulli prior MCMC runs and returns expected output", {
  result_bern <- mc_bernoulli_prior(initial_graph, bern_params)
  expect_edge_probs(result_bern, "Bernoulli")
})

test_that("Beta prior MCMC runs and returns expected output", {
  result_beta <- mc_betas_prior(initial_graph, beta_params)
  expect_edge_probs(result_beta, "Beta")
})

test_that("Hierarchical prior MCMC runs and returns expected output", {
  result_hier <- mc_hierarchical_prior(initial_graph, hier_params)
  expect_edge_probs(result_hier, "Hierarchical")
})
