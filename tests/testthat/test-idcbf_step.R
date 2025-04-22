test_that("idcbf_step returns valid graph and precision matrix", {
  set.seed(123)
  p <- 4
  G <- diag(p)
  diag(G) <- 0
  D <- diag(p)
  S <- diag(p) * 0.5
  K <- diag(p)
  delta <- 3
  deltan <- 8
  nu <- 1
  theta <- 1
  l <- 1

  prior <- mc_unif_prior
  prior_params <- list(iter = 1000, p = p)
  prior_ratio <- unif_prior_ratio # flat ratio

  result <- idcbf_step(Gs = G, D = D, S = S, K = K, p = p,
                           delta = delta, deltan = deltan, nu = nu,
                           theta = theta, l = l, prior = prior,
                           prior_params = prior_params,
                           prior_ratio = prior_ratio,
                           identifiable = TRUE)

  expect_type(result, "list")
  expect_named(result, c("G", "K", "alphas"))
  expect_true(is.matrix(result$G))
  expect_true(is.matrix(result$K))
  expect_equal(dim(result$G), c(p, p))
  expect_equal(dim(result$K), c(p, p))
  expect_true(all(diag(result$G) == 0))
  expect_true(is_identifiable(result$G))
})


test_that("permute() and inv_permute() work correctly", {
  set.seed(123)

  # Generate a symmetric test matrix
  p <- 5
  G <- matrix(sample(0:1, p^2, replace = TRUE), p, p)
  G[lower.tri(G)] <- 0
  G <- (G + t(G))
  diag(G) <- 0

  # Pick two nodes to permute
  i <- 2
  j <- 4

  # Apply permutation
  permuted <- permute(G, i, j)
  G_perm <- permuted$G
  perm <- permuted$perm

  # Check: i and j are now in positions (p-1, p)
  expect_equal(G_perm[p - 1, p], G[i, j])
  expect_equal(G_perm[p, p - 1], G[j, i])

  # Invert permutation
  G_recovered <- inv_permute(G_perm, perm)

  # Check if the original matrix is recovered
  expect_equal(G_recovered, G)
})

