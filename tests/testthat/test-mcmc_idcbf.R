library(stats)
library(testthat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

test_that("mcmc_idcbf sampler runs and returns correct structure", {
  set.seed(123)
  p <- 5
  L <- 3
  ns <- c(10, 15, 50)
  Sigmas <- stats::rWishart(length(ns), p, diag(p))
  sample_cov_list <- lapply(seq(dim(Sigmas)[3]), function(x) Sigmas[ , , x])

  prior <- GSFM::mc_hierarchical_prior
  prior_params <- list(iter = 1000, p = p)
  prior_params$m <- choose(p, 2)
  prior_ratio <- GSFM::hierarchy_prior_ratio

  results <- mcmc_idcbf(sample_cov_list,
                                ns = ns,
                                iter = 10,
                                burnin = 5,
                                prior = prior,
                                prior_params = prior_params,
                                prior_ratio = prior_ratio,
                                df_prior = 3,
                                identifiable = TRUE,
                                nu = rep(1, L),
                                theta = 1)

  expect_type(results, "list")
  expect_length(results, 5)
  expect_type(results[[1]], "list")
  expect_length(results[[1]], L)

  for (group_result in results[[1]]) {
    expect_named(group_result, c("G", "K", "alphas"))
    expect_true(is.matrix(group_result$G))
    expect_true(is.matrix(group_result$K))
  }
})



test_that("Posterior edge probabilities are consistent across multiple chains", {
  set.seed(123)

  p <- 5
  L <- 3
  ns <- c(50, 100, 500)
  g1 <- matrix(0, p, p)
  g2 <- g1
  g2[p-1, p] <- g2[p, p-1] <- 1
  g3 <- g2
  g3[p - 4, 2] <- g3[2, p-4] <- 1
  D <- diag(p)
  gs <- list(g1, g2, g3)
  Precisions <- lapply(gs, function(g) rgwish(1, g, D = D))

  x <- vector("list", L)
  sample_cov_list <- vector("list", L)
  for(l in 1:L) {
    x[[l]] <- matrix(nrow = ns[l], ncol = p)
    for(i in 1:ns[l]) {
      x[[l]][i,] <- mvrnorm(1, mu = rep(0, p), Sigma = solve(Precisions[[l]]) )
    }
    sample_cov_list[[l]] <- (ns[l]-1)*cov(x[[l]])
  }

  n_chains <- 3
  iter <- 1000
  burnin <- iter/2
  prior <- GSFM::mc_hierarchical_prior
  prior_params <- list(iter = 1000, p = p)
  prior_params$m <- choose(p, 2)
  prior_ratio <- GSFM::hierarchy_prior_ratio
  nu <- 1
  theta <- 1

  # Run multiple chains
  all_chains <- vector("list", n_chains)
  for (chain_id in 1:n_chains) {
    set.seed(1000 + chain_id)  # Different seed per chain
    all_chains[[chain_id]] <- mcmc_idcbf(
      sample_cov_list,
      ns = ns,
      iter = iter,
      burnin = burnin,
      prior = prior,
      prior_params = prior_params,
      prior_ratio = prior_ratio,
      df_prior = 3,
      identifiable = FALSE,
      nu = nu,
      theta = theta
    )
  }

  # Compute cumulative edge probabilities for each chain and group
  edge_probs_long <- list()
  for (chain_id in 1:n_chains) {
    results <- all_chains[[chain_id]]

    for (l in 1:L) {
      edge_counts <- matrix(0, p, p)
      cum_list <- list()

      for (t in seq_along(results)) {
        G <- results[[t]][[l]]$G
        edge_counts <- edge_counts + (G != 0)
        probs <- edge_counts / t
        probs[lower.tri(probs, diag = TRUE)] <- NA
        df <- as.data.frame(as.table(probs))
        colnames(df) <- c("i", "j", "prob")
        df <- df[!is.na(df$prob), ]
        df$group <- paste0("Group ", l)
        df$iter <- t
        df$chain <- paste0("Chain ", chain_id)
        cum_list[[t]] <- df
      }

      edge_probs_long[[paste0("chain", chain_id, "_group", l)]] <- do.call(rbind, cum_list)
    }
  }

  df_all <- bind_rows(edge_probs_long)

  # Pick top 3 most probable edges across all chains and groups
  top_edges <- df_all %>%
    filter(iter == iter) %>%
    group_by(group, i, j) %>%
    summarize(avg_prob = mean(prob), .groups = "drop") %>%
    arrange(desc(avg_prob)) %>%
    slice_head(n = 4) %>%
    mutate(edge = paste0("(", i, ",", j, ")")) %>%
    pull(edge)

  df_all <- df_all %>%
    mutate(edge = paste0("(", i, ",", j, ")")) %>%
    filter(edge %in% top_edges)

  # Plot
  plot_path <- tempfile(fileext = ".png")
  p_plot <- ggplot(df_all, aes(x = iter, y = prob, color = chain, linetype = chain)) +
    geom_line(size = 0.9) +
    facet_grid(group ~ edge) +
    labs(
      title = "Cumulative Edge Probability Convergence (Multiple Chains)",
      x = "Iteration", y = "Cumulative Probability",
      color = "Chain", linetype = "Chain"
    ) +
    theme_minimal()
  ggsave(plot_path, p_plot, width = 12, height = 6)
  cat("Multi-chain convergence plot saved to:", plot_path, "\n")

  # Optional check: make sure final probabilities are consistent across chains
  final_probs <- df_all %>%
    filter(iter == iter) %>%
    group_by(group, edge) %>%
    summarize(sd_prob = sd(prob), .groups = "drop")

  expect_lt(max(final_probs$sd_prob), 0.1)
})




