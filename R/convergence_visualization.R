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
# Samples is an iteration by number of parameters matrix
calculate_estimates <- function(samples, burnin = 0) {
  #burnin <- 0#nrow(samples)/2
  samples <- as.matrix(samples)
  estimates <- matrix(0, nrow(samples), ncol(samples))
  if(burnin > 0) {
    estimates <- matrix(0, burnin, ncol(samples))
  }
  for(i in 1:ncol(samples)) {
    samp <- samples[,i][(burnin+1):nrow(samples)]
    estimates[,i] <- cumsum(samp)/seq_along(samp)
  }
  return(estimates)
}
# Testing against other way to calculate cumulative posterior
###### Need convergence plots on the number of edges for the various runs
experiment_reps <- 3
steps <- length(all_chains[[1]])
r <- all_chains
number_edges <- matrix(0, steps, experiment_reps)
for(i in 1:length(r)) {
  run <- r[[i]]
  for(s in 1:length(run)) {
    g <- run[[s]][[l]]$G
    number_edges[s,i] <- sum(g[upper.tri(g)])
  }
}
expected_edges <- calculate_estimates(number_edges)

df <- data.frame(expected_edges, iters = 1:length(expected_edges[,1]))
df_long <- reshape2::melt(df, id.vars = "iters")

ggplot(df_long, aes(x = log(iters), y = value, run = variable, color = I("black"))) +
  geom_line(show.legend = FALSE, size = 1) +
  ylim(min(df_long$value) - 2, max(df_long$value) + 2) +
  labs(
    x = "Log Iteration",
    y = "Graph Size"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_text(size = 18, face = "plain", margin = margin(t = 20)),
    axis.title.y = element_text(size = 18, face = "plain", margin = margin(r = 20)),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 18, hjust = 0.5)
  )

chainid <- 3
group <- 3
numedges <- vector("numeric", length(all_chains[[chainid]]))
for(s in 1:length(all_chains[[chainid]])) {
  g <- all_chains[[chainid]][[s]][[group]]$G
  numedges[s] <- sum(g[upper.tri(g)])
}
plot(1:length(all_chains[[chainid]]), cumsum(numedges)/seq_along(numedges))

