# Markov Chain Sampler for Graph Priors from Dobra et al. (2011)

# ---- Helper Functions ----

# Compute the number of edges in an undirected graph
size <- function(g) {
  sum(g[upper.tri(g)])
}

# Uniform prior ratio: Pr(G1) / Pr(G2) = 1
unif_prior_ratio <- function(g1, g2, prior_params) {
  1
}

# Bernoulli prior ratio
# Pr(G1)/Pr(G2) = phi^(|G1|-|G2|) * (1-phi)^(|G2|-|G1|)
bern_prior_ratio <- function(g1, g2, prior_params) {
  phi <- prior_params$phi
  diff <- size(g1) - size(g2)
  phi^diff * (1 - phi)^(-diff)
}

# Beta prior ratio
beta_prior_ratio <- function(g1, g2, prior_params) {
  a <- prior_params$a
  b <- prior_params$b
  m <- prior_params$m
  beta(a + size(g1), b + m - size(g1)) / beta(a + size(g2), b + m - size(g2))
}

# Hierarchical prior ratio
hierarchy_prior_ratio <- function(g1, g2, prior_params) {
  m <- prior_params$m
  choose(m, size(g2)) / choose(m, size(g1))
}

# ---- Sampling Step ----

# Propose new graph by toggling a random edge, accept with Metropolis ratio
sample_prior_graph <- function(graph, calculate_prior_ratio, prior_params) {
  p <- prior_params$p
  candidate_edge <- sample(1:p, 2)
  i <- candidate_edge[1]
  j <- candidate_edge[2]

  candidate_graph <- graph
  candidate_graph[i, j] <- 1 - graph[i, j]
  candidate_graph[j, i] <- candidate_graph[i, j]

  ratio <- calculate_prior_ratio(candidate_graph, graph, prior_params)
  
  if (rbinom(1, 1, min(1, ratio))) {
    graph <- candidate_graph
  }

  graph
}

# ---- MCMC Wrappers ----

# Generic MCMC sampler
run_mc_sampler <- function(graph, prior_params, prior_ratio_fn) {
  p <- prior_params$p
  iter <- prior_params$iter
  prior_params$m <- choose(p, 2)  # total number of edges in a complete graph

  graphs_sum <- graph
  for (i in 1:iter) {
    graph <- sample_prior_graph(graph, prior_ratio_fn, prior_params)
    graphs_sum <- graphs_sum + graph
  }

  list(graph = graph, edge_probs = graphs_sum / iter)
}

# Samplers with specific priors
mc_unif_prior         <- function(graph, prior_params) run_mc_sampler(graph, prior_params, unif_prior_ratio)
mc_bernoulli_prior    <- function(graph, prior_params) run_mc_sampler(graph, prior_params, bern_prior_ratio)
mc_betas_prior        <- function(graph, prior_params) run_mc_sampler(graph, prior_params, beta_prior_ratio)
mc_hierarchical_prior <- function(graph, prior_params) run_mc_sampler(graph, prior_params, hierarchy_prior_ratio)

# ---- Direct Sampling from Priors ----

# Sample graph from uniform prior
sample_unif_prior <- function(graph, prior_params) {
  p <- prior_params$p
  upper <- rbinom((p - 1) * p / 2, 1, 0.5)
  g <- matrix(0, p, p)
  g[upper.tri(g)] <- upper
  g + t(g)
}

# Sample graph from Bernoulli prior
sample_bernoulli_prior <- function(graph, prior_params) {
  p <- prior_params$p
  phi <- prior_params$phi
  upper <- rbinom((p - 1) * p / 2, 1, phi)
  g <- matrix(0, p, p)
  g[upper.tri(g)] <- upper
  g + t(g)
}
