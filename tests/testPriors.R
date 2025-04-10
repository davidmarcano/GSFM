# ---- Testing Setup ----
source("../priors.R")

set.seed(123)  # for reproducibility

# Number of nodes in the graph
p <- 5

# Initial empty graph (p x p matrix of zeros)
initial_graph <- matrix(0, p, p)

# Number of MCMC iterations
iter <- 1000

# ---- Define Prior Parameters ----

# Uniform prior: no parameters needed beyond p and iter
unif_params <- list(p = p, iter = iter)

# Bernoulli prior
bern_params <- list(p = p, iter = iter, phi = 0.3)

# Beta prior
beta_params <- list(p = p, iter = iter, a = 1, b = 1)

# Hierarchical prior
hier_params <- list(p = p, iter = iter)

# ---- Run Samplers ----

cat("Running uniform prior MCMC...\n")
result_unif <- mc_unif_prior(initial_graph, unif_params)

cat("Running Bernoulli prior MCMC...\n")
result_bern <- mc_bernoulli_prior(initial_graph, bern_params)

cat("Running Beta prior MCMC...\n")
result_beta <- mc_betas_prior(initial_graph, beta_params)

cat("Running hierarchical prior MCMC...\n")
result_hier <- mc_hierarchical_prior(initial_graph, hier_params)

# ---- Print Results ----

print_edge_probs <- function(result, name) {
  cat("\n---", name, "edge probabilities ---\n")
  print(round(result$edge_probs, 2))
}

print_edge_probs(result_unif, "Uniform")
print_edge_probs(result_bern, "Bernoulli")
print_edge_probs(result_beta, "Beta")
print_edge_probs(result_hier, "Hierarchical")

# ---- Visualize as heatmap ----

visualize_edge_probs <- function(edge_probs, title) {
  image(t(edge_probs[nrow(edge_probs):1, ]), main = title, col = heat.colors(10), axes = FALSE)
}

# Uncomment below to visualize
par(mfrow = c(2, 2))
visualize_edge_probs(result_unif$edge_probs, "Uniform Prior")
visualize_edge_probs(result_bern$edge_probs, "Bernoulli Prior")
visualize_edge_probs(result_beta$edge_probs, "Beta Prior")
visualize_edge_probs(result_hier$edge_probs, "Hierarchical Prior")
par(mfrow = c(1, 1))
