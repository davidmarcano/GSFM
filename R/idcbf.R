library(BDgraph)


#' @title Permute Nodes to Target Edge Position
#' @description Moves nodes \code{i} and \code{j} to the \code{(p-1, p)} and \code{(p, p-1)} positions in a matrix.
#' @param G A square matrix (e.g., adjacency or precision matrix).
#' @param i Index of the first node to permute.
#' @param j Index of the second node to permute.
#' @return A list with:
#' \describe{
#'   \item{G}{The permuted matrix.}
#'   \item{perm}{The permutation vector applied to rows and columns.}
#' }
#' @keywords internal
permute <- function(G, i, j) {
  stopifnot(i != j, is.matrix(G), nrow(G) == ncol(G))
  p <- nrow(G)

  # Build target permutation
  nodes <- seq_len(p)
  rest <- setdiff(nodes, c(i, j))
  perm <- c(rest, i, j)

  G_perm <- G[perm, perm]
  return(list(G = G_perm, perm = perm))
}

#' @title Inverse Permutation of Permuted Matrix
#' @description Applies the inverse of a previously stored permutation to return matrix to original node order.
#' @param G A permuted square matrix.
#' @param perm The permutation vector used to permute the matrix.
#' @return The matrix with inverse permutation applied.
#' @keywords internal
inv_permute <- function(G, perm) {
  stopifnot(is.matrix(G), nrow(G) == ncol(G), length(perm) == nrow(G))
  G_inv <- matrix(0, nrow(G), ncol(G))
  G_inv[perm, perm] <- G
  return(G_inv)
}


#' @title Normalization Constant for G-Wishart Ratio
#' @description Computes a normalization term used in computing acceptance probabilities.
#' @param K A precision matrix.
#' @param D A symmetric matrix.
#' @return A scalar normalization constant.
#' @keywords internal
N_KU <- function(K, D) {
  p <- nrow(K)
  r <- chol(K)
  term1 <- r[p - 1, p - 1] * sqrt(2 * pi / D[p, p])
  lhs <- (r[p - 1, p - 1] * D[p - 1, p]) / D[p, p]
  rhs <- sum(r[-c(p, p - 1), p - 1] * r[-c(p, p - 1), p])/r[p-1, p-1]
  term2 <- exp(0.5 * D[p, p] * ((lhs - rhs)^2))
  return(term1 * term2)
}

#' @title Identifiable Double Conditional Bayes Factor Sampler for Multiple Gaussian Graphical Models
#' @description Performs one iteration of the identifiable double conditional Bayes factor (IDCBF) algorithm for a multi-group setting in Gaussian graphical models.
#'
#' @param Gs Either a 3D array of current graphs for all groups or a single graph matrix (if only one group).
#' @param D A symmetric positive-definite matrix used in the Wishart prior.
#' @param S Sample covariance matrix for the current group.
#' @param K Initial precision matrix (optional).
#' @param p Number of nodes (variables) in the graph.
#' @param delta Degrees of freedom for the Wishart prior on the proposal graph.
#' @param deltan Degrees of freedom for the posterior update.
#' @param nu Not used in this function but may be passed for compatibility with wrappers.
#' @param theta Not used in this function but may be passed for compatibility with wrappers.
#' @param l Group index.
#' @param prior Function to compute the prior probability of a graph (not used here directly).
#' @param prior_params A list of parameters passed to the `prior_ratio` function.
#' @param prior_ratio A function that returns the prior ratio between two graphs.
#' @param identifiable Logical. If TRUE, only accepts moves that preserve identifiability of the model.
#'
#' @return A list with:
#' \describe{
#'   \item{G}{The updated graph matrix for group \code{l}.}
#'   \item{K}{The corresponding updated precision matrix.}
#'   \item{alphas}{Vector of acceptance probabilities for each edge update attempt.}
#' }
#' @import BDgraph
#' @export
idcbf_step <- function(Gs, D, S, K, p, delta, deltan, nu, theta,
                       l, prior, prior_params, prior_ratio, identifiable) {
  if(is.list(Gs)) {
    G_current <- Gs[[l]]
  } else if(is.array(Gs) && !is.matrix(Gs)) {
    G_current <- Gs[,,l]
  } else if(is.matrix(Gs)) {
    G_current <- Gs
  } else {
    stop("Gs must be either a list or an array or a single matrix")
  }
  D_plus_S <- D + S
  K_current <- BDgraph::rgwish(1, G_current, deltan, D_plus_S)

  alphas <- matrix(0, p, p)


  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      permuted <- permute(G_current, i, j)
      perm <- permuted$perm
      #should be the same deterministic permutation for all
      G_permuted <- permuted$G
      K_permuted <- permute(K_current, i, j)$G
      D_permuted <- permute(D, i, j)$G
      S_permuted <- permute(S, i, j)$G

      G_tilde <- G_permuted
      edge_removed <- FALSE

      if (G_permuted[p - 1, p] == 1) {
        G_tilde[p - 1, p] <- G_tilde[p, p - 1] <- 0
        edge_removed <- TRUE
      } else {
        G_tilde[p - 1, p] <- G_tilde[p, p - 1] <- 1
      }

      # Adding/removal of an edge has an inverse relationship
      # with the complement of the graph. Thus, we need to only check
      # identifiability when we ADD an edge in G. Adding an edge in G
      # removes an edge in the complement of G.
      if (identifiable && !edge_removed && !is_identifiable(G_tilde)) next

      pratio <- prior_ratio(G_tilde, G_permuted, prior_params)
      K_tilde <- BDgraph::rgwish(1, G_tilde, delta, D_permuted)

      N_ratio <- if (edge_removed) {
        N_KU(K_tilde, D_permuted) / N_KU(K_permuted, D_permuted + S_permuted)
      } else {
        N_KU(K_permuted, D_permuted + S_permuted) / N_KU(K_tilde, D_permuted)
      }

      ratio <- if (edge_removed) 1 / pratio else pratio
      alpha <- N_ratio * ratio

      move <- FALSE
      if (rbinom(1, 1, min(alpha, 1))) {
        move <- TRUE
        G_permuted <- G_tilde
      }

      G_current <- inv_permute(G_permuted, perm)

      if (move) {
        K_current <- BDgraph::rgwish(1, G_current, deltan, D_plus_S)
        alphas[i, j] <- alphas[i, j] + 1
      }
    }
  }

  return(list(G = G_current, K = K_current, alphas = alphas[upper.tri(alphas)]))
}

