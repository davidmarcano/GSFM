library(igraph)
library(graphics)

#' Convert Number to Binary Vector
#'
#' This function converts a given number into its binary vector representation.
#' It returns the least significant bits as specified by 'noBits'.
#'
#' @param number The number to be converted to binary.
#' @param noBits The number of bits to return. If missing, the full binary representation is returned.
#' @return A numeric vector representing the binary representation of the given number.
#' @examples
#' GSFM:::number_to_binary(5, 3)  # returns c(1, 0, 1)
number_to_binary <- function(number, noBits) {
  binary_vector <- rev(as.numeric(intToBits(number)))

  if (!missing(noBits)) {
    binary_vector <- binary_vector[(length(binary_vector) - noBits + 1):length(binary_vector)]
  }

  return(binary_vector)
}

#' Bipartite Check Utility via BFS Coloring
#'
#' This utility function checks whether a graph component is bipartite.
#' It uses BFS to attempt to color the graph with two colors.
#'
#' @param graph An adjacency matrix representing the graph.
#' @param color_arr A numeric vector used to store the color of each vertex (-1 means uncolored).
#' @param src The source vertex from where BFS starts.
#' @return A list with two elements:
#'   - color_arr: A numeric vector with the color assignments.
#'   - bipartite: A logical indicating whether the component is bipartite.
#' @examples
#' graph <- matrix(c(
#'   0, 1, 0,
#'   1, 0, 1,
#'   0, 1, 0
#' ), 3, 3, byrow = TRUE)
#' GSFM:::is_bipartite_component(graph, rep(-1, 3), 1)
is_bipartite_component <- function(graph, color_arr, src) {
  p <- ncol(graph)
  color_arr[src] <- 1
  queue <- list(src)

  while (length(queue) > 0) {
    u <- queue[[1]]
    queue <- queue[-1]

    if (graph[u, u] == 1) {
      return(list(color_arr = color_arr, bipartite = FALSE))
    }

    for (v in 1:p) {
      if (graph[u, v] == 1 && color_arr[v] == -1) {
        color_arr[v] <- 1 - color_arr[u]
        queue <- c(queue, v)
      } else if (graph[u, v] == 1 && color_arr[v] == color_arr[u]) {
        return(list(color_arr = color_arr, bipartite = FALSE))
      }
    }
  }

  return(list(color_arr = color_arr, bipartite = TRUE))
}

#' Add Edges to a Connected Component in a Graph
#'
#' Traverses a graph starting from a source node and marks all reachable nodes as part of
#' the same connected component, updating the component matrix and visited node tracking.
#'
#' @param graph A binary adjacency matrix (p x p) representing the graph, where `graph[i, j] = 1`
#'   indicates an edge from node `i` to node `j`.
#' @param cc A binary matrix (p x p) representing the connected component being constructed.
#'   Initially contains 0s and is updated to 1s where connections are found.
#' @param connected An integer vector of length `p` where `-1` indicates unvisited nodes and
#'   `1` indicates visited nodes.
#' @param src An integer index indicating the source node to begin traversal from.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{cc}{The updated component matrix with added edges.}
#'     \item{connected}{The updated vector of visited nodes.}
#'   }
#' @examples
#' graph <- matrix(c(
#'   0, 1, 0,
#'   1, 0, 1,
#'   0, 1, 0
#' ), 3, 3, byrow = TRUE)
#' cc <- matrix(0, 3, 3)
#' connected <- rep(-1, 3)
#' GSFM:::add_edges_to_component(graph, cc, connected, 1)
add_edges_to_component <- function(graph, cc, connected, src) {
  queue <- list()
  queue <- c(queue, src)
  p <- ncol(graph)
  while(length(queue) != 0) {
    u <- queue[[1]]
    queue <- queue[-1]
    for(j in 1:p) {
      if(graph[u,j] == 1) {
        if (connected[j] == -1) {
          queue <- c(queue, j)
          connected[j] <- 1
        }
        cc[u,j] <- 1
        cc[j,u] <- 1
      }
    }
  }
  return(list(cc = cc, connected = connected))
}

#' Connected Component Extraction
#'
#' This function performs a BFS traversal to find all the connected components in the graph.
#' It returns a list of connected components where each component is a matrix of edges.
#'
#' @param graph An adjacency matrix representing the graph.
#' @return A list of connected components, each represented by an adjacency matrix.
#' @examples
#' graph <- matrix(c(
#'   0, 1, 0, 0,
#'   1, 0, 0, 0,
#'   0, 0, 0, 1,
#'   0, 0, 1, 0
#' ), 4, 4, byrow = TRUE)
#' GSFM:::find_connected_components(graph)
find_connected_components <- function(graph) {
  p <- ncol(graph)
  connected <- rep(-1, p)
  cc_list <- list()

  for (i in 1:p) {
    if (connected[i] == -1) {
      result <- add_edges_to_component(graph, matrix(0, p, p), connected, i)
      cc_list <- append(cc_list, list(list(result$cc, i)))
      connected <- result$connected
    }
  }

  return(cc_list)
}

#' Identifiability Check for Graphs
#'
#' This function checks if the graph is identifiable by examining its complement graph.
#' It looks for bipartite components in the complement graph (which indicates non-identifiability).
#'
#' @param graph An adjacency matrix representing the graph.
#' @param debugFlag A logical flag to plot the original and complement graphs for debugging.
#' @return A logical indicating whether the graph is identifiable (TRUE) or not (FALSE).
#' @export
#' @examples
#' graph <- matrix(c(
#'   0, 1, 1,
#'   1, 0, 1,
#'   1, 1, 0
#' ), 3, 3, byrow = TRUE)
#' is_identifiable(graph)
is_identifiable <- function(graph, debugFlag = FALSE) {
  p <- ncol(graph)
  complement_g <- 1 - graph
  diag(complement_g) <- 0

  if (debugFlag) {
    graphics::par(mfrow = c(1, 2))
    plot(igraph::graph_from_adjacency_matrix(graph, mode = "undirected"), main = "Original Graph")
    plot(igraph::graph_from_adjacency_matrix(complement_g, mode = "undirected"), main = "Complement Graph")
    graphics::par(mfrow = c(1, 1))
  }

  components <- find_connected_components(complement_g)

  for (comp in components) {
    cc_matrix <- comp[[1]]
    src <- comp[[2]]
    color_arr <- rep(-1, p)
    result <- is_bipartite_component(cc_matrix, color_arr, src)

    if (result$bipartite) {
      return(FALSE)
    }
  }

  return(TRUE)
}
