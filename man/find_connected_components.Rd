% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/identifiability.R
\name{find_connected_components}
\alias{find_connected_components}
\title{Connected Component Extraction}
\usage{
find_connected_components(graph)
}
\arguments{
\item{graph}{An adjacency matrix representing the graph.}
}
\value{
A list of connected components, each represented by an adjacency matrix.
}
\description{
This function performs a BFS traversal to find all the connected components in the graph.
It returns a list of connected components where each component is a matrix of edges.
}
\examples{
graph <- matrix(c(
  0, 1, 0, 0,
  1, 0, 0, 0,
  0, 0, 0, 1,
  0, 0, 1, 0
), 4, 4, byrow = TRUE)
GSFM:::find_connected_components(graph)
}
