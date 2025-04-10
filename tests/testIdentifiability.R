library(igraph)
source("../identifiability.R")

# Utility function to nicely display test results
run_identifiability_test <- function(graph_matrix, expected, test_name, debugFlagFlag = FALSE) {
  result <- is_identifiable(graph_matrix, debugFlag)
  status <- if (result == expected) "PASS" else "FAIL"
  cat(sprintf("[%s] %s: expected=%s, got=%s\n", status, test_name, expected, result))
}
debugFlag <- TRUE
# -----------------------------
# Test 1: Fully Connected Graph
# Complement is empty (no edges) => all components trivially bipartite
# Expect: FALSE (not identifiable)
# -----------------------------
p1 <- 4
graph1 <- matrix(1, p1, p1)
diag(graph1) <- 0
run_identifiability_test(graph1, FALSE, "Fully Connected Graph", debugFlag)

# -----------------------------
# Test 2: Empty Graph
# Complement is a complete graph => one big non-bipartite component
# Expect: TRUE (identifiable)
# -----------------------------
p2 <- 5
graph2 <- matrix(0, p2, p2)
run_identifiability_test(graph2, TRUE, "Empty Graph", debugFlag)

# -----------------------------
# Test 3: Bipartite Graph (4-cycle)
# Complement has 2 edges between non-adjacent nodes
# Those edges form disconnected components that are trivially bipartite
# Expect: FALSE (not identifiable)
# -----------------------------
graph3 <- matrix(c(
  0,1,0,1,
  1,0,1,0,
  0,1,0,1,
  1,0,1,0
), 4, 4, byrow=TRUE)
run_identifiability_test(graph3, FALSE, "4-Cycle Bipartite Graph", debugFlag)

# -----------------------------
# Test 4: Triangle Graph (3-cycle)
# Complement has disconnected vertices + an edge between one pair
# Expect: FALSE (non-identifiable) since single nodes are not identifiable
# -----------------------------
graph4 <- matrix(c(
  0,1,1,
  1,0,1,
  1,1,0
), 3, 3, byrow=TRUE)
run_identifiability_test(graph4, FALSE, "Triangle Graph", debugFlag)

# -----------------------------
# Test 5: Star Graph (center connected to all others)
# Complement includes a single node.
# Expect: FALSE (non-identifiable)
# -----------------------------
p5 <- 5
graph5 <- matrix(0, p5, p5)
graph5[1, 2:5] <- 1
graph5[2:5, 1] <- 1
run_identifiability_test(graph5, FALSE, "Star Graph", debugFlag)

# -----------------------------
# Test 6: Two disjoint edges
# Complement contains a 4-cycle (bipartite)
# Expect: FALSE (not identifiable)
# -----------------------------
graph6 <- matrix(0, 4, 4)
graph6[1,2] <- 1; graph6[2,1] <- 1
graph6[3,4] <- 1; graph6[4,3] <- 1
run_identifiability_test(graph6, FALSE, "Two Disjoint Edges", debugFlag)
