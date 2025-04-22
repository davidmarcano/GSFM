library(testthat)

# Define helper to wrap expected output test
test_identifiability <- function(graph_matrix, expected, label, debug = FALSE) {
  test_that(label, {
    result <- is_identifiable(graph_matrix, debug)
    expect_equal(result, expected)
  })
}

# -----------------------------
# Test 1: Fully Connected Graph
# Complement is empty (no edges) => all components trivially bipartite
# Expect: FALSE (not identifiable)
# -----------------------------
p1 <- 4
graph1 <- matrix(1, p1, p1)
diag(graph1) <- 0
test_identifiability(graph1, FALSE, "Fully Connected Graph")

# -----------------------------
# Test 2: Empty Graph
# Complement is a complete graph => one big non-bipartite component
# Expect: TRUE (identifiable)
# -----------------------------
p2 <- 5
graph2 <- matrix(0, p2, p2)
test_identifiability(graph2, TRUE, "Empty Graph")

# -----------------------------
# Test 3: Bipartite Graph (4-cycle)
# Complement has 2 edges between non-adjacent nodes
# Expect: FALSE (not identifiable)
# -----------------------------
graph3 <- matrix(c(
  0,1,0,1,
  1,0,1,0,
  0,1,0,1,
  1,0,1,0
), 4, 4, byrow = TRUE)
test_identifiability(graph3, FALSE, "4-Cycle Bipartite Graph")

# -----------------------------
# Test 4: Triangle Graph (3-cycle)
# Complement has disconnected vertices + an edge between one pair
# Expect: FALSE (not identifiable)
# -----------------------------
graph4 <- matrix(c(
  0,1,1,
  1,0,1,
  1,1,0
), 3, 3, byrow = TRUE)
test_identifiability(graph4, FALSE, "Triangle Graph")

# -----------------------------
# Test 5: Star Graph
# One central node connected to all others
# Expect: FALSE (not identifiable)
# -----------------------------
p5 <- 5
graph5 <- matrix(0, p5, p5)
graph5[1, 2:5] <- 1
graph5[2:5, 1] <- 1
test_identifiability(graph5, FALSE, "Star Graph")

# -----------------------------
# Test 6: Two Disjoint Edges
# Complement contains a 4-cycle (bipartite)
# Expect: FALSE (not identifiable)
# -----------------------------
graph6 <- matrix(0, 4, 4)
graph6[1,2] <- 1; graph6[2,1] <- 1
graph6[3,4] <- 1; graph6[4,3] <- 1
test_identifiability(graph6, FALSE, "Two Disjoint Edges")
