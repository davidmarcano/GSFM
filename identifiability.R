# We note that although we are checking for bipartiteness here, non-bipartite-ness
# of the complement graph is not enough for identifiability. For example, 
# components with only one node are non-bipartite AND non-identifiable. These
# edge cases are handled below since we allow single nodes to be colored
# by the two-color algoirthm.

# ---- Convert Number to Binary Vector ----
# This function converts a given number into a binary vector representation.
# It returns the least significant bits as specified by 'noBits'.
number_to_binary <- function(number, noBits) {
    # Convert the number to its binary representation (in reverse order)
    binary_vector <- rev(as.numeric(intToBits(number)))
    
    # If 'noBits' is provided, truncate the vector to the required number of bits
    if (!missing(noBits)) {
        binary_vector <- binary_vector[(length(binary_vector) - noBits + 1):length(binary_vector)]
    }
    
    # Return the resulting binary vector
    return(binary_vector)
}

# ---- Bipartite Check Utility via BFS Coloring ----
# This utility function checks whether a graph component is bipartite.
# It uses BFS to attempt to color the graph with two colors.
is_bipartite_component <- function(graph, color_arr, src) {
    p <- ncol(graph)  # Get the number of vertices in the graph
    
    # Assign the first color to the source vertex
    color_arr[src] <- 1
    
    # Initialize a queue for BFS traversal
    queue <- list(src)
    
    # Perform BFS to check if the component can be colored with two colors
    while (length(queue) > 0) {
        u <- queue[[1]]  # Dequeue the first vertex
        queue <- queue[-1]  # Remove the vertex from the queue
        
        # If there's a self-loop, it's not bipartite
        if (graph[u, u] == 1) {
            return(list(color_arr = color_arr, bipartite = FALSE))
        }
        
        # Check all adjacent vertices of u
        for (v in 1:p) {
            # If there's an edge and v is uncolored, assign the alternate color to v
            if (graph[u, v] == 1 && color_arr[v] == -1) {
                color_arr[v] <- 1 - color_arr[u]
                queue <- c(queue, v)  # Enqueue v for further processing
            } 
            # If there's an edge and v has the same color as u, it's not bipartite
            else if (graph[u, v] == 1 && color_arr[v] == color_arr[u]) {
                return(list(color_arr = color_arr, bipartite = FALSE))
            }
        }
    }

    # If we get here, the component is bipartite
    return(list(color_arr = color_arr, bipartite = TRUE))
}

# ---- Connected Component Extraction ----
# This function performs a BFS traversal to find all the connected components in the graph.
# It returns a list of connected components where each component is a matrix of edges.
add_edges_to_component <- function(graph, cc, connected, src) {
    queue <- list(src)
    p <- ncol(graph)  # Number of vertices in the graph
    
    # Perform BFS to add all edges in the connected component
    while (length(queue) > 0) {
        u <- queue[[1]]  # Dequeue the first vertex
        queue <- queue[-1]  # Remove the vertex from the queue
        
        for (v in 1:p) {
            # If there's an edge between u and v, add it to the connected component
            if (graph[u, v] == 1) {
                if (connected[v] == -1) {
                    queue <- c(queue, v)  # Enqueue v for further processing
                    connected[v] <- 1  # Mark v as connected
                }
                cc[u, v] <- 1  # Add the edge to the component
                cc[v, u] <- 1  # Ensure the component is undirected
            }
        }
    }
    
    # Return the updated connected component matrix and connected vertices array
    return(list(cc = cc, connected = connected))
}

# Function to find all connected components in the graph
find_connected_components <- function(graph) {
    p <- ncol(graph)  # Number of vertices in the graph
    connected <- rep(-1, p)  # Initialize an array to track connected vertices
    cc_list <- list()  # Initialize a list to store connected components
    
    # Traverse through each vertex
    for (i in 1:p) {
        if (connected[i] == -1) {
            # If the vertex is not yet visited, find its connected component
            result <- add_edges_to_component(graph, matrix(0, p, p), connected, i)
            # Store the connected component and the source node
            cc_list <- append(cc_list, list(list(result$cc, i)))
            connected <- result$connected  # Update the connected array
        }
    }
    
    # Return the list of connected components
    return(cc_list)
}

# ---- Identifiability Check ----
# This function checks if the graph is identifiable by examining its complement graph.
# It looks for bipartite components in the complement graph (which indicates non-identifiability).
is_identifiable <- function(graph, debugFlag = FALSE) {
    p <- ncol(graph)  # Get the number of vertices in the graph
    
    # Step 1: Compute the complement graph (reverse the edges)
    complement_g <- 1 - graph
    diag(complement_g) <- 0  # Remove self-loops in the complement graph
    
    # Step 2: Debugging - Plot the original and complement graphs
    if (debugFlag) {
        library(igraph)
        par(mfrow = c(1, 2))  # Set up side-by-side plotting
        plot(graph_from_adjacency_matrix(graph, mode = "undirected"), main = "Original Graph")
        plot(graph_from_adjacency_matrix(complement_g, mode = "undirected"), main = "Complement Graph")
        par(mfrow = c(1, 1))  # Reset plotting layout
    }
    
    # Step 3: Find all connected components in the complement graph
    components <- find_connected_components(complement_g)

    # Step 4: Check each connected component for bipartiteness
    for (comp in components) {
        cc_matrix <- comp[[1]]  # The matrix representing the connected component
        src <- comp[[2]]  # The source node used for BFS
        color_arr <- rep(-1, p)  # Initialize the color array
        result <- is_bipartite_component(cc_matrix, color_arr, src)
        
        # If any component is bipartite, return FALSE (the graph is not identifiable)
        if (result$bipartite) {
            return(FALSE)
        }
    }
    
    # Step 5: If no bipartite component is found, the graph is identifiable
    return(TRUE)
}
