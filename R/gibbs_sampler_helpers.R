

#' Calculate Shifted Observations
#'
#' Centers the observed data by subtracting the intercept (alpha) and factor components (lambda*f).
#' This creates "shifted" observations used in conditional calculations.
#'
#' @param x Numeric matrix of observed data (n x p)
#' @param alpha Numeric vector of intercept parameters (length p)
#' @param lambda Numeric vector of factor loadings (length p)
#' @param f Numeric vector of factor scores (length n)
#' @return Numeric matrix of shifted observations (n x p)
#'
#' @examples
#' # Sample data
#' x <- matrix(rnorm(100), 20, 5)
#' alpha <- rnorm(5)
#' lambda <- abs(rnorm(5))
#' f <- rnorm(20)
#'
#' # Calculate shifted values
#' x_shifted <- calculate_shifted_x(x, alpha, lambda, f)
#'
#' @export
calculate_shifted_x <- function(x, alpha, lambda, f) {
  # Input validation
  if (ncol(x) != length(alpha) || ncol(x) != length(lambda)) {
    stop("Dimensions of x, alpha, and lambda do not match")
  }
  if (nrow(x) != length(f)) {
    stop("Number of observations in x doesn't match length of f")
  }

  # Vectorized calculation
  x - outer(f, lambda) - matrix(alpha, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
}

#' Initialize MCMC chains for lattice model
#'
#' Creates and initializes all the necessary chains for the Gibbs sampler.
#'
#' @param steps Number of MCMC steps
#' @param L Number of groups
#' @param p Number of variables
#' @param n Vector of sample sizes for each group
#' @param has_start Logical indicating if starting values are provided
#' @param start_lambda Starting values for lambda parameters
#' @param start_alpha Starting values for alpha parameters
#' @param start_mus Starting values for mu parameters
#' @param start_kw Starting values for Kw matrix
#' @param start_f Starting values for factor scores
#' @param start_delta Starting values for delta parameters
#' @param start_gl Starting values for group graphs
#' @param start_kl Starting values for group precision matrices
#' @param start_nu Starting values for nu parameters
#' @param start_theta Starting values for theta matrix
#' @param start_gamma Starting values for gamma matrix
#' @param start_corr_gamma Starting values for correlation matrices
#' @param non_continuous Logical indicating non-continuous data
#' @param z Starting values for latent variables (if non-continuous)
#' @return List containing initialized chains for all parameters
initialize_chains <- function(steps, L, p, n, has_start, start_lambda, start_alpha,
                              start_mus, start_kw, start_f, start_delta, start_gl,
                              start_kl, start_nu, start_theta, start_gamma,
                              start_corr_gamma, non_continuous, z, params) {
  g_chain = initialize_g_chain(steps, p, L, start_gl, has_start) # need to initialize first
  list(
    alpha_chain = initialize_chain(steps, p, L, start_vals = start_alpha, has_start = has_start),
    lambda_chain = initialize_chain(steps, p, L, start_vals = start_lambda, has_start = has_start),
    f_temp = initialize_f_temp(L, n, start_f, has_start),
    g_chain = g_chain,
    k_chain = initialize_k_chain(steps, p, L, start_kl, has_start, g_chain[1,,,]),
    acceptprob_chain = array(0, c(steps, p*(p-1)/2, L)),
    delta_chain = initialize_delta_chain(steps, L, start_delta, has_start, params = params),
    mu_chain = initialize_mu_chain(steps, L, start_mus, has_start),
    kw_chain = initialize_kw_chain(steps, L, start_kw, has_start, params = params),
    #nu_chain = initialize_nu_chain(steps, p, start_nu, has_start, params = params),
    #theta_chain = initialize_theta_chain(steps, L, start_theta, has_start, params = params),
    #gamma_chain = initialize_gamma_chain(steps, L, start_gamma, has_start, params = params),
    corr_gamma_chain = initialize_corr_gamma_chain(steps, p, L, start_corr_gamma, has_start),
    z_chain = if(non_continuous) list(z) else NULL,
    xpred = array(0, c(steps, p, L))
  )
}

#' Initialize a generic MCMC chain
#'
#' Helper function to initialize a single MCMC chain.
#'
#' @param steps Number of MCMC steps
#' @param dim1 First dimension size
#' @param dim2 Second dimension size (optional)
#' @param dim3 Third dimension size (optional)
#' @param start_vals Starting values
#' @param has_start Logical indicating if starting values are provided
#' @return Initialized chain array/matrix
#' @keywords internal
initialize_chain <- function(steps, dim1, dim2 = NULL, dim3 = NULL, start_vals, has_start) {
  if(is.null(dim2) && is.null(dim3)) {
    if(has_start) matrix(start_vals, steps, dim1) else matrix(0, steps, dim1)
  } else if(is.null(dim3)) {
    if(has_start) array(start_vals, c(steps, dim1, dim2)) else array(0, c(steps, dim1, dim2))
  } else {
    if(has_start) array(start_vals, c(steps, dim1, dim2, dim3)) else array(0, c(steps, dim1, dim2, dim3))
  }
}

#' Set hyperparameters for lattice model
#'
#' Defines the default hyperparameters used in the Gibbs sampler.
#'
#' @param p Number of variables
#' @param L Number of groups
#' @param W Spatial adjacency matrix
#' @return List containing all hyperparameters
set_hyperparameters <- function(p, L, W) {
  list(
    rho = 0.07,
    Dl = diag(1, p, p),
    deltal = 3,
    deltaw = 3,
    n0 = 1,
    c = 2,
    d = 1,
    w = 0.9,
    alpha = 2,
    beta = 5,
    a = 1,
    b = 4,
    I_L = diag(1, L, L),
    W = W
  )
}

#' Initialize copula structures for non-continuous data
#'
#' Creates data structures needed for copula modeling of non-continuous variables.
#'
#' @param x List of data matrices (one per group)
#' @param L Number of groups
#' @return List containing:
#' \itemize{
#'   \item lessthan - Array indicating values less than breakpoints
#'   \item greaterthan - Array indicating values greater than breakpoints
#'   \item obs_buckets - Matrix indicating bucket membership for each observation
#' }
initialize_copula_structures <- function(x, L) {
  copula_lessthan <- vector("list", L)
  copula_greaterthan <- vector("list", L)
  observation_buckets <- vector("list", L)

  for(l in 1:L) {
    subs <- create_copula_buckets(x[[l]])
    copula_lessthan[[l]] <- subs$lessthan
    copula_greaterthan[[l]] <- subs$greaterthan
    observation_buckets[[l]] <- subs$obs_buckets
  }

  list(
    lessthan = copula_lessthan,
    greaterthan = copula_greaterthan,
    obs_buckets = observation_buckets
  )
}

#' Set up parallel processing
#'
#' Initializes a parallel cluster for distributed computation.
#'
#' @return A parallel cluster object
setup_parallel_processing <- function() {
  numCores <- 6
  cl <- makeCluster(numCores, outfile = "TEST_2_5")
  clusterEvalQ(cl, .libPaths())
  return(cl)
}

#' Update group parameters in parallel
#'
#' Distributes the group parameter updates across multiple cores.
#'
#' @param cl Parallel cluster object
#' @param L Number of groups
#' @param chains List containing current MCMC chains
#' @param x List of data matrices
#' @param n Vector of sample sizes
#' @param p Number of variables
#' @param params List of hyperparameters
#' @param non_continuous Logical indicating non-continuous data
#' @param copula Logical indicating copula model
#' @param copula_structs List of copula structures
#' @param s Current MCMC iteration
#' @return List of updated parameters for each group
update_groups_parallel <- function(cl, L, chains, x, n, p, params,
                                   non_continuous, copula, copula_structs, s) {
  if(non_continuous) {
    # Prepare arguments for parallel processing
    args <- list(
      alpha_chain = chains$alpha_chain[s-1,,],
      lambda_chain = chains$lambda_chain[s-1,,],
      gl = chains$g_chain[s-1,,,],
      kl = chains$k_chain[s-1,,,],
      nu = chains$nu_chain[s-1,,],
      theta = chains$theta_chain[s-1,,],
      mul = chains$mu_chain[s-1,],
      delta_group = chains$delta_chain[s-1,],
      fl = chains$f_temp,
      x = x,
      n = n,
      p = p,
      Dl = params$Dl,
      deltal = params$deltal,
      identifiable = TRUE,
      non_continuous = non_continuous,
      prior = prior,
      prior_params = prior_params,
      prior_ratio = prior_ratio,
      alpha = params$alpha,
      beta = params$beta,
      w = params$w,
      a = params$a,
      b = params$b,
      c = params$c,
      d = params$d,
      n0 = params$n0,
      s = s,
      copy_x = copy_x,
      cumulativeMargProb = cumulativeMargProb,
      copula = copula,
      copula_lessthan_buckets = copula_structs$lessthan,
      copula_greaterthan_buckets = copula_structs$greaterthan,
      obs_buckets = copula_structs$obs_buckets
    )

    parLapply(cl, 1:L, update_group_parameters, args)
  } else {
    # Simplified arguments for continuous case
    args <- list(
      alpha_chain = chains$alpha_chain[s-1,,],
      lambda_chain = chains$lambda_chain[s-1,,],
      gl = chains$g_chain[s-1,,,],
      kl = chains$k_chain[s-1,,,],
      nu = chains$nu_chain[s-1,,],
      theta = chains$theta_chain[s-1,,],
      mul = chains$mu_chain[s-1,],
      delta_group = chains$delta_chain[s-1,],
      fl = chains$f_temp,
      x = x,
      n = n,
      p = p,
      Dl = params$Dl,
      deltal = params$deltal,
      identifiable = TRUE,
      non_continuous = FALSE,
      prior = prior,
      prior_params = prior_params,
      prior_ratio = prior_ratio,
      alpha = params$alpha,
      beta = params$beta,
      w = params$w,
      a = params$a,
      b = params$b,
      c = params$c,
      d = params$d,
      n0 = params$n0,
      s = s
    )

    parLapply(cl, 1:L, update_group_parameters, args)
  }
}

#' Update MCMC chains with new results
#'
#' Incorporates the results from parallel group updates into the main chains.
#'
#' @param chains List containing current MCMC chains
#' @param results List of results from parallel updates
#' @param L Number of groups
#' @param s Current MCMC iteration
#' @return Updated chains list
update_chains <- function(chains, results, L, s) {
  for(l in 1:L) {
    chains$alpha_chain[s,,l] <- results[[l]]$alpha_chain
    chains$lambda_chain[s,,l] <- results[[l]]$lambda_chain
    chains$f_temp[[l]] <- results[[l]]$f_temp
    chains$g_chain[s,,,l] <- results[[l]]$graph_chain
    chains$k_chain[s,,,l] <- results[[l]]$K_chain
    chains$corr_gamma_chain[s,,,l] <- results[[l]]$corr_gamma_chain
    chains$delta_chain[s,l] <- results[[l]]$tau_chain

    if(!is.null(results[[l]]$latent_zs)) {
      chains$z_chain[[l]] <- results[[l]]$latent_zs
    }

    chains$xpred[s,,l] <- if(is.null(results[[l]]$xpred)) {
      rep(NA, p)
    } else {
      results[[l]]$xpred
    }
  }
  return(chains)
}

#' Update shared parameters
#'
#' Updates parameters that are shared across groups (theta, gamma, nu).
#'
#' @param chains List containing current MCMC chains
#' @param params List of hyperparameters
#' @param s Current MCMC iteration
#' @return Updated chains list
# update_shared_parameters <- function(chains, params, s) {
#   mh_updates <- mh_theta_gamma(
#     chains$g_chain[s,,,],
#     chains$theta_chain[s-1,,],
#     chains$gamma_chain[s-1,,],
#     chains$nu_chain[s-1,,],
#     params$alpha,
#     params$beta,
#     params$w
#   )
#
#   chains$theta_chain[s,,] <- mh_updates$theta
#   chains$gamma_chain[s,,] <- mh_updates$gamma
#   chains$nu_chain[s,,] <- mh_nu(
#     chains$g_chain[s,,,],
#     chains$nu_chain[s-1,,],
#     chains$theta_chain[s,,],
#     params$a,
#     params$b
#   )
#
#   return(chains)
# }

#' Save progress during MCMC
#'
#' Periodically saves the state of the MCMC chains to disk.
#'
#' @param chains List containing current MCMC chains
#' @param s Current MCMC iteration
save_progress <- function(chains, s) {
  filename <- ifelse(s %% 200 < 100, "13yearsFirst", "13yearsLast")
  save(list = setdiff(ls(all.names = TRUE), "results"),
       file = paste0(filename, s %% 200, ".RData"))
}


#' Initialize factor scores chain
#'
#' Initializes the chain for factor scores.
#'
#' @param L Number of groups
#' @param n Vector of sample sizes
#' @param start_f Starting values for factors
#' @param has_start Logical indicating if starting values are provided
#' @return Initialized factor scores list
#' @keywords internal
initialize_f_temp <- function(L, n, start_f, has_start) {
  if(has_start) {
    start_f
  } else {
    lapply(n, function(nl) rnorm(nl))
  }
}

#' Initialize precision matrices chain
#'
#' Initializes the chain for group precision matrices.
#'
#' @param steps Number of MCMC steps
#' @param p Number of variables
#' @param L Number of groups
#' @param start_kl Starting values for precision matrices
#' @param has_start Logical indicating if starting values are provided
#' @return Initialized precision matrices array
#' @keywords internal
initialize_k_chain <- function(steps, p, L, start_kl, has_start, g1) {
  if(has_start) {
    k <- array(start_kl, c(steps, p, p, L))
  } else {
    k <- array(0, c(steps, p, p, L))
    for(l in 1:L) {
      k[1,,,l] <- BDgraph::rgwish(1, g1[,,l], D = diag(1, p, p))
    }
  }
  return(k)
}

#' Initialize graph structures chain
#'
#' Initializes the chain for group graph structures.
#'
#' @param steps Number of MCMC steps
#' @param p Number of variables
#' @param L Number of groups
#' @param start_gl Starting values for graphs
#' @param has_start Logical indicating if starting values are provided
#' @return Initialized graph structures array
#' @keywords internal
initialize_g_chain <- function(steps, p, L, start_gl, has_start, identifiable = TRUE) {
  if(has_start) {
    g <- array(start_gl, c(steps, p, p, L))
  } else {
    graph <- matrix(0, p, p)
    gl <- lapply(1:L, function(l) {
      repeat {
        G <- prior(graph, prior_params)$graph
        if(!identifiable || is_identifiable(G)) return(G)
      }
    })
    g <- array(0, c(steps, p, p, L))
    g[1,,,] <- simplify2array(gl)
  }
  return(g)
}

#' Initialize delta parameters chain
#'
#' Initializes the chain for delta parameters.
#'
#' @param steps Number of MCMC steps
#' @param L Number of groups
#' @param start_delta Starting values for delta
#' @param has_start Logical indicating if starting values are provided
#' @return Initialized delta parameters matrix
#' @keywords internal
initialize_delta_chain <- function(steps, L, start_delta, has_start, params = NULL) {
  if(has_start) {
    matrix(start_delta, steps, L)
  } else {
    matrix(rinvgamma(steps*L, shape = params$c, scale = (params$c*params$d)/2), steps, L)
  }
}

#' Initialize mu parameters chain
#'
#' Initializes the chain for mu parameters.
#'
#' @param steps Number of MCMC steps
#' @param L Number of groups
#' @param start_mus Starting values for mu
#' @param has_start Logical indicating if starting values are provided
#' @return Initialized mu parameters matrix
#' @keywords internal
initialize_mu_chain <- function(steps, L, start_mus, has_start) {
  if(has_start) {
    matrix(start_mus, steps, L)
  } else {
    matrix(0, steps, L)
  }
}

#' Initialize Kw matrix chain
#'
#' Initializes the chain for the spatial precision matrix Kw.
#'
#' @param steps Number of MCMC steps
#' @param L Number of groups
#' @param start_kw Starting values for Kw
#' @param has_start Logical indicating if starting values are provided
#' @return Initialized Kw matrix array
#' @keywords internal
initialize_kw_chain <- function(steps, L, start_kw, has_start, params = NULL) {
  if(has_start) {
    array(start_kw, c(steps, L, L))
  } else {
    array(rgwish(1, params$W, params$deltaw,
                 (params$deltaw - 2)*solve(params$I_L - params$rho*params$W)), c(steps, L, L))
  }
}

#' Initialize nu parameters chain
#'
#' Initializes the chain for nu parameters.
#'
#' @param steps Number of MCMC steps
#' @param p Number of variables
#' @param start_nu Starting values for nu
#' @param has_start Logical indicating if starting values are provided
#' @return Initialized nu parameters array
#' @keywords internal
# initialize_nu_chain <- function(steps, p, start_nu, has_start, params = NULL) {
#   if(has_start) {
#     array(start_nu, c(steps, p, p))
#   } else {
#     array(init_nu(p, params$a, params$b), c(steps, p, p))
#   }
# }

#' Initialize theta matrix chain
#'
#' Initializes the chain for theta matrix.
#'
#' @param steps Number of MCMC steps
#' @param L Number of groups
#' @param start_theta Starting values for theta
#' @param has_start Logical indicating if starting values are provided
#' @return Initialized theta matrix array
#' @keywords internal
# initialize_theta_chain <- function(steps, L, start_theta, has_start, params = NULL) {
#   if(has_start) {
#     array(start_theta, c(steps, L, L))
#   } else {
#     theta_gamma <- make_theta_gamma(L, params$alpha, params$beta, params$w)
#     array(theta_gamma$theta, c(steps, L, L))
#   }
# }

#' Initialize gamma matrix chain
#'
#' Initializes the chain for gamma matrix.
#'
#' @param steps Number of MCMC steps
#' @param L Number of groups
#' @param start_gamma Starting values for gamma
#' @param has_start Logical indicating if starting values are provided
#' @return Initialized gamma matrix array
#' @keywords internal
# initialize_gamma_chain <- function(steps, L, start_gamma, has_start, params = NULL) {
#   if(has_start) {
#     array(start_gamma, c(steps, L, L))
#   } else {
#     theta_gamma <- make_theta_gamma(L, params$alpha, params$beta, params$w)
#     array(theta_gamma$gamma, c(steps, L, L))
#   }
# }

#' Initialize correlation matrices chain
#'
#' Initializes the chain for group correlation matrices.
#'
#' @param steps Number of MCMC steps
#' @param p Number of variables
#' @param L Number of groups
#' @param start_corr_gamma Starting values for correlation matrices
#' @param has_start Logical indicating if starting values are provided
#' @return Initialized correlation matrices array
#' @keywords internal
initialize_corr_gamma_chain <- function(steps, p, L, start_corr_gamma, has_start) {
  if(has_start) {
    array(start_corr_gamma, c(steps, p, p, L))
  } else {
    array(0, c(steps, p, p, L))
  }
}
