source("../priors.R")
source("../identifiability.R")
#source("mggm_dcbf.R")
#source("../ggm_dcbf.R")
library(MASS)
library(invgamma)
library(parallel)
library(tmvtnorm)
library(truncnorm)
library(statmod)
#library(TruncatedNormal)

update_factors <- function(K, x, lambda, alpha) {
  lamlam <- t(lambda) %*% lambda
  #print(str(K))
  lamk <- (t(lambda) %*% solve(K) %*% lambda)/lamlam
  sigma <- sqrt(lamk/(lamlam + lamk))
  
  updated_f <- vector("numeric", nrow(x))
  for(j in 1:nrow(x)) {
    mu <- (t(lambda) %*% (x[j,] - alpha))/(lamlam + lamk)
    updated_f[j] <- rnorm(1, mean = mu, sd = sigma)
  }
  updated_f <- updated_f - mean(updated_f)
  #(updated_f - mean(updated_f))#/sd(updated_f)
  #updated_f <- scale(updated_f, TRUE, TRUE)
  return(updated_f)
  
}

update_alphas <- function(K, x, n, n0, f, lambda) {
  p <- ncol(x)
  covariance <- solve(n*K + n0*diag(rep(1, p)))
  x_hat <- matrix(0, n, p)
  
  for(j in 1:n) {
    x_hat[j,] <- x[j,] - lambda*f[j]
  }
  new_alphas <- mvrnorm(1, mu = covariance %*% K %*% colSums(x_hat), Sigma = covariance)
  #new_alphas <- (new_alphas - mean(new_alphas))/sd(new_alphas)
  return(new_alphas)
}

update_lambdas <- function(K, x, f, delta, alpha, lambda) {
  n <- nrow(x)
  p <- ncol(x)
  x_hat <- matrix(0, n, ncol(x))
  for(j in 1:n) {
      x_hat[j,] <- x[j,] - alpha
  }
  
  deltafsum <- 1/delta + sum(f^2)
  mu <- 1/deltafsum * colSums(f*x_hat)
  #mu[1] <- abs(mu[1])
  covariance <- solve(deltafsum*K)
  # Need to force numerical symmetry...because of rounding errors?
  covariance[lower.tri(covariance)] <- t(covariance)[lower.tri(covariance)]
  Hprecision <- deltafsum*K
  Hprecision[lower.tri(Hprecision)] <- t(Hprecision)[lower.tri(Hprecision)]
  samp <- mvrnorm(1, mu = mu, Sigma = covariance)
  
  # print("truncated sampling")
  # # print(sum(f^2))
  # # print(1/deltafsum)
  # # print(colSums(f*x_hat))
  # print("Mu")
  # print(mu)
  # print(covariance)
  # print(lambda)
  # print(K)
  # # print(sum(f^2))
  # # print(delta)
  # # print(deltafsum)
  # # print(deltafsum*K)
  # # print(mu)
  # # print(deltafsum)
  # print(Hprecision)
  
  # samp <- tmvtnsim::rtmvnorm(mu, covariance, lower = c(0, rep(-Inf, p-1)), 
  #                             upper=rep(Inf, p))
  # samp <- tmvtnorm::rtmvnorm(1, mean = mu, sigma = covariance, lower = c(0, rep(-Inf, p-1)),
  #                                     algorithm = "gibbs", start.value = lambda)
  #samp <- rtmvnorm(1, mean = mu, H = Hprecision, lower = c(0, rep(-Inf, p-1)),
  #                  algorithm = "gibbs")
  #print(samp)
  # #samp <- rtmvnorm(n, mu, sigma, lb, ub)
  # 
  # resamples <- 100
  # #print(samp)
  # while((samp[1] < 0) | (is.nan(samp[1])) | is.infinite(samp[1])) {
  #     #samp <- mvrnorm(1, mu = mu, Sigma = covariance)
  #     print("truncated sampling failed")
  #     samp <- rtmvnorm(1, mean = mu, H = Hprecision, lower = c(0, rep(-Inf, p-1)),
  #                      algorithm = "gibbs")
  #     resamples <- resamples - 1
  #     if(resamples == 0) {
  #         print("reached 100 resamples, breaking out")
  #         samp <- mvrnorm(1, mu = mu, Sigma = covariance)
  #         break
  #     }
  # }
  # print("made it out of truncated sampling")
  # print(samp)
  # constrained <- TRUE
  # if(constrained) {
  #     lbound <- c(0, rep(-Inf, length(lambda) - 1))
  #     ubound <- c(rep(Inf, length(lambda)))
  #     #samp <- mvrnorm(1, mu = mu, Sigma = sigma)
  #     samp <- rtmvnorm(1, mean = mu, sigma = covariance, lower = lbound, upper = ubound,
  #                      algorithm = "gibbs",
  #                      burn.in.samples = 1000,
  #                      thinning = 5)
  #     #print(samp)
  #     # count <- 0
  #     # while(is.na(samp[1]) || samp[1] <= 0) {
  #     #     samp <- mvrnorm(1, mu = mu, Sigma = covariance)
  #     #     if(is.na(samp[1])) {
  #     #         mu = rep(0, p)
  #     #     }
  #     #     samp <- rtmvnorm(1, mean = mu, sigma = covariance, lower = lbound, upper = ubound,
  #     #                      algorithm = "gibbs",
  #     #                      burn.in.samples = 1000,
  #     #                      thinning = 5)
  #     #     if(count %% 1000 == 0 ) {
  #     #         stop(1)
  #     #     }
  #     #     count <- count + 1
  #     # }
  # }
  #samp <- (samp - mean(samp))/sd(samp)
  return(matrix(samp,1))
}

update_lambdas_spikeslab <- function(K, x, f, delta, alpha, lambda, pi_lambda) {
    n <- nrow(x)
    p <- ncol(x)
    x_hat <- matrix(0, n, ncol(x))
    for(j in 1:n) {
        x_hat[j,] <- x[j,] - alpha
    }
    
    deltafsum <- 1/delta + sum(f^2)
    mu <- 1/deltafsum * colSums(f*x_hat)
    covariance <- solve(deltafsum*K)
    Hprecision <- deltafsum*K
    # Need to force numerical symmetry...because of rounding errors?
    Hprecision[lower.tri(Hprecision)] <- t(Hprecision)[lower.tri(Hprecision)]

    samp <- rtmvnorm(1, mean = mu, H = Hprecision, lower = c(0, rep(-Inf, p-1)),
                     algorithm = "gibbs")
    #samp <- rtmvnorm(n, mu, sigma, lb, ub)
    
    resamples <- 100
    #print(samp)
    while((samp[1] < 0) | (is.nan(samp[1])) | is.infinite(samp[1])) {
        #samp <- mvrnorm(1, mu = mu, Sigma = covariance)
        print("truncated sampling failed")
        samp <- rtmvnorm(1, mean = mu, H = Hprecision, lower = c(0, rep(-Inf, p-1)),
                         algorithm = "gibbs")
        resamples <- resamples - 1
        if(resamples == 0) {
            print("reached 100 resamples, breaking out")
            samp <- mvrnorm(1, mu = mu, Sigma = covariance)
            break
        }
    }
    # With probability \pi we zero some of these out
    zeroes_lam <- rbinom(p, 1, pi_lambda)
    samp[which(zeroes_lam == 0)] <- 0

    return(samp)
}

update_pi_prob <- function(lambda) {
    shape1 <- 1 + sum(lambda != 0)
    shape2 <- 1 + length(lambda) - sum(lambda != 0)
    return(rbeta(1, shape1, shape2))
}

update_delta <- function(c, d, lambda, K) {
    p <- length(lambda)
    B <- chol(K)
    lambda_hat <- B %*% t(lambda)
    shape <- c + 0.5*p
    scale <- 0.5*(c*d + (t(lambda_hat)%*%lambda_hat))
    return(rinvgamma(1, shape = shape, rate = scale))
}

calculate_shifted_x <- function(x, alpha, lambda, f) {
    n <- nrow(x)
    x_hat <- matrix(0, n, ncol(x))
    for(j in 1:n) {
        x_hat[j,] <- x[j,] - alpha - lambda*f[j]
    }
    return(x_hat)
}

# copula_less/greaterthan_buckets - array of n x p x (max number of buckets)
# obs_buckets - matrix of n x p
## WORKING HERE - CURRENT PROBLEM: the subsetting below
# is as if each observation is "projected" onto its corresponding
# bucket. As a result, it is as if each observation has been
# replaced with the lower bound of the bucket (for lessthan_buckets)
# or the upper bound of the bucket (for greatthan_buckets). I believe
# this is fine, but this does not maintain the ordering WITHIN
# buckets, instead only maintains the ordering ACROSS buckets
# in the latent space.
# TODO: the subsetting can be done better since we only have 6 buckets
# the lb/ub will be the same within a bucket. Thus, we only need
# to calculate 6 lbs and 6 ubs. Should move this calculation to be outside.
# SHOULD THE UPPER BOUND TO BUCKET 6 BE INFINITY???? Right now
# bucket 6 contains the zis for anything that is 45+
###min(zi[buckets[,1] == 6])
# > max(zi[buckets[,1] == 6])
# [1] 46.55988
# > max(zi[buckets[,1] == 5])
# [1] -0.9553825
# > max(zi[buckets[,1] == 4])
# [1] -3.671026
# > max(zi[buckets[,1] == 3])
# [1] -7.418606
# > max(zi[buckets[,1] == 2])
# [1] -9.801843
# > max(zi[buckets[,1] == 1])
update_zs <- function(K, G, x, lambda, alpha, f, z, copula = FALSE,
                      lessthan_buckets, greaterthan_buckets, buckets) {
    n <- nrow(z)  # Number of observations
    p <- ncol(z)  # Number of variables
    
    # Preallocate new z matrix
    new_z <- matrix(0, n, p)
    
    if (copula) {
        # Copula-based update
        new_z <- z  # Start with the current z values
        
        for (i in 1:p) {
            zi <- new_z[, i]  # Extract current column
            sigma <- 1 / K[i, i]  # Compute standard deviation
            
            # Get unique bucket indices for the variable
            unique_buckets <- sort(unique(buckets[, i]))
            num_buckets <- length(unique_buckets)
            
            # Preallocate lower and upper bounds
            lbs <- rep(-Inf, num_buckets)
            ubs <- rep(Inf, num_buckets)
            
            # Compute bounds for each unique bucket
            for (b in 1:num_buckets) {
                bucket <- unique_buckets[b]
                
                # Find max value for the lower bound
                lessthan_current <- zi[lessthan_buckets[, i, bucket]]
                if (length(lessthan_current) > 0) {
                    lbs[b] <- max(lessthan_current)
                }
                
                # Find min value for the upper bound
                greaterthan_current <- zi[greaterthan_buckets[, i, bucket]]
                if (length(greaterthan_current) > 0) {
                    ubs[b] <- min(greaterthan_current)
                }
            }
            
            # Sample new values for z
            for (j in 1:n) {
                current_bucket <- buckets[j, i]
                mu <- alpha[i] + lambda[i] * f[j] - 
                    sum((K[i, ] / K[i, i] * (z[j, ] - alpha - lambda * f[j]))[G[i, ] == 1])
                
                new_z[j, i] <- rtruncnorm(1, a = lbs[current_bucket], 
                                          b = ubs[current_bucket], 
                                          mean = mu, sd = sigma)
            }
        }
    } else {
        # Standard (non-copula) update
        for (j in 1:n) {
            for (i in 1:p) {
                mu <- alpha[i] + lambda[i] * f[j] - 
                    sum((K[i, ] / K[i, i] * (z[j, ] - alpha - lambda * f[j]))[G[i, ] == 1])
                sigma <- sqrt(1 / K[i, i])
                
                # Determine truncation bounds based on x[j, i]
                lower_bound <- ifelse(x[j, i], 0, -Inf)
                upper_bound <- ifelse(x[j, i], Inf, 0)
                
                new_z[j, i] <- rtruncnorm(1, a = lower_bound, 
                                          b = upper_bound, 
                                          mean = mu, sd = sigma)
            }
        }
    }
    
    return(new_z)
}

# update_zs <- function(K, G, x, lambda, alpha, f, z, copula = FALSE,
#                       lessthan_buckets, greaterthan_buckets,
#                       buckets) {
#     n <- nrow(z)
#     p <- ncol(z)
#     print(copula)
#     if(copula) {
#         new_z <- z
#         for(i in 1:p) {
#             zi <- new_z[,i]
#             sigmas <- 1/K[i,i]
#             
#             uniquebuckets <- sort(unique(buckets[,i]))
#             lbs <- vector("numeric", length = length(uniquebuckets))
#             ubs <- vector("numeric", length = length(uniquebuckets))
#             
#             for(b in 1:length(uniquebuckets)) {
#                 bucket <- uniquebuckets[b]
#                 lessthan_current <- zi[lessthan_buckets[,i,bucket]]
#                 greaterthan_current <- zi[greaterthan_buckets[,i,bucket]]
#                 if(length(lessthan_current) > 0) {
#                     lbs[b] <- max(lessthan_current)
#                     # currently lb does not match exactly with the max of the previous
#                     # groups, this is due to the intervals for groups being v <= x <- v+1
#                     # but the copula formulation not considering the equality.
#                     # Possible solution: make the breakpoints +/- 0.5.
#                 } else {
#                     lbs[b] <- -Inf
#                 }
#                 if(length(greaterthan_current) > 0) {
#                     ubs[b] <- min(greaterthan_current)
#                 } else {
#                     ubs[b] <- Inf
#                 }
#             }
#             
#             for(j in 1:n) {
#                 # copula with only binary
#                 # if(x[j,i]) {
#                 #     ub <- min(ub, zi[j])
#                 # }
#                 # else {
#                 #     lb <- max(lb, zi[j])
#                 # }
#                 # lessthan_current <- zi[x[,i] < x[j,i]]
#                 # greaterthan_current <- zi[x[,i] > x[j,i]]
#                 current_bucket <- buckets[j,i]
#                 mu <- alpha[i] + lambda[i]*f[j] - sum((K[i,]/K[i,i] * (z[j,] - alpha - lambda*f[j]))[G[i,]==1])
#                 new_z[j,i] <- rtruncnorm(1, a = lbs[current_bucket], 
#                                          b = ubs[current_bucket], 
#                                          mean = mu, sd = sigmas)
#             }
#         }
#     } else {
#         new_z <- matrix(0, n, p)
#         for(j in 1:n) {
#             for(i in 1:p) {
#                 mu <- alpha[i] + lambda[i]*f[j] - sum((K[i,]/K[i,i] * (z[j,] - alpha - lambda*f[j]))[G[i,]==1])
#                 sigmas <- 1/K[i,i]
#                 if(x[j,i]) {
#                     new_z[j,i] <- rtruncnorm(1, a = 0, mean = mu, sd = sqrt(sigmas))
#                 }
#                 else {
#                     new_z[j,i] <- rtruncnorm(1, b = 0, mean = mu, sd = sqrt(sigmas))
# 
#                 }
#             }
#         }
#     }
#     #new_z <- scale(new_z, TRUE, FALSE)
#     return(new_z)
# }

################################################################################
########## Updates and helper functions specific to the GDP(3, 1) prior

calculate_lambda_hat_gdp <- function(psi_hat, f, x_hat_col) {
    s <- 0
    for(j in 1:length(x_hat_col)) {
        s <- s + f[j] * x_hat_col[j]
    }
    return(psi_hat * s)
}

update_lambdas_gdp <- function(K, x, f, psi, alpha) {
    new_lambdas <- vector("numeric", length = ncol(x))
    p <- ncol(x)
    x_hat <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    T_psi <- diag(psi)
    for(j in 1:nrow(x)) {
        x_hat[j,] <- x[j,] - alpha
    }
    covar <- solve(sum(f^2)*K + solve(T_psi))
    mu <- as.numeric(covar %*% (K %*% colSums(f*x_hat)))
    # print(covar)
    # print(mu)
    
    new_lambdas <- rtmvnorm(1, mean = mu, sigma = covar, lower = c(0, rep(-Inf, p-1)),
                     algorithm = "gibbs")
    #samp <- rtmvnorm(n, mu, sigma, lb, ub)
    
    resamples <- 100
    #print(samp)
    while((new_lambdas[1] < 0) | (is.nan(new_lambdas[1])) | is.infinite(new_lambdas[1])) {
        #new_lambdas <- mvrnorm(1, mu = mu, Sigma = covar)
        print("truncated sampling failed")
        new_lambdas <- rtmvnorm(1, mean = mu, sigma = covar, lower = c(0, rep(-Inf, p-1)),
                         algorithm = "gibbs")
        resamples <- resamples - 1
        if(resamples == 0) {
            print("reached 100 resamples, breaking out")
            new_lambdas <- mvrnorm(1, mu = mu, Sigma = covar)
            break
        }
    }
    #new_lambdas <- mvtnorm::rmvnorm(1, mean = mu, sigma = covar)
    return(new_lambdas)
}

update_psi <- function(lambda, xi) {
    new_psi <- vector("numeric", length = length(lambda))
    for(i in 1:length(lambda)) {
        inv_psi <- rinvgauss(1, mean = abs(xi[i]/lambda[i]), shape = xi[i]^2)
        new_psi[i] <- 1/inv_psi
    }
    return(new_psi)
}

update_xi <- function(lambda, a, b) {
    new_xi <- vector("numeric", length = length(lambda))
    for(i in 1:length(lambda)) {
        new_xi[i] <- rgamma(1, a + 1, b + abs(lambda[i]))
    }
    return(new_xi)
}
################################################################################

# Here we are only using the first row of z (the continuous latent data). Should
# we be using more?
# latenttoobserved <- function(K, G, x, lambda, alpha, f, z) {
#     xpred <- matrix(0, nrow(z), ncol(z))
#     for(i in 1:nrow(z)) {
#         # we draw f_j from N(0,1) for a new observation
#         fj <- rnorm(1)
#         for(j in 1:ncol(z)) {
#             mu <- alpha[j] + lambda[j]*fj - sum((K[j,]/K[j,j] * (z[i,j] - alpha - lambda*fj))[G[j,]==1])
#             sigmas <- 1/K[i,i]
#             
#             # Convert latent to observed space without truncation.
#             xpred[i,j] <- as.numeric(rnorm(1, mu, sqrt(sigmas)) > 0)
#             
#         }
#     }
#     
#     return(xpred)
# }


UpdateExpectedCellsSampler <- function(NumberOfVariables,
                                mysamples,
                                CumulativeMarginalProbs)
{
    ExpectedCellsVector <- vector("numeric", length = NumberOfVariables)
    for(i in 1:NumberOfVariables)
    {
        s = pnorm(mysamples[i])#gsl_cdf_ugaussian_P(gsl_matrix_get(mysamples,i,0));
        q1 = 0;
        q2 = 3; # This number is 1 + dimensions for the i'th variable. Here we hardcode
        # for binary variables.
        q <- 2
        if(CumulativeMarginalProbs[i,q]<s)
        {
            ExpectedCellsVector[i] = 1;
        }
        else
        {
            ExpectedCellsVector[i] = 0;
        }
        
    }	
    # s = ExpectedCells->Get();
    # ExpectedCells->Set(1+s);
    
    return(ExpectedCellsVector)
}
UpdateExpectedCellsPosNegSampler <- function(NumberOfVariables,
                                      mysamples,
                                      CumulativeMarginalProbs)
{
    ExpectedCellsVector <- vector("numeric", length = NumberOfVariables)
    for(i in 1:NumberOfVariables)
    {
        if(mysamples[i] > 0)
        {
            ExpectedCellsVector[i] = 1;
        }
        else
        {
            ExpectedCellsVector[i] = 0;
        }
        
    }	
    # s = ExpectedCells->Get();
    # ExpectedCells->Set(1+s);
    
    return(ExpectedCellsVector)
}

randomMVNSampler <- function(nSamplesToGenerate,
                      NumberOfVariables,
                      corrData,
                      CumulativeMarginalProbs,
                      alpha, lambda, f)
{
    ExpectedCells <- matrix(0, nSamplesToGenerate, NumberOfVariables)
    for(asample in 1:nSamplesToGenerate)
    {
        #yMatrix <- mvrnorm(1, mu = alpha + lambda*f[asample], Sigma = corrData)
        yMatrix <- mvrnorm(1, mu = alpha + lambda*rnorm(1), Sigma = corrData)
        
        ExpectedCells[asample,] <- UpdateExpectedCellsSampler(NumberOfVariables,
                                                       yMatrix,
                                                       CumulativeMarginalProbs)
    }
    return(ExpectedCells)
    
}

# Set the precision matrix to actual zeros based on the zeros in G
setZeros <- function(G,K) {
    K[G == 0] <- 0
    return(K)
}



# # Gibbs sampler - First assume a slab and spike prior

gibbs_sampler_ss <- function(steps, x, alpha_start, lambda_start, f_start, n, prior, 
                             prior_params, prior_ratio, identifiable, G1, K1, 
                             delta_start, pi_prob, c, d, n0, burnin = steps/2, normalmodel = FALSE,
                             probit = FALSE, z = NULL) {
    p <- ncol(x)
    n <- nrow(x)
    print(p)
    print(n)
    samplecor <- cor(x)
    
    f_chain <- matrix(0, steps, n)
    #f_chain[1,] <- f_start
    f_chain[1,] <- rnorm(n, 0, 1)
    
    alpha_chain <- matrix(0, steps, p)
    #alpha_chain[1,] <- alpha_start
    alpha_chain[1,] <- rmvnorm(1, rep(0, p), 1/n0*diag(1, p, p))
    
    delta_chain <- vector("numeric", length = steps)
    #delta_chain[1] <- delta_start
    delta_chain[1] <- rinvgamma(1, shape = c, scale = (c*d)/2)
    
    acceptprob_chain <- matrix(0, steps, p*(p-1)/2)
    
    graph_chain <- vector("list", steps)
    # Initialize the Graph. sample from the given prior.
    graph <- matrix(0, p, p)
    G <- prior(graph, prior_params)$graph
    if(identifiable) {
        while(!is_identifiable(G)) {
            G <- prior(graph, prior_params)$graph
        }
    }
    
    graph_chain[[1]] <- G
    #G <- G1
    del <- 3
    D <- diag(p)
    deltan <- del + n
    K_chain <- vector("list", steps)
    K <- my_rgwish(1, G, del, D)
    K_diag <- diag(K)
    K[G == 0] <- 0
    diag(K) <- K_diag
    #K <- K1
    K_chain[[1]] <- K
    Gamma_chain <- vector("list", steps)
    Gamma_chain[[1]] <- cov2cor(solve(K))
    print(K)
    print(solve(K))
    #Gamma_chain[[1]] <- diag(1/sqrt(diag(solve(K))))%*%solve(K)%*%diag(1/sqrt(diag(solve(K))))
    
    lbound <- c(0, rep(-Inf, p - 1))
    ubound <- c(rep(Inf, p))
    lambda_chain <- matrix(0, steps, p)
    #lambda_chain[1,] <- lambda_start
    lambda_chain[1,] <- rtmvnorm(1, mean = rep(0, p), sigma = delta_chain[1]*solve(K), 
             lower = lbound, upper = ubound)
    
    print("alpha, lambda:")
    print(alpha_chain[1,])
    print(lambda_chain[1,])
    count <- 1
    
    
    # For conforming to the data types
    pi_prob_chain <- vector("numeric", length = steps)
    #pi_prob_chain[1] <- pi_prob
    pi_prob_chain[1] <- sample_pi_prob(1,1,1)
    #pi_prob_chain[1] <- rbeta(1, 1, 1)
    psi_chain <- matrix(0, steps, p)
    xi_chain <- matrix(0, steps, p)
    
    xi_chain[1,] <- rgamma(p, 3, 1)
    for(i in 1:p) {
        psi_chain[1,i] <- rexp(1, (xi_chain[1,i]^2)/2)
    }
    kls <- vector("numeric", steps)
    kls[1] <- 0
    
    
    ### Probit Model
    L <- 1
    zs <- vector("list", steps)
    zs[[1]] <- z
    # going to sample 1 random mvns at each step
    xpred <- array(0, c(steps, p, L))
    cumulativeMargProb <- matrix(0, p, 3)
    if(probit) {
        #calculate marginal probability of observed
        cumulativeMargProb[,2] <- 1-(colSums(x)/n)
        cumulativeMargProb[,3] <- rep(1, p)
        copy_x <- x
        x <- z
    }
    
    for(s in 2:(steps)) {
        shifted <- calculate_shifted_x(x, alpha_chain[s-1,], lambda_chain[s-1,], f_chain[s-1,])
        # # sample scatter matrix:
        #S <- t(shifted) %*% shifted
        S <- (n-1)*cov(shifted)
        #S <- n*solve(K_chain[[s-1]])
        dcbf <- ggm_dcbf(G, D, S, K, p,
                         del, deltan, prior, prior_params, prior_ratio,
                         identifiable)
        G <- dcbf$G
        K <- dcbf$K
        # print(K)
        # print(G)
        K_diag <- diag(K)
        K[G == 0] <- 0
        diag(K) <- K_diag
        
        graph_chain[[s]] <- G
        K_chain[[s]] <- K
        #acceptprob_chain[s,] <- dcbf$alphas
        Gamma_chain[[s]] <- cov2cor(solve(K))
        #Gamma_chain[[s]] <- diag(1/sqrt(diag(solve(K))))%*%solve(K)%*%diag(1/sqrt(diag(solve(K))))
        # We try scaling the other way:
        #Gamma_chain[[s]] <- cov2cor(K)
        
        if(probit) {
            K <- solve(Gamma_chain[[s]])
            #K <- Gamma_chain[[s]]
        }
        
        f_chain[s,] <- update_factors(K, x, lambda_chain[s-1,], alpha_chain[s-1,])
        
        
        alpha_chain[s,] <- update_alphas(K, x, n, n0, f_chain[s,], lambda_chain[s-1,])
        
        l_count <- 0
        while(sum(lambda_chain[s,]^2) == 0) {
            # lambda_chain[s,] <- update_lambdas(K, x, f_chain[s,], delta_chain[s-1],
            #                                         alpha_chain[s,])
            # lambda_chain[s,] <- update_lambdas_spikeslab(K, x, f_chain[s,], delta_chain[s-1],
            #                                    alpha_chain[s,], lambda_chain[s-1,], pi_prob_chain[s-1])
            #     
            lambda_chain[s,] <- update_lambdas_gdp(K, x, f_chain[s,],
                                                   psi_chain[s-1,], alpha_chain[s,])
                
            l_count <- l_count + 1
            if(l_count > 2) {
                print("norm is 0")
                print(l_count)
            }
                
        }
        
        delta_chain[s] <- update_delta(c, d, lambda_chain[s,], K)
        
        xi_chain[s,] <- update_xi(lambda_chain[s,], 1, 3)
        psi_chain[s,] <- update_psi(lambda_chain[s,], xi_chain[s,])
        
        pi_prob_chain[s] <- update_pi_prob(lambda_chain[s,])
        
        ## Probit Model update
        if(probit) {
            zs[[s]] <- update_zs(K, G, copy_x, lambda_chain[s,], 
                                 alpha_chain[s,], f_chain[s,], x)
            x <- zs[[s]]
            ######################################
            # Trying model averaging for the expected cell counts
            # xpred[((s-1)*n):((s*n)-1),,1] <- randomMVNSampler(n, p, solve(K), cumulativeMargProb,
            #                                              alpha_chain[s,],lambda_chain[s,],
            #                                              f_chain[s,])
            xpred[s,,1] <- randomMVNSampler(1, p, solve(K), cumulativeMargProb,
                                                              alpha_chain[s,],lambda_chain[s,],
                                                              f_chain[s,])
            ######################################
            
        }
        
        count <- count + 1
        if(count %% 1000 == 0 ) {
            print(count)
        }
        print(count)
        
        
        
    }
    return(list(K = apply(simplify2array(K_chain[-c(1:burnin)]), c(1:2), mean),
                G =  apply(simplify2array(graph_chain[-c(1:burnin)]), c(1:2), mean),
                f = colMeans(f_chain[-c(1:burnin),]), 
                alpha = colMeans(alpha_chain[-c(1:burnin),]), 
                lambda = colMeans(lambda_chain[-c(1:burnin),]), 
                delta = mean(delta_chain[-c(1:burnin)]),
                lambda_chain = lambda_chain,
                f_chain = f_chain,
                graph_chain = graph_chain,
                K_chain = K_chain,
                alpha_chain = alpha_chain,
                tau_chain = delta_chain,
                pi_chain = pi_prob_chain,
                acceptprob_chain = acceptprob_chain,
                z_chain = zs,
                Gamma_chain = Gamma_chain,
                xpred = xpred))
}


