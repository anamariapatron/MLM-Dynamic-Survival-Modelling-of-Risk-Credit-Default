# Function to simulate from a landmark proportional hazards model
# seed: seed to be used in the simulation
# DES: a list of design matrices
# Betas: a matrix containing the regression coefficients at each landmark point
# Thetas: a matrix containing the parameters of the baseline hazard at each landmark point
# LMs: landmark points, they must contain 0 as the first point.
# dist: baseline hazard distribution (W, LN, G)

simLMPH <- function(seed, DES, Betas, Thetas, LMs, dist) {
  set.seed(seed)
  
  # Number of landmark intervals
  K <- length(LMs)-1
  n <- nrow(DES[[1]])  # infer number of individuals from first design matrix
  
  # Add time 0 to start the timeline
  tpoints <- as.vector(LMs)
  time <- rep(0, n)  # time since origin
  alive <- rep(TRUE, n)  # track who is still alive
  out_time <- rep(NA, n)  # to store final failure times
  
  for (k in 1:K) {
    des_k <- DES[[k]]
    beta_k <- Betas[k, ]
    theta_k <- Thetas[k, ]
    
    # Linear predictor
    linpred <- des_k %*% beta_k
    mlambda <- exp(-linpred)
    
    # Get baseline survival and quantile functions
    if (dist == "W") {
      S0f <- function(t) pweibull(t, shape = theta_k[1], scale = theta_k[2], lower.tail = FALSE)
      quantf <- function(p) qweibull(p, shape = theta_k[1], scale = theta_k[2])
    } else if (dist == "LN") {
      S0f <- function(t) plnorm(t, meanlog = theta_k[1], sdlog = theta_k[2], lower.tail = FALSE)
      quantf <- function(p) qlnorm(p, meanlog = theta_k[1], sdlog = theta_k[2])
    } else if (dist == "G") {
      S0f <- function(t) pgamma(t, shape = theta_k[1], rate = theta_k[2], lower.tail = FALSE)
      quantf <- function(p) qgamma(p, shape = theta_k[1], rate = theta_k[2])
    } 
    else {
      stop("Unsupported distribution.")
    }
    
    for (i in 1:n) {
      if (!alive[i]) next
      
      u <- runif(1)
      mlambda_i <- mlambda[i]
      t_prev <- tpoints[k]
      
      # Survival up to previous time point
      S_L <- S0f(t_prev)
      p_target <- 1 - exp( log(u*S_L)*mlambda_i) 
      
      # If p_target < 0, this means survival probability underflows. Cap at small epsilon.
      p_target <- max(p_target, .Machine$double.eps)
      t_rel <- quantf(p_target)
      
      # If failure occurs before next landmark, store and mark dead
      if (t_rel <= LMs[k]) {
        out_time[i] <- t_rel
        alive[i] <- FALSE
      }
      # else continue to next landmark (do nothing)
    }
  }
  
  # For those who survive all landmark intervals, simulate one final time
  if (any(alive)) {
    k <- K  # final interval (extrapolation)
    des_k <- DES[[k]]
    beta_k <- Betas[k, ]
    theta_k <- Thetas[k, ]
    
    linpred <- des_k %*% beta_k
    lambda <- exp(linpred)
    
    if (dist == "W") {
      S0f <- function(t) pweibull(t, shape = theta_k[1], scale = theta_k[2], lower.tail = FALSE)
      quantf <- function(p) qweibull(p, shape = theta_k[1], scale = theta_k[2])
    } else if (dist == "LN") {
      S0f <- function(t) plnorm(t, meanlog = theta_k[1], sdlog = theta_k[2], lower.tail = FALSE)
      quantf <- function(p) qlnorm(p, meanlog = theta_k[1], sdlog = theta_k[2])
    } else if (dist == "G") {
      S0f <- function(t) pgamma(t, shape = theta_k[1], rate = theta_k[2], lower.tail = FALSE)
      quantf <- function(p) qgamma(p, shape = theta_k[1], rate = theta_k[2])
    } 
    
    for (i in which(alive)) {
      u <- runif(1)
      t_prev <- tpoints[K+1]
      S_L <- S0f(t_prev)
      mlambda_i <- mlambda[i]
      p_target <- 1 - exp( log(u*S_L)*mlambda_i)
      p_target <- max(p_target, .Machine$double.eps)
      t_rel <- quantf(p_target)
      out_time[i] <- t_rel
    }
  }
  
  return(out_time)
}




library(HazReg)

plotLMPH <- function(individuals, DES, Betas, Thetas, LMs, dist = "W", ymax, lout) {
  
  # Time grid for smooth plotting
  time_grid <- seq(0, max(LMs), length.out = lout)
  
  # Add 0 to LMs to define interval edges
  tpoints <- as.vector(LMs)
  K <- length(tpoints)
  
  # Hazard function by distribution
  get_basehaz <- function(t, theta, dist) {
    if (dist == "W") {
      # Weibull hazard
      
      return(hweibull(t, theta[1], theta[2]))
    }
    else if (dist == "LN") {
      
      return(hlnorm(t, theta[1], theta[2]))
    }
    else if (dist == "G") {
      
      return(hgamma(t, theta[1], theta[2]))
    }
    else {
      stop("Distribution not implemented.")
    }
  }
  
  # Colors for plotting
  colors <- rainbow(length(individuals))
  
  # Begin plotting
  plot(NULL, xlim = c(0, max(time_grid)), ylim = c(0, ymax),
       xlab = "Time", ylab = "Hazard", main = "Individual-Specific Hazard Functions")
  
  for (idx in seq_along(individuals)) {
    i <- individuals[idx]
    hazards <- numeric(length(time_grid))
    
    for (k in 1:K) {
      t_start <- tpoints[k]
      t_end <- ifelse(k < K, tpoints[k + 1], max(time_grid))
      tk_grid <- time_grid[time_grid >= t_start & time_grid < t_end]
      
      if (length(tk_grid) == 0) next
      
      des_k <- DES[[k]]
      beta_k <- Betas[k, ]
      theta_k <- Thetas[k, ]
      
      linpred <- as.numeric(des_k[i, ] %*% beta_k)
      base_haz <- get_basehaz(tk_grid, theta_k, dist)
      haz_k <- base_haz * exp(linpred)
      
      hazards[time_grid %in% tk_grid] <- haz_k
    }
    
    lines(time_grid, hazards, col = colors[idx], lwd = 2)
  }
  
  legend("topright", legend = paste("Individual", individuals), col = colors,
         lty = 1, lwd = 2)
}


