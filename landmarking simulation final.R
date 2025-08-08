# CLEAR WORKSPACE
rm(list = ls())

# LOAD LIBRARIES
library(HazReg)
library(dplyr)
library(lubridate)

set.seed(1234)

# SIMULATION PARAMETERS
n <- 1000                   # number of individuals
n_iterations <- 5          # number of landmarking iterations
base_date <- as.Date("2020-01-01")  # reference date for observation windows

# BASELINE COVARIATES (2 variables)
x <- as.matrix(cbind(rnorm(n), rnorm(n)))

# LANDMARK TIMES
landmarks <- seq(0, n_iterations * 1, by = 1)  # 0, 12, 24, ..., 60
K <- length(landmarks) - 1  # number of intervals


# COEFFICIENTS (BETAS) PER INTERVAL, b pequeños
betas <- matrix(c(-0.2, 0.3,
                  -0.2, 0.3,
                  -0.2, 0.3,
                  -0.2, 0.3,
                  -0.2, 0.3), nrow = K, byrow = TRUE)


thetas <- matrix(c(0.5, 100,
                   0.5, 90,
                   0.5, 80,
                   0.5, 70,
                   0.5, 60), nrow = K, byrow = TRUE)

betas = rbind(c(0.1,-0.2),c(0.1,-0.2),c(0.1,-0.2),c(0.15,-0.25),c(0.15,-0.25),c(0.15,-0.25))/2

thetas = rbind(c(10,1.5),c(11,1.7),c(12,1.9),c(13,2),c(14,3),c(15,4))*1.5
# MAIN DATAFRAME SETUP
df_wide <- data.frame(individual = 1:n)
df_wide$initial_date <- as_date("2010-01-01") + days(sample(0:364, n, replace = TRUE))

# CREATE COLUMNS FOR EACH ITERATION
for (i in 1:n_iterations) {
  df_wide[[paste0("event_time_iteration_", i)]] <- NA_real_
  df_wide[[paste0("observation_date_iteration_", i)]] <- as.Date(NA)
  df_wide[[paste0("status_iteration_", i)]] <- NA_integer_
  df_wide[[paste0("x1_iteration_", i)]] <- rnorm(10) #covariate 1
  df_wide[[paste0("x2_iteration_", i)]] <- rnorm(10) #covariate 2
}

# TRACK WHO IS STILL ALIVE
alive <- rep(TRUE, n)

# LANDMARK-BASED SIMULATION FUNCTION (UPDATED)
simLMPH <- function(seed, x, Betas, Thetas, LMs, dist, alive_init) {
  set.seed(seed)
  K <- length(LMs) - 1
  n <- sum(alive_init)
  tpoints <- as.vector(LMs)
  out_time <- rep(NA, n)
  alive <- alive_init
  
  for (k in 1:K) {
    survivors <- which(alive)
    if (length(survivors) == 0) break
    
    # Subset design matrix for survivors in this interval
    des_k <- x[survivors, , drop = FALSE]
    beta_k <- Betas[k, ]
    theta_k <- Thetas[k, ]
    linpred <- des_k %*% beta_k
    mlambda <- exp(-linpred)
    
    # Baseline functions
    if (dist == "W") {
      S0f <- function(t) pweibull(t, shape = theta_k[1], scale = theta_k[2], lower.tail = FALSE)
      quantf <- function(p) qweibull(p, shape = theta_k[1], scale = theta_k[2])
    } else {
      stop("Only Weibull distribution is implemented.")
    }
    
    for (j in seq_along(survivors)) {
      i <- survivors[j]
      u <- runif(1)
      t_prev <- tpoints[k]
      S_L <- S0f(t_prev)
      p_target <- 1 - exp(log(u * S_L) * mlambda[j])
      p_target <- max(p_target, .Machine$double.eps)
      t_rel <- quantf(p_target)
      if (t_rel <= LMs[k + 1]) {
        out_time[i] <- t_rel
        alive[i] <- FALSE
      }
    }
  }
  
  # Final extrapolation for remaining survivors
  if (any(alive)) {
    survivors <- which(alive)
    des_k <- x[survivors, , drop = FALSE]
    beta_k <- Betas[K, ]
    theta_k <- Thetas[K, ]
    linpred <- des_k %*% beta_k
    mlambda <- exp(-linpred)
    
    S0f <- function(t) pweibull(t, shape = theta_k[1], scale = theta_k[2], lower.tail = FALSE)
    quantf <- function(p) qweibull(p, shape = theta_k[1], scale = theta_k[2])
    
    for (j in seq_along(survivors)) {
      i <- survivors[j]
      u <- runif(1)
      t_prev <- tpoints[K + 1]
      S_L <- S0f(t_prev)
      p_target <- 1 - exp(log(u * S_L) * mlambda[j])
      p_target <- max(p_target, .Machine$double.eps)
      t_rel <- quantf(p_target)
      out_time[i] <- t_rel
    }
  }
  
  return(out_time)
}
# MAIN SIMULATION LOOP (CORREGIDO)
alive <- rep(TRUE, nrow(df_wide))  # Inicialización

for (i in seq_len(n_iterations)) {
  survivors <- which(alive)
  if (length(survivors) == 0) break
  
  cat("Iteration", i, "- survivors:", length(survivors), "\n")
  
  # Simulate event times (la función ya usa submatrices internamente)
  sim_times <- simLMPH(
    seed = 1000 + i * 10,
    x = x[survivors, , drop = FALSE],
    Betas = betas,
    Thetas = thetas,
    LMs = landmarks,
    dist = "W",
    alive_init = rep(TRUE, length(survivors))
  )
  
  # Generate observation dates randomly within the 1-month window
  month_start <- base_date %m+% months((i - 1) * 12)
  month_end <- month_start %m+% months(1) - days(1)
  
  obs_dates <- as.Date(runif(length(survivors),
                             min = as.numeric(month_start),
                             max = as.numeric(month_end)),
                       origin = "1970-01-01")
  
  # Months since initial date
  months_elapsed <- interval(df_wide$initial_date[survivors], obs_dates) %/% months(1)
  
  # Event status (1 = event, 0 = censored)
  status_vals <- as.integer(months_elapsed >= sim_times)
  
  # Save results to dataframe
  df_wide[[paste0("event_time_iteration_", i)]][survivors] <- sim_times
  df_wide[[paste0("observation_date_iteration_", i)]][survivors] <- obs_dates
  df_wide[[paste0("status_iteration_", i)]][survivors] <- status_vals
  
  # Update alive vector
  alive[survivors] <- status_vals == 0
}


for (i in 1:n_iterations) {
  status_col <- df_wide[[paste0("status_iteration_", i)]]
  if (all(is.na(status_col))) {
    cat("Iteración", i, ": sin datos (NA)\n")
  } else {
    cat("Iteración", i, ":",
        round(mean(status_col, na.rm = TRUE), 3),
        "muertos\n")
  }
}


# -------------------------------------------------------------
# PLOT INDIVIDUAL-SPECIFIC HAZARDS
# -------------------------------------------------------------

DES <- list()
for (i in 1:n_iterations) {
  combined_matrix <- cbind(df_wide[[paste0("x1_iteration_", i)]], df_wide[[paste0("x2_iteration_", i)]])
  DES[[i]] <- combined_matrix
}


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


plotLMPH(individuals = c(1,5), DES = DES, Betas = betas, Thetas = thetas, LMs = landmarks, dist = "W", ymax  = 1, lout = 1000)

plotLMPH(individuals = c(1,5,6,7), DES = DES, Betas = betas, Thetas = thetas, LMs = landmarks, dist = "W", ymax  = 0.01, lout = 1000)


