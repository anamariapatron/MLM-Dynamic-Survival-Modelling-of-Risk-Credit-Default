# CLEAR WORKSPACE
rm(list = ls())

# LOAD LIBRARIES
library(HazReg)
library(dplyr)
library(lubridate)
setwd("~/GitHub/MLM-On-The-Applications-Of-Survival-Modelling-In-The-Mortgage-Sector")

set.seed(1234)

# SIMULATION PARAMETERS
n <- 1000                   # number of individuals
n_iterations <- 5           # number of landmarking iterations
base_date <- as.Date("2020-01-01")  # reference date for observation windows

# BASELINE COVARIATES (2 variables)
x <- as.matrix(cbind(rnorm(n), rnorm(n)))

# LANDMARK TIMES
landmarks <- seq(0, n_iterations * 1, by = 1)  # e.g., 0,1,2,3,4,5
K <- length(landmarks) - 1  # number of intervals

# Adjusted coefficients (betas) - smaller magnitude
betas <- matrix(c(-0.1, 0.15,
                  -0.1, 0.15,
                  -0.1, 0.15,
                  -0.1, 0.15,
                  -0.1, 0.15),
                nrow = K, byrow = TRUE)

# Adjusted Weibull parameters - larger scales (longer survival)
thetas <- matrix(c(0.5, 400,
                   0.5, 280,
                   0.5, 360,
                   0.5, 400,
                   0.5, 450), nrow = K, byrow = TRUE)


# SETUP MAIN DATAFRAME
df_wide <- data.frame(individual = 1:n)
df_wide$initial_date <- as_date("2010-01-01") + days(sample(0:364, n, replace = TRUE))

# CREATE COLUMNS FOR EACH ITERATION
for (i in 1:n_iterations) {
  df_wide[[paste0("event_time_iteration_", i)]] <- NA_real_
  df_wide[[paste0("observation_date_iteration_", i)]] <- as.Date(NA)
  df_wide[[paste0("status_iteration_", i)]] <- NA_integer_
  # Covariates for each iteration (random normal values)
  df_wide[[paste0("x1_iteration_", i)]] <- rnorm(n)
  df_wide[[paste0("x2_iteration_", i)]] <- rnorm(n)
}

# TRACK WHO IS STILL ALIVE
alive <- rep(TRUE, n)

# LANDMARK-BASED SIMULATION FUNCTION (UPDATED WITH TRUNCATED UNIFORM)
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
    
    des_k <- x[survivors, , drop = FALSE]
    beta_k <- Betas[k, ]
    theta_k <- Thetas[k, ]
    linpred <- des_k %*% beta_k
    mlambda <- exp(-linpred)
    
    if (dist == "W") {
      # Baseline survival and quantile functions for Weibull
      S0f <- function(t) pweibull(t, shape = theta_k[1], scale = theta_k[2], lower.tail = FALSE)
      quantf <- function(p) qweibull(p, shape = theta_k[1], scale = theta_k[2])
    } else {
      stop("Only Weibull distribution is implemented.")
    }
    
    for (j in seq_along(survivors)) {
      i <- survivors[j]
      # Calculate baseline survival at landmark time L_k (truncation point)
      S_L <- S0f(tpoints[k])
      # Generate truncated uniform variable between 0 and S(L_k)
      u <- runif(1, min = 0, max = S_L)
      # Adjust survival probability for individual with covariates
      p_target <- u^(1 / mlambda[j])  # invert the conditional survival
      p_target <- max(p_target, .Machine$double.eps)  # avoid zero
      
      # Calculate event time via inverse baseline survival
      t_rel <- quantf(p_target)
      
      # Check if event time falls within current interval
      if (t_rel <= tpoints[k + 1]) {
        out_time[i] <- t_rel
        alive[i] <- FALSE
      }
    }
  }
  
  # Final extrapolation for survivors who haven't failed yet
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
      # At final step, truncate at max time (e.g., last landmark)
      S_L <- S0f(tpoints[K + 1])
      u <- runif(1)
      p_target <- u^(1 / mlambda[j])
      p_target <- max(p_target, .Machine$double.eps)
      t_rel <- quantf(p_target)
      out_time[i] <- t_rel
    }
  }
  
  return(out_time)
}

# MAIN SIMULATION LOOP WITH UPDATED FUNCTION
alive <- rep(TRUE, nrow(df_wide))

for (i in seq_len(n_iterations)) {
  survivors <- which(alive)
  if (length(survivors) == 0) break
  
  cat("Iteration", i, "- survivors:", length(survivors), "\n")
  
  # Use covariates at iteration i
  covariates_i <- cbind(df_wide[[paste0("x1_iteration_", i)]],
                        df_wide[[paste0("x2_iteration_", i)]])
  
  # Subset covariates to survivors
  x_survivors <- covariates_i[survivors, , drop = FALSE]
  
  # Simulate event times using truncated uniform method
  sim_times <- simLMPH(
    seed = 1234,
    x = x_survivors,
    Betas = betas,
    Thetas = thetas,
    LMs = landmarks,
    dist = "W",
    alive_init = rep(TRUE, length(survivors))
  )
  
  # Generate observation dates randomly within the 1-year window for this iteration
  month_start <- base_date %m+% months((i - 1) * 12)
  month_end <- month_start %m+% months(1) - days(1)
  
  obs_dates <- as.Date(runif(length(survivors),
                             min = as.numeric(month_start),
                             max = as.numeric(month_end)),
                       origin = "1970-01-01")
  
  # Calculate months elapsed since initial date to observation date
  months_elapsed <- interval(df_wide$initial_date[survivors], obs_dates) %/% months(1)
  
  # Event status: 1 = event occurred, 0 = censored
  status_vals <- as.integer(months_elapsed >= sim_times)
  
  # Save results to dataframe
  df_wide[[paste0("event_time_iteration_", i)]][survivors] <- sim_times
  df_wide[[paste0("observation_date_iteration_", i)]][survivors] <- obs_dates
  df_wide[[paste0("status_iteration_", i)]][survivors] <- status_vals
  
  # Update who is still alive for next iteration
  alive[survivors] <- status_vals == 0
}

# Print summary of deaths per iteration
for (i in 1:n_iterations) {
  status_col <- df_wide[[paste0("status_iteration_", i)]]
  if (all(is.na(status_col))) {
    cat("Iteration", i, ": no data (NA)\n")
  } else {
    cat("Iteration", i, ":",
        round(mean(status_col, na.rm = TRUE), 3),
        "deaths\n")
  }
}

# -------------------------------------------------------------
# PLOT INDIVIDUAL-SPECIFIC  HAZARDS
# -------------------------------------------------------------

DES <- list()
for (i in 1:n_iterations) {
  combined_matrix <- cbind(df_wide[[paste0("x1_iteration_", i)]], df_wide[[paste0("x2_iteration_", i)]])
  DES[[i]] <- combined_matrix
}

#TO DO:

# -------------------------------------------------------------
# PLOT INDVIDUAL CUMULATIVE HAZARDS
# -------------------------------------------------------------


plotLMPH_cumhaz <- function(individuals, DES, Betas, Thetas, LMs, dist = "W", ymax, lout) {
  time_grid <- seq(0, max(LMs), length.out = lout)
  tpoints <- as.vector(LMs)
  K <- length(tpoints)
  
  get_basehaz <- function(t, theta, dist) {
    if (dist == "W") {
      shape <- theta[1]
      scale <- theta[2]
      hz <- (shape / scale) * (t / scale)^(shape - 1)
      hz[t <= 0] <- 0
      return(hz)
    } else stop("Distribution not implemented.")
  }
  
  colors <- rainbow(length(individuals))
  
  plot(NULL, xlim = c(0, max(time_grid)), ylim = c(0, ymax),
       xlab = "Time", ylab = "Cumulative Hazard", main = "Individual Cumulative Hazard Functions")
  
  for (idx in seq_along(individuals)) {
    i <- individuals[idx]
    cumhaz <- numeric(length(time_grid))
    
    for (k in 1:(K - 1)) {
      t_start <- tpoints[k]
      t_end <- tpoints[k + 1]
      tk_grid <- time_grid[time_grid >= t_start & time_grid < t_end]
      if (length(tk_grid) < 2) next
      
      des_k <- DES[[k]]
      beta_k <- Betas[k, ]
      theta_k <- Thetas[k, ]
      
      linpred <- as.numeric(des_k[i, ] %*% beta_k)
      base_haz <- get_basehaz(tk_grid, theta_k, dist)
      haz_k <- base_haz * exp(linpred)
      
      dh <- diff(tk_grid)
      haz_mid <- (head(haz_k, -1) + tail(haz_k, -1)) / 2
      cumint <- cumsum(haz_mid * dh)
      
      start_idx <- which.min(abs(time_grid - tk_grid[1]))
      
      if (k == 1) {
        cumhaz[start_idx:(start_idx + length(cumint) - 1)] <- cumint
      } else {
        prev_end <- max(which(cumhaz != 0))
        cumhaz[(prev_end + 1):(prev_end + length(cumint))] <- cumint + cumhaz[prev_end]
      }
    }
    
    # detect the last non-zero point
    last_nonzero <- max(which(cumhaz != 0))
    
    if (last_nonzero > 1) {
      lines(time_grid[1:last_nonzero], cumhaz[1:last_nonzero], col = colors[idx], lwd = 2)
    }
  }
  
  legend("topright", legend = paste("Individual", individuals), col = colors,
         lty = 1, lwd = 2)
}

png("graficos/individual_cumulative_hazard.png", width = 800, height = 600)



plotLMPH_cumhaz(individuals = c(2,3,85), DES = DES, Betas = betas, Thetas = thetas, 
                LMs = landmarks, dist = "W", ymax = 0.2, lout = 1000)

dev.off()

# -------------------------------------------------------------
# PLOT POPULATION CUMULATIVE HAZARDS
# -------------------------------------------------------------



get_population_cumhaz_data <- function(DES, Betas, Thetas, LMs, dist = "W", lout = 1000) {
  
  # Build a time grid (avoid t = 0 to prevent Inf hazards)
  time_grid <- seq(1e-6, max(LMs), length.out = lout)
  tpoints <- as.vector(LMs)
  K <- length(tpoints) - 1
  
  # Hazard functions
  get_basehaz <- function(t, theta, dist) {
    if (dist == "W") {
      return((theta[1] / theta[2]) * (t / theta[2])^(theta[1] - 1))
    } else if (dist == "LN") {
      stop("Lognormal not implemented for pop hazard here")
    } else if (dist == "G") {
      stop("Gamma not implemented for pop hazard here")
    } else {
      stop("Unknown distribution")
    }
  }
  
  # Initialize outputs
  hazard <- numeric(length(time_grid))
  cumhaz <- numeric(length(time_grid))
  
  # Loop intervals
  cumhaz_so_far <- 0
  for (k in 1:K) {
    t_start <- tpoints[k]
    t_end   <- tpoints[k + 1]
    tk_idx  <- which(time_grid >= t_start & time_grid < t_end)
    if (length(tk_idx) == 0) next
    
    des_k   <- DES[[k]]
    beta_k  <- Betas[k, ]
    theta_k <- Thetas[k, ]
    
    mean_linpred <- mean(des_k %*% beta_k)
    haz_k <- get_basehaz(time_grid[tk_idx], theta_k, dist) * exp(mean_linpred)
    
    hazard[tk_idx] <- haz_k
    cumhaz[tk_idx] <- cumhaz_so_far + cumsum(haz_k * (max(time_grid) - min(time_grid)) / lout)
    
    cumhaz_so_far <- tail(cumhaz[tk_idx], 1)
  }
  
  survival <- exp(-cumhaz)
  
  data.frame(
    time = time_grid,
    hazard = hazard,
    cumulative_hazard = cumhaz,
    survival = survival
  )
}


# Compute data
cumhaz_data <- get_population_cumhaz_data(
  DES = DES,
  Betas = betas,
  Thetas = thetas,
  LMs = landmarks,
  dist = "W",
  lout = 1000
)

head(cumhaz_data, 10)

# Plot cumulative hazard
png("graficos/population_cumulative_hazard.png", width = 800, height = 600)

valid_idx <- which(cumhaz_data$cumulative_hazard[-nrow(cumhaz_data)] != 0)

if (length(valid_idx) > 0) {
  plot(cumhaz_data$time[-nrow(cumhaz_data)][valid_idx],
       cumhaz_data$cumulative_hazard[-nrow(cumhaz_data)][valid_idx],
       type = "l", col = "blue", lwd = 2,
       xlab = "Time", ylab = "Cumulative Hazard",
       main = "Population Cumulative Hazard")
  grid()
} else {
  message("No hay valores de cumulative_hazard distintos de 0 para graficar.")
}

dev.off()



