# ===============================
#  SIMULATION
# ===============================

rm(list = ls())
library(lubridate)
set.seed(1234)

# --- Parameters
n <- 1000
n_iterations <- 5
base_date <- as.Date("2020-01-01")

landmarks <- seq(0, n_iterations, by = 1)  
K <- length(landmarks) - 1  

betas <- matrix(c(-0.1, 0.15,
                  -0.1, 0.15,
                  -0.1, 0.15,
                  -0.1, 0.15,
                  -0.1, 0.15), nrow = K, byrow = TRUE)

thetas <- matrix(c(0.5, 400,
                   0.5, 280,
                   0.5, 360,
                   0.5, 400,
                   0.5, 450), nrow = K, byrow = TRUE)

# --- Data
initial_date <- as_date("2010-01-01") + days(sample(0:364, n, replace = TRUE))

x1_matrix <- matrix(rnorm(n * n_iterations), nrow = n, ncol = n_iterations)
x2_matrix <- matrix(rnorm(n * n_iterations), nrow = n, ncol = n_iterations)

time_matrix   <- matrix(NA_real_, nrow = n, ncol = n_iterations)
status_matrix <- matrix(NA_integer_, nrow = n, ncol = n_iterations)
obsdate_matrix <- matrix(as.Date(NA), nrow = n, ncol = n_iterations)

# --- Simulation function (Landmark PH Weibull)
simLMPH <- function(seed, x, Betas, Thetas, LMs, dist) {
  set.seed(seed)
  K <- length(LMs) - 1
  n <- nrow(x)
  tpoints <- as.vector(LMs)
  out_time <- rep(NA, n)
  
  for (k in 1:K) {
    beta_k <- Betas[k, ]
    theta_k <- Thetas[k, ]
    linpred <- x %*% beta_k
    mlambda <- exp(-linpred)
    
    if (dist == "W") {
      S0f <- function(t) pweibull(t, shape = theta_k[1], scale = theta_k[2], lower.tail = FALSE)
      quantf <- function(p) qweibull(p, shape = theta_k[1], scale = theta_k[2])
    } else {
      stop("Only Weibull implemented")
    }
    
    for (i in 1:n) {
      S_L <- S0f(tpoints[k])  
      u <- runif(1, min = 0, max = S_L)
      p_target <- max(u^(1 / mlambda[i]), .Machine$double.eps)
      t_rel <- quantf(p_target)
      
      if (is.na(out_time[i]) && t_rel <= tpoints[k + 1]) {
        out_time[i] <- t_rel
      }
    }
  }
  
  # Final extrapolation for those without event
  no_event <- which(is.na(out_time))
  if (length(no_event) > 0) {
    beta_k <- Betas[K, ]
    theta_k <- Thetas[K, ]
    linpred <- x[no_event, , drop = FALSE] %*% beta_k
    mlambda <- exp(-linpred)
    
    S0f <- function(t) pweibull(t, shape = theta_k[1], scale = theta_k[2], lower.tail = FALSE)
    quantf <- function(p) qweibull(p, shape = theta_k[1], scale = theta_k[2])
    
    for (j in seq_along(no_event)) {
      i <- no_event[j]
      u <- runif(1)
      p_target <- max(u^(1 / mlambda[j]), .Machine$double.eps)
      t_rel <- quantf(p_target)
      out_time[i] <- t_rel
    }
  }
  
  return(out_time)
}

# --- Main loop (everyone gets a time at each iteration)
alive <- rep(TRUE, n)

for (i in seq_len(n_iterations)) {
  cat("Iteration", i, "\n")
  
  covariates_i <- cbind(x1_matrix[, i], x2_matrix[, i])
  
  # Simulate event times for ALL individuals
  sim_times <- simLMPH(
    seed = 1234 + i,  # different seed per iteration
    x = covariates_i,
    Betas = betas,
    Thetas = thetas,
    LMs = landmarks,
    dist = "W"
  )
  
  # Save all times (even if already dead)
  time_matrix[, i] <- sim_times
  
  # Observation dates only make sense for those alive
  month_start <- base_date %m+% months((i - 1) * 12)
  month_end <- month_start %m+% months(1) - days(1)
  obs_dates <- as.Date(runif(n,
                             min = as.numeric(month_start),
                             max = as.numeric(month_end)),
                       origin = "1970-01-01")
  obsdate_matrix[, i] <- obs_dates
  
  # Status: event if elapsed time >= simulated time, else censored
  months_elapsed <- interval(initial_date, obs_dates) %/% months(1)
  status_vals <- as.integer(months_elapsed >= sim_times)
  
  # If already dead in a previous iteration → keep status = NA
  status_vals[!alive] <- NA
  
  status_matrix[, i] <- status_vals
  
  # Update alive vector
  alive <- ifelse(is.na(status_vals), FALSE, status_vals == 0)
}

# --- Summary
for (i in 1:n_iterations) {
  status_col <- status_matrix[, i]
  cat("Iteration", i, ":",
      sum(status_col == 1, na.rm = TRUE), "events\n")
}






