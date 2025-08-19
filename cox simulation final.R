# set directory
getwd()
#setwd("~/Documents/GitHub/MLM-On-The-Applications-Of-Survival-Modelling-In-The-Mortgage-Sector")

# ---- Libraries ----
rm(list = ls())
library(HazReg)
library(dplyr)
library(lubridate)

set.seed(1234)

# ---- Parameters ----
n <- 1000
n_iterations <- 5
landmarks <- seq(0, n_iterations * 1, by = 1)  # e.g., 0,1,2,3,4,5

x <- as.matrix(cbind(rnorm(n), rnorm(n)))
theta0 <- c(0.1, 2)
betah0 <- c(-0.5, 0.75) # e.g. income, interest rate
base_date <- as.Date("2020-01-01")

# Betas and thetas constant across iterations (for now)
betas <- matrix(rep(betah0, n_iterations), nrow = n_iterations, byrow = TRUE)
thetas <- matrix(rep(theta0, n_iterations), nrow = n_iterations, byrow = TRUE)

# ---- Initialize wide dataset ----
df_wide <- data.frame(individual = 1:n)
df_wide$initial_date <- as_date("2010-01-01") + days(sample(0:364, n, replace = TRUE))

alive <- rep(TRUE, n)

# Pre-create result columns
for (i in 1:n_iterations) {
  df_wide[[paste0("event_time_iteration_", i)]] <- NA_real_
  df_wide[[paste0("observation_date_iteration_", i)]] <- as.Date(NA)
  df_wide[[paste0("status_iteration_", i)]] <- NA_integer_
}

# ---- Simulation loop ----
for (i in seq_len(n_iterations)) {
  survivors <- which(alive)
  if (length(survivors) == 0) break
  
  cat("Iteration", i, "- survivors:", length(survivors), "\n")
  
  # Simulate event times
  sim_times <- simGH(
    seed = 1234,
    n = length(survivors),
    des = x[survivors, , drop = FALSE],
    theta = theta0,
    beta = betah0,
    hstr = "PH",
    baseline = "W"
  )
  
  # Observation window
  month_start <- base_date %m+% months((i - 1) * 12)
  month_end <- month_start %m+% months(1) - days(1)
  obs_dates <- as.Date(runif(length(survivors),
                             min = as.numeric(month_start),
                             max = as.numeric(month_end)),
                       origin = "1970-01-01")
  
  months_elapsed <- interval(df_wide$initial_date[survivors], obs_dates) %/% months(1)
  status_vals <- as.integer(months_elapsed >= sim_times)
  
  # Save
  df_wide[[paste0("event_time_iteration_", i)]][survivors] <- sim_times
  df_wide[[paste0("observation_date_iteration_", i)]][survivors] <- obs_dates
  df_wide[[paste0("status_iteration_", i)]][survivors] <- status_vals
  
  alive[survivors] <- status_vals == 0
}

# ---- Summary of survivors ----
survivors_summary <- data.frame(
  iteration = 1:n_iterations,
  survivors = sapply(1:n_iterations, function(i) {
    sum(df_wide[[paste0("status_iteration_", i)]] == 0, na.rm = TRUE)
  })
)

print(survivors_summary)
write.csv(survivors_summary, file = "survivors_summary.csv", row.names = FALSE)
write.csv(df_wide, file = "df_wide.csv", row.names = FALSE)

# ---- Hazard Functions ----
get_basehaz <- function(t, theta, dist="W") {
  if (dist == "W") {
    shape <- theta[1]; scale <- theta[2]
    hz <- (shape/scale) * (t/scale)^(shape - 1)
    hz[t <= 0] <- 0
    return(hz)
  } else stop("Distribution not implemented.")
}

getPopulationCumHaz <- function(DES, Betas, Thetas, LMs, dist="W", lout=1000) {
  time_grid <- seq(0, max(LMs), length.out=lout)
  K <- length(LMs)-1
  
  cumhaz <- numeric(length(time_grid))
  cumhaz_so_far <- 0
  
  for (k in 1:K) {
    t_start <- LMs[k]; t_end <- LMs[k+1]
    tk_idx <- which(time_grid >= t_start & time_grid < t_end)
    if (length(tk_idx) == 0) next
    
    des_k <- DES[[k]]; beta_k <- Betas[k,]; theta_k <- Thetas[k,]
    mean_linpred <- mean(des_k %*% beta_k)
    
    haz_k <- get_basehaz(time_grid[tk_idx] - t_start, theta_k, dist) * exp(mean_linpred)
    dt <- c(diff(time_grid[tk_idx]), 0)
    cumhaz[tk_idx] <- cumhaz_so_far + cumsum(haz_k * dt)
    cumhaz_so_far <- tail(cumhaz[tk_idx], 1)
  }
  
  out <- data.frame(
    time = time_grid,
    cumulative_hazard = cumhaz,
    survival = exp(-cumhaz)
  )
  # drop flat zero part
  out <- out[out$cumulative_hazard > 0, ]
  return(out)
}

plotCoxCumHazIndividual <- function(individuals, DES, Betas, Thetas, LMs, dist="W", lout=1000, ymax=NULL) {
  time_grid <- seq(0, max(LMs), length.out=lout)
  K <- length(LMs) - 1
  colors <- rainbow(length(individuals))
  
  plot(NULL, xlim=c(0, max(time_grid)), ylim=c(0, ifelse(is.null(ymax), 3, ymax)),
       xlab="Time", ylab="Cumulative Hazard",
       main="Individual Cumulative Hazard")
  
  for (idx in seq_along(individuals)) {
    i <- individuals[idx]
    cumhaz <- numeric(length(time_grid))
    cumhaz_so_far <- 0
    
    for (k in 1:K) {
      t_start <- LMs[k]; t_end <- LMs[k+1]
      tk_idx <- which(time_grid >= t_start & time_grid < t_end)
      if (length(tk_idx) == 0) next
      
      des_k <- DES[[k]]; beta_k <- Betas[k,]; theta_k <- Thetas[k,]
      linpred <- as.numeric(des_k[i,] %*% beta_k)
      
      haz_k <- get_basehaz(time_grid[tk_idx] - t_start, theta_k, dist) * exp(linpred)
      dt <- c(diff(time_grid[tk_idx]), 0)
      cumhaz[tk_idx] <- cumhaz_so_far + cumsum(haz_k * dt)
      cumhaz_so_far <- tail(cumhaz[tk_idx], 1)
    }
    
    valid_idx <- which(cumhaz > 0)
    lines(time_grid[valid_idx], cumhaz[valid_idx], col=colors[idx], lwd=2)
  }
  legend("topleft", legend=paste("Ind", individuals), col=colors, lty=1, lwd=2)
}

# ---- Build design list ----
DES <- list()
for (i in 1:n_iterations) {
  DES[[i]] <- x   # covariates same across iterations
}

# ---- Population Plot ----
pop_cumhaz <- getPopulationCumHaz(DES = DES, Betas = betas, Thetas = thetas, LMs = landmarks, dist="W")

png("graficos/population_cumhazCox.png", width=800, height=600)
plot(pop_cumhaz$time, pop_cumhaz$cumulative_hazard, type='l', lwd=2, col="blue",
     xlab="Time", ylab="Cumulative Hazard", main="Population Cumulative Hazard")
dev.off()

# ---- Individual Plot ----
png("graficos/individual_cumhazCox.png", width=800, height=600)
plotCoxCumHazIndividual(
  individuals = c(2, 3, 85),
  DES = DES,
  Betas = betas,
  Thetas = thetas,
  LMs = landmarks,
  dist="W"
)
dev.off()

