rm(list = ls())

library(HazReg)
library(dplyr)
library(lubridate)

set.seed(1234)

# Parameters
n <- 1000
n_months <- 12
n_iterations <- 5

x <- as.matrix(cbind(rnorm(n), rnorm(n)))
theta0 <- c(0.1, 2)
betah0 <- c(-0.5, 0.75)
base_date <- as.Date("2020-01-01")
landmarks <- seq(0, by = n_months, length.out = n_iterations)

# Initialize wide dataframe
df_wide <- data.frame(individual = 1:n)
df_wide$initial_date <- as_date("2010-01-01") + days(sample(0:364, n, replace = TRUE))

# Initialize vectors to track cumulative event times and survival status
event_time_total <- rep(0, n)    # cumulative event times per individual
alive <- rep(TRUE, n)            # who is alive (not had event yet)

total_months <- n_months * n_iterations

# Pre-create columns for observation dates, predicted event times, and status per month
for (m in 1:total_months) {
  df_wide[[paste0("observation_date_month_", m)]] <- as.Date(NA)
  df_wide[[paste0("predicted_event_time_month_", m)]] <- NA_real_
  df_wide[[paste0("status_month_", m)]] <- NA_integer_
}

# Create columns to hold cumulative predicted event time at each iteration
for (i in 1:n_iterations) {
  df_wide[[paste0("event_time_iteration_", i)]] <- NA_real_
}

# Month counter to track global month index across iterations
month_counter <- 1

for (i in seq_len(n_iterations)) {
  survivors <- which(alive)
  if (length(survivors) == 0) break   # stop if no survivors remain
  
  cat("Iteration", i, "- survivors:", length(survivors), "\n")
  
  # Simulate *incremental* event times ONLY for survivors
  # This simulates how long each survivor will survive *starting now* 
  sim_incremental <- simGH(seed = 1234 + i * 100,
                           n = length(survivors),
                           des = x[survivors, , drop = FALSE],
                           theta = theta0,
                           beta = betah0,
                           hstr = "PH", baseline = "W")
  
  # Update cumulative event times: add the new incremental times to previous totals
  event_time_total[survivors] <-  sim_incremental
  
  # Save cumulative predicted event times at current iteration for survivors
  df_wide[[paste0("event_time_iteration_", i)]][survivors] <- event_time_total[survivors]
  
  # For each month within the iteration, generate observation dates, predicted event times, and status
  for (m in 1:n_months) {
    global_month <- landmarks[i] + m      # global month index across all iterations
    
    month_start <- base_date %m+% months(global_month - 1)
    month_end <- month_start %m+% months(1) - days(1)
    
    # Generate random observation dates uniformly in the month for survivors
    obs_dates <- rep(as.Date(NA), n)
    obs_dates[survivors] <- as.Date(runif(length(survivors),
                                          min = as.numeric(month_start),
                                          max = as.numeric(month_end)),
                                    origin = "1970-01-01")
    
    # Predicted event times = cumulative event time for survivors, NA otherwise
    predicted_vals <- rep(NA_real_, n)
    predicted_vals[survivors] <- event_time_total[survivors]
    
    # Calculate months elapsed from initial_date to observation_date for survivors
    months_elapsed <- rep(NA_integer_, n)
    months_elapsed[survivors] <- interval(df_wide$initial_date[survivors], obs_dates[survivors]) %/% months(1)
    
    # Status = 1 if months elapsed >= predicted event time (event occurred), else 0
    status_vals <- rep(NA_integer_, n)
    status_vals[survivors] <- as.integer(months_elapsed[survivors] >= predicted_vals[survivors])
    
    # Save results into wide dataframe
    df_wide[[paste0("observation_date_month_", global_month)]] <- obs_dates
    df_wide[[paste0("predicted_event_time_month_", global_month)]] <- predicted_vals
    df_wide[[paste0("status_month_", global_month)]] <- status_vals
  }
  
  # Update alive vector: survivors are those with status=0 at last month of iteration
  last_month <- landmarks[i] + n_months
  alive <- df_wide[[paste0("status_month_", last_month)]] == 0
  alive[is.na(alive)] <- FALSE
}

# print number of survivors after each iteration for verification
for (i in 1:n_iterations) {
  last_month <- landmarks[i] + n_months
  alive_count <- sum(df_wide[[paste0("status_month_", last_month)]] == 0, na.rm = TRUE)
  cat(sprintf("Survivors after iteration %d: %d\n", i, alive_count))
}

