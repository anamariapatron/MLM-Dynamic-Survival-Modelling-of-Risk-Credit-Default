

# set directory
getwd()
setwd("~/Documents/GitHub/MLM-On-The-Applications-Of-Survival-Modelling-In-The-Mortgage-Sector")

#libraries
rm(list = ls())
library(HazReg)
library(dplyr)
library(lubridate)

set.seed(1234)

# Parameters
n <- 1000
n_iterations <- 5

x <- as.matrix(cbind(rnorm(n), rnorm(n)))
theta0 <- c(0.1, 2)
betah0 <- c(-0.5, 0.75)
base_date <- as.Date("2020-01-01")

# Initialize main data frame
df_wide <- data.frame(individual = 1:n)
df_wide$initial_date <- as_date("2010-01-01") + days(sample(0:364, n, replace = TRUE))

# Track survival status
alive <- rep(TRUE, n)

# Pre-create result columns for each iteration
for (i in 1:n_iterations) {
  df_wide[[paste0("event_time_iteration_", i)]] <- NA_real_
  df_wide[[paste0("observation_date_iteration_", i)]] <- as.Date(NA)
  df_wide[[paste0("status_iteration_", i)]] <- NA_integer_
}

# Iterate over each time point
for (i in seq_len(n_iterations)) {
  survivors <- which(alive)
  if (length(survivors) == 0) break
  
  cat("Iteration", i, "- survivors:", length(survivors), "\n")
  
  # Simulate event times for survivors with unique seed per iteration
  sim_times <- simGH(
    seed = 1000 + i * 10,  # different seed each iteration
    n = length(survivors),
    des = x[survivors, , drop = FALSE],
    theta = theta0,
    beta = betah0,
    hstr = "PH",
    baseline = "W"
  )
  
  # Generate one random observation date for each survivor within the iteration window
  month_start <- base_date %m+% months((i - 1) * 12)
  month_end <- month_start %m+% months(1) - days(1)
  obs_dates <- as.Date(runif(length(survivors),
                             min = as.numeric(month_start),
                             max = as.numeric(month_end)),
                       origin = "1970-01-01")
  
  # Calculate months elapsed between initial date and observation date
  months_elapsed <- interval(df_wide$initial_date[survivors], obs_dates) %/% months(1)
  
  # Determine status: 1 if event occurred, 0 if still at risk
  status_vals <- as.integer(months_elapsed >= sim_times)
  
  # Store outputs in the wide dataframe
  df_wide[[paste0("event_time_iteration_", i)]][survivors] <- sim_times
  df_wide[[paste0("observation_date_iteration_", i)]][survivors] <- obs_dates
  df_wide[[paste0("status_iteration_", i)]][survivors] <- status_vals
  
  # Update alive vector
  alive[survivors] <- status_vals == 0
}

# Print summary of survivors after each iteration
for (i in 1:n_iterations) {
  count <- sum(df_wide[[paste0("status_iteration_", i)]] == 0, na.rm = TRUE)
  cat(sprintf("Survivors after iteration %d: %d\n", i, count))
}

write.csv(df_wide, file = "df_wide.csv", row.names = FALSE)

