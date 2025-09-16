# ===============================
# Libraries
# ===============================

library(lubridate)
library(truncnorm)
library(ggplot2)
set.seed(1234)

# ===============================
# Parameters
# ===============================

# Function to compute lambda given k, an interval, and the target probability p
k <- 1.5   # shape parameter (can tune)
p_target <- 0.2  # 20% default per interval


LM_start_date <- as.Date("2015-01-01")
LMs <-  c(60, 62, 64, 66, 68, 70, 72)

find_lambda <- function(t0, t1, k, p_target = 0.2) {
  f <- function(lambda) {
    ratio <- exp(-((t1 / lambda)^k - (t0 / lambda)^k))
    ratio - (1 - p_target)
  }
  uniroot(f, interval = c(1e-6, 1e6))$root
}

lambdas <- sapply(1:(length(LMs) - 1), function(j) {
  find_lambda(LMs[j], LMs[j + 1], k, p_target)
})



# ===============================
#  SIMULATION
# ===============================

rm(list = ls())
# Example: n individuals with different credit start dates
n <- 1000
n_iterations <- 7
base_date <- as.Date("2020-01-01")  
credit_start_dates <- as_date("2015-01-01") + days(sample(0:364, n, replace = TRUE))

# Function to generate individual landmarks
generate_landmarks <- function(credit_start) {
  lm1 <- credit_start + years(5)           # first landmark: +5 years
  lm2 <- lm1 + weeks(6)                    # second: +6 weeks
  lm3 <- lm2 + weeks(6)                    # third
  lm4 <- lm3 + weeks(6)                    # fourth
  lm5 <- lm4 + weeks(6)                    # fifth
  lm6 <- lm4 + weeks(6)
  lm7 <- lm4 + weeks(6)
  return(as_date(c(lm1, lm2, lm3, lm4, lm5, lm6, lm7)))
}
# Generate the list of individual landmarks
landmarks_list <- lapply(credit_start_dates, generate_landmarks)


landmarks <- seq(0, n_iterations, by = 1) #assumption cohort  
K <- length(landmarks) - 1  


betas <- matrix(c(
  50,  0.9,  0.8, -0.2,
  51, 0.92, 0.81, -0.21,
  49, 0.88, 0.79, -0.19,
  52, 0.95, 0.82, -0.22,
  48, 0.89, 0.78, -0.18,
  5,  0.9,  0.8, -0.2,
  51, 0.91, 0.81, -0.21
), nrow = 7, byrow = TRUE)
# Interest rate (+), # Inflation (+),# LTI (+), # Age (-)

#parameters of weibull each landmark: Shape (k), Scale (λ)
#thetas <- matrix(c(0.5, 400, 0.5, 280, 0.5, 360, 0.5, 400, 0.5, 450), nrow = K, byrow = TRUE)
thetas <- matrix(c(0.5, 22.25, 0.5, 22.49, 0.5, 22.49756 , 0.5, 22.73315 , 0.5, 22.96397,0.5, 23.19023,0.5, 23.41217), nrow = K, byrow = TRUE)


# --- Data
# credits are originated during 2010, throughout the year
#Covariates

# 1. Interest rate (UK Bank Rate)
# Truncated normal between 0.001 and 0.06, mean 0.03, sd 0.01
x1_matrix <- matrix(rtruncnorm(n * n_iterations, a = 0.001, b = 0.06, mean = 0.03, sd = 0.01),
                    nrow = n, ncol = n_iterations)

# 2. Inflation (CPI)
# Truncated normal between 0 and 0.10, mean 0.02, sd 0.015
x2_matrix <- matrix(rtruncnorm(n * n_iterations, a = 0, b = 0.10, mean = 0.02, sd = 0.015),
                    nrow = n, ncol = n_iterations)

# 3. Loan-to-income (LTI)
# Truncated normal between 1 and 6, mean 3.5, sd 0.5
x3_matrix <- matrix(rtruncnorm(n * n_iterations, a = 1, b = 6, mean = 3.5, sd = 2),
                    nrow = n, ncol = n_iterations)

# 4. Age of borrower
# Truncated normal between 20 and 70, mean 40, sd 10
x4_matrix <- matrix(rtruncnorm(n * n_iterations, a = 20, b = 70, mean = 40, sd = 10),
                    nrow = n, ncol = n_iterations)



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
    mlambda <- exp(linpred)
    
    if (dist == "W") {
      S0f <- function(t) pweibull(t, shape = theta_k[1], scale = theta_k[2], lower.tail = FALSE)
      quantf <- function(p) qweibull(p, shape = theta_k[1], scale = theta_k[2])
    } else {
      stop("Only Weibull implemented")
    }
    
    for (i in 1:n) {
      S_L <- S0f(tpoints[k])  
      u <- runif(1)
      p_target <- max(1-u^(1 / mlambda[i]), .Machine$double.eps)
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
    mlambda <- exp(linpred)
    
    S0f <- function(t) pweibull(t, shape = theta_k[1], scale = theta_k[2], lower.tail = FALSE)
    quantf <- function(p) qweibull(p, shape = theta_k[1], scale = theta_k[2])
    
    for (j in seq_along(no_event)) {
      i <- no_event[j]
      u <- runif(1)
      p_target <- max(1-u^(1 / mlambda[j]), .Machine$double.eps)
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
  
  # Select the covariates for this iteration
  covariates_i <- cbind(x1_matrix[, i], x2_matrix[, i],x3_matrix[, i],x4_matrix[, i] )
  
  # Initialize vectors for simulated times and observation dates
  sim_times <- numeric(n)
  obs_dates <- as.Date(rep(NA, n))
  
  # --- Simulation per individual
  for (j in 1:n) {
    # Simulate event time using each individual's own landmarks
    sim_times[j] <- simLMPH(
      seed = 1234,               # different seed per iteration and individual
      x = matrix(covariates_i[j, ], nrow = 1),
      Betas = betas,
      Thetas = thetas,
      LMs = landmarks_list[[j]],         # individual landmarks
      dist = "W"
    )
    
    # Assign observation date based on the corresponding landmark
    # Optionally, add a small random noise (0-3 days) for variability
    obs_dates[j] <- landmarks_list[[j]][i] + sample(0:3, 1)
  }
  
  # Save results
  time_matrix[, i] <- sim_times
  obsdate_matrix[, i] <- obs_dates
  
  # Determine status: event if elapsed time >= simulated time, else censored
  months_elapsed <- interval(credit_start_dates, obs_dates) %/% months(1)
  status_vals <- as.integer(months_elapsed >= sim_times)
  
  # If the individual was already dead in previous iterations → set NA
  status_vals[!alive] <- NA
  status_matrix[, i] <- status_vals
  
  # Update the alive vector
  alive <- ifelse(is.na(status_vals), FALSE, status_vals == 0)
}
initial_dates_matrix <- credit_start_dates


# --- Summary
for (i in 1:n_iterations) {
  status_col <- status_matrix[, i]
  cat("Iteration", i, ":",
      sum(status_col == 1, na.rm = TRUE), "events\n")
}


# ===============================
#  plots: deaths
# ===============================


library(ggplot2)
LMs <-  c(60, 62, 64, 66, 68, 70, 72)
cum_deaths <- cumsum(colSums(status_matrix, na.rm = TRUE))
df_cum <- data.frame(time = LMs, cum_deaths = cum_deaths)

ggplot(df_cum, aes(x = time, y = cum_deaths)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "darkblue") +
  labs(title = "Cumulative Deaths Over Time", x = "landmarking - not date", y = "Cumulative Deaths") +
  theme_minimal()


# ===============================
#  plots: deaths times
# ===============================
#cumulative deaths
library(ggplot2)
library(dplyr)

# Convert matrices to long format
df_events <- data.frame(
  individual = rep(1:nrow(status_matrix), times = ncol(status_matrix)),
  time = as.vector(time_matrix),
  status = as.vector(status_matrix)
)

# Keep only death events
df_deaths <- subset(df_events, status == 1 & !is.na(time))

# Create data frame for shaded "quarters" (every 3 time units)
time_min <- floor(min(df_deaths$time, na.rm = TRUE))
time_max <- ceiling(max(df_deaths$time, na.rm = TRUE))
quarter_starts <- seq(time_min, time_max, by = 3)
quarters <- data.frame(
  start = quarter_starts,
  end = quarter_starts + 3
)

# Plot
ggplot() +
  # Shade each 3-time-unit "quarter"
  geom_rect(data = quarters, aes(xmin = start, xmax = end, ymin = 0, ymax = max(df_deaths$individual)),
            fill = "grey90", alpha = 0.3) +
  # Plot death points
  geom_point(data = df_deaths, aes(x = time, y = individual), color = "red", size = 2, alpha = 0.7) +
  # Labels
  labs(
    title = "Death times per individual",
    subtitle = "Each point represents the time an individual died",
    x = "Time elapsed since t0",
    y = "Individual"
  ) +
  # Theme adjustments
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

 
# ===============================
#  plots: deaths dates
# ===============================


# Convert matrices to long format
df_events <- data.frame(
  individual = rep(1:nrow(status_matrix), times = ncol(status_matrix)),
  date = as.Date(as.vector(obsdate_matrix)),
  status = as.vector(status_matrix)
)

# Keep only death events
df_deaths <- subset(df_events, status == 1 & !is.na(date))

# Create data frame for shaded quarters
quarter_starts <- seq(floor_date(min(df_deaths$date), unit = "quarter"),
                      ceiling_date(max(df_deaths$date), unit = "quarter"),
                      by = "3 months")
quarters <- data.frame(
  start = quarter_starts,
  end = quarter_starts + months(3)
)

# Plot
ggplot() +
  geom_rect(data = quarters, aes(xmin = start, xmax = end, ymin = 0, ymax = max(df_deaths$individual)),
            fill = "grey90", alpha = 0.3) +
  geom_point(data = df_deaths, aes(x = date, y = individual), color = "red", size = 2, alpha = 0.7) +
  geom_vline(xintercept = as.Date("2020-12-01"), color = "black", linetype = "solid", size = 1) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  labs(
    title = "Death Dates per Individual",
    subtitle = "Each point represents the date an individual died",
    x = "Quarter start",
    y = "Individual"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


####### points for modelling in where model is going t be re-calibrated #######
mean_dates <- apply(obsdate_matrix, 2, function(col) {
  as.Date(mean(as.numeric(col)), origin = "1970-01-01")
})
mean_dates <- as.Date(mean_dates)

LM_start_date <- as.Date("2015-01-01")
LMs <-  c(60, 62, 64, 66, 68, 70, 72)

labels_df <- data.frame(
  mean_dates = LM_start_date %m+% months(LMs),
  label = paste0("t", seq_along(mean_dates))
)

ggplot() +
  geom_rect(data = quarters, aes(xmin = as.Date(start), xmax = as.Date(end), ymin = 0, ymax = max(df_deaths$individual)),
            fill = "grey90", alpha = 0.3) +
  geom_point(data = df_deaths, aes(x = as.Date(date), y = individual), color = "red", size = 2, alpha = 0.7) +
  geom_vline(data = labels_df, aes(xintercept = mean_dates, color = label), linetype = "solid", size = 1) +
  geom_text(data = labels_df, aes(x = mean_dates, 
                                  y = max(df_deaths$individual) , 
                                  label = label, color = label),
            vjust = 0.2, hjust = 1, size = 5, fontface = "bold") +
  scale_x_date(date_breaks = "3 month", date_labels = "%b %Y") +
  labs(
    title = "Default Dates per Individual",
    subtitle = "Each point shows the date an individual defaulted",
    x = "Month",
    y = "Individual"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  guides(color = "none")


 

