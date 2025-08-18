###############################################
####### Survival modelling for the case #######
####### of Mortgage sector in UK        #######
##############################################
library(survival)

set.seed(123)
n <- 500
beta <- 0.5
lambda <- 0.1
v <- 1.5
change_time <- 3

x0 <- rnorm(n)
x1 <- x0 + rnorm(n, mean = 0.2)

u <- runif(n)
T1 <- (-log(u) / (lambda * exp(beta * x0)))^(1 / v)
survived_change <- T1 > change_time
u2 <- runif(n)
T2 <- rep(Inf, n)
T2[survived_change] <- change_time + (-log(u2[survived_change]) / (lambda * exp(beta * x1[survived_change])))^(1 / v)

T_total <- pmin(T1, T2)
C <- rexp(n, rate = 0.05)
time <- pmin(T_total, C)
event <- as.numeric(T_total <= C)

df <- data.frame(id = 1:n, time = time, event = event, x0 = x0, x1 = x1)

landmarks <- seq(1, 6, by = 1)

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

for (t_star in landmarks) {
  df_lm <- df[df$time > t_star, ]
  if (nrow(df_lm) == 0) next
  
  df_lm$time_post <- df_lm$time - t_star
  df_lm$event_post <- df_lm$event
  
  # Covariable según el tiempo landmark
  covariate <- if (t_star < change_time) "x0" else "x1"
  df_lm$group <- ifelse(df_lm[[covariate]] > 0, "High", "Low")
  
  fit_km <- survfit(Surv(time_post, event_post) ~ group, data = df_lm)
  
  plot(fit_km, col = c("blue", "red"), lty = 1:2,
       xlab = paste("Time since landmark t* =", t_star),
       ylab = "Survival probability",
       main = paste("Landmark at t* =", t_star))
  legend("bottomleft", legend = c("Low", "High"), col = c("blue", "red"), lty = 1:2)
}

write.csv(df, "simulated_survival_data.csv", row.names = FALSE)