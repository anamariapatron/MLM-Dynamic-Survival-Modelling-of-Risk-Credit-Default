# ===============================
#  Modelling 
# ===============================
library(Landmarking)
library(dplyr)
library(survival)
library(ggplot2)
library(eha)
library(flexsurv)
library(GHSurv)

n <- nrow(time_matrix)
m <- ncol(time_matrix)


# ===============================
#  Approach 1 LANDMARKING
# ===============================
df_long <- data.frame()
for (i in 1:m) {
  df_long <- rbind(df_long,
                   data.frame(
                     id = 1:n,
                     iteration = i,
                     time = time_matrix[, i],
                     status = status_matrix[, i],
                     x1 = x1_matrix[, i],
                     x2 = x2_matrix[, i],
                     x3 = x3_matrix[, i],
                     x4 = x4_matrix[, i],
                     obs_date = obsdate_matrix[j, i]
                   ))
}



# Tus landmarks
#LM_start_date <- as.Date("2015-01-01")
#LM_dates <- as.Date(c("2020-01-01", "2020-03-01", "2020-05-01", "2021-07-01", "2021-09-01")) #a los 60 meses, 80 meses, 82 meses, etc
#landmarks <- round(interval(LM_start_date, LM_dates) %/% months(1))
landmarks <-  c(60, 62, 64, 66, 68, 70, 72)
df_long$obs_date <- as.Date(df_long$obs_date)
df_long$obs_time <- as.numeric(difftime(df_long$obs_date,initial_dates_matrix,  units = "days")) / 30.44
df_long$month_final <- landmarks[df_long$iteration]

# data_model_landmark_LOCF <- fit_LOCF_landmark(
#   data_long = df_long,
#   x_L = c(60, 62, 64, 66, 68, 70, 72),   # tiempos de landmark
#   x_hor = c(60, 62, 64, 66, 68, 70, 72)+12, # horizontes de predicción
#   covariates = c("x1", "x2"),
#   covariates_time = c(NULL, NULL),
#   k = 10,  # número de puntos para aproximar la integral de la función de riesgo
#   individual_id = "id",
#   event_time = "time",
#   event_status = "status",
#   survival_submodel = "cause_specific"  # modelo de supervivencia causa específica
# )
# t_vals <- seq(0, max(df_aft$time), length.out = 100)
# newdata <- df_aft %>% summarise(iteration = mean(iteration), month_final = mean(month_final))
# H_vals <- cumulative_hazard(aft_model, newdata, t_vals)
# 
# df_plot <- data.frame(time = t_vals, cumulative_hazard = H_vals)
# 
# landmarks <-  c(60, 62, 64, 66, 68, 70, 72)
# # Landmark times, etiquetas y colores
# labels <- c("t1", "t2", "t3", "t4", "t5", "t6", "t7")
# colors <- c("#FF9999", "#FFCC33", "#2EB67D","#57C9BE", "#7FDBFF", "#9999FF", "#FF99B8") 
# # Crear data.frame de landmarks
# df_lm <- data.frame(
#   time = landmarks,
#   label = labels,
#   color = colors
# )


# ===============================
#  Approach 2 AFT
# ===============================
library(dplyr)


# dataframe
df_aft <- df_long %>%
  group_by(id) %>%
  summarise(
    event_index = which(status == 1)[1],  # NA si no hay evento
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    status = ifelse(!is.na(event_index), 1, 0),
    time = ifelse(!is.na(event_index),
                  df_long$obs_time[df_long$id == id][event_index],
                  72),
    iteration = ifelse(!is.na(event_index),
                       df_long$iteration[df_long$id == id][event_index],
                       7),
    x1 = ifelse(!is.na(event_index),
                df_long$x1[df_long$id == id][event_index],
                df_long$x1[df_long$id == id & df_long$iteration == 7]),
    
    x2 = ifelse(!is.na(event_index),
                df_long$x2[df_long$id == id][event_index],
                df_long$x2[df_long$id == id & df_long$iteration == 7]),
    x3 = ifelse(!is.na(event_index),
                df_long$x3[df_long$id == id][event_index],
                df_long$x3[df_long$id == id & df_long$iteration == 7]),
    x4 = ifelse(!is.na(event_index),
                df_long$x4[df_long$id == id][event_index],
                df_long$x4[df_long$id == id & df_long$iteration == 7]),
    month_final = ifelse(!is.na(event_index),
                         df_long$month_final[df_long$id == id][event_index],
                         72)
  ) %>%
  ungroup()






# ===============================
#  Approach 3 general hazards
# ===============================

# =====================================
# Modellation: PH and AFT
# =====================================



# Cox Proportional Hazards
cox_fit <- coxph(Surv(time, status) ~ x1 + x2+x3+x4, data = df_aft)

# Weibull PH (parametric proportional hazards)
fit_weibull_ph <- flexsurvreg(
  Surv(time, status) ~ x1 + x2+x3+x4,
  data = df_aft,
  dist = "weibullPH"
)

# Weibull AFT
fit_weibull_aft <- flexsurvreg(
  Surv(time, status) ~ x1 + x2+x3+x4,
  data = df_aft,
  dist = "weibull"
)

 
# 3. Compare cumulative hazard functions
# Time grid

times <- seq(50, 72, length.out = 200)

# --- Weibull PH cumulative hazard
plot(fit_weibull_ph, type = "cumhaz", col = "blue", lwd = 2,
     xlab = "Time", ylab = "Cumulative Hazard",
     main = "Cumulative Hazard Function - Static Modelling",
     xlim = c(50, 72))

# --- Weibull AFT cumulative hazard
plot(fit_weibull_aft, type = "cumhaz", col = "green", lwd = 2,
     add = TRUE)

# --- Add legend
legend("topleft", legend = c("KM", "PH", "AFT"),
       col = c("black", "blue", "green"), lwd = 2, bty = "n")


# Landmark lines
landmarks <- c(60, 62, 64, 66, 68, 70, 72)
labels <- c("t1", "t2", "t3", "t4", "t5", "t6", "t7")
colors <- c("#FF9999", "#FFCC33", "#2EB67D","#57C9BE", "#7FDBFF", "#9999FF", "#FF99B8")

# for(i in seq_along(landmarks)) {
#   abline(v = landmarks[i], col = colors[i], lwd = 2, lty = 1)  # dashed vertical line
#   text(x = landmarks[i], y = par("usr")[4], labels[i], pos = 3, col = colors[i], cex = 0.8)
# }


