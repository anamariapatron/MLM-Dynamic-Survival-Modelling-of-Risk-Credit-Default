# ===============================
#  Modelling 
# ===============================
library(Landmarking)
library(dplyr)
library(survival)
library(ggplot2)
library(eha)

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
    month_final = ifelse(!is.na(event_index),
                         df_long$month_final[df_long$id == id][event_index],
                         72)
  ) %>%
  select(id, status, time, iteration, month_final, x1, x2) %>%
  ungroup()




# AFT model using a weibull baseline distribution #TO DO: weibull?
fit.aft <- aftreg(Surv(time, status) ~ x1 + x2,
                  data = df_aft,
                  dist = "weibull")

summary(fit.aft)

# PH model using a weibull baseline distribution
fit.ph <- phreg(Surv(time, status) ~ x1 + x2,
                data = df_aft,
                dist = "weibull")

summary(fit.ph)

# Comparison of the three models using AIC
AIC.aft <- -2*(fit.aft$loglik[2]) + 2*length(fit.aft$coefficients)
AIC.ph <- -2*(fit.ph$loglik[2]) + 2*length(fit.ph$coefficients)

AIC.aft

AIC.ph

plot(fit.aft)
# Filtrar el dataframe para el rango deseado
df_plot_trunc <- df_plot %>% 
  filter(time >= 60 & time <= 70)

df_lm_trunc <- df_lm %>% 
  filter(time >= 60 & time <= 70)

# Graficar solo en el rango 60-70
ggplot(df_plot_trunc, aes(x = time, y = cumulative_hazard)) +
  geom_line(color = "blue", size = 1.2) +
  geom_vline(data = df_lm_trunc, aes(xintercept = time, color = label),
             linetype = "dashed", show.legend = FALSE) +
  geom_text(data = df_lm_trunc,
            aes(x = time,
                y = max(df_plot_trunc$cumulative_hazard) * 0.95,
                label = label,
                color = label),
            angle = 90, vjust = -0.5, hjust = 1, show.legend = FALSE) +
  scale_color_manual(values = setNames(colors, labels)) +
  labs(
    title = "Cumulative Hazard from AFT Weibull Model (Times 60-70)",
    x = "Time",
    y = "Cumulative Hazard"
  ) +
  theme_minimal()







# Función para cumulative hazard
cumulative_hazard <- function(model, newdata, t) {
  # t = vector de tiempos
  lp <- predict(model, newdata = newdata, type = "lp")  # linear predictor
  k <- 1 / model$scale
  lambda <- exp(-lp * model$scale)
  H <- (t / lambda)^k
  return(H)
}


t_vals <- seq(0, max(df_aft$time), length.out = 100)
newdata <- df_aft %>% summarise(iteration = mean(iteration), month_final = mean(month_final))
H_vals <- cumulative_hazard(aft_model, newdata, t_vals)

df_plot <- data.frame(time = t_vals, cumulative_hazard = H_vals)

landmarks <-  c(60, 62, 64, 66, 68, 70, 72)
# Landmark times, etiquetas y colores
labels <- c("t1", "t2", "t3", "t4", "t5", "t6", "t7")
colors <- c("#FF9999", "#FFCC33", "#2EB67D","#57C9BE", "#7FDBFF", "#9999FF", "#FF99B8") 
# Crear data.frame de landmarks
df_lm <- data.frame(
  time = landmarks,
  label = labels,
  color = colors
)

# Graficar
ggplot(df_plot, aes(x = time, y = cumulative_hazard)) +
  geom_line(color = "blue", size = 1.2) +
  
  # Líneas verticales en los landmarks
  geom_vline(data = df_lm, aes(xintercept = time, color = label),
             linetype = "dashed", show.legend = FALSE) +
  
  # Etiquetas de landmarks
  geom_text(data = df_lm,
            aes(x = time,
                y = max(df_plot$cumulative_hazard) * 0.95,
                label = label,
                color = label),
            angle = 90, vjust = -0.5, hjust = 1, show.legend = FALSE) +
  
  scale_color_manual(values = setNames(colors, labels)) +
  
  labs(
    title = "Cumulative Hazard from AFT Weibull Model",
    x = "Time",
    y = "Cumulative Hazard"
  ) +
  theme_minimal()

# ===============================
#  Approach 3 general hazards
# ===============================


