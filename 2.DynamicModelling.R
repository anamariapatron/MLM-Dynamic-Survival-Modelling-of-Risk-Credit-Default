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

#quiat lo oreo
# ===============================
#  Approach 2 landmarking
# ===============================

df_long <- data.frame()
for (i in 1:m) {
  df_long <- rbind(df_long,
                   data.frame(
                     id = 1:n,
                     iteration = i,
                     time_raw = time_matrix[, i],
                     status_raw = status_matrix[, i],
                     x1 = x1_matrix[, i],
                     x2 = x2_matrix[, i],
                     x3 = x3_matrix[, i],
                     x4 = x4_matrix[, i],
                     obs_date= obsdate_matrix[j, i]
                   ))
}


landmarks <-  c(56, 58, 60, 62, 64)
df_long$obs_date <- as.Date(df_long$obs_date)
df_long$obs_time <- as.numeric(difftime(df_long$obs_date,initial_dates_matrix,  units = "days")) / 30.44
df_long$month_final <- landmarks[df_long$iteration]
df_long <- df_long %>%
  filter(!is.na(status_raw))

df_long <- df_long %>%
  group_by(id) %>%
  mutate(
    event_time = if (any(status_raw == 1)) {
      time_raw[status_raw == 1][1]   # fecha del primer evento
    } else {
      max(time_raw, na.rm = TRUE)    # si nunca hubo evento, el máximo
    }
  ) %>%
  ungroup()



df_long <- df_long %>%
  group_by(id) %>%
  mutate(event_status = if_else(any(status_raw == 1), 1, 0)) %>%
  ungroup()


df_long2 <-
  return_ids_with_LOCF(
    data_long = df_long,
    individual_id = "id",
    covariates = c("x1", "x2","x3","x4"),
    covariates_time = c(rep("obs_time", 4)),
    x_L = c(56, 58, 60, 62, 64)
  )

 
data_model_landmark_LOCF <-
  fit_LOCF_landmark(
    data_long = df_long2,
    x_L = c(56, 58, 60, 62, 64),
    x_hor = c(56, 58, 60, 62, 64)+12,
    covariates =  c("x1", "x2","x3","x4"),
    covariates_time = c(rep("obs_time", 4)),
    k = 1,
    individual_id = "id",
    event_time = "event_time",
    event_status = "event_status",
    survival_submodel = "standard_cox"
  )

plot(
  density(100 * data_model_landmark_LOCF[["60"]]$data$event_prediction),
  xlab = "Predicted risk of CVD event (%)",
  main = "Landmark age 60"
)

data_model_landmark_LOCF[["60"]]$data$event_prediction
data_model_landmark_LOCF[["60"]]$data$event_time



data_model_landmark_LOCF$model_survival



# ===============================
#  plots
# ===============================

# Elegimos la persona que nos interesa
library(ggplot2)
library(dplyr)
library(tidyr)

# Supongamos que queremos la persona con id = 2
person_id <- 3

# Seleccionamos los landmarks que tenemos
LMs <- names(data_model_landmark_LOCF)

# Crear dataframe con cumulative hazard para esta persona en cada landmark
cumhaz_person <- lapply(LMs, function(LM) {
  df_lm <- data_model_landmark_LOCF[[LM]]$data
  df_person <- df_lm %>% filter(id == person_id)
  if(nrow(df_person) == 0) return(NULL)
  
  data.frame(
    LM = as.numeric(LM),
    cumhaz = -log(1 - df_person$event_prediction)
  )
}) %>% bind_rows()

# Graficar
ggplot(cumhaz_person, aes(x = LM, y = cumhaz)) +
  geom_point(size = 3, color = "blue") +
  geom_line(color = "blue", lwd = 1) +
  labs(
    x = "Landmark time",
    y = "Cumulative hazard",
    title = paste("Cumulative hazard for person", person_id)
  ) +
  theme_minimal()



# ===============================
#  predictions
# ===============================

par(mfrow = c(2, 3))

for (a in landmarks) {
  preds <- 100 * data_model_landmark_LOCF[[as.character(a)]]$data$event_prediction
  
  plot(
    density(preds, na.rm = TRUE),
    xlab = "Predicted risk of CVD event (%)",
    main =  paste("Landmark age", a),
    lwd = 2,
    col = "red"
  )
}
# ===============================
#  dynamic cumulative hazard
# ===============================
 
pdf("landmark_densities.pdf", width = 12, height = 8)

 
par(mfrow = c(2,3), mar = c(4,4,2,1))  # ajustar márgenes

 
col_lm  <- rgb(1, 0.5, 0.5, 0.5)   
col_ph  <- rgb(0.5, 0.5, 1, 0.5)  
col_aft <- rgb(0.5, 1, 0.5, 0.5)   

for (a in landmarks) {
  preds <- 100 * data_model_landmark_LOCF[[as.character(a)]]$data$event_prediction
  
  dens_aft <- density(risk_list_aft[[as.character(a)]], na.rm = TRUE)
  dens_ph  <- density(risk_list_ph[[as.character(a)]], na.rm = TRUE)
  dens_base <- density(preds, na.rm = TRUE)
  
  
  plot(
    dens_base,
    xlab = "Predicted risk of default event (%)",
    main = paste("Estimation at point", a),
    lwd = 2,
    ylim = c(0, max(dens_base$y, dens_ph$y, dens_aft$y)),
    type = "n"
  )
  
  # rellenar densidades con colores pastel
  polygon(dens_base, col = col_lm, border = "red", lwd = 2)
  polygon(dens_ph, col = col_ph, border = "blue", lwd = 2, lty = 2)
  polygon(dens_aft, col = col_aft, border = "darkgreen", lwd = 2, lty = 3)
  
  # leyenda
  legend("topright",
         legend = c("LM", "PH", "AFT"),
         col = c("red", "blue", "darkgreen"),
         lwd = 2,
         lty = c(1,2,3),
         bty = "n")
}

dev.off()
