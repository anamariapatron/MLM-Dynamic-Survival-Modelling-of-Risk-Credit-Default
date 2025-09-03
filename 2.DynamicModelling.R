# ===============================
#  Modelling 
# ===============================
library(Landmarking)
library(dplyr)
library(survival)
library(ggplot2)
library(eha)
library(tidyr)

n <- nrow(time_matrix)
m <- ncol(time_matrix)

 
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


landmarks <-  c(54,56, 58, 60, 62, 64)
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
    x_L = c(54,56, 58, 60, 62, 64)
  )

 
data_model_landmark_LOCF <-
  fit_LOCF_landmark(
    data_long = df_long2,
    x_L = c(54,56, 58, 60, 62, 64),
    x_hor = c(54,56, 58, 60, 62, 64)+12,
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

# Supongamos que queremos la persona con id = 3
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




