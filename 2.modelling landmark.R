library(ggplot2)
library(dplyr)
library(Landmarking)
set.seed(1)

### base de datos
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
                     obs_date= obsdate_matrix[j, i]
                   ))
}


landmarks <-  c(60, 62, 64)
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
    covariates = c("x1", "x2"),
    covariates_time = c(rep("obs_time", 2)),
    x_L = c(60, 62, 64)
  )

 
data_model_landmark_LOCF <-
  fit_LOCF_landmark(
    data_long = df_long2,
    x_L = c(60, 62, 64),
    x_hor = c(60, 62, 64)+6,
    covariates = c("x1", "x2"),
    covariates_time = c(rep("obs_time", 2)),
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


data_model_landmark_LOCF[["61"]]$data$event_prediction


df_long2[["60"]]$status
df_long2[["60"]]$time








# Elegimos la persona que nos interesa
person_id <- 2

# Recopilamos los riesgos por cada landmark
risks_person <- lapply(names(data_model_landmark_LOCF), function(LM) {
  df <- data_model_landmark_LOCF[[LM]]$data
  df_person <- df %>% filter(id == person_id)
  if(nrow(df_person) == 0) return(NULL)
  data.frame(
    LM = as.numeric(LM),
    risk = df_person$event_prediction
  )
}) %>% bind_rows()

# Graficar riesgo instantáneo vs landmark
ggplot(risks_person, aes(x = LM, y = risk)) +
  geom_point(size = 3, color = "red") +
  geom_line(color = "red", lwd = 1) +
  labs(
    x = "Landmark time",
    y = "Predicted risk (instantaneous)",
    title = paste("Individual risk for person", person_id)
  ) +
  theme_minimal()
