library(Landmarking)
set.seed(1)
data(data_repeat_outcomes)
head(data_repeat_outcomes)


data_repeat_outcomes <-
  return_ids_with_LOCF(
    data_long = data_repeat_outcomes,
    individual_id = "id",
    covariates = c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd"),
    covariates_time = c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd"),
    x_L = c(60, 61)
  )


data_model_landmark_LOCF <-
  fit_LOCF_landmark(
    data_long = data_repeat_outcomes,
    x_L = c(60, 61, 62),
    x_hor = c( 60, 61,62)+1,
    covariates = c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd"),
    covariates_time = c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd"),
    k = 2,
    individual_id = "id",
    event_time = "event_time",
    event_status = "event_status",
    survival_submodel = "cause_specific"
  )

plot(
  density(100 * data_model_landmark_LOCF[["60"]]$data$event_prediction),
  xlab = "Predicted risk of CVD event (%)",
  main = "Landmark age 60"
)

data_model_landmark_LOCF[["60"]]$data$event_prediction
data_model_landmark_LOCF[["60"]]$data$event_time


data_model_landmark_LOCF[["61"]]$data$event_prediction


#######################################
library(ggplot2)
library(dplyr)
library(tidyr)

# Supongamos que queremos la persona con id = 2
person_id <- 2

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

#######################################



data_model_landmark_LME <-
  fit_LME_landmark(
    data_long = data_repeat_outcomes[["60"]],
    x_L = c(60),
    x_hor = c(65),
    cross_validation_df =
      cross_validation_list,
    fixed_effects = c("ethnicity", "smoking", "diabetes"),
    fixed_effects_time =
      "response_time_sbp_stnd",
    random_effects = c("sbp_stnd", "tchdl_stnd"),
    random_effects_time = c("response_time_sbp_stnd", "response_time_tchdl_stnd"),
    individual_id = "id",
    standardise_time = TRUE,
    lme_control = nlme::lmeControl(maxIter =
                                     100, msMaxIter = 100),
    event_time = "event_time",
    event_status = "event_status",
    survival_submodel = "cause_specific"
  )