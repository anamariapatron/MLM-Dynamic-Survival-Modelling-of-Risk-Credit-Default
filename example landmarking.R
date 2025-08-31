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
    x_L = c(60, 61),
    x_hor = c(65, 66),
    covariates = c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd"),
    covariates_time = c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd"),
    k = 10,
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

