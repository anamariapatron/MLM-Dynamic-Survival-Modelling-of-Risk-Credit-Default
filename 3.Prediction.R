# ===============================
#  Prediction
# ===============================
library(Landmarking)
library(dplyr)
library(survival)
library(ggplot2)
library(eha)

landmarks <- c(60, 62, 64, 66, 68, 70, 72)
# ===============================
#  aft and ph
# ===============================

risk_list_ph <- lapply(landmarks, function(t_star) {
  surv_pred <- summary(
    fit_weibull_ph,
    type = "survival",
    t = t_star,
    newdata = df_aft
  )
  sapply(surv_pred, function(x) 1 - x$est) * 100
})


names(risk_list_ph) <- landmarks

risk_list_aft <- lapply(landmarks, function(t_star) {
  surv_pred <- summary(
    fit_weibull_aft,
    type = "survival",
    t = t_star,
    newdata = df_aft
  )
  sapply(surv_pred, function(x) 1 - x$est) * 100
})

names(risk_list_aft) <- landmarks




# ===============================
#  landmarking
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
#  plot density curves
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
