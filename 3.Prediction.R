# ===============================
#  Prediction
# ===============================
library(Landmarking)
library(dplyr)
library(survival)
library(ggplot2)
library(eha)

landmarks <- c(54,56,58,60, 62, 64)
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
  
  d <- density(preds, na.rm = TRUE)
  
  plot(
    d,
    xlab = "Predicted risk of default event (%)",
    main = paste("Landmark month", a),
    lwd = 2,
    col = "deeppink3",    
    ylim = c(0, max(d$y) * 1.2)  
  )
  
  polygon(d, col = adjustcolor("pink", alpha.f = 0.5), border = NA)
  
  # Línea encima
  lines(d, col = "deeppink3", lwd = 2)
}


# ===============================
#  aft
# ===============================
par(mfrow = c(2, 3))

for (a in landmarks) {
  
  d <- density(risk_list_aft[[as.character(a)]], na.rm = TRUE)
  
  plot(
    d,
    xlab = "Predicted risk of default event (%)",
    main = paste("Landmark month", a),
    lwd = 2,
    col = "green3",    
    ylim = c(0, max(d$y) * 1.2)  
  )
  
  
  
  # Línea encima
  lines(d, col = "green3", lwd = 2)
}



# ===============================
#  plot density curves 
# ===============================

pdf("landmark_densities_horizontal.pdf", width = 12, height = 8)
# Ajustamos mfrow, márgenes internos (mar) y externos (oma) para subplots más grandes
par(mfrow = c(2, 3), 
    mar = c(4, 4, 2, 2),   # márgenes dentro de cada subplot
    oma = c(6, 0, 0, 0),   # margen externo inferior para etiqueta X y leyenda
    xpd = TRUE, 
    cex.main = 1.5, 
    cex.lab = 1.3, 
    cex.axis = 1.2)

for (a in landmarks) {
  
  preds <- 100 * data_model_landmark_LOCF[[as.character(a)]]$data$event_prediction
  preds <- preds[preds <= 4]
  
  d1 <- density(risk_list_aft[[as.character(a)]], na.rm = TRUE, from = 0, to = 4)
  d2 <- density(preds, na.rm = TRUE, from = 0, to = 4)
  
  ymax <- max(d1$y, d2$y, na.rm = TRUE) * 1.2
  
  plot(
    0, 0, type = "n",
    xlab = "",
    ylab = "Density",
    main = paste("Landmark month", a),
    xlim = c(0, 4),
    ylim = c(0, ymax)
  )
  
  polygon(c(d1$x, rev(d1$x)), c(rep(0, length(d1$y)), rev(d1$y)), 
          col = adjustcolor("green", alpha.f = 0.3), border = "green3")
  lines(d1, col = "green3", lwd = 2)
  
  polygon(c(d2$x, rev(d2$x)), c(rep(0, length(d2$y)), rev(d2$y)), 
          col = adjustcolor("pink", alpha.f = 0.3), border = "deeppink3")
  lines(d2, col = "deeppink3", lwd = 2)
}


mtext("Predicted risk of default event (%)", side = 1, outer = TRUE, line = 4, cex = 1.3)


par(xpd = NA)
legend(x = mean(par("usr")[1:2]) - 2,  # resta 0.2 para moverla a la izquierda
       y = -0.3 * max(par("usr")[3:4]), 
       legend = c("AFT", "Landmarking"),
       col = c("green3", "deeppink3"),
       lwd = 2,
       pt.bg = c(adjustcolor("green", alpha.f = 0.3), adjustcolor("pink", alpha.f = 0.3)),
       pch = 15,
       bty = "n",
       cex = 1.3,
       horiz = TRUE)



dev.off()


