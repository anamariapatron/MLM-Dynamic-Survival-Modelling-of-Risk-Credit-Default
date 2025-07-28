rm(list=ls())

# Required packages
library(HazReg)

# Sample size
n <- 1000

# Simulated design matrix
set.seed(1234)
x <- as.matrix(cbind(rnorm(n), rnorm(n)))

#----------------------------
# PGW-GH simulation
#----------------------------

# True parameters
theta0 <- c(0.1,2)
betah0 <- c(-0.5,0.75)

# Data simulation (first iteration, before first landmark)
simdat <- simGH(seed = 1234, n = n, des = x, 
                theta = theta0, beta = betah0, 
                hstr = "PH", baseline = "W")

# Crear un data frame con los tiempos de evento y el estado (suponiendo que todos los eventos ocurren)
simdat_df <- data.frame(time = simdat, status = rep(1, length(simdat)))  # Estado 1: evento ocurrido

# Función para realizar simulación y análisis de landmarking
run_landmarking_simulation <- function(landmark_times, n, x, theta0, betah0) {
  results <- list()
  
  for (t in landmark_times) {
    # Realizamos la simulación
    simdat <- simGH(seed = 1234, n = n, des = x, 
                    theta = theta0, beta = betah0, 
                    hstr = "PH", baseline = "W")
    
    # Crear un data frame con los tiempos de evento
    simdat_df <- data.frame(time = simdat, status = rep(1, length(simdat)))  # Estado 1 para todos los eventos
    
    # Filtrar los datos hasta el "landmark"
    simdat_landmark <- simdat_df[simdat_df$time <= t, ]
    
    
    # Guardamos los resultados
    results[[as.character(t)]] <- simdat_landmark
  }
  
  return(results)
}

# Definir los diferentes puntos de "landmark"
landmark_times <- c(5, 10, 15, 20)

# Ejecutar la simulación para diferentes landmarks
simulation_results <- run_landmarking_simulation(landmark_times, n, x, theta0, betah0)

# Ver resultados
simulation_results


