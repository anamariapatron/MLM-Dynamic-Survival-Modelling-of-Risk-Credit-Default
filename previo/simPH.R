rm(list=ls())

# Required packages
#library(devtools)
#install_github("FJRubio67/HazReg")
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

