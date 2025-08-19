source("SimLM.R")

set.seed(1234)
des1 <- cbind(rnorm(10),rnorm(10))
des2 <- cbind(rnorm(10),rnorm(10))
des3 <- cbind(rnorm(10),rnorm(10))
des4 <- cbind(rnorm(10),rnorm(10))
des5 <- cbind(rnorm(10),rnorm(10))
des6 <- cbind(rnorm(10),rnorm(10))

DES = list(des1,des2,des3,des4,des5,des6)

Betas = rbind(c(0.1,-0.2),c(0.1,-0.2),c(0.1,-0.2),c(0.15,-0.25),c(0.15,-0.25),c(0.15,-0.25))/2

Thetas = rbind(c(10,1.5),c(11,1.7),c(12,1.9),c(13,2),c(14,3),c(15,4))*1.5

LMs = c(0,1,2,3,4,5)

simLMPH(seed = 123, DES = DES, Betas = Betas, Thetas = Thetas, LMs = LMs, dist = "W")


plotLMPH(individuals = c(1,5), DES = DES, Betas = Betas, Thetas = Thetas, LMs = LMs, dist = "W", ymax  = 0.01, lout = 1000)

