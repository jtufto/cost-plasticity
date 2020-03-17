library(MASS)
library(TMB)
library(glmmTMB)
compile("cost_plasticity.cpp")
dyn.load(dynlib("cost_plasticity"))

# Experimental design
m <- 100 # number of families
i <- rep(1:m, each=10) 
n <- length(i)

# True parameter values
G <- diag(c(1,.1)) # G-matrix
sigma_E <- .1
abbarbar <- c(0,.5) # mean elevation and slope in the population

# Simulated data
set.seed(1)
ab <- mvrnorm(m, abbarbar, G/2)[i,] + mvrnorm(n, c(0,0), Sigma = G/2) # individual elevations and slopes
epsilon <- rnorm(n) # environment of development
z <- ab[,1] + ab[,2]*epsilon + rnorm(n, sd=sigma_E) # phenotypes


data <- list(z=z, 
             epsilon=epsilon, 
             i=i)
parameters <- list(a=abbarbar[1],
                   b=abbarbar[2],
                   logSigma=0)
obj <- MakeADFun(data, parameters, DLL="cost_plasticity")
fit <- nlminb(obj$par, obj$fn, obj$gr)
fit$par
lm(z~epsilon)
