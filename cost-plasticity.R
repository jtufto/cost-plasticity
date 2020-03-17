library(MASS)
library(TMB)
compile("cost_plasticity.cpp")
dyn.load(dynlib("cost_plasticity"))

# Experimental design
m <- 1000 # number of families
i <- rep(1:m, each=10) 
n <- length(i)

# True parameter values
G <- diag(c(1,.2)) # G-matrix
sigma_E <- 1
abbarbar <- c(0,.5) # mean elevation and slope in the population

# Simulated data (will add fitness data later)
set.seed(1)
ab <- mvrnorm(m, abbarbar, G/2)[i,] + mvrnorm(n, c(0,0), Sigma = G/2) # individual elevations and slopes
epsilon <- rnorm(n) # environment of development
z <- ab[,1] + ab[,2]*epsilon + rnorm(n, sd=sigma_E) # phenotypes

# Compute marginal likelihood
data <- list(z=z, 
             epsilon=epsilon, 
             i=i-1) 
parameters <- list(abbarbar=abbarbar,
                   log_sigma_E=0,
                   log_sigma_ab=c(0,0),
                   rho=0,
                   abbar=matrix(0,m,2),
                   ab=matrix(0,n,2))
obj <- MakeADFun(data, parameters, DLL="cost_plasticity", random = c("abbar","ab"))

# Fit the model
fit <- nlminb(obj$par, obj$fn, obj$gr)
fit$par
rep <- sdreport(obj)
summary(rep,"report")


