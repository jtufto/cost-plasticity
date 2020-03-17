library(MASS)
library(TMB)
compile("cost_plasticity.cpp")
dyn.load(dynlib("cost_plasticity"))

# Experimental design
m <- 1000 # number of families
i <- rep(1:m, each=10) 
n <- length(i)

# True parameter values
G <- diag(c(1,.1)) # G-matrix
sigma_E <- 1
abbarbar <- c(0,.5) # mean elevation and slope in the population
beta <- c(0,1,-.2,.5,0)

# Simulated data phenotypes and fitnesses
set.seed(1)
ab <- mvrnorm(m, abbarbar, G/2)[i,] + mvrnorm(n, c(0,0), Sigma = G/2) # individual elevations and slopes
epsilon <- rnorm(n) # environment of development
z <- ab[,1] + ab[,2]*epsilon + rnorm(n, sd=sigma_E) # phenotypes
eta <- beta[1] + beta[2]*z + beta[3]*z^2 + beta[4]*ab[,2] + beta[5]*ab[,2]^2
y <- rpois(n, exp(eta))

# Compute marginal likelihood
data <- list(z = z, 
             epsilon = epsilon, 
             y = y,
             i = i-1) 
parameters <- list(abbarbar = abbarbar,
                   beta = beta,
                   log_sigma_E = 0,
                   log_sigma_ab = c(0,0),
                   rho = 0,
                   abbar = matrix(0,m,2),
                   ab = matrix(0,n,2))
map <- list(
  beta = factor(c(1,2,3,4,NA))
)
obj <- MakeADFun(data, parameters, DLL="cost_plasticity", random = c("abbar","ab"), map=map)

# Fit the model
fit <- nlminb(obj$par, obj$fn, obj$gr)
fit$par
rep <- sdreport(obj)
summary(rep,c("fixed","report"))


