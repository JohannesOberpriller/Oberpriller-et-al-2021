### Sampler Validation 
## This scirpt implements the validation of the sampler 
## Data is generated in the same way as model fits 
set.seed(123)
require(BayesianTools)
require(BASFERROR)
source("./helper_functions.R")
source("./analysis_functions.R")


forecastingtimestep = as.integer(100)
runtime = as.integer(85)
simulation_start = as.integer(1920)

# Get reference data
source("./GetReferenceData_add_noise.R")

# make reference data divideable by 100 which is the forecasting timestep 

GPPdata <- refGPP[1:(NDAYS-81)] 
ETdata <- refET[1:(NDAYS-81)]


## Sensitivity analysis and parameter selection
parSel <- c(48, 47, 11, 34,51, 52, 54, 55, 56)


## Get helper functions
source("./SSM/PMCMCfun_func.R")

## Make BayesianSetup
lower <- c(params[parSel]*0.3, 0.,0.,0.,0.,0.)
upper <- c(params[parSel]*1.7, .5,.5,1.,1.,1.)
best <- c(params[parSel], 0.2,0.2,0.2,0.2,0.2)

source("./SSM/createPrior.R")

# set up the bayesian setup 

BS <- createBayesianSetup(likelihood = LL, prior = prior)

#intialize storage of errors

counterERR <- 0
NERR <- numeric()
LUERR <- numeric()


# Run MCMC

settings <- list(iterations = 100000, consoleUpdates  = 1) 



out <- runMCMC(BS, sampler = "DEzs", settings = settings)

save(out, file = "./Results/MCMCchain_func_val.RData")
save(WAERR, file = "./Results/WA_func_val.RData")
save(RESPERR, file = "./Results/RESP_func_val.RData")
save(NPPERR, file = "./Results/NPP_func_val.RData")

save(predET, file = "./Results/PredET_func_val.RData")
save(predGPP, file = "./Results/PredGPP_func_val.RData")
