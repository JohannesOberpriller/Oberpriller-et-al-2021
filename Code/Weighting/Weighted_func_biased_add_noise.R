#!/usr/bin/env Rscript
## Functional biased model -- normal calibration 

library(DEoptim)
require(BayesianTools)
library(BASFERROR)

#### Setting up everything for Calibration ####

params <- df_params[,1]
forecastingtimestep = as.integer(100)
runtime = as.integer(85)
simulation_start = as.integer(1920)
## generate the reference data 
source("./GetReferenceData_add_noise.R")

indLL <- seq(1,(NDAYS-(forecastingtimestep -1)), by = forecastingtimestep)   # Marks the starting point of LL calculation
indYears <- matrix_weather[indLL, 1]
indDays <- matrix_weather[indLL, 2]

GPPdata <- refGPP[(indLL+(forecastingtimestep -1))]  # At the end of every 14 days data is used.
ETdata <- refET[(indLL+(forecastingtimestep -1))]

fullGPPdata <- refGPP[1:(NDAYS)]
fullETdata <- refET[1:(NDAYS)]

ET_timestep = as.integer(10)
ET_obs_times <- seq(1,(length(fullETdata)-(ET_timestep -1)), by = ET_timestep)


## Parameter values and selection
parSel <- c(48, 11, 34,51, 52, 56)

## Create BayesianSetup
lower <- c(params[parSel]*0.7, 0.,0.)
upper <- c(params[parSel]*1.3, 1.,1.)
best <- c(params[parSel], 0.2, 0.2)




# Priors

## set up density for the parameters to calculate their density

density <- function(par){
  
  d1 <- sum(dunif(par[1:6], min = lower[1:6], max = upper[1:6], log = T))
  
  d2 <- sum(dgamma(par[7:8], shape = 2, scale = 0.1, log = T))
  
  return(d1 + d2)
  
}

## set up sampler for sampling parameters

sampler <- function(n = 1){
  ## flat priors fot the model parameters
  d1 <- runif(n, lower[1], upper[1])
  d2 <- runif(n, lower[2], upper[2])
  d3 <- runif(n, lower[3], upper[3])
  d4 <- runif(n, lower[4], upper[4])
  d5 <- runif(n, lower[5], upper[5])
  d6 <- runif(n, lower[6], upper[6])
  
  # gamma priors for the standard deviations 
  
  d8 <- rgamma(n, shape = 2, scale = 0.1)
  d9 <- rgamma(n, shape = 2, scale = 0.1)
  
  
  return(cbind(d1,d2,d3,d4,d5,d6,d8,d9))
  
}

## prior object for the use in the MCMC sampler BayesianTools 

prior <- createPrior(density = density, sampler = sampler,
                     lower = lower, upper = upper, best = best)

## Likelihood function to be used in statistical inference 

LL <- function(pars){
  ## overwrite default model parameters 
  params[parSel] <- pars[1:6]

  # prepare weather data for basfor   
  weather_data <- weather_BASFOR(as.integer(indYears[1]), as.integer(indDays[1]), NDAYS,clim)
  
  # preapre inital state variables for basfor 
  STATEVARS = rep(0,14)
  # run the model 
  out <- run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
                       bias = as.integer(1), randerr = as.integer(0),
                       ft = as.integer(1), p = params, w = weather_data,
                       calf = calendar_fert, calN = calendar_Ndep,
                       calpT = calendar_prunT, caltT = calendar_thinT, as.integer(NDAYS),
                       NOUT = 24, sv = STATEVARS, stateers = c(1,1,1), procerr = c(1,1,1,1,1))
  
  # calculate residuals 
  diffGPP <- abs(out[,19] - fullGPPdata)
  diffET <- abs(out[ET_obs_times,21] - fullETdata[ET_obs_times])
  
  
  sdGPP <- sdET <- numeric(length(fullGPPdata))
  
  # calculate the likelihood 
  lik <- sum(dnorm(diffGPP, sd = pars[7], log = T)) +
    sum(dnorm(diffET, sd = pars[8], log = T))*(length(diffGPP)/length(diffET))
  
  return(lik)
  
}

## Function handling the export of packages to the cluster 

packageFun <- function(packages = NULL, dlls = NULL) {
  if(!is.null(packages)){
    for(i in packages) library(i, character.only = TRUE)
  }
  if(!is.null(dlls)){
    for(i in dlls) try(dyn.load(i), silent = T)
  }
}

## negative log-likelihood for the optimizer

neg_LL <- function(parms){
  ll <- tryCatch(-LL(parms),error = function(e){return(Inf)})
  return(ll)
}


#### First Calibration ####

cl2 <- parallel::makeCluster(24)

objects2 <- ls(envir = .GlobalEnv)

packages2 <- (.packages())

parallel::clusterCall(cl2, packageFun, packages2)

parallel::clusterExport(cl2, varlist = objects2)


out_optim <- DEoptim(neg_LL, lower = lower, upper = upper, control = list(itermax = 200, trace = F, cluster = cl2, NP = 120))


stopCluster(cl2)

bayesianSetup <- createBayesianSetup(LL, prior = prior, parallel = 24)

out_MCMC <- runMCMC(bayesianSetup, settings = list(iterations = 500000, nrChains = 1, Z = out_optim$member$pop, startValue = out_optim$member$pop))

save(out_MCMC, file = "./Results/MCMC_unbal_func_biased_add_noise_weighted.RData")

rm(cl2,out_optim, out_MCMC, bayesianSetup)

#### Second Calibration ####


## generate the data for calibration 

source("./GetReferenceData_add_noise2.R")

GPPdata <- refGPP[(indLL+(forecastingtimestep -1))]  
ETdata <- refET[(indLL+(forecastingtimestep -1))]

fullGPPdata <- refGPP[1:(NDAYS)]
fullETdata <- refET[1:(NDAYS)]

## perform calibration 

cl2 <- parallel::makeCluster(24)

objects2 <- ls(envir = .GlobalEnv)

packages2 <- (.packages())

parallel::clusterCall(cl2, packageFun, packages2)

parallel::clusterExport(cl2, varlist = objects2)


out_optim <- DEoptim(neg_LL, lower = lower, upper = upper, control = list(itermax = 200, trace = F, cluster = cl2, NP = 120))


stopCluster(cl2)

bayesianSetup <- createBayesianSetup(LL, prior = prior, parallel = 24)

out_MCMC <- runMCMC(bayesianSetup, settings = list(iterations = 500000, nrChains = 1, Z = out_optim$member$pop, startValue = out_optim$member$pop))

save(out_MCMC, file = "./Results/MCMC_unbal_func_biased_add_noise_weighted2.RData")

rm(cl2,out_optim, out_MCMC, bayesianSetup)

#### Third Calibration ####


## generate the data for calibration

source("./GetReferenceData_add_noise3.R")

GPPdata <- refGPP[(indLL+(forecastingtimestep -1))]  
ETdata <- refET[(indLL+(forecastingtimestep -1))]

fullGPPdata <- refGPP[1:(NDAYS)]
fullETdata <- refET[1:(NDAYS)]

## perform calibration 

cl2 <- parallel::makeCluster(24)

objects2 <- ls(envir = .GlobalEnv)

packages2 <- (.packages())

parallel::clusterCall(cl2, packageFun, packages2)

parallel::clusterExport(cl2, varlist = objects2)


out_optim <- DEoptim(neg_LL, lower = lower, upper = upper, control = list(itermax = 200, trace = F, cluster = cl2, NP = 120))


stopCluster(cl2)

bayesianSetup <- createBayesianSetup(LL, prior = prior, parallel = 24)

out_MCMC <- runMCMC(bayesianSetup, settings = list(iterations = 500000, nrChains = 1, Z = out_optim$member$pop, startValue = out_optim$member$pop))

save(out_MCMC, file = "./Results/MCMC_unbal_func_biased_add_noise_weighted3.RData")

rm(cl2,out_optim, out_MCMC, bayesianSetup)

#### Fourth Calibration ####


## generate the data for calibration

source("./GetReferenceData_add_noise4.R")

GPPdata <- refGPP[(indLL+(forecastingtimestep -1))]  
ETdata <- refET[(indLL+(forecastingtimestep -1))]

fullGPPdata <- refGPP[1:(NDAYS)]
fullETdata <- refET[1:(NDAYS)]

## perform calibration 

cl2 <- parallel::makeCluster(24)

objects2 <- ls(envir = .GlobalEnv)

packages2 <- (.packages())

parallel::clusterCall(cl2, packageFun, packages2)

parallel::clusterExport(cl2, varlist = objects2)


out_optim <- DEoptim(neg_LL, lower = lower, upper = upper, control = list(itermax = 200, trace = F, cluster = cl2, NP = 120))


stopCluster(cl2)

bayesianSetup <- createBayesianSetup(LL, prior = prior, parallel = 24)

out_MCMC <- runMCMC(bayesianSetup, settings = list(iterations = 500000, nrChains = 1, Z = out_optim$member$pop, startValue = out_optim$member$pop))

save(out_MCMC, file = "./Results/MCMC_unbal_func_biased_add_noise_weighted4.RData")

rm(cl2,out_optim, out_MCMC, bayesianSetup)

#### Fifth Calibration ####


## generate the data for calibration

source("./GetReferenceData_add_noise5.R")

GPPdata <- refGPP[(indLL+(forecastingtimestep -1))]  
ETdata <- refET[(indLL+(forecastingtimestep -1))]

fullGPPdata <- refGPP[1:(NDAYS)]
fullETdata <- refET[1:(NDAYS)]

## perform calibration 

cl2 <- parallel::makeCluster(24)

objects2 <- ls(envir = .GlobalEnv)

packages2 <- (.packages())

parallel::clusterCall(cl2, packageFun, packages2)

parallel::clusterExport(cl2, varlist = objects2)


out_optim <- DEoptim(neg_LL, lower = lower, upper = upper, control = list(itermax = 200, trace = F, cluster = cl2, NP = 120))


stopCluster(cl2)

bayesianSetup <- createBayesianSetup(LL, prior = prior, parallel = 24)

out_MCMC <- runMCMC(bayesianSetup, settings = list(iterations = 500000, nrChains = 1, Z = out_optim$member$pop, startValue = out_optim$member$pop))

save(out_MCMC, file = "./Results/MCMC_unbal_func_biased_add_noise_weighted5.RData")