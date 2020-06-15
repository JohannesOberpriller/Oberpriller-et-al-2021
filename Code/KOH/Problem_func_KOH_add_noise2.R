#!/usr/bin/env Rscript
## Functional biased model -- normal calibration 
require(MASS)
require(BayesianTools)
library(kernlab)
library(fitdistrplus)
library(actuar)
library(DEoptim)

forecastingtimestep = as.integer(100)
runtime = as.integer(85)
simulation_start = as.integer(1920)

source("./GetReferenceData_add_noise2.R")

indLL <- seq(1,(NDAYS-(forecastingtimestep -1)), by = forecastingtimestep)   # Marks the starting point of LL calculation
indYears <- matrix_weather[indLL, 1]
indDays <- matrix_weather[indLL, 2]

GPPdata <- refGPP[(indLL+(forecastingtimestep -1))]  # At the end of every 14 days data is used.
ETdata <- refET[(indLL+(forecastingtimestep -1))]

fullGPPdata <- refGPP[1:(NDAYS)]
fullETdata <- refET[1:(NDAYS)]

GPP_timestep = as.integer(5)
GPP_obs_times <- seq(1,(length(fullGPPdata)-(GPP_timestep -1)), by = GPP_timestep)


## Parameter values and selection
parSel <- c(48, 11, 34,51, 52, 56)

## Create BayesianSetup
lower <- c(params[parSel]*0.7, 0.,0,0,0)
upper <- c(params[parSel]*1.3, 1.,1.,1,1)
best <- c(params[parSel], 0.2, 0.2,0.2,0.2)


# Priors



density <- function(par){
  
  d1 <- sum(dunif(par[1:6], min = lower[1:6], max = upper[1:6], log = T))
  
  d2 <- sum(dgamma(par[7:10], shape = 2, scale = 0.1, log = T))
  
  return(d1 + d2)
  
}


sampler <- function(n = 1){
  
  d1 <- runif(n, lower[1], upper[1])
  d2 <- runif(n, lower[2], upper[2])
  d3 <- runif(n, lower[3], upper[3])
  d4 <- runif(n, lower[4], upper[4])
  d5 <- runif(n, lower[5], upper[5])
  d6 <- runif(n, lower[6], upper[6])
  
  
  
  d8 <- rgamma(n, shape = 2, scale = 0.1)
  d9 <- rgamma(n, shape = 2, scale = 0.1)
  d10 <- rgamma(n, shape = 2, scale = 0.1)
  d11 <- rgamma(n, shape = 2, scale = 0.1)
  
  
  return(cbind(d1,d2,d3,d4,d5,d6,d8,d9, d10, d11))
  
}



prior <- createPrior(density = density, sampler = sampler,
                     lower = lower, upper = upper, best = best)

get_sampling_points <- function(timeseries,indLL,NDAYS){
  tryCatch({
    weibulldist_fit <- fitdist(timeseries, distr = "weibull", method = "mle", lower = c(0, 0), start = list(scale = 1, shape = 1))
    weights_fit <- pinvweibull(timeseries,shape = weibulldist_fit$estimate[["shape"]], scale = weibulldist_fit$estimate[["scale"]])
    rand_points_GPP <- sample(x = 1:(NDAYS-81) ,2*length(indLL), prob = weights_fit)
    return(rand_points_GPP)
  },
  error = function(e){
    rand_points_GPP <- sample(x = 1:(NDAYS-81) ,2*length(indLL))
    return(rand_points_GPP)
  })
}



LL <- function(pars){
  
  params[parSel] <- pars[1:6]
  
  
  weather_data <- weather_BASFOR(as.integer(indYears[1]), as.integer(indDays[1]), NDAYS,clim)
  
  
  out <- run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
                       bias = as.integer(1), randerr = as.integer(0),
                       ft = as.integer(1), p = params, w = weather_data,
                       calf = calendar_fert, calN = calendar_Ndep,
                       calpT = calendar_prunT, caltT = calendar_thinT, NDAYS,
                       NOUT = 24, sv = STATEVARS, stateers = c(1,1,1), procerr = c(1,1,1,1,1))
  
  
  
  diffGPP <-  fullGPPdata - out[,19]
  diffET <- fullETdata - out[,21]
  
  
  GPP_sampled = cbind(out[(length(out[,19])-2000):length(out[,19]),19], weather_data[(length(out[,19])-2000):length(out[,19]),3:7])
  GPP_complet = cbind(out[,19], weather_data[1:length(out[,19]),3:7])
  
  ET_sampled = cbind(out[(length(out[,19])-2000):length(out[,21]),21], weather_data[(length(out[,19])-2000):length(out[,21]),3:7])
  ET_complet = cbind(out[,21], weather_data[1:length(out[,21]),3:7])
  
  fit_GPP = gausspr(x = GPP_sampled , y = diffGPP[(length(out[,19])-2000):length(out[,19])])
  fit_ET = gausspr(x = ET_sampled,  y = diffET[(length(out[,21])-2000):length(out[,21])])
  
  predict_GPP = predict(fit_GPP, GPP_complet)
  predict_ET = predict(fit_ET, ET_complet)
  
  
  diffGPP <-  abs(fullGPPdata - out[,19] - predict_GPP)
  diffET <- abs(fullETdata - out[,21] - predict_ET)
  
  
  lik <- sum(dnorm(diffGPP, sd = pars[7], log = T)) +
    sum(dnorm(diffET, sd = pars[8], log = T)) + 
    sum(dnorm(predict_GPP, sd = pars[9], log = T)) +
    sum(dnorm(predict_ET, sd = pars[10], log = T)) 
  gc()
  return(lik)
  
}


packageFun <- function(packages = NULL, dlls = NULL) {
  if(!is.null(packages)){
    for(i in packages) library(i, character.only = TRUE)
  }
  if(!is.null(dlls)){
    for(i in dlls) try(dyn.load(i), silent = T)
  }
}



neg_LL <- function(parms){
  ll <- tryCatch(-LL(parms),
                 error = function(e){return(Inf)})
  return(ll)
}


#bayesianSetup <- createBayesianSetup(LL, prior = prior, parallel = 'F')

cl2 <- parallel::makeCluster(24)

objects2 <- ls(envir = .GlobalEnv)

packages2 <- (.packages())

parallel::clusterCall(cl2, packageFun, packages2)

parallel::clusterExport(cl2, varlist = objects2)



out_optim <- DEoptim(neg_LL, lower = lower, upper = upper, control = list(itermax = 200, trace = F, cluster = cl2, NP = 120))

stopCluster(cl2)



bayesianSetup <- createBayesianSetup(LL, prior = prior, parallel = 24)


out_MCMC <- runMCMC(bayesianSetup, settings = list(iterations = 500000, nrChains = 1, Z = out_optim$member$pop, startValue = out_optim$member$pop))


save(out_MCMC, file = "./Results/MCMC_func_KOH_add_noise2.RData")