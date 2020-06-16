# define the forecasting timestep 

forecastingtime = as.integer(100)

#function to run the model 

model <- function(particle, par, weather_data, RESTART = as.integer(0)){
  
  stateers <- rnorm(3, mean = 1, sd = par[12:14])
  
  procerr <- rnorm(5, mean = 1, sd = 0.2)
  
  pars <- params
  pars[parSel] <- par[1:9]
  
  
  out <- run_mod_model(rs = RESTART, statespace = as.integer(1), 
                       bias = as.integer(0), randerr = as.integer(0),
                       ft = as.integer(1), p = pars, w = weather_data,
                       calf = calendar_fert, calN = calendar_Ndep,
                       calpT = calendar_prunT, caltT = calendar_thinT, forecastingtime,
                       NOUT = 24, sv = particle[4:17], stateers = stateers, procerr = procerr)

  return(out)
  
}

## Calculate approximate Likelihood
calculateApproxLikelihood <- function(weights){
  np <- length(weights[1,]) # number of particles
  
  w_hat <- numeric()
  
  for(i in 1:nrow(weights)){
    w_hat[i] <- logSumExp(weights[i,], mean = T)
  }
  
  margLL <- sum(w_hat)
  
  return(margLL)
}




## Point wise likelihood function

likelihood <- function(predictedGPP, predictedET, 
                       observedGPP, observedET, sdpar){
  
  diffGPP <-  abs(predictedGPP - observedGPP)
  diffET <- abs(predictedET - observedET)
  
  lik <- dnorm(diffGPP, sd = max(0.01, sdpar[1]*observedGPP), log = T) + 
    dnorm(diffET, sd = max(0.01, sdpar[2]*observedET), log = T)
  
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


## Make cluster
require(parallel)
nrCores <- 11

## Create cluster


tmpdlls <- getLoadedDLLs()
dlls <- vector(mode = "character", length = length(tmpdlls))
counter <- 0
for(i in tmpdlls){
  counter <- counter+1
  dlls[counter] <- i[[2]]
}

cl <- makeCluster(nrCores, useXDR = F)




## Add objects, dlls and packages
objects <- ls(envir = .GlobalEnv)
# #objects <- objects[objects != "matrix_weather"]
packages <- (.packages())
# #print(dlls)
parallel::clusterCall(cl, packageFun, packages, dlls)

# Export packages and objects to cluster

#parallel::clusterCall(cl, packageFun, packages = packages)
parallel::clusterExport(cl, varlist = objects)



## Likelihood function
LL <- function(parms){
  
  # rearrange the data 
  
  dataGPP = GPPdata
  dataET = ETdata
  
  # define the number of particles 
  
  numParticles = 220
  # subset the weather data 
  weather_data <- weather_BASFOR(as.integer(indYears[1]), as.integer(indDays[1]), forecastingtime, 
                                 clim[indLL[1]:(indLL[1]+forecastingtime-1),])
  # set up empty matrices to store results 
  
  resampling_weights = matrix(NA, nrow = length(indLL), ncol = numParticles)
  weights <- matrix(NA, nrow = length(dataGPP), ncol = numParticles)
  acceptederrorsWA <- matrix(NA, nrow = length(indLL), ncol = numParticles)
  acceptederrorsNPP <- matrix(NA, nrow = length(indLL), ncol = numParticles)
  acceptederrorsRESP <- matrix(NA, nrow = length(indLL), ncol = numParticles)
  predictionsGPP <- matrix(NA, nrow = length(indLL), ncol = numParticles) 
  predictionsET <- matrix(NA, nrow = length(indLL), ncol = numParticles)  
  
  # draw inital distribution 
  particles <-f0(numParticles)  ## generate start particles

  # make first step   
  particles <- matrix(parallel::parRapply(cl, x = particles, FUN = model,
                                          parms, weather_data, RESTART = as.integer(0)),
                      nrow = numParticles, byrow = TRUE )

  predictedGPP = particles[,(forecastingtime*19+1):(forecastingtime*20)]
  predictedET = particles[,(forecastingtime*21+1):(forecastingtime*22)]
  
  # calculate the weights 
  
  weights[1:forecastingtime,] <- likelihood(predictedGPP = predictedGPP, predictedET = predictedET,
                                            observedGPP = GPPdata[1:forecastingtime], 
                                            observedET = ETdata[1:forecastingtime], sdpar = parms[10:11]) # calculate weights for all
  
  
  #### Resample particles based on weights
  
  resampling_weights[1,] <- exp(weights[forecastingtime,] - logSumExp(weights[forecastingtime,]))
  
  indX <- sample.int(n = numParticles, size = numParticles, prob = resampling_weights[1,], replace = TRUE)
  
  orderseq = seq(from = 1, to = (forecastingtime*24), by = forecastingtime) + forecastingtime - 1
  
  
  particles <- particles[,orderseq]
  
  predictionsGPP[1,] <- particles[indX, 19]
  predictionsET[1,] <- particles[indX, 21]
  
  acceptederrorsWA[1,] <- particles[indX,22]
  acceptederrorsRESP[1,] <- particles[indX,23]
  acceptederrorsNPP[1,] <- particles[indX,24]
  
  particles <- particles[indX,]
  
  
  #iterate the previous steps to the end of the time 
  
  for(k in 2:(length(indLL))){ 
    # subsample the weather data 
    weather_data <- weather_BASFOR(as.integer(indYears[k]), as.integer(indDays[k]), forecastingtime, 
                                   clim[indLL[k]:(indLL[k]+forecastingtime-1),])
    
    .GlobalEnv$abbruchparticle <- particles
    
    # run the model for all particles 
    
    particles <- matrix(parallel::parRapply(cl, x = particles, FUN = model,
                                             parms, weather_data, RESTART = as.integer(1)),
                        nrow = numParticles, byrow = TRUE)
    
    
    
    
    predictedGPP = particles[,(forecastingtime*19+1):(forecastingtime*20)]
    predictedET = particles[,(forecastingtime*21+1):(forecastingtime*22)]
    
    #calculate the weights 
    
    weights[(forecastingtime*(k-1)+1):(forecastingtime*k),] <- likelihood(predictedGPP = predictedGPP,
                                                                          predictedET = predictedET,
                                                                          observedGPP = GPPdata[(forecastingtime*(k-1)+1):(forecastingtime*(k))], 
                                                                          observedET = ETdata[(forecastingtime*(k-1)+1):(forecastingtime*(k))],
                                                                          sdpar = parms[10:11])
    # calculate teh resampling weights 
    
    resampling_weights[k,] <- exp(weights[(forecastingtime*k),] - logSumExp(weights[(forecastingtime*k),]))
    
    
    
    #### Resample particles based on resampling_weights
    
    indX <- sample.int(n= numParticles,size = numParticles, prob = resampling_weights[k,], replace = TRUE)
    
    orderseq = seq(from = 1, to = (forecastingtime*24), by = forecastingtime) + forecastingtime - 1


    particles <- particles[,orderseq]
    
    predictionsGPP[k,] <- particles[indX, 19]
    predictionsET[k,] <- particles[indX, 21]
    
    acceptederrorsNPP[k,] <- particles[indX, 22]
    acceptederrorsWA[k,] <- particles[indX, 23]
    acceptederrorsRESP[k,] <- particles[indX, 24]
    
    particles <- particles[indX,]
  }
  
  
  .GlobalEnv$counterERR <- .GlobalEnv$counterERR + 1
  
  .GlobalEnv$WAERR[.GlobalEnv$counterERR] <- mean(acceptederrorsWA)
  
  .GlobalEnv$RESPERR[.GlobalEnv$counterERR] <- mean(acceptederrorsRESP)
  .GlobalEnv$NPPERR[.GlobalEnv$counterERR] <- mean(acceptederrorsNPP)
  
  .GlobalEnv$predET[.GlobalEnv$counterERR] <- mean(predictionsET)
  .GlobalEnv$predGPP[.GlobalEnv$counterERR] <- mean(predictionsGPP)
  
  
  currentll <- calculateApproxLikelihood(weights[indLL+forecastingtime-1,]) # get approximate ll values
  
  
  return(currentll)
}


inits <- output[nrow(output),]

numParticles = 220

f0 <- function(numParticles){
  out <- matrix(0, ncol = 24, nrow = numParticles)
  return(out)
}

logSumExp<- function(x, mean = F) {
  
  nObs = length(x)
  
  if(any(x == Inf)) stop("positive infinity values in log probabilities")
  if(any(x == -Inf )){
    message("encountered -Inf in logSumExp - value was removed")    
    x = x[x != -Inf] 
  }
  
  offset <- max(x)
  if (mean == T) out = log(sum(exp(x - offset))/nObs) + offset
  else out = log(sum(exp(x - offset))) + offset
  return(out)
  }
