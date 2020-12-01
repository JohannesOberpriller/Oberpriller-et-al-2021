set.seed(123)
library(BayesianTools)
library(kernlab)
library(scales)
library(mvnfast)

forecastingtimestep = as.integer(100)
runtime = as.integer(85)
simulation_start = as.integer(1920)

source("./GetReferenceData_add_noise.R")


indLL <- seq(1,(NDAYS-(forecastingtimestep -1)), by = forecastingtimestep)   # Marks the starting point of LL calculation
indYears <- matrix_weather[indLL, 1]
indDays <- matrix_weather[indLL, 2]

GPPdata <- refGPP[(indLL+(forecastingtimestep -1))]  # At the end of every 14 days data is used.
ETdata <- refET[(indLL+(forecastingtimestep -1))]

fullGPPdata <- refGPP[1:(NDAYS)]
fullETdata <- refET[1:(NDAYS)]

length(fullGPPdata)

GPP_timestep = as.integer(5)
GPP_obs_times <- seq(1,(length(fullGPPdata)-(GPP_timestep -1)), by = GPP_timestep)


## Parameter values and selection
parSel <- c(48, 11, 34,51, 52, 56)

## Create BayesianSetup
lower <- c(params[parSel]*0.7, 0.,0,0,0)
upper <- c(params[parSel]*1.3, 1.,1.,1,1)
best <- c(params[parSel], 0.2, 0.2,0.2,0.2)


forecastingtimestep = as.integer(100)
runtime = as.integer(95)
simulation_start = as.integer(1920)

# 2. get the climate data, which are model inputs, in order to run the model 

climate <- readRDS("./data/climate_complete.rds")
calendar_Ndep <- as.matrix(readRDS("./data/Ndep_complete.rds"))
calendar_Ndep <- calendar_Ndep[which(calendar_Ndep[,1] >= 1919 ),]

# 3.  Set the dates when fertilization, pruning and thinning was 
# We set this to -1, indicating it did not happen 

calendar_fert  <- matrix( -1, nrow=100, ncol=3 )
calendar_prunT <- matrix( -1, nrow=100, ncol=3 )
calendar_thinT <- matrix( -1, nrow=100, ncol=3 )

# 4. getting the default parameters of BASFOR which are saved under df_params 
# we do this for coniferous trees

params <- df_params[,1]



# 5. bring the climate data into and appropriate shape 
# and exclude data which is before the time and do the same 
# for the nitrogen deposition 

clim <- data.frame("year" = climate[,1],
                   "doy" = climate[,2],
                   "GR" = climate[,3],
                   "T" = climate[,4],
                   "RAIN" = climate[,5],
                   "WN" = climate[,6],
                   "VP" = climate[,7])

clim <- clim[which(clim[,1] >= 1920),]

fixed_data <- prepare_data(calendar_Ndep, clim, simulation_start, runtime)
clim <- fixed_data[[1]]
calendar_Ndep <- fixed_data[[2]]


# 5. genearate starting values of the state variables and get weather values ready 

matrix_weather <- weather_BASFOR(as.integer(1920), as.integer(1), NDAYS, clim)
STATEVARS <- rep(0,14)


## Run model to set new default parameters 

output <- run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
                        bias = as.integer(0), randerr = as.integer(0), 
                        ft = as.integer(1), p = params, w = matrix_weather,
                        calf = calendar_fert, calN = calendar_Ndep,
                        calpT = calendar_prunT, caltT = calendar_thinT, 
                        n = NDAYS, NOUT = 24, sv = STATEVARS, 
                        stateers = c(1,1,1), procerr = c(1,1,1,1,1))

## Change some initial values
params[57] <- TREEDENS0 <- output[30681,4]
params[4] <- CRtree0 <- output[30681,5] / TREEDENS0
params[2] <- CBtree0 <- output[30681,7] / TREEDENS0
params[3] <- CLtree0 <- output[30681,8] / TREEDENS0
params[5] <- CStree0 <- output[30681,6] / TREEDENS0
params[33] <- CLITT0 <- output[30681,10]
params[43] <- NMIN0 <- output[30681,16]

indLL <- seq(1,(NDAYS-(forecastingtimestep -1)), by = forecastingtimestep)   # Marks the starting point of LL calculation
indYears <- matrix_weather[indLL, 1]
indDays <- matrix_weather[indLL, 2]

weather_data <- weather_BASFOR(as.integer(indYears[1]), as.integer(indDays[1]), NDAYS,
                               clim)

LL_approx <- function(pars){
  
  params[parSel] <- pars[1:6]
  
  weather_data <- weather_BASFOR(as.integer(indYears[1]), as.integer(indDays[1]), NDAYS,clim)
  
  out <- run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
                       bias = as.integer(0), randerr = as.integer(0),
                       ft = as.integer(1), p = params, w = weather_data,
                       calf = calendar_fert, calN = calendar_Ndep,
                       calpT = calendar_prunT, caltT = calendar_thinT, NDAYS,
                       NOUT = 24, sv = STATEVARS, stateers = c(1,1,1), procerr = c(1,1,1,1,1))

  diffGPP = fullGPPdata-out[,19]
  diffET = fullETdata-out[,21]

  
  
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

load("./Results/MCMC_func_KOH_add_noise.RData")

posterior_sample = getSample(out_MCMC, start = 200, numSamples = 240)

posterior_sample = posterior_sample[1:240,]

print(posterior_sample)

sample_one = apply(FUN = LL_approx, X = posterior_sample, MARGIN = 1)

forecastingtimestep = as.integer(100)
runtime = as.integer(85)
simulation_start = as.integer(1920)

source("./GetReferenceData_add_noise.R")


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


forecastingtimestep = as.integer(100)
runtime = as.integer(95)
simulation_start = as.integer(1920)

# 2. get the climate data, which are model inputs, in order to run the model 

climate <- readRDS("./data/climate_complete.rds")
calendar_Ndep <- as.matrix(readRDS("./data/Ndep_complete.rds"))
calendar_Ndep <- calendar_Ndep[which(calendar_Ndep[,1] >= 1919 ),]

# 3.  Set the dates when fertilization, pruning and thinning was 
# We set this to -1, indicating it did not happen 

calendar_fert  <- matrix( -1, nrow=100, ncol=3 )
calendar_prunT <- matrix( -1, nrow=100, ncol=3 )
calendar_thinT <- matrix( -1, nrow=100, ncol=3 )

# 4. getting the default parameters of BASFOR which are saved under df_params 
# we do this for coniferous trees

params <- df_params[,1]



# 5. bring the climate data into and appropriate shape 
# and exclude data which is before the time and do the same 
# for the nitrogen deposition 

clim <- data.frame("year" = climate[,1],
                   "doy" = climate[,2],
                   "GR" = climate[,3],
                   "T" = climate[,4],
                   "RAIN" = climate[,5],
                   "WN" = climate[,6],
                   "VP" = climate[,7])

clim <- clim[which(clim[,1] >= 1920),]

fixed_data <- prepare_data(calendar_Ndep, clim, simulation_start, runtime)
clim <- fixed_data[[1]]
calendar_Ndep <- fixed_data[[2]]


# 5. genearate starting values of the state variables and get weather values ready 

matrix_weather <- weather_BASFOR(as.integer(1920), as.integer(1), NDAYS, clim)
STATEVARS <- rep(0,14)


## Run model to set new default parameters 

output <- run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
                        bias = as.integer(0), randerr = as.integer(0), 
                        ft = as.integer(1), p = params, w = matrix_weather,
                        calf = calendar_fert, calN = calendar_Ndep,
                        calpT = calendar_prunT, caltT = calendar_thinT, 
                        n = NDAYS, NOUT = 24, sv = STATEVARS, 
                        stateers = c(1,1,1), procerr = c(1,1,1,1,1))

## Change some initial values
params[57] <- TREEDENS0 <- output[30681,4]
params[4] <- CRtree0 <- output[30681,5] / TREEDENS0
params[2] <- CBtree0 <- output[30681,7] / TREEDENS0
params[3] <- CLtree0 <- output[30681,8] / TREEDENS0
params[5] <- CStree0 <- output[30681,6] / TREEDENS0
params[33] <- CLITT0 <- output[30681,10]
params[43] <- NMIN0 <- output[30681,16]

indLL <- seq(1,(NDAYS-(forecastingtimestep -1)), by = forecastingtimestep)   # Marks the starting point of LL calculation
indYears <- matrix_weather[indLL, 1]
indDays <- matrix_weather[indLL, 2]

weather_data <- weather_BASFOR(as.integer(indYears[1]), as.integer(indDays[1]), NDAYS,
                               clim)


# LL_full <- function(pars){
#   
#   print(pars)
#   params[parSel] <- pars[1:6]
#   
#   print(params)
#   
#   
#   weather_data <- weather_BASFOR(as.integer(indYears[1]), as.integer(indDays[1]), NDAYS,clim)
#   
#   
#   out <- run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
#                        bias = as.integer(1), randerr = as.integer(0),
#                        ft = as.integer(1), p = params, w = weather_data,
#                        calf = calendar_fert, calN = calendar_Ndep,
#                        calpT = calendar_prunT, caltT = calendar_thinT, NDAYS,
#                        NOUT = 24, sv = STATEVARS, stateers = c(1,1,1), procerr = c(1,1,1,1,1))
#   
# 
#   diffGPP <-  fullGPPdata - out[,19]
#   diffET <- fullETdata - out[,21]
#   
#   
#   GPP_sampled = cbind(out[(length(out[,19])-2000):length(out[,19]),19], weather_data[(length(out[,19])-2000):length(out[,19]),3:7])
#   GPP_complet = cbind(out[,19], weather_data[1:length(out[,19]),3:7])
#   
#   ET_sampled = cbind(out[(length(out[,19])-2000):length(out[,21]),21], weather_data[(length(out[,19])-2000):length(out[,21]),3:7])
#   ET_complet = cbind(out[,21], weather_data[1:length(out[,21]),3:7])
#   
#   fit_GPP = gausspr(x = GPP_sampled , y = diffGPP[(length(out[,19])-2000):length(out[,19])], variance.model =T)
#   fit_ET = gausspr(x = ET_sampled,  y = diffET[(length(out[,21])-2000):length(out[,21])])
#   
#   predict_GPP = predict(fit_GPP, GPP_complet)
#   predict_ET = predict(fit_ET, ET_complet)
#   
#   kernel_GPP = rbfdot(sigma = fit_GPP@kernelf@kpar$sigma)
#   kernel = kernelMatrix(kernel_GPP, GPP_complet)
#   
#   inverse_kernel = chol2inv(chol(kernel))
# 
#   lik_GPP = - 1/2*sum(predict_GPP*(inverse_kernel%*%predict_GPP)) # + log(1/sqrt(2*pi*det(inverse_kernel))) 
#   print(str(inverse_kernel))
#   print(lik_GPP)
#   rm(inverse_kernel)
#   
#   
# 
#   kernel_ET = rbfdot(sigma = fit_ET@kernelf@kpar$sigma)
#   kernel = kernelMatrix(kernel_ET, ET_complet)
# 
#   inverse_kernel = chol2inv(chol(kernel))
# 
#   lik_ET =  -1/2*sum(predict_ET*(inverse_kernel%*%predict_ET)) # + log(1/sqrt(2*pi*det(inverse_kernel)))
# 
#   rm(inverse_kernel)
#   diffGPP <-  abs(fullGPPdata - out[,19] - predict_GPP)
#   diffET <- abs(fullETdata - out[,21] - predict_ET)
# 
# 
#   lik <- sum(dnorm(diffGPP, sd = pars[7], log = T)) +
#     sum(dnorm(diffET, sd = pars[8], log = T)) + 
#     lik_GPP + lik_ET
#   gc()
#   return(lik)
#   
# }

LL_full2 <- function(pars){

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

  fit_GPP = gausspr(x = GPP_sampled , y = diffGPP[(length(out[,19])-2000):length(out[,19])], variance.model =T)
  fit_ET = gausspr(x = ET_sampled,  y = diffET[(length(out[,21])-2000):length(out[,21])])
  
  
  
  predict_GPP = predict(fit_GPP, GPP_complet)
  predict_ET = predict(fit_ET, ET_complet)
  

  kernel_GPP = rbfdot(sigma = fit_GPP@kernelf@kpar$sigma)
  kernel = kernelMatrix(kernel_GPP, GPP_sampled)


  lik_GPP = dmvn(c(fit_GPP@fitted), mu =rep(0,length(fit_GPP@fitted)), sigma = kernel,log = T)

  kernel_ET = rbfdot(sigma = fit_ET@kernelf@kpar$sigma)
  kernel = kernelMatrix(kernel_ET, ET_sampled)


  lik_ET = dmvn(c(fit_ET@fitted), mu =rep(0,length(fit_ET@fitted)), sigma = kernel,log = T)

   # rm(inverse_kernel)
  diffGPP <-  abs(fullGPPdata - out[,19] - predict_GPP)
  diffET <- abs(fullETdata - out[,21] - predict_ET)


  lik <- sum(dnorm(diffGPP, sd = pars[7], log = T)) +
    sum(dnorm(diffET, sd = pars[8], log = T)) +
    lik_GPP + lik_ET
  gc()
  return(lik)

}


posterior_sample = posterior_sample[sample(1:240,1), ]

print(posterior_sample)

sample_two = LL_full2(posterior_sample)


saveRDS(sample_two, "./Results/likelihood_full.Rds")
saveRDS(sample_one, "./Results/likelihood_of_approx.Rds")


full = readRDS("./Results/likelihood_full.Rds")
approx = readRDS("./Results/likelihood_of_approx.Rds")

plot(density(approx), main = "Differences approximation- full likelihood", ylab = "Density",
     xlab = "Log-posterior", col = scales::alpha("steelblue", 0.8))
polygon(density(approx),col = scales::alpha("steelblue", 0.8 ))
abline(v = full, col = scales::alpha("darkred", 0.8), lwd = 2)
legend("topright", fill = c(scales::alpha("steelblue", 0.8 ),scales::alpha("darkred", 0.8)), 
       legend = c("Posterior distribution with approximation",
                  "Sampled log-posterior full"))

