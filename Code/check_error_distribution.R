set.seed(42)
library(BayesianTools)
library(DHARMa)
library(BASFERROR)
library(scales)

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




### here comes the real check of error distributions ####


### Analysis from correct model ###

load("./Results/Model_validation_chain_less_parameter_add_noise.RData")

posterior_samples = getSample(out_MCMC, start = 200, numSamples = 360)

GPP_simulations = matrix(nrow = nrow(posterior_samples),
                             ncol = NDAYS)
ET_simulations = matrix(nrow = nrow(posterior_samples),
                            ncol = NDAYS)

residuals_GPP = matrix(nrow = nrow(posterior_samples),
                         ncol = NDAYS)
residuals_ET = matrix(nrow = nrow(posterior_samples),
                        ncol = NDAYS)

for(i in 1:nrow(posterior_samples)){
  params[parSel] <- posterior_samples[i,1:6]
  out <- run_mod_model(rs = as.integer(0), statespace = as.integer(0),
                       bias = as.integer(0), randerr = as.integer(0),
                       ft = as.integer(1), p = params, w = weather_data,
                       calf = calendar_fert, calN = calendar_Ndep,
                       calpT = calendar_prunT, caltT = calendar_thinT, NDAYS,
                       NOUT = 24, sv = STATEVARS, stateers = c(1,1,1), procerr = c(1,1,1,1,1))

  GPP_simulations[i,] = out[,19]
  ET_simulations[i,] = out[,21]

  residuals_GPP[i,] = out[,19]- fullGPPdata
  residuals_ET[i,] = out[,21] - fullETdata
}


load("./Results/MCMC_func_biased_add_noise.RData")

posterior_samples = getSample(out_MCMC, start = 200, numSamples = 360)

GPP_simulations_biased = matrix(nrow = nrow(posterior_samples),
                         ncol = NDAYS)
ET_simulations_biased = matrix(nrow = nrow(posterior_samples),
                        ncol = NDAYS)

residuals_GPP_biased = matrix(nrow = nrow(posterior_samples),
                       ncol = NDAYS)
residuals_ET_biased = matrix(nrow = nrow(posterior_samples),
                      ncol = NDAYS)

for(i in 1:nrow(posterior_samples)){
  params[parSel] <- posterior_samples[i,1:6]
  out <- run_mod_model(rs = as.integer(0), statespace = as.integer(0),
                       bias = as.integer(1), randerr = as.integer(0),
                       ft = as.integer(1), p = params, w = weather_data,
                       calf = calendar_fert, calN = calendar_Ndep,
                       calpT = calendar_prunT, caltT = calendar_thinT, NDAYS,
                       NOUT = 24, sv = STATEVARS, stateers = c(1,1,1), procerr = c(1,1,1,1,1))
  
  GPP_simulations_biased[i,] = out[,19]
  ET_simulations_biased[i,] = out[,21]
  
  residuals_GPP_biased[i,] = out[,19]- fullGPPdata
  residuals_ET_biased[i,] = out[,21] - fullETdata
}

residuals_GPP_biased_corrected = matrix(nrow = nrow(posterior_samples),
                              ncol = NDAYS)
residuals_ET_biased_corrected = matrix(nrow = nrow(posterior_samples),
                             ncol = NDAYS)

for(i in 1:nrow(posterior_samples)){
  params[parSel] <- posterior_samples[i,1:6]
  out <- run_mod_model(rs = as.integer(0), statespace = as.integer(0),
                       bias = as.integer(0), randerr = as.integer(0),
                       ft = as.integer(1), p = params, w = weather_data,
                       calf = calendar_fert, calN = calendar_Ndep,
                       calpT = calendar_prunT, caltT = calendar_thinT, NDAYS,
                       NOUT = 24, sv = STATEVARS, stateers = c(1,1,1), procerr = c(1,1,1,1,1))
  

  residuals_GPP_biased_corrected[i,] = out[,19]- fullGPPdata
  residuals_ET_biased_corrected[i,] = out[,21] - fullETdata
}

density_GPP_biased = density(as.vector(residuals_GPP_biased), bw = 0.3)
plot(density_GPP_biased, col = "red", xlim = c(-2,2), xlab = "Residuals", " Density of residuals", ylim = c(0,2))
polygon(density_GPP_biased, col = alpha("darkred", 0.4))
abline(v = median(as.vector(residuals_GPP_biased)), col = alpha("darkred",0.8), lwd = 2)
abline(v = median(as.vector(residuals_GPP)), col = alpha("steelblue",0.8) ,lwd = 2)
polygon(density(as.vector(residuals_GPP), bw = 0.1),  col = alpha("steelblue",0.4))
legend("topright", fill = c(alpha("steelblue", 0.8 ),alpha("darkred", 0.8)), 
       legend = c("Median residuals true model",
                  "Median residuals biased models"))







