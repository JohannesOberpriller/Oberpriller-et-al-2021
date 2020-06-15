## Getting data for sampler validation 
set.seed(12345)
require(BASFERROR)
source("./helper_functions.R")

# set the forecasting timestep, the runtime and the start of the simulation 

forecastingtimestep = as.integer(100)
runtime = as.integer(85)
simulation_start = as.integer(1920)

# get the climate data, which are model inputs, in order to run the model 

climate <- readRDS("./data/climate_complete.rds")
calendar_Ndep <- as.matrix(readRDS("./data/Ndep_complete.rds"))
calendar_Ndep <- calendar_Ndep[which(calendar_Ndep[,1] >= 1919 ),]

# Set the dates when fertilization, pruning and thinning was 
# We set this to -1, indicating it did not happen 

calendar_fert  <- matrix( -1, nrow=100, ncol=3 )
calendar_prunT <- matrix( -1, nrow=100, ncol=3 )
calendar_thinT <- matrix( -1, nrow=100, ncol=3 )

# getting the default parameters of BASFOR which are saved under df_params 
# we do this for coniferous trees


params <- df_params[,1]


# bring the climate data into and appropriate shape 
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

NDAYS = nrow(clim)

# genearate starting values of the state variables and get weather values ready 

matrix_weather <- weather_BASFOR(as.integer(1920), as.integer(1), NDAYS, clim)
STATEVARS <- rep(0,7)


## Run model to set new default parameters 


output <-  run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
                         bias = as.integer(0), randerr = as.integer(0), 
                         ft = as.integer(1), p = params, w = matrix_weather,
                         calf = calendar_fert, calN = calendar_Ndep,
                         calpT = calendar_prunT, caltT = calendar_thinT, 
                         n = NDAYS, NOUT = 24, sv = STATEVARS, 
                         stateers = c(1,1,1), procerr = c(1,1,1,1,1))



## Change some initial values
params[57] <- TREEDENS0 <- output[NDAYS,4]
params[4] <- CRtree0 <- output[NDAYS,5] / TREEDENS0
params[2] <- CBtree0 <- output[NDAYS,7] / TREEDENS0
params[3] <- CLtree0 <- output[NDAYS,8] / TREEDENS0
params[5] <- CStree0 <- output[NDAYS,6] / TREEDENS0
params[33] <- CLITT0 <- output[NDAYS,10]
params[43] <- NMIN0 <- output[NDAYS,16]



## Thin data from the runs and prepare next run 

indLL <- seq(1,(NDAYS-forecastingtimestep+1), by = 100)   # Marks the starting point of LL calculation
indYears <- matrix_weather[indLL, 1]
indDays <- matrix_weather[indLL, 2]



weather_data <- weather_BASFOR(as.integer(indYears[1]), as.integer(indDays[1]), NDAYS, clim)


# actually run the model with the parameters 


output <- run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
                        bias = as.integer(0), randerr = as.integer(0), 
                        ft = as.integer(1), p = params, w = weather_data,
                        calf = calendar_fert, calN = calendar_Ndep,
                        calpT = calendar_prunT, caltT = calendar_thinT, n = NDAYS,
                        NOUT = 24, sv = STATEVARS, stateers = c(1,1,1), procerr = c(1,1,1,1,1))



# extract the reference to fit the model

refGPP <- output[,19] + rnorm(length(output[,19]), mean = 0, sd  = 0.2)



refET <- output[,21] + rnorm(length(output[,19]), mean = 0, sd = 0.2)