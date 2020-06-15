require(sensitivity)
require(BASFERROR)
source("analysis_functions.R")

# getting the default parameters of BASFOR which are saved under df_params 
# we do this for coniferous trees

params <- df_params[,1]

#### Define lower and upper bounds
lower <- params*0.7
upper <- params*1.3

# defining the timestep, runtime and year of simulation start 

forecastingtimestep = as.integer(100)
runtime = as.integer(95)
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

output <- run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
                        bias = as.integer(0), randerr = as.integer(0), 
                        ft = as.integer(1), p = params, w = matrix_weather,
                        calf = calendar_fert, calN = calendar_Ndep,
                        calpT = calendar_prunT, caltT = calendar_thinT, 
                        n = NDAYS, NOUT = 24, sv = STATEVARS, 
                        stateers = c(1,1,1), procerr = c(1,1,1,1,1))

str(output)

refGPP <- output[,19]
refET <- output[,21]

##### Sensitivity analysis based on sum of squares for GPP and ET

## define the LL which is the sum of squares of GPP and ET
morrisLL <- function(X){
  
  out <- numeric()
  for(i in 1:nrow(X)){
    params <- X[i,]
    model_output <- run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
                                  bias = as.integer(0), randerr = as.integer(0), 
                                  ft = as.integer(1), p = params, w = matrix_weather,
                                  calf = calendar_fert, calN = calendar_Ndep,
                                  calpT = calendar_prunT, caltT = calendar_thinT, 
                                  n = NDAYS, NOUT = 24, sv = STATEVARS, 
                                  stateers = c(1,1,1), procerr = c(1,1,1,1,1))
    
    out[i] <- sum((refGPP - model_output[,19])^2 + (refET - model_output[,21])^2)
    
    
  }
  
  ind <- which(is.na(out))
  out[ind] <- mean(out, na.rm = T)
  return(out)
}

# Define the oat sensitivity analysis 

one_at_a_time <-function(parsel, lower, upper,ETorGPP){
  step = (upper -lower)/10.
  
  pallette = rainbow(10)
  
  pdf(NULL)
  dev.control(displaylist="enable")
  
  
  params[parsel] = lower + step
  model_output <- run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
                                bias = as.integer(0), randerr = as.integer(0), 
                                ft = as.integer(1), p = params, w = matrix_weather,
                                calf = calendar_fert, calN = calendar_Ndep,
                                calpT = calendar_prunT, caltT = calendar_thinT, 
                                n = NDAYS, NOUT = 24, sv = STATEVARS, 
                                stateers = c(1,1,1), procerr = c(1,1,1,1,1))
  
  plot(model_output[(length(model_output[,ETorGPP])-300):length(model_output[,ETorGPP]),ETorGPP] ,type = "l", col = pallette[1], ylab = "", xlab = "")
  
  for(i in 2:10){
    params[parsel] = lower + step*i
    model_output <- run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
                                  bias = as.integer(0), randerr = as.integer(0), 
                                  ft = as.integer(1), p = params, w = matrix_weather,
                                  calf = calendar_fert, calN = calendar_Ndep,
                                  calpT = calendar_prunT, caltT = calendar_thinT, 
                                  n = NDAYS, NOUT = 24, sv = STATEVARS, 
                                  stateers = c(1,1,1), procerr = c(1,1,1,1,1))
    
    lines(model_output[(length(model_output[,ETorGPP])-300):length(model_output[,ETorGPP]),ETorGPP], add = T, col = pallette[i])
  
  }

  p1.base <- recordPlot()
  invisible(dev.off())

  return(p1.base)
  
}




#### Run the sensitivity analysis

sens <- morris(model = morrisLL, factors = length(params), r = 5000,
               design = list(type = "oat", levels = 5, grid.jump = 3),
               binf = lower, bsup = upper, scale = TRUE)
 

save(sens, file = "morrisSensitivity.RData")

mu.star <- apply(sens$ee, 2, function(x) mean(abs(x)))
sigma <- apply(sens$ee, 2, sd)

which.max(mu.star)
plot(mu.star[-parSel], sigma[-parSel])
points(mu.star[parSel], sigma[parSel], col = "red")

## Parameter selection 
parSel <- c(48, 47, 11, 34,51, 52, 54, 55, 56, 23, 34)

