#setwd("/home/johannes/Documents/PhD/Code/state-space-error/Code/Results")
set.seed(123)
library(scales)
library(BayesianTools)
library(BASFERROR)
library(miceadds)
library(kernlab)
library(fitdistrplus)
library(actuar)
library(kimisc)
source("./helper_functions.R")



isEmpty <- function(x) {
  # function to check if the value is zero 
  return(identical(x, numeric(0)))
}

post_predictiv_uncertainty_parallel <- function(filename, KOH = F, biastype = NULL, parallel = F){
  
  ## load the posterior samples 
  load(filename)
  
  # get sample from the posteriors, from each posterior get 120 samples 

  
  sam <- getSample(out_MCMC, numSamples = 120, start = 400)
  for(i in 2:5){
    main = substring(filename, 1, nchar(filename)-6)
    file = paste0(main,as.character(i),".RData")
  load(file)
  sami <- getSample(out_MCMC, numSamples = 120, start = 400)
  sam = rbind(sam, sami)
  }
  
  ## get the same set up as data for runs :
  

  # 1. set the forecasting timestep, the runtime and the start of the simulation 
  
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
  
  
  
  NDAYS = nrow(clim)
  
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
  
  ## Thin data from the runs and prepare next run 
  
  indLL <- seq(1,(NDAYS-(forecastingtimestep -1)), by = forecastingtimestep)   # Marks the starting point of LL calculation
  indYears <- matrix_weather[indLL, 1]
  indDays <- matrix_weather[indLL, 2]
  
  weather_data <- weather_BASFOR(as.integer(indYears[1]), as.integer(indDays[1]), NDAYS,
                                 clim)
  
  # Get empty matrices for the results 
  
  ET_matrix <- matrix(nrow = 7304, ncol = nrow(sam)+1)
  GPP_matrix <- matrix(nrow = 7304, ncol = nrow(sam)+1)
  ET_matrix_with_noise <- matrix(nrow = 7304, ncol = nrow(sam)+1)
  GPP_matrix_with_noise <- matrix(nrow = 7304, ncol = nrow(sam)+1)
  
  # Make the first run which is basically the truth for the rest 
  
  correct_run <- run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
                               bias = as.integer(0), randerr = as.integer(0), 
                               ft = as.integer(1), p = params, w = weather_data,
                               calf = calendar_fert, calN = calendar_Ndep,
                               calpT = calendar_prunT, caltT = calendar_thinT, 
                               n = NDAYS , NOUT = 24, sv = STATEVARS, 
                               stateers = c(1,1,1), procerr = c(1,1,1,1,1))
  
  
  ET_matrix[,1] =  correct_run[(nrow(correct_run)-7303):nrow(correct_run),21]
  GPP_matrix[,1] = correct_run[(nrow(correct_run)-7303):nrow(correct_run),19]
  
  ET_matrix_with_noise[,1] = correct_run[(nrow(correct_run)-7303):nrow(correct_run),21] + rnorm(n = 7304, mean = 0, sd = 0.2)
  GPP_matrix_with_noise[,1] = correct_run[(nrow(correct_run)-7303):nrow(correct_run),19]+ + rnorm(n = 7304, mean = 0, sd = 0.2)
  
  
  if(parallel == T){
    
    # set up cluster and export everything important 
    # in order to do the runs 
    
    cl <- parallel::makeCluster(48)
    packages2 <- (.packages())
    
    parallel::clusterCall(cl, packageFun, packages2)
    objects = c("sam", "params","biastype","indYears", 
                 "indDays", "NDAYS", "clim", "calendar_Ndep",
                "calendar_fert", "calendar_prunT", "calendar_thinT", "params")
    functions = c("make_single_run","get_n_unif_samples","get_unif_sample", "isEmpty")
    parallel::clusterExport(cl, varlist = functions)
    parallel::clusterExport(cl, varlist = objects, envir=environment())
    
  for(k in 0:4){
    # run data script
    if( k == 0){
      source("./GetReferenceData_add_noise.R")
    }
    else{
      filename = paste0("./GetReferenceData_add_noise",as.character(k+1),".R")
      source(filename)
    }
    # export stuff 
    objects_indiv = c("refET","refGPP")
    parallel::clusterExport(cl, varlist = objects_indiv, envir=environment())
    # run analysis 
    indiv_samples = (k*120+1):((k+1)*120)
    runs = parallel::parApply(cl = cl, FUN = make_single_run, X = sam[indiv_samples,], MARGIN = 1, KOH = KOH, biastype = biastype,
                 indYears = indYears, indDays = indDays, NDAYS = NDAYS, clim = clim,  calendar_Ndep = calendar_Ndep,
                 calendar_fert, calendar_prunT, calendar_thinT, params = params)
    #get it into correct shape
    counter = 1
    for(j in indiv_samples){
      ET_matrix[,j+1] = runs[[counter]][[1]]
      GPP_matrix[,j+1] = runs[[counter]][[2]]
      ET_matrix_with_noise[,j+1] = runs[[counter]][[3]]
      GPP_matrix_with_noise[,j+1] = runs[[counter]][[4]]
      counter = counter +1
    }
  }
    
  }
  else{
    # make the runs 
  
    runs = apply(FUN = make_single_run, X = sam, MARGIN = 1, KOH = KOH, biastype = biastype, 
               indYears = indYears, indDays = indDays, NDAYS = NDAYS, clim = clim,  calendar_Ndep = calendar_Ndep,
               calendar_fert = calendar_fert, calendar_prunT = calendar_prunT, calendar_thinT = calendar_thinT, params = params)
    #extract the results
    for(i in 1:length(runs)){
      ET_matrix[,i+1] = runs[[i]][[1]]
      GPP_matrix[,i+1] = runs[[i]][[2]]
      ET_matrix_with_noise[,i+1] = runs[[i]][[3]]
      GPP_matrix_with_noise[,i+1] = runs[[i]][[4]]
    }
  }
  
  # calculate minimum and maximum for results and return them 
  
  minET <- apply(FUN = min, X = ET_matrix_with_noise[,1:ncol(ET_matrix_with_noise)], MARGIN = 1)
  maxET <- apply(FUN = max, X = ET_matrix_with_noise[,1:ncol(ET_matrix_with_noise)], MARGIN = 1)
  minGPP <- apply(FUN = min, X = GPP_matrix_with_noise[,1:ncol(GPP_matrix_with_noise)], MARGIN = 1)
  maxGPP <- apply(FUN = max, X = GPP_matrix_with_noise[,1:ncol(GPP_matrix_with_noise)], MARGIN = 1)

  GPP_bounds = matrix(c(minGPP,maxGPP), nrow = 2, byrow = T, ncol = 7304)
  ET_bounds = matrix(c(minET,maxET), nrow = 2, byrow = T, ncol = 7304)

  
  return(list(ET_matrix,GPP_matrix, ET_bounds, GPP_bounds, ET_matrix_with_noise, GPP_matrix_with_noise, sam))
  
}


make_single_run <- function(sample, KOH = F, biastype = NULL, indYears, indDays, NDAYS, clim, calendar_Ndep, calendar_fert, calendar_prunT, calendar_thinT, params){

if(biastype == "func"){
  
    # set up the statevariables 
    STATEVARS <- rep(0,14)
    
    # parameters we calibrated 
    parSel <- c(48, 11, 34, 51, 52, 56)
    
    # swape them with the posterior samples 
    params[parSel] <- sample[1:6]
    
    # set up the weather in BASFOR form 
    
    weather_data <- weather_BASFOR(as.integer(indYears[1]), as.integer(indDays[1]), NDAYS, clim)
    
    # run the model 
    
    sample_runs <- run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
                                 bias = as.integer(1), randerr = as.integer(0), 
                                 ft = as.integer(1), p = params, w = weather_data,
                                 calf = calendar_fert, calN = calendar_Ndep,
                                 calpT = calendar_prunT, caltT = calendar_thinT, 
                                 n = NDAYS , NOUT = 24, sv = STATEVARS, 
                                 stateers = c(1,1,1), procerr = c(1,1,1,1,1))
    
    # if the KOH is true fit a posteriori an GP to it 
    if(KOH == T){
    
      # get the reference data 
      GPPdata <- refGPP
      ETdata <- refET
      
      # calculate the deviation of data and model 
      
      diffGPP <- GPPdata - sample_runs[1:length(GPPdata),19]
      diffET <- ETdata - sample_runs[1:length(ETdata),21]
      
      # divide into a calibration range and a prediction range 
      calib_range = (nrow(sample_runs)-7304):(length(GPPdata))
      predict_range = (nrow(sample_runs)-7304):(nrow(sample_runs))
      
      # get the prediction data for GPP from the calib range  
      GPP_sampled = cbind(sample_runs[calib_range,19], weather_data[calib_range,3:7])
      GPP_complet = cbind(sample_runs[,19], weather_data[1:length(sample_runs[,21]),3:7])
      
      # get the prediction data for ET from the calib range
      ET_sampled = cbind(sample_runs[calib_range,21], weather_data[calib_range,3:7])
      ET_complet = cbind(sample_runs[,21], weather_data[1:length(sample_runs[,21]),3:7])
      
      # fit the gaussian processes 
      fit_GPP = gausspr(x = GPP_sampled , y = diffGPP[calib_range], kernel = "rbfdot", variance.model = T)
      fit_ET = gausspr(x = ET_sampled,  y = diffET[calib_range], kernel = "rbfdot",variance.model = T)
      
      # predict with the for all data points
      predict_GPP = predict(fit_GPP, GPP_complet)
      predict_ET = predict(fit_ET, ET_complet)
      
      # get the predicted results and return them
      sample_runs[,19] = predict_GPP + sample_runs[,19]
      sample_runs[,21] = predict_ET + sample_runs[,21]
      
      
    }
    
    # put it together with the other results with noise and without noise 
    
    GPP_matrix = sample_runs[(nrow(sample_runs)-7303):nrow(sample_runs),19] 
    ET_matrix = sample_runs[(nrow(sample_runs)-7303):nrow(sample_runs),21]
    
    ET_matrix_with_noise = ET_matrix + rnorm(length(ET_matrix), mean = 0, sd = sample[8])
    GPP_matrix_with_noise = GPP_matrix + rnorm(length(GPP_matrix), mean = 0, sd = sample[7])
}
  

else{
  # set up the statevariables 
    STATEVARS <- rep(0,14)
    
    # parameters we calibrated 
    parSel <- c(48, 11, 34, 51, 52, 56)
    
    # swape them with the posterior samples 
    params[parSel] <- sample[1:6]
    
    # set up the weather in BASFOR form 
    
    weather_data <- weather_BASFOR(as.integer(indYears[1]), as.integer(indDays[1]), NDAYS, clim)
    
    # run the model 
    sample_runs <- run_mod_model(rs = as.integer(0), statespace = as.integer(0), 
                                 bias = as.integer(0), randerr = as.integer(0), 
                                 ft = as.integer(1), p = params, w = weather_data,
                                 calf = calendar_fert, calN = calendar_Ndep,
                                 calpT = calendar_prunT, caltT = calendar_thinT, 
                                 n = NDAYS , NOUT = 24, sv = STATEVARS, 
                                 stateers = c(1,1,1), procerr = c(1,1,1,1,1))
    
    # if the KOH is true fit a posteriori an GP to it 
    if(KOH == T){
      
      # get the reference data 
      GPPdata <- refGPP
      ETdata <- refET
      
      # calculate the deviation of data and model 
      diffGPP <- GPPdata - sample_runs[1:length(GPPdata),19]
      diffET <- ETdata - sample_runs[1:length(ETdata),21]
  
      # divide into a calibration range and a prediction range 
      calib_range = (nrow(sample_runs)-7304):(length(GPPdata))
      predict_range = (nrow(sample_runs)-7304):(nrow(sample_runs))
      
      # get the prediction data for GPP from the calib range  
      GPP_sampled = cbind(sample_runs[calib_range,19], weather_data[calib_range,1:7])
      GPP_complet = cbind(sample_runs[,19], weather_data[1:length(sample_runs[,21]),1:7])
      
      # get the prediction data for ET from the calib range
      ET_sampled = cbind(sample_runs[calib_range,21], weather_data[calib_range,1:7])
      ET_complet = cbind(sample_runs[,21], weather_data[1:length(sample_runs[,21]),1:7])

      # fit the gaussian processes
      fit_GPP = gausspr(x = GPP_sampled , y = diffGPP[calib_range],kernel = "rbfdot")
      fit_ET = gausspr(x = ET_sampled,  y = diffET[calib_range], kernel = "rbfdot")

      # predict with the for all data points
      predict_GPP = predict(fit_GPP, GPP_complet)
      predict_ET = predict(fit_ET, ET_complet)

      # get the predicted results and return them
      sample_runs[,19] = predict_GPP + sample_runs[,19]
      sample_runs[,21] = predict_ET + sample_runs[,21]
    }
    
    # put it together with the other results with noise and without noise     
    GPP_matrix = sample_runs[(nrow(sample_runs)-7303):nrow(sample_runs),19] 
    ET_matrix = sample_runs[(nrow(sample_runs)-7303):nrow(sample_runs),21]
    
    ET_matrix_with_noise = ET_matrix + rnorm(length(ET_matrix), mean = 0, sd = sample[8])
    GPP_matrix_with_noise = GPP_matrix + rnorm(length(GPP_matrix), mean = 0, sd = sample[7])
  
}

  return(list(ET_matrix, GPP_matrix, ET_matrix_with_noise, GPP_matrix_with_noise))

}



sd_deviation_function <- function(samples){
  
  ## Calculates mean and standard deviation of calibrated and predicted per time step
  
  sd_GPP_cal = apply(FUN = sd, X = samples[[2]][1:3652,2:ncol(samples[[2]])], 
                                   MARGIN = 1, na.rm = T)
  
  mean_GPP_cal = apply(FUN = mean, X = samples[[2]][1:3652,2:ncol(samples[[2]])], 
                     MARGIN = 1, na.rm = T)
  
  sd_GPP_uncal = apply(FUN = sd, X = samples[[2]][3653:nrow(samples[[2]]),2:ncol(samples[[2]])], 
                                     MARGIN = 1, na.rm = T)
  
  mean_GPP_uncal = apply(FUN = mean, X = samples[[2]][3653:nrow(samples[[2]]),2:ncol(samples[[2]])], 
                       MARGIN = 1, na.rm = T)
  
  sd_ET_cal = apply(FUN = sd, X = samples[[1]][1:3652,2:ncol(samples[[1]])],MARGIN = 1,  na.rm = T)
  
  mean_ET_cal = apply(FUN = mean, X = samples[[1]][1:3652,2:ncol(samples[[1]])],MARGIN = 1,  na.rm = T)
  
  sd_ET_uncal = apply(FUN = sd, X = samples[[1]][3653:nrow(samples[[1]]),2:ncol(samples[[1]])], MARGIN = 1, na.rm = T)
  
  mean_ET_uncal = apply(FUN = mean, X = samples[[1]][3653:nrow(samples[[1]]),2:ncol(samples[[1]])], MARGIN = 1, na.rm = T)
  
  ## sorts out all values, where the standard deviation is very low ## 
  
  sd_0_GPP = which(sd_GPP_cal!=0)
  sd_0_GPP_un = which(sd_GPP_uncal!=0)
  sd_0_ET = which(sd_ET_cal!=0)
  sd_0_ET_un = which(sd_ET_uncal!=0)
  
  
  ## calculates the difference in units of standard deviation 
  
  diffGPP = mean(abs(mean_GPP_cal[sd_0_GPP] - samples[[2]][1:3652,1][sd_0_GPP])/sd_GPP_cal[sd_0_GPP])
  diffGPPUn = mean(abs(mean_GPP_uncal[sd_0_GPP_un] - samples[[2]][3653:nrow(samples[[2]]),1][sd_0_GPP_un])/sd_GPP_uncal[sd_0_GPP_un])
  diffET = mean(abs(mean_ET_cal[sd_0_ET] - samples[[1]][1:3652,1][sd_0_ET])/sd_ET_cal[sd_0_ET])
  diffETUn = mean(abs(mean_ET_uncal[sd_0_ET_un] - samples[[1]][3653:nrow(samples[[1]]),1][sd_0_ET_un] )/sd_ET_uncal[sd_0_ET_un] )
  
  
  
  scores_complete = cbind(diffGPP,diffGPPUn,diffET, diffETUn)
  
  return(scores_complete)
  
}



calculate_errors <- function(samples){
  # function to calculate the errors of the time sereis predictions 
  # calculates the mean of the absolute deviations 
  
  # set up matrix 
  errors <- matrix( nrow = 2, ncol = 2)
  rownames(errors) <- c("GPP","ET")
  colnames(errors) <- c("calibrated","uncalibrated")
  uncalibrated_GPP_error = vector(mode = "numeric", length = 3652)
  calibrated_GPP_error = vector(mode = "numeric", length = 3652)
  uncalibrated_ET_error = vector(mode = "numeric", length = 3652)
  calibrated_ET_error = vector(mode = "numeric", length = 3652)
  
  # get the values 
  
  calibrated_GPP_error = apply(FUN = mean, X = (samples[[2]][1:3652,2:ncol(samples[[2]])] -  samples[[2]][1:3652,1]), MARGIN = 1) 
  calibrated_ET_error =  apply(FUN = mean, X = (samples[[1]][1:3652,2:ncol(samples[[2]])] -  samples[[1]][1:3652,1]), MARGIN = 1) 
  uncalibrated_GPP_error = apply(FUN = mean, X = (samples[[2]][3653:nrow(samples[[2]]),2:ncol(samples[[2]])] - samples[[2]][3653:nrow(samples[[2]]),1]), MARGIN = 1) 
  uncalibrated_ET_error = apply(FUN = mean, X = (samples[[1]][3653:nrow(samples[[2]]),2:ncol(samples[[2]])] - samples[[1]][3653:nrow(samples[[2]]),1]), MARGIN = 1) 
  
  # build the mean abs of the time series 
  
  errors[1,1] = (mean(abs(calibrated_GPP_error)))
  errors[1,2] = (mean(abs(uncalibrated_GPP_error)))
  errors[2,2] = (mean(abs(uncalibrated_ET_error)))
  errors[2,1] = (mean(abs(calibrated_ET_error)))
  return(errors)
}



calculate_parameter_error <- function(filename){
  # load data 
  load(filename)
  
  # get data from the posterior 
  sam <- getSample(out_MCMC, numSamples = 10000)
  for(i in 2:5){
    main = substring(filename, 1, nchar(filename)-6)
    file = paste0(main,as.character(i),".RData")
    load(file)
    sami <- getSample(out_MCMC, numSamples = 10000)
    sam = rbind(sam, sami)
  }
  
  # get default parameter 
  params <- df_params[,1]
  # get parameters we calibrated 
  parSel <- c(48, 11, 34,51, 52, 56)
  best <- c(params[parSel], 0.2, 0.2)
  parameter_error = vector( mode = "numeric", length = ncol(sam))
  for(i in 1:6){
    parameter_error[i] = mean((sam[,i] - best[i])/best[i])
  }
  # calculate the absolute mean of the weighted parameter deviations s 
  return(abs(mean(parameter_error)))
}


packageFun <- function(packages = NULL, dlls = NULL) {
  if(!is.null(packages)){
    for(i in packages) library(i, character.only = TRUE)
  }
  if(!is.null(dlls)){
    for(i in dlls) try(dyn.load(i), silent = T)
  }
}