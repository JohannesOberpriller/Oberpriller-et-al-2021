## The make calendar Ndep function takes as input a calendar with the climate data, 
## the simulation start and the runtime and trasnforms this to be used in BASFERROR

make_calendar_Ndep <- function(calendar_Ndep, climate, simulation_start, runtime){
  
  # check if there are less rows than 98, which is somehow the maximum of basfor  
  nrows <- nrow(calendar_Ndep)
  if(nrows > 98){
    message("Can not use all the N_dep information because of time limit in BASFOR")
  }
  
  # get first year and last year from runtime and simulation_start 
  
  firstrow <- which(calendar_Ndep[,1] <= simulation_start)[[length(which(calendar_Ndep[,1] < simulation_start))]]
  lastrow <- which(calendar_Ndep[,1] >= simulation_start + runtime)[[1]]
  calendar_Ndep <- calendar_Ndep[firstrow:lastrow,]
  nrows <- nrow(calendar_Ndep)
  
  # checking what is in the climate data  
  
  new_calendar_Ndep = matrix(nrow = 100,ncol = 3)
  
  maxyear_climate = max(climate[,1])
  maxyear_Ndep = max(calendar_Ndep[,1])
  
  minyear_climate = min(climate[,1])
  minyear_Ndep = min(calendar_Ndep[,1])

  # arrange the calenders given the content of the packages 
  
  if((maxyear_climate >= maxyear_Ndep) & (minyear_climate <= minyear_Ndep)){
    message("N_deposition data has to end one year later and start one year earlier than climate data!! Automatic adjustment")
    new_calendar_Ndep[1,] = calendar_Ndep[1,]
    new_calendar_Ndep[1,1] = -1
    new_calendar_Ndep[2:(nrows+1),] = calendar_Ndep
    new_calendar_Ndep[nrows + 2,] = new_calendar_Ndep[nrows,]
    new_calendar_Ndep[nrows + 2,1] = new_calendar_Ndep[nrows + 2,1] + 1
    for(i in (nrows+3):100){
      new_calendar_Ndep[i,] = c(-1,-1,-1)
    }
  }
  
  else if(minyear_climate <= minyear_Ndep){
    message("N_deposition data has to start one year earlier than climate data! Automatic adjustment!")
    new_calendar_Ndep[1,] = calendar_Ndep[1,]
    new_calendar_Ndep[1,1] = -1
    new_calendar_Ndep[2:(nrows+1),] = calendar_Ndep
    for(i in (nrows+2):100){
      new_calendar_Ndep[i,] = c(-1,-1,-1)
    }
  }
  
  
  else if(maxyear_climate >= maxyear_Ndep){
    message("N_deposition data has to end one year later than climate data! Automatic adjustment!")
    new_calendar_Ndep[1:(nrows),] = calendar_Ndep
    new_calendar_Ndep[nrows + 1, ] = new_calendar_Ndep[nrows,]
    new_calendar_Ndep[nrows + 1,1] = new_calendar_Ndep[nrows + 1,1] + 1
    for(i in (nrows+2):100){
      new_calendar_Ndep[i,] = c(-1,-1,-1)
    }
  }
  
  else{
    new_calendar_Ndep[1:(nrows),] = calendar_Ndep
    for(i in (nrows+1):100){
      new_calendar_Ndep[i,] = c(-1,-1,-1)
    }
  }
  return(new_calendar_Ndep)
}

## The prepare_data function wraps everything up in order to 
## make the data useable for BASFERROR


prepare_data <- function(calendar_Ndep, climate, simulation_start, runtime){
  # check if there are less rows than 100, which is somehow the maximum of basfor  
  if(runtime > 100){
    stop("Actual BASFOR implementation does not support a simulation with more than 100 years")
  }
  new_climate <- climate[which((climate[,1] > simulation_start) & (climate[,1] < simulation_start + runtime)),]
  calendar_Ndep <- make_calendar_Ndep(calendar_Ndep, new_climate, simulation_start, runtime)
  return(list(new_climate, calendar_Ndep))
}

