# Building and Installing the library

To install the library, you have to open a terminal within this folder (Crtl + T in this folder) and run the following in the Terminal: 
1. R CMD build BASFERROR 
2. R CMD install BASFERROR_0.0.1.tar.gz

After this the library is installed to R home folder. 

# Model Description 

The BASic FORest model, BASFOR, is a process-based model for forest biogeochemistry which simulates the C-, N- and water-cycles.
The original version can be found [here](https://github.com/davcam/BASFOR/blob/master/vignettes/Quickstart.Rmd)

## Differences to the original version 

There are some main differences to the original version: 
1. We implemented bias and stat. erorr into the model 
2. We implemented a restart option into the model 
3. We made the model able to be run as a state-space-model

This had the following consequences for the run_model function:
1. It is now called run_mod_model (for run modified model).  
2. There are additional arguments to the model call 
  +  rs: A Boolean 0 or 1, which says if the model should be restarted                       
  + statespace A Boolean 0 or 1, which says if the model should be run as statespace model
  + bias: A Boolean 0 or 1, which says if the model should be run with structural error 
  + randerr: A Boolean 0 or 1, which says if the model should be run with random process errors 
  + stateers: A vector, indicating the statistical errors ( c(1,1,1) would mean no error)
  + procerr: A vector, indicating the process errors ( c(1,1,1,1,1) would mean no error)
 
 The original model is thus a submodel of the modified model with every boolean set to zero.
