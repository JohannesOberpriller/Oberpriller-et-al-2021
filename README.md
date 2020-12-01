# Robust inference for complex system models 

This repository contains the code for the robust inference with complex models. <br/>

## First step is to compile the BASFERROR model. <br/>
For this you have to go to the models folder and type following commands into your terminal

1. R CMD build BASFERROR 
2. R CMD INSTALL BASFERROR_0.1.tar.gz (or if deviating another tar.gz called BASFERROR_x.x.tar.gz) 

## First step is to compile the BASFERROR model. <br/>
For this you have to go to the models folder and type following commands into your terminal

1. R CMD build BASFERROR (builds the model library for you)
2. R CMD INSTALL BASFERROR_0.1.tar.gz (installs the library to your R-path)

## Second step is to run the experiments 

The Code/ folder contains the following folders: 
1. data/ provides the data from the Profound data base
2. KOH containing the KOH calibration scripts 
3. SSM containg the state-space model calibration scripts
4. Weighting containting the calibration with weighting 
5. Standard_Calibration containing the standard calibration scripts 



Prior to the calibration you can run the SensitivityAnalysis.R to find the most 
sensitive parameters. 

After running all scripts the Problem_analysis.R script gives you the 
plots and the summary statistics contained in the paper.




## More documentation can be found in the individual folders for the individual methods 


