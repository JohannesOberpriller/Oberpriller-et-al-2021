# Folder for running the simulations in Oberpriller et al 2020

This contains the following subfolders: 
+ KOH : Simulations for the KOH method 
+ SSM : Simulations for the State-space model 
+ Standard_Calibration: Simulations for the standard calibration
+ Weighting  : Simulations for weighting of data 
+ data : Env. data to run the model
+ KOH_original: Code to run the original Kennedy and O'Hagan method 

Additionally there are scripts used to prepare and analyse the results:

- GetReferenceData_add_noise.R: Produces reference data 
- Problem_analysis.R: Produces the results 
- SensitivityAnalysis.R: Identifies the most sensitive parameters for fitting 
- analysis_functions/helper_functions.R: functions used in Problem_analysis.R 
- justify_approximation.R: a simulation scirpt, for justification of the approximation for the likelihood. 
- check_error_distribution.R: sample script to check if statistical or structural error are dominating

More detailed information can be found in the subfolders 
*Scripts should be run via the source button in Rstudio and can not be run within jobs, as working directory changes*

