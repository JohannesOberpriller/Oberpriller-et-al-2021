# Applying the KOH to BASFOR

This folder contains the code to apply the KOH method with our modifications to the BAFOR model: 
+ Problem_KOH_add_noise.R : The KOH method applied to the "true" model.
+ Problem_func_KOH_add_noise.R : The KOH method applied to the model with structural error. 

In the code we fit a Gaussian process using the kernlab package to the resdiuals of the data and model in every MCMC step.
Because of the huge amount of data, we subsample the data to the most informative data points ( around 10 % of the data) for faster fitting. 
Moreover to avoid the inversion of the covariance matrix, we do only fit one parameter 
Let *GP* be the Gaussian process and m<sub>ET</sub> and m<sub>GPP</sub> be the model predictions for GPP and ET, d<sub>ET</sub> and d<sub>GPP</sub> the data for ET and GPP. 
Then in every MCMC step : 
+ residuals<sub>ET</sub> = d<sub>ET</sub> - m<sub>ET</sub>
+ subsampled_predictors = c(d<sub>ET</sub>, climate)[informative_data_points]
+ GP<sub>ET</sub> = fitGP(y =  residuals<sub>ET</sub>, x = subsampled predictors)

+ residuals<sub>GPP</sub> = d<sub>GPP</sub> - m<sub>GPP</sub>
+ subsampled_predictors = c(d<sub>GPP</sub>, climate)[informative_data_points]
+ GP<sub>GPP</sub> = fitGP(y =  residuals<sub>GPP</sub>, x = subsampled predictors)

Then we extraplote to the full model calibration data 
+ GP_full<sub>ET</sub> = predict(GP<sub>ET</sub>, c(d<sub>ET</sub>, climate))
+ GP_full<sub>GPP</sub> = predict(GP<sub>GPP</sub>, c(d<sub>GPP</sub>, climate))

Then, we calculate the likelihood <br>
*Scripts should be run via the source button in Rstudio and can not be run within jobs, as working directory changes*
