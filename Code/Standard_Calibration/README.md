# Standard calibration of the model 
In this folder, we perform standard calibration for the model. 
This means, we optimize the likelihood given by 
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;L(\theta)&space;=&space;\frac{1}{\sqrt{2\pi&space;\sigma^2}}&space;exp[-\frac{1}{2\sigma^2}(d_{ET}-m(\theta)_{ET})^2]&plus;&space;\frac{1}{\sqrt{2\pi&space;\sigma^2}}&space;exp[-\frac{1}{2\sigma^2}(d_{GPP}-m(\theta)_{GPP})^2]" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;L(\theta)&space;=&space;\frac{1}{\sqrt{2\pi&space;\sigma^2}}&space;exp[-\frac{1}{2\sigma^2}(d_{ET}-m(\theta)_{ET})^2]&plus;&space;\frac{1}{\sqrt{2\pi&space;\sigma^2}}&space;exp[-\frac{1}{2\sigma^2}(d_{GPP}-m(\theta)_{GPP})^2]" title="L(\theta) = \frac{1}{\sqrt{2\pi \sigma^2}} exp[-\frac{1}{2\sigma^2}(d_{ET}-m(\theta)_{ET})^2]+ \frac{1}{\sqrt{2\pi \sigma^2}} exp[-\frac{1}{2\sigma^2}(d_{GPP}-m(\theta)_{GPP})^2]" /></a>
to get posterior estimates of the parameters <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\theta" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\theta" title="\theta" /></a> 
given the data d and the model predictions m.

# Optimisation 

We first sample values from the priors and optimize the likelihood function with them 
Then we use these optimized values to feed them into the DEZs-MCMC sampler as start values as well as Z-matrix. 
The sampling is done with the BayesianTools package.

This general pattern is applied to four different cases: 

+ Model_validation_add_noise.R : The "true" model with balanced data 
+ Problem_func_biased_add_noise.R : The model with structural error with balanced data 
+ Problem_unbal_biased_add_noise.R : The "true" model with unbalanced data 
+ Problem_unbal_func_biased_add_noise.R : The model with structural error with unbalanced data <br>
*Scripts should be run via the source button in Rstudio and can not be run within jobs, as working directory changes*
