## Anaylsis of the Problematic functions 
set.seed(123)
source("./analysis_functions.R")
require(BayesianTools)
require(BASFERROR)
library(grid)
library(gridExtra)
library(scales)
library(ggplot2)
library(cowplot)

### setting the default parameters ###

params = df_params[,1]

#### Consequences of model errror ####

### Correct model, balanced data ###

post_model_val <- post_predictiv_uncertainty_parallel("./Results/Model_validation_chain_less_parameter_add_noise.RData", biastype = "none")

predictive_error_Model_val <- calculate_errors(post_model_val)

uncert_Model_val <- sd_deviation_function(post_model_val)

parameter_error_Model_val <- calculate_parameter_error("./Results/Model_validation_chain_less_parameter_add_noise.RData")


### Structural model error, balanced data ###

post_func_bias <- post_predictiv_uncertainty_parallel("./Results/MCMC_func_biased_add_noise.RData", biastype = "func")

predictive_error_func_bias <- calculate_errors(post_func_bias)

parameter_error_func_bias <- calculate_parameter_error("./Results/MCMC_func_biased_add_noise.RData")

uncert_func_bias <- sd_deviation_function(post_func_bias)


### Correct model, unbalanced data ### 

post_unbal_bias <- post_predictiv_uncertainty_parallel("./Results/MCMC_unbal_biased_add_noise.RData", biastype = "none")

predictive_error_unbal_bias <- calculate_errors(post_unbal_bias)

parameter_error_unbal_bias <- calculate_parameter_error("./Results/MCMC_unbal_biased_add_noise.RData")

uncert_unbal_bias <- sd_deviation_function(post_unbal_bias)


### Structural model error, balanced data ###

post_unbal_func_bias <- post_predictiv_uncertainty_parallel("./Results/MCMC_unbal_func_biased_add_noise.RData", biastype = "func")

predictive_error_unbal_func_bias <- calculate_errors(post_unbal_func_bias)

parameter_error_unbal_func_bias <- calculate_parameter_error("./Results/MCMC_unbal_func_biased_add_noise.RData")

uncert_unbal_func_bias <- sd_deviation_function(post_unbal_func_bias)


### Setting up the matrix for storage ###

results_overview <- matrix(ncol = 5, nrow = 4)

results_overview_uncert <- matrix(ncol = 4, nrow = 4)



colnames(results_overview) <- c("Parameter error", "Calibrated error GPP", "Uncalibrated error GPP", "Calibrated error ET", "Uncalibrated error ET")
rownames(results_overview) <- rep(c("Model Validation", "Structural Bias"),2)  


results_overview[1,1] <- parameter_error_Model_val 
results_overview[1,2] <- predictive_error_Model_val[1,1]
results_overview[1,3] <- predictive_error_Model_val[1,2]
results_overview[1,4] <- predictive_error_Model_val[2,1]
results_overview[1,5] <- predictive_error_Model_val[2,2]

results_overview_uncert[1,] <- uncert_Model_val

results_overview[2,1] <- parameter_error_func_bias
results_overview[2,2] <- predictive_error_func_bias[1,1]
results_overview[2,3] <- predictive_error_func_bias[1,2]
results_overview[2,4] <- predictive_error_func_bias[2,1]
results_overview[2,5] <- predictive_error_func_bias[2,2]


results_overview_uncert[2,] <- uncert_func_bias

results_overview[3,1] <- parameter_error_unbal_bias 
results_overview[3,2] <- predictive_error_unbal_bias[1,1]
results_overview[3,3] <- predictive_error_unbal_bias[1,2]
results_overview[3,4] <- predictive_error_unbal_bias[2,1]
results_overview[3,5] <- predictive_error_unbal_bias[2,2]

results_overview_uncert[3,] <- uncert_unbal_bias



results_overview[4,1] <- parameter_error_unbal_func_bias
results_overview[4,2] <- predictive_error_unbal_func_bias[1,1]
results_overview[4,3] <- predictive_error_unbal_func_bias[1,2]
results_overview[4,4] <- predictive_error_unbal_func_bias[2,1]
results_overview[4,5] <- predictive_error_unbal_func_bias[2,2]

results_overview_uncert[4,] <- uncert_unbal_func_bias

balanced_error = results_overview[c(1,2),1]
unbalanced_error = results_overview[c(3,4),1]
problem = c(rep("Correct Model",2), rep("Systematic Model Error",2))
Data = rep(c("Balanced Data","Unbalanced Data"),2)
ParameterError = results_overview[c(1,3,2,4),1]


data_problem = data.frame(problem, Data, ParameterError)

parameters_plot_paper <- ggplot() + 
  geom_bar(data = data_problem, aes(fill = Data, y=ParameterError*100, x=problem), position="dodge", stat="identity", color = "black") +
  labs(y = element_blank()) +
  theme(text = element_text(), legend.text=element_text()) +
  theme(axis.text.x = element_text(angle = 90 ), legend.position = c(0.9,0.7), axis.title.y=element_blank(),
        legend.direction = "vertical", legend.title = element_blank()) +  coord_flip(expand = c(0, 0)) + theme(aspect.ratio = 0.1) +
  theme(axis.text.x = element_text(angle = 0)) + ylim(c(0, max(ParameterError*100)+5)) + 
  scale_fill_manual(values = rep(c("Black", "White"),2),guide=guide_legend(reverse=T)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

## Getting data in shape for plots

GPPCalibrated = results_overview[c(1,3,2,4),2]
GPPuncert = round(results_overview_uncert[c(1,3,2,4),1],1)
GPPCalibrated_problem = data.frame(problem, Data, GPPCalibrated, GPPuncert)


GPPUnCalibrated = results_overview[c(1,3,2,4),3]
GPPUnuncert = round(results_overview_uncert[c(1,3,2,4),2],1)
GPPUnCalibrated_problem = data.frame(problem, Data, GPPUnCalibrated, GPPUnuncert)


ETCalibrated = results_overview[c(1,3,2,4),4]
ETUncert = round(results_overview_uncert[c(1,3,2,4),3],1)
ETCalibrated_problem = data.frame(problem, Data, ETCalibrated, ETUncert)

ETUnCalibrated = results_overview[c(1,3,2,4),5]
ETUnuncert = round(results_overview_uncert[c(1,3,2,4),4],1)
ETUnCalibrated_problem = data.frame(problem, Data, ETUnCalibrated,ETUnuncert)

## Plotting GPP

GPPplot <- ggplot(GPPCalibrated_problem, aes(fill = Data, y=GPPCalibrated, x=problem), stat = "identity") + 
  geom_bar(position="dodge", stat="identity", color = "black") + labs(y = element_blank()) +
  theme( axis.title.y=element_blank(), legend.position = "none", legend.title = element_blank()) +
  coord_flip(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(0.,0.2,0.4,0.6,0.8), limits = c(0., 0.82), 
                     labels = c(0,0.2,0.4,0.6,0.8)) +
  scale_fill_manual(values = rep(c("Black", "White"),2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(aes(label=GPPuncert), position=position_dodge(width=0.9), hjust= - .5)



GPPUnplot <- ggplot(GPPUnCalibrated_problem, aes(fill = Data, y=GPPUnCalibrated, x=problem)) + 
  geom_bar(position="dodge", stat="identity", color ="black")  + labs(y = element_blank()) +
  theme( legend.position = "none", axis.title.y=element_blank()) +  coord_flip(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(0.,0.2,0.4,0.6,0.8), limits = c(0., 0.82), 
                     labels = c(0,0.2,0.4,0.6,0.8)) +
  scale_fill_manual(values = rep(c("Black", "White"),2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text(aes(label=GPPUnuncert), position=position_dodge(width=0.9), hjust= - 0.5)

## Plotting ET

ETplot <- ggplot(ETCalibrated_problem, aes(fill = Data, y=ETCalibrated, x=problem)) + 
  geom_bar(position="dodge", stat="identity", color = "black")  + labs(y = element_blank()) +
  theme(axis.text.y = element_blank(), legend.position = "none", axis.title.y=element_blank(), legend.title = element_blank()) +  coord_flip(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(0.,0.05,0.10, 0.15,0.20), limits = c(-0., 0.21), 
                     labels =c(0.,0.05,0.10, 0.15,0.20)) +
  scale_fill_manual(values = rep(c("Black", "White"),2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(aes(label=ETUncert), position=position_dodge(width=0.9), hjust= - .5)


ETUnplot <- ggplot(ETUnCalibrated_problem, aes(fill = Data, y=ETUnCalibrated, x=problem)) + 
  geom_bar(position="dodge", stat="identity", color = "black") + labs(y = element_blank()) +
  theme(axis.text.y = element_blank(), legend.position = "None",
        legend.direction = "vertical", axis.title.y=element_blank()) +  coord_flip(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(0.,0.05,0.10, 0.15,0.20), limits = c(-0., 0.21), 
                     labels =c(0.,0.05,0.10, 0.15,0.20)) +
  scale_fill_manual(values = rep(c("Black", "White"),3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text(aes(label=ETUnuncert), position=position_dodge(width=0.9), hjust= - .5)


## Combining the plots 

Complete_overview_paper <- grid.arrange(arrangeGrob(arrangeGrob(GPPplot, top = "                                             GPP", left = "Calibration" ) ,
                                              arrangeGrob(GPPUnplot, left = "Prediction" ), heights = c(0.5,0.5)),
                                  arrangeGrob(arrangeGrob(ETplot, top = "ET"), arrangeGrob(ETUnplot)), ncol = 2,
                                  widths = c(1.3,1))


Completest_overview_paper <- plot_grid(arrangeGrob(parameters_plot_paper, left = ""),
                                          arrangeGrob(Complete_overview_paper),
                                          rel_heights = c(0.3,0.7), labels = c("a)","b)"), align = "h", ncol = 1)


#### Balancing the dataweights  ####

## Correct model,  unbalanced data 

unbal_model_weighted <- post_predictiv_uncertainty_parallel("./Results/MCMC_unbal_biased_add_noise_weighted.RData", biastype = "none")

predictive_unbal_model_weighted2 <- calculate_errors(unbal_model_weighted)

uncert_unbal_model_weighted <- sd_deviation_function(unbal_model_weighted)

parameter_unbal_model_weighted <- calculate_parameter_error("./Results/MCMC_unbal_biased_add_noise_weighted.RData")

## Structural model error, unbalanced data 

unbal_func_bias_weighted <- post_predictiv_uncertainty_parallel("./Results/MCMC_unbal_func_biased_add_noise_weighted.RData", biastype = "func")

predictive_unbal_func_bias_weighted2 <- calculate_errors(unbal_func_bias_weighted)

uncert_unbal_func_bias_weighted <- sd_deviation_function(unbal_func_bias_weighted)

parameter_unbal_func_bias_weighted <- calculate_parameter_error("./Results/MCMC_unbal_func_biased_add_noise_weighted.RData")

### Filling the matrix to store results ###

results_weighted_overview <- matrix(ncol = 5, nrow = 2)

results_weighted_overview_uncert <- matrix(ncol = 4, nrow = 2)

colnames(results_weighted_overview) <- c("Parameter error", "Calibrated error GPP", "Uncalibrated error GPP", "Calibrated error ET", "Uncalibrated error ET")
rownames(results_weighted_overview) <- c("Model Validation", "Structural Bias")

results_weighted_overview[1,1] <- parameter_unbal_model_weighted 
results_weighted_overview[1,2] <- predictive_unbal_model_weighted2[1,1]
results_weighted_overview[1,3] <- predictive_unbal_model_weighted2[1,2]
results_weighted_overview[1,4] <- predictive_unbal_model_weighted2[2,1]
results_weighted_overview[1,5] <- predictive_unbal_model_weighted2[2,2]

results_weighted_overview_uncert[1,] <- uncert_unbal_model_weighted


results_weighted_overview[2,1] <- parameter_unbal_func_bias_weighted
results_weighted_overview[2,2] <- predictive_unbal_func_bias_weighted2[1,1]
results_weighted_overview[2,3] <- predictive_unbal_func_bias_weighted2[1,2]
results_weighted_overview[2,4] <- predictive_unbal_func_bias_weighted2[2,1]
results_weighted_overview[2,5] <- predictive_unbal_func_bias_weighted2[2,2]

results_weighted_overview_uncert[2,] <- uncert_unbal_func_bias_weighted


### Generating plots for the weighted cases ###


unbalanced_error_unweighted = results_overview[c(3,4),1]
unbalanced_error_weighted = results_weighted_overview[c(1,2),1]
problem = c(rep("Correct Model",2), rep("Systematic Model Error",2))
Data_weighted = rep(c("Unbalanced Data","Weighted Data"),2)
ParameterError_weighted = vector(length = 4)
for(i in 1:2){
  ParameterError_weighted[(i-1)*2+1] = unbalanced_error_unweighted[i]
  ParameterError_weighted[(i)*2] = unbalanced_error_weighted[i]
}



data_problem_weighted = data.frame(problem, Data_weighted, ParameterError_weighted)


parameters_plot_weighted_paper <- ggplot() + 
  geom_bar(data = data_problem_weighted, aes(fill = Data_weighted, y=ParameterError_weighted*100, x=problem),
           position="dodge", stat="identity", color = "black") +
  labs(y = element_blank()) +
  theme(text = element_text(), legend.text=element_text()) +
  theme(axis.text.x = element_text(angle = 90 ), legend.position = c(0.9,0.7), axis.title.y=element_blank(), line = element_blank(),
        legend.direction = "vertical", legend.title = element_blank()) +  coord_flip(expand = c(0, 0)) + theme(aspect.ratio = 0.1) +
  theme(axis.text.x = element_text(angle = 0)) + ylim(c(0, max(ParameterError_weighted*100)+5))+
  scale_fill_manual(values = rep(c("White", "Darkgrey"),3),guide=guide_legend(reverse=T)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

## Reshaping of GPP Data for Plots


GPP_error_unweighted = results_overview[c(3,4),2]
GPP_error_weighted = results_weighted_overview[,2]
GPPCalibrated_weighted = vector(length = 4)
for(i in 1:2){
  GPPCalibrated_weighted[(i-1)*2+1] = GPP_error_unweighted[i]
  GPPCalibrated_weighted[(i)*2] = GPP_error_weighted[i]
}


GPPUn_error_unweighted = results_overview[c(3,4),3]
GPPUn_error_weighted = results_weighted_overview[,3]
GPPUnCalibrated_weighted = vector(length =4)
for(i in 1:2){
  GPPUnCalibrated_weighted[(i-1)*2+1] = GPPUn_error_unweighted[i]
  GPPUnCalibrated_weighted[(i)*2] = GPPUn_error_weighted[i]
}

GPP_error_unweighted_uncert = results_overview_uncert[c(3,4),1]
GPP_error_weighted_uncert = results_weighted_overview_uncert[,1]
GPPCalibrated_weighted_uncert = vector(length = 4)
for(i in 1:2){
  GPPCalibrated_weighted_uncert[(i-1)*2+1] = GPP_error_unweighted_uncert[i]
  GPPCalibrated_weighted_uncert[(i)*2] = GPP_error_weighted_uncert[i]
}


GPPUn_error_unweighted_uncert = results_overview_uncert[c(3,4),2]
GPPUn_error_weighted_uncert = results_weighted_overview_uncert[,2]
GPPUnCalibrated_weighted_uncert = vector(length = 4)
for(i in 1:2){
  GPPUnCalibrated_weighted_uncert[(i-1)*2+1] = GPPUn_error_unweighted_uncert[i]
  GPPUnCalibrated_weighted_uncert[(i)*2] = GPPUn_error_weighted_uncert[i]
}



GPPuncert_weighted = round(GPPCalibrated_weighted_uncert,1)
GPPCalibrated_problem_weighted = data.frame(problem, Data_weighted, GPPCalibrated_weighted, GPPuncert_weighted)

GPPUnuncert_weighted = round(GPPUnCalibrated_weighted_uncert,1)
GPPUnCalibrated_problem_weighted = data.frame(problem, Data_weighted, GPPUnCalibrated_weighted, GPPUnuncert_weighted)

## Reshaping of ET Data for Plots

ET_error_unweighted = results_overview[c(3,4),4]
ET_error_weighted = results_weighted_overview[,4]
ETCalibrated_weighted = vector(length = 4)
for(i in 1:2){
  ETCalibrated_weighted[(i-1)*2+1] = ET_error_unweighted[i]
  ETCalibrated_weighted[(i)*2] = ET_error_weighted[i]
}


ETUn_error_unweighted = results_overview[c(3,4),5]
ETUn_error_weighted = results_weighted_overview[,5]
ETUnCalibrated_weighted = vector(length = 4)
for(i in 1:2){
  ETUnCalibrated_weighted[(i-1)*2+1] = ETUn_error_unweighted[i]
  ETUnCalibrated_weighted[(i)*2] = ETUn_error_weighted[i]
}


ET_error_unweighted_uncert = results_overview_uncert[c(3,4),3]
ET_error_weighted_uncert = results_weighted_overview_uncert[,3]
ETCalibrated_weighted_uncert = vector(length =4)
for(i in 1:2){
  ETCalibrated_weighted_uncert[(i-1)*2+1] = ET_error_unweighted_uncert[i]
  ETCalibrated_weighted_uncert[(i)*2] = ET_error_weighted_uncert[i]
}


ETUn_error_unweighted_uncert = results_overview_uncert[c(3,4),4]
ETUn_error_weighted_uncert = results_weighted_overview_uncert[,4]
ETUnCalibrated_weighted_uncert = vector(length = 4)
for(i in 1:2){
  ETUnCalibrated_weighted_uncert[(i-1)*2+1] = ETUn_error_unweighted_uncert[i]
  ETUnCalibrated_weighted_uncert[(i)*2] = ETUn_error_weighted_uncert[i]
}

ETCalibrated_weighted_uncert = round(ETCalibrated_weighted_uncert,1)
ETUnCalibrated_weighted_uncert = round(ETUnCalibrated_weighted_uncert,1)



ETCalibrated_problem_weighted = data.frame(problem, Data_weighted, ETCalibrated_weighted,ETCalibrated_weighted_uncert)
ETUnCalibrated_problem_weighted = data.frame(problem, Data_weighted, ETUnCalibrated_weighted,ETUnCalibrated_weighted_uncert)

### Making the plots 

GPPplot_weighted<- ggplot(GPPCalibrated_problem_weighted, aes(fill = Data_weighted, y=GPPCalibrated_weighted, x=problem)) + 
  geom_bar(position="dodge", stat="identity", color = "black") + labs(y = element_blank()) +
  theme( axis.title.y=element_blank(), legend.position = "none", legend.title = element_blank()) + coord_flip(expand = c(0, 0))+
  scale_y_continuous(breaks=c(0.,0.2,0.4,0.6,0.8), limits = c(0., 0.82), 
                     labels = c(0,0.2,0.4,0.6,0.8)) +
  scale_fill_manual(values = rep(c("White", "Darkgrey"),3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(aes(label=GPPuncert_weighted), position=position_dodge(width=0.9), hjust= - .5)



GPPUnplot_weighted <- ggplot(GPPUnCalibrated_problem_weighted, aes(fill = Data_weighted, y=GPPUnCalibrated_weighted, x=problem)) + 
  geom_bar(position="dodge", stat="identity", color = "black")  + labs(y = element_blank()) +
  theme( legend.position = "none", axis.title.y=element_blank()) +  coord_flip(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(0.,0.2,0.4,0.6,0.8), limits = c(0., 0.82), 
                     labels = c(0,0.2,0.4,0.6,0.8)) +
  scale_fill_manual(values = rep(c("White", "Darkgrey"),3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(aes(label=GPPUnuncert_weighted), position=position_dodge(width=0.9), hjust= c(-.5,-.5,-.5,-.2))


ETplot_weighted <- ggplot(ETCalibrated_problem_weighted, aes(fill = Data_weighted, y=ETCalibrated_weighted, x=problem)) + 
  geom_bar(position="dodge", stat="identity", color = "black")  + labs(y = element_blank()) +
  theme(axis.text.y = element_blank(), legend.position = "none", axis.title.y=element_blank(), legend.title = element_blank()) +
  coord_flip(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(0.,0.05,0.10, 0.15,0.20), limits = c(-0., 0.21), 
                     labels =c(0.,0.05,0.10, 0.15,0.20)) +
  scale_fill_manual(values = rep(c("White", "Darkgrey"),2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text(aes(label=ETCalibrated_weighted_uncert), position=position_dodge(width=0.9), hjust= - .5)


ETUnplot_weighted <- ggplot(ETUnCalibrated_problem_weighted, aes(fill = Data_weighted, y=ETUnCalibrated_weighted, x=problem)) + 
  geom_bar(position="dodge", stat="identity", color = "black") + labs(y = element_blank()) +
  theme(axis.text.y = element_blank(), legend.position = "none",
        legend.direction = "vertical", axis.title.y=element_blank()) +  coord_flip(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(0.,0.05,0.10, 0.15,0.20), limits = c(-0., 0.21), 
                     labels =c(0.,0.05,0.10, 0.15,0.20)) +
  scale_fill_manual(values = rep(c("White", "Darkgrey"),2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text(aes(label=ETUnCalibrated_weighted_uncert), position=position_dodge(width=0.9), hjust= - .5)


Complete_overview_weighted_paper <- grid.arrange(arrangeGrob(arrangeGrob(GPPplot_weighted, top = "                                             GPP", left = "Calibration" ) ,
                                                       arrangeGrob(GPPUnplot_weighted, left = "Prediction" ), heights = c(0.5,0.5)),
                                           arrangeGrob(arrangeGrob(ETplot_weighted,  top = "ET"), arrangeGrob(ETUnplot_weighted)), ncol = 2
                                           , widths = c(1.3,1))



Completest_weighted_overview <-  plot_grid(arrangeGrob(parameters_plot_weighted_paper, left = ""),
                                              arrangeGrob(Complete_overview_weighted_paper),
                                              ncol = 1,
                                              align = "h", rel_heights = c(0.3,0.7), labels = c("a)", "b)"))

#### KOH solutions after calibration #####################


#### Correct model, balanced data, KOH 

unbal_no_bias_KOH <- post_predictiv_uncertainty_parallel(filename = "./Results/Model_validation_chain_less_parameter_add_noise.RData",
                                                         biastype = "none", KOH = T, parallel = T)


predictive_unbal_no_bias_KOH <- calculate_errors(unbal_no_bias_KOH)

uncert_unbal_no_bias_KOH <- sd_deviation_function(unbal_no_bias_KOH)


### Structural wrong model, balanced data, KOH 

unbal_func_bias_KOH <- post_predictiv_uncertainty_parallel(filename = "./Results/MCMC_func_biased_add_noise.RData"
                                                           , biastype = "func", KOH = T, parallel = T)

predictive_unbal_func_bias_KOH <- calculate_errors(unbal_func_bias_KOH)

uncert_unbal_func_bias_KOH <- sd_deviation_function(unbal_func_bias_KOH)


results_KOH_overview <- matrix(ncol = 4, nrow = 2)

results_KOH_overview_uncert <- matrix(ncol = 4, nrow = 2)

colnames(results_KOH_overview) <- c("Calibrated error GPP", "Uncalibrated error GPP", "Calibrated error ET", "Uncalibrated error ET")
rownames(results_KOH_overview) <- c("Correct Model", "Structural Bias")



results_KOH_overview[1,1] <- predictive_unbal_no_bias_KOH[1,1]
results_KOH_overview[1,2] <- predictive_unbal_no_bias_KOH[1,2]
results_KOH_overview[1,3] <- predictive_unbal_no_bias_KOH[2,1]
results_KOH_overview[1,4] <- predictive_unbal_no_bias_KOH[2,2]

results_KOH_overview_uncert[1,] = uncert_unbal_no_bias_KOH[1,]


results_KOH_overview[2,1] <- predictive_unbal_func_bias_KOH[1,1]
results_KOH_overview[2,2] <- predictive_unbal_func_bias_KOH[1,2]
results_KOH_overview[2,3] <- predictive_unbal_func_bias_KOH[2,1]
results_KOH_overview[2,4] <- predictive_unbal_func_bias_KOH[2,2]

results_KOH_overview_uncert[2,] = uncert_unbal_func_bias_KOH[1,] 


## Preparing GPP Data for Plot 

GPP_error_unKOH_uncert = results_overview_uncert[c(1,2),1]
GPP_error_KOH_uncert = results_KOH_overview_uncert[c(1,2),1]
GPPCalibrated_KOH_uncert = vector(length = 4)
for(i in 1:2){
  GPPCalibrated_KOH_uncert[(i-1)*2+1] = round(GPP_error_unKOH_uncert[i],1)
  GPPCalibrated_KOH_uncert[(i)*2] = round(GPP_error_KOH_uncert[i],1)
}

GPPUn_error_unKOH_uncert = results_overview_uncert[c(1,2),2]
GPPUn_error_KOH_uncert = results_KOH_overview_uncert[c(1,2),2]
GPPUnCalibrated_KOH_uncert = vector(length = 4)
for(i in 1:2){
  GPPUnCalibrated_KOH_uncert[(i-1)*2+1] = round(GPPUn_error_unKOH_uncert[i],1)
  GPPUnCalibrated_KOH_uncert[(i)*2] = round(GPPUn_error_KOH_uncert[i],1)
}


GPP_error_unKOH = results_overview[c(1,2),2]
GPP_error_KOH = results_KOH_overview[c(1,2),1]
GPPCalibrated_KOH = vector(length = 4)
for(i in 1:2){
  GPPCalibrated_KOH[(i-1)*2+1] = GPP_error_unKOH[i]
  GPPCalibrated_KOH[(i)*2] = GPP_error_KOH[i]
}


GPPUn_error_unKOH = results_overview[c(1,2),3]
GPPUn_error_KOH = results_KOH_overview[c(1,2),2]
GPPUnCalibrated_KOH = vector(length = 4)
for(i in 1:2){
  GPPUnCalibrated_KOH[(i-1)*2+1] = GPPUn_error_unKOH[i]
  GPPUnCalibrated_KOH[(i)*2] = GPPUn_error_KOH[i]
}

Data_KOH = rep(c("Uncorrected Error","Corrected Error"),2)

GPPCalibrated_problem_KOH = data.frame(problem, Data_KOH, GPPCalibrated_KOH,GPPCalibrated_KOH_uncert)

GPPUnCalibrated_problem_KOH = data.frame(problem, Data_KOH, GPPUnCalibrated_KOH,GPPUnCalibrated_KOH_uncert)

### Preparing ET data for plots 

ET_error_unKOH = results_overview[c(1,2),4]
ET_error_KOH = results_KOH_overview[c(1,2),3]
ETCalibrated_KOH = vector(length = 4)
for(i in 1:2){
  ETCalibrated_KOH[(i-1)*2+1] = ET_error_unKOH[i]
  ETCalibrated_KOH[(i)*2] = ET_error_KOH[i]
}


ETUn_error_unKOH = results_overview[c(1,2),5]
ETUn_error_KOH = results_KOH_overview[c(1,2),4]
ETUnCalibrated_KOH = vector(length = 4)
for(i in 1:2){
  ETUnCalibrated_KOH[(i-1)*2+1] = ETUn_error_unKOH[i]
  ETUnCalibrated_KOH[(i)*2] = ETUn_error_KOH[i]
}

ET_error_unKOH_uncert = results_overview_uncert[c(1,2),3]
ET_error_KOH_uncert = results_KOH_overview_uncert[c(1,2),3]
ETCalibrated_KOH_uncert = vector(length = 4)
for(i in 1:2){
  ETCalibrated_KOH_uncert[(i-1)*2+1] = round(ET_error_unKOH_uncert[i],1)
  ETCalibrated_KOH_uncert[(i)*2] = round(ET_error_KOH_uncert[i],1)
}




ETUn_error_unKOH_uncert = results_overview_uncert[c(1,2),4]
ETUn_error_KOH_uncert = results_KOH_overview_uncert[c(1,2),4]
ETUnCalibrated_KOH_uncert = vector(length = 4)
for(i in 1:2){
  ETUnCalibrated_KOH_uncert[(i-1)*2+1] = round(ETUn_error_unKOH_uncert[i],1)
  ETUnCalibrated_KOH_uncert[(i)*2] = round(ETUn_error_KOH_uncert[i],1)
}


ETCalibrated_problem_KOH = data.frame(problem, Data_KOH, ETCalibrated_KOH, ETCalibrated_KOH_uncert)

ETUnCalibrated_problem_KOH = data.frame(problem, Data_KOH, ETUnCalibrated_KOH, ETUnCalibrated_KOH_uncert)



## Making the GPP plots

GPPplot_KOH <- ggplot(GPPCalibrated_problem_KOH,aes(fill = Data_KOH, y=GPPCalibrated_KOH, x=problem)) + 
  geom_bar(position="dodge", stat="identity", color = "black") + labs(y = element_blank()) +
  theme( axis.title.y=element_blank(), legend.position = "none", legend.title = element_blank()) +
  coord_flip(expand = c(0, 0)) + 
  scale_y_continuous(breaks=c(0.,0.2,0.4,0.6,0.8), limits = c(0., 0.82), 
                     labels = c(0,0.2,0.4,0.6,0.8)) +
  scale_fill_manual(values = rep(c( "LightGrey", "Black"),2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text(aes(label=GPPCalibrated_KOH_uncert), position=position_dodge(width=0.9), hjust= - .5)


GPPUnplot_KOH <- ggplot(GPPUnCalibrated_problem_KOH, aes(fill = Data_KOH, y=GPPUnCalibrated_KOH, x=problem)) +
  geom_bar(position="dodge", stat="identity", color = "black")  + labs(y = element_blank()) +
  theme( legend.position = "none", axis.title.y=element_blank()) +  coord_flip(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(0.,0.2,0.4,0.6,0.8), limits = c(0., 0.82), 
                     labels = c(0,0.2,0.4,0.6,0.8)) +
  scale_fill_manual(values = rep(c( "LightGrey", "Black"),2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text(aes(label=GPPUnCalibrated_KOH_uncert), position=position_dodge(width=0.9), hjust= - .5)


## Making the ET plots

ETplot_KOH <- ggplot(ETCalibrated_problem_KOH, aes(fill = Data_KOH, y=ETCalibrated_KOH, x=problem)) + 
  geom_bar(position="dodge", stat="identity", color = "black")  + labs(y = element_blank()) +
  theme(axis.text.y = element_blank(), axis.title.y=element_blank(), legend.position = c(0.8,0.7),
        legend.title = element_blank()) +  coord_flip(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(0.,0.05,0.10, 0.15,0.20), limits = c(-0., 0.21), 
                     labels =c(0.,0.05,0.10, 0.15,0.20)) +
  scale_fill_manual(values = rep(c( "LightGrey","Black"),2),guide=guide_legend(reverse=T)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(aes(label=ETCalibrated_KOH_uncert), position=position_dodge(width=0.9), hjust= - .5)



ETUnplot_KOH <- ggplot(ETUnCalibrated_problem_KOH, aes(fill = Data_KOH, y=ETUnCalibrated_KOH, x=problem)) + 
  geom_bar(position="dodge", stat="identity", color = "black") + labs(y = element_blank()) +
  theme(axis.text.y = element_blank(), legend.position = "none",
        legend.direction = "vertical", axis.title.y=element_blank()) +  coord_flip(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(0.,0.05,0.10, 0.15,0.20), limits = c(-0., 0.21), 
                     labels =c(0.,0.05,0.10, 0.15,0.20)) +
  scale_fill_manual(values = rep(c( "LightGrey","Black"),2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(aes(label=ETUnCalibrated_KOH_uncert), position=position_dodge(width=0.9), hjust= - .5)



Complete_overview_KOH_paper <- grid.arrange(arrangeGrob(arrangeGrob(GPPplot_KOH, top = "                                             GPP", left = "Calibration" ) ,
                                                        arrangeGrob(GPPUnplot_KOH, left = "Prediction" ), heights = c(0.5,0.5)),
                                            arrangeGrob(arrangeGrob(ETplot_KOH,top = "ET"), arrangeGrob(ETUnplot_KOH) ), ncol = 2, widths = c(1.3,1.))





#### KOH solutions while calib ####




### Correct Model, balanced data, KOH while calibration

unbal_no_bias_KOH_while <- post_predictiv_uncertainty_parallel(filename = "./Results/MCMC_KOH_add_noise.RData",
                                                          biastype = "none", KOH = T, parallel = T)


predictive_unbal_no_bias_KOH_while <- calculate_errors(unbal_no_bias_KOH_while)

uncert_unbal_no_bias_KOH_while <- sd_deviation_function(unbal_no_bias_KOH_while)

parameter_unbal_no_bias_KOH_while <- calculate_parameter_error("./Results/MCMC_KOH_add_noise.RData")



### Structural model error, balanced data, KOH while calibration

unbal_func_bias_KOH_while <- post_predictiv_uncertainty_parallel(filename = "./Results/MCMC_func_KOH_add_noise.RData", 
                                                            biastype = "func", KOH = T, parallel = T)

predictive_unbal_func_bias_KOH_while <- calculate_errors(unbal_func_bias_KOH_while)

uncert_unbal_func_bias_KOH_while <- sd_deviation_function(unbal_func_bias_KOH_while)

parameter_unbal_func_bias_KOH_while <- calculate_parameter_error("./Results/MCMC_func_KOH_add_noise.RData")



### Reanranging Data for Plots

results_KOH_while_overview <- matrix(ncol = 5, nrow = 2)

results_KOH_while_overview_uncert <- matrix(ncol = 4, nrow = 2)

colnames(results_KOH_while_overview) <- c("Parameter error", "Calibrated error GPP", "Uncalibrated error GPP", "Calibrated error ET", "Uncalibrated error ET")
rownames(results_KOH_while_overview) <- c("Correct Model", "Structural Bias")


results_KOH_while_overview[1,1] <- parameter_unbal_no_bias_KOH_while
results_KOH_while_overview[1,2] <- predictive_unbal_no_bias_KOH_while[1,1]
results_KOH_while_overview[1,3] <- predictive_unbal_no_bias_KOH_while[1,2]
results_KOH_while_overview[1,4] <- predictive_unbal_no_bias_KOH_while[2,1]
results_KOH_while_overview[1,5] <- predictive_unbal_no_bias_KOH_while[2,2]

results_KOH_while_overview_uncert[1,] = uncert_unbal_no_bias_KOH_while[1,]

results_KOH_while_overview[2,1] <- parameter_unbal_func_bias_KOH_while
results_KOH_while_overview[2,2] <- predictive_unbal_func_bias_KOH_while[1,1]
results_KOH_while_overview[2,3] <- predictive_unbal_func_bias_KOH_while[1,2]
results_KOH_while_overview[2,4] <- predictive_unbal_func_bias_KOH_while[2,1]
results_KOH_while_overview[2,5] <- predictive_unbal_func_bias_KOH_while[2,2]

results_KOH_while_overview_uncert[2,] = uncert_unbal_func_bias_KOH_while[1,] 


### Parameter plot 

error_uncorrected = results_overview[c(1,2),1]
error_KOH_while = results_KOH_while_overview[c(1,2),1]
problem =  c(rep("Correct Model",2), rep("Systematic Model Error",2))
Data_KOH = rep(c("Uncorrected Error","Corrected Error"),2)



ParameterError_KOH_while = vector(length = 4)
for(i in 1:2){
  ParameterError_KOH_while[(i-1)*2+1] = error_uncorrected[i]
  ParameterError_KOH_while[(i)*2] = error_KOH_while[i]
}


data_problem_KOH_while =  data.frame(problem, Data_KOH, ParameterError_KOH_while)

parameters_plot_KOH_while_paper <- ggplot() + 
  geom_bar(data = data_problem_KOH_while, aes(fill = Data_KOH, y=ParameterError_KOH_while*100, x=problem)
           , position="dodge", stat="identity",color = "black") +
  labs(y = element_blank()) +
  theme(text = element_text(), legend.text=element_text()) +
  theme(axis.text.x = element_text(angle = 90),  legend.position = c(0.9,0.7), axis.title.y = element_blank(), line = element_blank(),
        legend.direction = "vertical", legend.title = element_blank()) + coord_flip(expand = c(0, 0)) + theme(aspect.ratio = 0.1) +
  theme(axis.text.x = element_text(angle = 0)) + ylim(c(0, max(ParameterError_KOH_while*100)+2)) +
  scale_fill_manual(values = rep(c("Lightgrey","Black"),2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


## Preparing GPP data  

GPP_error_unKOH_while_uncert = results_overview_uncert[c(1,2),1]
GPP_error_KOH_while_uncert = results_KOH_while_overview_uncert[c(1,2),1]
GPPCalibrated_KOH_while_uncert = vector(length = 4)
for(i in 1:2){
  GPPCalibrated_KOH_while_uncert[(i-1)*2+1] = round(GPP_error_unKOH_while_uncert[i],1)
  GPPCalibrated_KOH_while_uncert[(i)*2] = round(GPP_error_KOH_while_uncert[i],1)
}

GPPUn_error_unKOH_while_uncert = results_overview_uncert[c(1,2),2]
GPPUn_error_KOH_while_uncert = results_KOH_while_overview_uncert[c(1,2),2]
GPPUnCalibrated_KOH_while_uncert = vector(length = 4)
for(i in 1:2){
  GPPUnCalibrated_KOH_while_uncert[(i-1)*2+1] = round(GPPUn_error_unKOH_while_uncert[i],1)
  GPPUnCalibrated_KOH_while_uncert[(i)*2] = round(GPPUn_error_KOH_while_uncert[i],1)
}


GPP_error_unKOH_while = results_overview[c(1,2),2]
GPP_error_KOH_while = results_KOH_while_overview[c(1,2),2]
GPPCalibrated_KOH_while = vector(length = 4)
for(i in 1:2){
  GPPCalibrated_KOH_while[(i-1)*2+1] = GPP_error_unKOH_while[i]
  GPPCalibrated_KOH_while[(i)*2] = GPP_error_KOH_while[i]
}


GPPUn_error_unKOH_while = results_overview[c(1,2),3]
GPPUn_error_KOH_while = results_KOH_while_overview[c(1,2),3]
GPPUnCalibrated_KOH_while = vector(length = 4)
for(i in 1:2){
  GPPUnCalibrated_KOH_while[(i-1)*2+1] = GPPUn_error_unKOH_while[i]
  GPPUnCalibrated_KOH_while[(i)*2] = GPPUn_error_KOH_while[i]
}


GPPCalibrated_problem_KOH_while  = data.frame(problem, Data_KOH, GPPCalibrated_KOH_while ,GPPCalibrated_KOH_while_uncert)
GPPUnCalibrated_problem_KOH_while = data.frame(problem, Data_KOH, GPPUnCalibrated_KOH_while,GPPUnCalibrated_KOH_while_uncert)


### Preparing ET data for Plots 


ET_error_unKOH_while = results_overview[c(1,2),4]
ET_error_KOH_while = results_KOH_while_overview[c(1,2),4]
ETCalibrated_KOH_while = vector(length = 4)
for(i in 1:2){
  ETCalibrated_KOH_while[(i-1)*2+1] = ET_error_unKOH_while[i]
  ETCalibrated_KOH_while[(i)*2] = ET_error_KOH_while[i]
}


ETUn_error_unKOH_while = results_overview[c(1,2),5]
ETUn_error_KOH_while = results_KOH_while_overview[c(1,2),5]
ETUnCalibrated_KOH_while = vector(length = 4)
for(i in 1:2){
  ETUnCalibrated_KOH_while[(i-1)*2+1] = ETUn_error_unKOH_while[i]
  ETUnCalibrated_KOH_while[(i)*2] = ETUn_error_KOH_while[i]
}

ET_error_unKOH_while_uncert = results_overview_uncert[c(1,2),3]
ET_error_KOH_while_uncert = results_KOH_while_overview_uncert[c(1,2),3]
ETCalibrated_KOH_while_uncert = vector(length = 4)
for(i in 1:2){
  ETCalibrated_KOH_while_uncert[(i-1)*2+1] = round(ET_error_unKOH_while_uncert[i],1)
  ETCalibrated_KOH_while_uncert[(i)*2] = round(ET_error_KOH_while_uncert[i],1)
}


ETUn_error_unKOH_while_uncert = results_overview_uncert[c(1,2),4]
ETUn_error_KOH_while_uncert = results_KOH_while_overview_uncert[c(1,2),4]
ETUnCalibrated_KOH_while_uncert = vector(length = 4)
for(i in 1:2){
  ETUnCalibrated_KOH_while_uncert[(i-1)*2+1] = round(ETUn_error_unKOH_while_uncert[i],1)
  ETUnCalibrated_KOH_while_uncert[(i)*2] = round(ETUn_error_KOH_while_uncert[i],1)
}



ETCalibrated_problem_KOH_while = data.frame(problem, Data_KOH, ETCalibrated_KOH_while, ETCalibrated_KOH_while_uncert)

ETUnCalibrated_problem_KOH_while = data.frame(problem, Data_KOH, ETUnCalibrated_KOH_while, ETUnCalibrated_KOH_while_uncert)


GPPplot_KOH <- ggplot(GPPCalibrated_problem_KOH_while, aes(fill = Data_KOH, y=GPPCalibrated_KOH_while, x=problem)) + 
  geom_bar(position="dodge", stat="identity", color = "black") + labs(y = element_blank()) +
  theme( axis.title.y=element_blank(), legend.position = "none", legend.title = element_blank()) +
  coord_flip(expand = c(0, 0)) + 
  scale_y_continuous(breaks=c(0.,0.2,0.4,0.6,0.8), limits = c(0., 0.82), 
                     labels = c(0,0.2,0.4,0.6,0.8)) +
  scale_fill_manual(values = rep(c( "LightGrey", "Black"),2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(aes(label=GPPCalibrated_KOH_while_uncert), position=position_dodge(width=0.9), hjust= - .5)



GPPUnplot_KOH <- ggplot(GPPUnCalibrated_problem_KOH_while, aes(fill = Data_KOH, y=GPPUnCalibrated_KOH_while, x=problem)) +
  geom_bar(position="dodge", stat="identity", color = "black")  + labs(y = element_blank()) +
  theme( legend.position = "none", axis.title.y=element_blank()) +  coord_flip(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(0.,0.2,0.4,0.6,0.8), limits = c(0., 0.82), 
                     labels = c(0,0.2,0.4,0.6,0.8)) +
  scale_fill_manual(values = rep(c( "LightGrey", "Black"),2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text(aes(label=GPPUnCalibrated_KOH_while_uncert), position=position_dodge(width=0.9), hjust= - .5)


ETplot_KOH <- ggplot(ETCalibrated_problem_KOH_while, aes(fill = Data_KOH, y=ETCalibrated_KOH_while, x=problem)) + 
  geom_bar(position="dodge", stat="identity", color = "black")  + labs(y = element_blank()) +
  theme(axis.text.y = element_blank(), axis.title.y=element_blank(), legend.position = "none",
        legend.title = element_blank()) +  coord_flip(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(0.,0.05,0.10, 0.15,0.20), limits = c(-0., 0.21), 
                     labels =c(0.,0.05,0.10, 0.15,0.20)) +
  scale_fill_manual(values = rep(c( "LightGrey","Black"),2),guide=guide_legend(reverse=T)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(aes(label=ETCalibrated_KOH_while_uncert), position=position_dodge(width=0.9), hjust= - .5)


ETUnplot_KOH <- ggplot(ETUnCalibrated_problem_KOH_while, aes(fill = Data_KOH, y=ETUnCalibrated_KOH_while, x=problem)) + 
  geom_bar(position="dodge", stat="identity", color = "black") + labs(y = element_blank()) +
  theme(axis.text.y = element_blank(), legend.position = "none",
        legend.direction = "vertical", axis.title.y=element_blank()) +  coord_flip(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(0.,0.05,0.10, 0.15,0.20), limits = c(-0., 0.21), 
                     labels =c(0.,0.05,0.10, 0.15,0.20)) +
  scale_fill_manual(values = rep(c( "LightGrey","Black"),2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(aes(label=ETUnCalibrated_KOH_while_uncert), position=position_dodge(width=0.9), hjust= - .5)



Complete_overview_paper <- grid.arrange(arrangeGrob(arrangeGrob(GPPplot, top = "                                             GPP", left = "Calibration" ) ,
                                                    arrangeGrob(GPPUnplot, left = "Prediction" ), heights = c(0.5,0.5)),
                                        arrangeGrob(arrangeGrob(ETplot, top = "ET"), arrangeGrob(ETUnplot)), ncol = 2,
                                        widths = c(1.3,1))


Complete_overview_KOH_while_paper <- grid.arrange(arrangeGrob(arrangeGrob(GPPplot_KOH, top = "                                             GPP", left = "Calibration" ) ,
                                                        arrangeGrob(GPPUnplot_KOH, left = "Prediction" ), heights = c(0.5,0.5)),
                                            arrangeGrob(arrangeGrob(ETplot_KOH,top = "ET"), arrangeGrob(ETUnplot_KOH) ), ncol = 2,
                                            widths = c(1.3,1.))

Completest_overview_paper_KOH_while <- plot_grid(arrangeGrob(parameters_plot_KOH_while_paper, left = ""),
                                        arrangeGrob(Complete_overview_KOH_while_paper),
                                        rel_heights = c(0.3,0.7), labels = c("a)","b)"), align = "h", ncol = 1)






