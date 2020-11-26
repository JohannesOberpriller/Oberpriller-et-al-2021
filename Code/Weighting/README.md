# Weighting of data streams 

This folder contains the code for weighting of data streams. 
+ Weighted_unbal_biased_add_noise.R: The "true" model with weighted likelihoods following the procedure described in the manuscript and below. 
+ Weighted_func_unbal_biased_add_noise.R: The model with structural errors with weighted likelihoods following the procedure described in the manuscript and below. 

In general, when we have **N** data streams with *n<sub>1</sub>* ... *n<sub>N</sub>* data points the model with structural error is going to overfit
to the data stream with the most data points. 

To provide the model to overfit we test the stratgey of weighting the data inversly to their number of data points. 
In detail, let *n<sub>i</sub> be the lowest number of datapoints of the different data streams: s.t. n = min{*n<sub>1</sub>* ... *n<sub>N</sub>*}, 
where inf is the minimum. 
Then for the individual likelihood contributions of the data streams *L<sub>1</sub>* ... *L<sub>N</sub>*, s.t. L<sub>complete</sub> = *L<sub>1</sub>* + ... + *L<sub>N</sub>*, 
we change the likelihood to *L<sub>j</sub>* = n<sub>i</sub>/n<sub>j</sub>, s.t. the complete likelihood is L<sub>complete_new</sub> = n<sub>i</sub>/n<sub>1</sub>*L<sub>1</sub>* + ... + n<sub>i/</sub>n<sub>N</sub>*L<sub>N</sub>*.

**More details are in the individual scripts.** 
