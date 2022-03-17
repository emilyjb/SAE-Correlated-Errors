# SAE-Correlated-Errors
This repository has code to implement the simulations in "Small area estimation for the Fay-Herriot model in which the measurement error is correlated with the sampling error in the response."  The code in this repository implements the simulations that compare alternative predictors for a small area model in which the sampling error in the response is correlated with the measurement error in the covariate. We do not use a fixed seed, so the results will be similar but not identical to those reported. 

* The file "equalpsiMLflipped.R" runs the simulation with equal psi. 
* The file "unequalpsiMLflipped.R" runs the simulations with unequal psi. 
* The file "unequalpsiMLflippedTMod.R" runs the simulations t or chi-square distributions.
* The file "unequalpsiMLflippedMultiCov.R" runs the simulations with unequal psi. 
* The file "getoutputloopabc.R" gathers the output. 
* The remaining files are sourced in the above programs. 
