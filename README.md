# SAE-Correlated-Errors

**Paper title:** "An application of a small area procedure with correlation between measurement error and sampling error to the Conservation Effects Assessment Project"
**Authors:** Emily Berg and Sepideh Mosaferi

This repository has code to implement the simulations for the above paper. The paper develops a small area estimation procedure for a model in which the measurement error in the covariate is correlated with the sampling error in the response. The code in this repository compares the method proposed in the paper to several alternative predictors. We do not use a fixed seed, so the results will be similar but not identical to those reported. 

* The file "equalpsiMLflipped.R" runs the simulation with equal psi. 
* The file "unequalpsiMLflipped.R" runs the simulations with unequal psi. 
* The file "unequalpsiMLflippedTMod.R" runs the simulations t or chi-square distributions.
* The file "unequalpsiMLflippedMultiCov.R" runs the simulations with unequal psi. 
* The file "getoutputloopabc.R" gathers the output. 
* The remaining files are sourced in the above programs. 
