# PerfEst
Toolbox for simulating tissue blood circulation, and estimating the perfusion. 

The main script is runTissueEstimation.m. To generate a prior mean, the script getPrior.m can be used. An ensemble of model realizations can be generated using generateInitialEnsemble.m. In order to run the main script, an initialization script must be present. See the Frog case for an example. 

Currently one example is given in the Frog-folder. Here blood flow through a frog tongue is simulated, and the perfusion is estimated. The folder contains 8 cases, each with different resolution for the simulation grid. Each case is run by first executing setup.m, and then executing runTissueEstimation.m 
