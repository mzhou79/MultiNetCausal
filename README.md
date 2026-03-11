# Causal-effect-estimation-with-multiple-network-induced-dependencies
## Overview
Evaluating causal effects with network data is often challenged by treatment interference, where the treatment assigned to one unit has a causal effect on others’ outcomes, and within-variable dependence or network dependence, where each variable is autocorrelated through network ties. Often, statistical associations between one’s treatment and others’ outcomes can be due to treatment interference, network dependence, or other shared variable factors. Much of the recent literature on causal inference with network data assumes that the network structures inducing these different types of dependence completely overlap. In this work, we propose a chain graph representation that reflects multiple networks capturing distinct types of causal and statistical dependence and apply auto-g-computation methods to estimate direct and indirect effects. 

This repository provides the code that can reproduce the simulation results in the manuscript for both settings, and the code for data application.

## Code

- `Setting1.R`: includes function and data generating process to reproduce setting 1.

- `Setting2.R`: includes function and data generating process to reproduce setting 2.

- `uconnect.R`: includes the code for data cleaning, network visulization, and the application to uConnect Study.
