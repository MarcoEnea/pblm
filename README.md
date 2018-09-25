# pblm
An R package to fit bivariate additive categorical regression models

This is an R package, written by Marco Enea, to fit bivariate additive categorical 
regression for moderate-to-small size datasets. The two responses can be  
nominal, ordinal or mixed nominal/ordinal. Partial proportional odds models with 
(non-)uniform association structure can be fitted, with the possibility to specify 
several logit types and parametrizations for the marginals and the association, 
including the Dale model. The marginal parameters and the association structure 
can also be smoothed, and/or the parameter space regularized, by using the penalty 
terms. The package also implements P-splines. 
