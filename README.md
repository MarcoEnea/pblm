# pblm
An R package to fit bivariate additive categorical regression models

This is an R package, by Marco Enea, to fit bivariate additive categorical regression for moderate size datasets. The two responses can be nominal, ordinal or mixed nominal/ordinal. Partial proportional odds models with (non-)uniform association structure can be fitted, with the possibility to specify several logit types and parametrizations for the marginals and the association, including the Dale model. The marginal parameters and the association structure can also be smoothed, and/or the parameter space regularized, by using penalty terms. The package also implements P-splines. 

The Supplementary Material folder contains R code for examples and simulations reported in the pubblication: 

Enea, M., Lovison, G. 2018. A penalized approach for the bivariate ordered logistic model with applications to social and medical data. Statistical Modelling, DOI: 10.1177/1471082X18782063.  
