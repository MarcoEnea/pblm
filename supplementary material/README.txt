*********************************************************************************
This file describes the supplementary material accompanying the code of the paper:

Enea, M., Lovison, G. 2018. A penalized approach for the bivariate ordered logistic model with applications to social and medical data. Statistical Modelling, DOI: 10.1177/1471082X18782063.


**********************************************************************************


1. Overview

The folder "supplementary material" contains:

   -  the compressed file "pblm_0.1-8.tar.gz"

   -  the folder "LR_P simulations"

   -  the folder "performance simulations"

   -  the folder "british males dataset analysis", 

   -  the folder "liver simulated data analysis"



2. The "pblm" package 

The file "pblm_0.1-8.tar.gz" contains the "pblm" R package to install.
This is an R package, written by Marco Enea, to fit bivariate additive categorical 
regression for moderate-to-small size datasets. The two responses can be  
nominal, ordinal or mixed nominal/ordinal. Partial proportional odds models with 
(non-)uniform association structure can be fitted, with the possibility to specify 
several logit types and parametrizations for the marginals and the association, 
including the Dale model. The marginal parameters and the association structure 
can also be smoothed, and/or the parameter space regularized, by using the penalty 
terms discussed ("ARC1" and "ARC2") or outlined ("ridge", "lasso" and "lassoV") in 
the paper, even though the implemetation of the latter two penalty terms needs to be 
better checked out and possibly subjected to enhancements. The penalty for mimicking 
inequality constraints is specified by default, when ordered responses are involved.  
The package also implements P-splines. Mikis Stasinopoulos and Robert Rigby 
contributed to that as authors of the gamlss functions pb() and ps() which the 
corresponding pblm functions pb() and pbs() are based on. They also helped the first 
author with their useful suggestions on the implementation of the backfitting 
algorithm within the pblm package. Finally, common "methods" such as summary, 
residuals and predict are available. At time of writing, the package is not on CRAN.



3.  LR_P simulations

This folder contains seven sub-folders with all simulations described in Section 4.
Six of them are related to the simulation schemes (sim1 or sim2) and vary with the 
number of response levels. The folder "sim_all" contains the code to build Figures 
1-4. However, this latter code cannot be run directly, as it requires the loading of 
the workspaces of all simulations. Accordingly, in order to reproduce all figures, 
each single simulation scheme must be previously run. The sub-folder sim1_ncat3 
contains three R codes, according to the sample size (200, 500, and 1000), for the 
simulationscheme 1 when D_1=D_2=3. The execution of file "sim_n200_ncat3_ncov2dic.R" 
will produce the file "sim_n200_ncat3_ncov2dic.pdf" with the histograms of the 
simulation varying with the penalty parameter lambda. Once the three codes have been 
run and their workspaces saved, the code "Figure_sim1.R" can be executed to produce a 
figure (sim1_ncat3_ncov2dic_all.pdf) which compares the three scenarios for some lamba 
values.

  


4. performance simulations

This contains three R files to run all simulations described in Section 5, and an 
additional R file "add_fun.R" containing additional R function not included into 
the pblm package.


5. british males dataset analysis

This folder contains all analyses and codes to run described in Section 6.1.


6. liver (simulated) data analysis

This folder contains data simulated from the original liver disease dataset described
in Section 6.2, as well as the code to perform all analyses and simulations. 
      


 
