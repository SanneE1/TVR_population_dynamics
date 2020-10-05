library(ipmr)
library(dplyr)
library(parallel)
library(pbapply)

source("analysis/simulations/ipmr_functions.R")

params_list <- read.csv("data/simulation/perennial1_HEQU_parameters.csv") 
params_list <- as.list(setNames(as.numeric(params_list$Value), params_list$Parameter))

clim_list <- read.csv("data/simulation/perennial1_HEQU_climate.csv") 
clim_list <- as.list(setNames(as.numeric(clim_list$Value), clim_list$Parameter))

clim_sd <- rep(seq(from = 1, to = 2, length.out = 2), 90)
clim_corr <- rep(rep(c(-0.9,0,0.9), each = 2), 30)


## make clusters
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

## export libraries to workers
clusterEvalQ(cl, c(library(ipmr), library(dplyr)))
## provide seed number to workers
clusterSetRNGStream(cl, iseed = 04103)
## export objects to workers
clusterExport(cl, c("params_list", "clim_list", "clim_sd", "clim_corr", "P_lambdas", "create_seq"))

lambdas <- pblapply(cl = cl,
          as.list(c(1:250)), 
          function(x) P_lambdas(n_it = 10000, 
                                clim_sd = clim_sd[x], 
                                clim_corr = clim_corr[x], 
                                params_list = params_list, 
                                clim_params = clim_list,
                                n_mesh = 100,
                                save_K = F))

stopCluster(cl)
