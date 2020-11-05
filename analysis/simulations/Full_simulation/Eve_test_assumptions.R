suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ipmr))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))



start <- Sys.time()
start
#  ----------------------------------------------------------------------------------------------------------------------------
# parsing arguments
#  ----------------------------------------------------------------------------------------------------------------------------

Parsoptions <- list (
  
  make_option(
    opt_str = c("-f", "--function"),
    dest    = "foption",
    help    = "Specify the ipmr function (based on HEQU or OPIM) to use",
    metavar = "HEQU|OPIM"),
  
  make_option(
    opt_str = c("-i", "--iterations"),
    dest    = "n_it",
    help    = "number of iterations to run the IPM",
    default = 10000),
  
  make_option(
    opt_str = c("-m", "--meshpoints"),
    dest    = "n_mesh",
    help    = "number of meshpoints when integrating the IPM",
    default = 100)
  
)

parser <- OptionParser(
  usage       = "Rscript %prog [options] output parameters climate_par",
  option_list = Parsoptions,
  description = "\nan Run lagged ipmr functions",
  epilogue    = ""
)

cli <- parse_args(parser, positional_arguments = 3)

foption <- cli$options$foption
n_mesh <- as.integer(cli$options$n_mesh)
n_it <- as.integer(cli$options$n_it)

taskID <- as.integer(Sys.getenv("SGE_TASK_ID"))

output <- cli$args[1]
parameters <- cli$args[2]
clim_par <- cli$args[3]


foption
taskID
n_mesh
n_it

if(foption == "HEQU"){
  source("/home/evers/lagged_buffering/analysis/simulations/Full_simulation/ipmr_functions_HEQU.R")
}
if(foption == "OPIM") {
  source("/home/evers/lagged_buffering/analysis/simulations/Full_simulation/ipmr_functions_OPIM.R")
}

params_list <- read.csv(parameters) 
params_list <- as.list(setNames(as.numeric(params_list$Value), params_list$Parameter))

str(params_list)

clim_list <- read.csv(clim_par) 
clim_list <- as.list(setNames(as.numeric(clim_list$Value), clim_list$Parameter))

str(clim_list)

## create parameter list -------------------------------------------------------------------------


clim_sd <- rep(seq(from = 0, to = 2, length.out = 10), 90)
clim_corr <- rep(rep(c(-0.9,0,0.9), each = 10), 30)

clim_sd[taskID]
clim_corr[taskID]

message("Within P kernel 1 year lagged")
lambda <- P_lambdas(n_it = n_it, 
                    clim_sd = clim_sd[taskID], 
                    clim_corr = clim_corr[taskID], 
                    params_list = params_list, 
                    clim_params = clim_list,
                    n_mesh = n_mesh,
                    save_K = F)






print("ipm done")  

lambda

n_lambda <- lambda$non_lagged
s_lambda <- lambda$lagged_s
g_lambda <- lambda$lagged_g

df <- tibble(clim_corr = clim_corr[taskID], 
             clim_sd = clim_sd[taskID],
             n_lambda, s_lambda, g_lambda)

str(df)
output

saveRDS(df, file = output)

Sys.time() - start
