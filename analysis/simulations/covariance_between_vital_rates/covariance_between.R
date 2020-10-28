suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ipmr))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(faux))


start <- Sys.time()
start
#  ----------------------------------------------------------------------------------------------------------------------------
# parsing arguments
#  ----------------------------------------------------------------------------------------------------------------------------

Parsoptions <- list (
  
  make_option(
    opt_str = c("-f", "--function"),
    dest    = "foption",
    help    = "Specify the ipmr function to use",
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

if(foption == "HEQU") {
  source("/home/evers/lagged_buffering/analysis/simulations/covariance_between_vital_rates/ipmr_functions_HEQU.R")
}
if(foption == "OPIM")
  source("/home/evers/lagged_buffering/analysis/simulations/covariance_between_vital_rates/ipmr_functions_OPIM.R")



params_list <- read.csv(parameters) 
params_list <- as.list(setNames(as.numeric(params_list$Value), params_list$Parameter))

str(params_list)

clim_list <- read.csv(clim_par) 
clim_list <- as.list(setNames(as.numeric(clim_list$Value), clim_list$Parameter))

str(clim_list)

clim_sd <- rep(seq(from = 0, to = 2, length.out = 10), 90)
clim_corr <- rep(rep(c(-0.9,0,0.9), each = 10), 30)

clim_sd[taskID]
clim_corr[taskID]

 
lambda <- P_lambdas(n_it = n_it, 
                    clim_sd = clim_sd[taskID], 
                    clim_corr = clim_corr[taskID], 
                    params_list = params_list, 
                    clim_params = clim_list,
                    n_mesh = n_mesh,
                    save_K = T)

print("ipm done")  
# b <- lapply(lambda$M_no_cov_ipm[[1]], function(x) as.vector(x)) %>% bind_rows %>% t
# print("start corr1")
# c <- corrr::correlate(b)
# print("correlation1 done")
# d <- as.matrix(c[,-1])
# n_corr_hist <- d
# n_corr_sum <- mean(d, na.rm = T)
# n_corr_sd <- sd(d, na.rm = T)
n_lambda <- lambda$no_cov_lambda

rm(b,c,d)

print("start corr2")
# e <- lapply(lambda$M_ipm_cov[[1]], function(x) as.vector(x)) %>% bind_rows %>% t
# print("correlation2 done")
# f <- corrr::correlate(e)
# g <- as.matrix(f[,-1])
# c_corr_hist <- g
# c_corr_sum <- mean(g, na.rm = T)
# c_corr_sd <- sd(g, na.rm = T)
c_lambda <- lambda$cov_lambda


df <- tibble(clim_corr[taskID], clim_sd[taskID],
             n_lambda, c_lambda) #,
             # n_corr_sum, c_corr_sum,
             # n_corr_sd, c_corr_sd,
             # list(n_corr_hist), list(c_corr_hist))
str(df)
output

saveRDS(df, file = output)

Sys.time() - start

