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
    help    = "Specify the ipmr function to use",
    metavar = "P_1yr|P_neg_1yr|P_2yr|PF_1yr|PF_neg_1yr"),
  
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

source("/home/evers/lagged_buffering/analysis/simulations/autocorrelation_single_vital_rate/ipmr_functions.R")
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

if(foption == "P_1yr") {
  message("Within P kernel 1 year lagged")
  lambda <- P_lambdas(n_it = n_it, 
                      clim_sd = clim_sd[taskID], 
                      clim_corr = clim_corr[taskID], 
                      params_list = params_list, 
                      clim_params = clim_list,
                      n_mesh = n_mesh,
                      save_K = T)
  
} 


print("ipm done")  
b <- lapply(lambda$M_no_auto_ipm[[1]], function(x) as.vector(x)) %>% bind_rows %>% t
print("start corr1")
c <- corrr::correlate(b)
print("correlation1 done")
d <- as.matrix(c[,-1])
n_corr_hist <- d
n_corr_sum <- mean(d, na.rm = T)
n_corr_sd <- sd(d, na.rm = T)
n_lambda <- lambda$no_auto_lambda

rm(b,c,d)

print("start corr2")
e <- lapply(lambda$M_ipm_g_auto[[1]], function(x) as.vector(x)) %>% bind_rows %>% t
print("correlation2 done")
f <- corrr::correlate(e)
g <- as.matrix(f[,-1])
a_corr_hist <- g
a_corr_sum <- mean(g, na.rm = T)
a_corr_sd <- sd(g, na.rm = T)
a_lambda <- lambda$ipm_g_auto_lambda


df <- tibble(clim_corr[taskID], clim_sd[taskID],
             n_lambda, a_lambda,
             n_corr_sum, a_corr_sum,
             n_corr_sd, a_corr_sd,
             list(n_corr_hist), list(a_corr_hist))
str(df)
output

saveRDS(df, file = output)

Sys.time() - start

