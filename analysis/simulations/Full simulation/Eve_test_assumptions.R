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

source("/home/evers/lagged_buffering/analysis/simulations/Full simulation/ipmr_functions.R")
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


if(foption == "PF_1yr") {
  lambda <- PF_lambdas(n_it = 10000, 
                       clim_sd = clim_sd[taskID], 
                       clim_corr = clim_corr[taskID], 
                       params_list = params_list, 
                       clim_params = list(s_temp = 1.233, 
                                          g_temp = -0.066,
                                          fpC_temp = -0.617,
                                          fnC_temp = -0.34))
}

if(foption == "P_2yr") {
  lambda <- P_2yr(n_it = 10000, 
                  clim_sd = clim_sd[taskID], 
                  clim_corr = clim_corr[taskID])
} 


print("ipm done")  
b <- lapply(lambda$M_non_lagged[[1]], function(x) as.vector(x)) %>% bind_rows %>% t
print("start corr1")
c <- corrr::correlate(b)
print("correlation1 done")
d <- as.matrix(c[,-1])
n_corr_hist <- d
n_corr_sum <- mean(d, na.rm = T)
n_corr_sd <- sd(d, na.rm = T)
n_lambda <- lambda$non_lagged

rm(b,c,d)

print("start corr2")
e <- lapply(lambda$M_s_lagged[[1]], function(x) as.vector(x)) %>% bind_rows %>% t
print("correlation2 done")
f <- corrr::correlate(e)
g <- as.matrix(f[,-1])
s_corr_hist <- g
s_corr_sum <- mean(g, na.rm = T)
s_corr_sd <- sd(g, na.rm = T)
s_lambda <- lambda$lagged_s

rm(e,f,g)

h <- lapply(lambda$M_g_lagged[[1]], function(x) as.vector(x)) %>% bind_rows %>% t
print("start corr3")
i <- corrr::correlate(h)
print("correlation3 done")
j <- as.matrix(i[,-1])
g_corr_hist <- j
g_corr_sum <- mean(j, na.rm = T)
g_corr_sd <- sd(j, na.rm = T)
g_lambda <- lambda$lagged_g

rm(h,i,j)

df <- tibble(clim_corr, clim_sd,
                 n_lambda, s_lambda, g_lambda,
                 n_corr_sum, s_corr_sum, g_corr_sum,
                 n_corr_sd, s_corr_sd, g_corr_sd,
             n_corr_hist, s_corr_hist, g_corr_hist)
str(df)
output

saveRDS(df, file = output)

Sys.time() - start
