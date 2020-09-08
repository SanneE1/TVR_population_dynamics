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
    metavar = "P_1yr|P_neg_1yr|P_2yr|PF_1yr|PF_neg_1yr")
  
)

parser <- OptionParser(
  usage       = "Rscript %prog [options] output meshpoints iterations",
  option_list = Parsoptions,
  description = "\nan Run lagged ipmr functions",
  epilogue    = ""
)

cli <- parse_args(parser, positional_arguments = 1)

foption <- cli$options$foption
taskID <- as.integer(Sys.getenv("SGE_TASK_ID"))
output <- cli$args[1]
n_mesh <- cli$args[2]
n_it <- cli$args[3]

foption
taskID

source("/home/evers/lagged_buffering/analysis/simulations/ipmr_functions.R")


## create parameter list -------------------------------------------------------------------------


params_list <- list(
  s_int = -0.229,
  s_slope = 1.077,
  # s_temp = 1.233,
  g_int = 0.424,
  g_slope = 0.846,
  # g_temp = -0.066,
  g_sd = 1.076,
  fp_int = -3.970,
  fp_slope = 1.719,
  fpC_int = -3.859385,
  fpC_slope = 1.719768,
  fpC_temp = -0.6169492,
  fn_int = -0.6652477,
  fn_slope = 0.7048809,
  fnC_int = -0.5661762,
  fnC_slope = 0.7048782,
  fnC_temp = -0.3398345,
  germ_mean = 0.1262968,
  germ_sd = 0.2725941,
  fd_mean = 1.178749,
  fd_sd = 0.8747764
)



clim_sd <- rep(seq(from = 0, to = 2, length.out = 16), 60)
clim_corr <- rep(rep(c(-0.9,0,0.9), each = 16), 20)

clim_sd[taskID]
clim_corr[taskID]

if(foption == "P_1yr") {
  message("Within P kernel 1 year lagged")
  lambda <- P_lambdas(n_it = n_it, 
                  clim_sd = clim_sd[taskID], 
                  clim_corr = clim_corr[taskID], 
                  params_list = params_list, 
                  clim_params = list(s_temp = 1.233, 
                                     g_temp = 0.066),
                  n_mesh = n_mesh)
  
  
  b <- lapply(lambda$M_non_lagged[[1]], function(x) as.vector(x)) %>% bind_rows %>% t
  c <- corrr::correlate(b)
  d <- as.matrix(c[,-1])
  n_corr_hist[[n]] <- d
  n_corr_sum[n] <- mean(d, na.rm = T)
  n_lambda[n] <- lambda$non_lagged
  
  e <- lapply(lambda$M_s_lagged[[1]], function(x) as.vector(x)) %>% bind_rows %>% t
  f <- corrr::correlate(e)
  g <- as.matrix(f[,-1])
  s_corr_hist[[n]] <- g
  s_corr_sum[n] <- mean(g, na.rm = T)
  s_lambda[n] <- lambda$lagged_s
  
  h <- lapply(lambda$M_g_lagged[[1]], function(x) as.vector(x)) %>% bind_rows %>% t
  i <- corrr::correlate(h)
  j <- as.matrix(i[,-1])
  g_corr_hist[[n]] <- j
  g_corr_sum[n] <- mean(j, na.rm = T)
  g_lambda[n] <- lambda$lagged_g
} 

if(foption == "P_neg_1yr"){
  lambda <- P_lambdas(n_it = n_it, 
                      clim_sd = clim_sd[taskID], 
                      clim_corr = clim_corr[taskID], 
                      params_list = params_list, 
                      clim_params = list(s_temp = -1.233, 
                                         g_temp = -0.066),
                      n_mesh = n_mesh)
  
  b <- lapply(lambda$M_non_lagged[[1]], function(x) as.vector(x)) %>% bind_rows %>% t
  c <- corrr::correlate(b)
  d <- as.matrix(c[,-1])
  n_corr_hist[[n]] <- d
  n_corr_sum[n] <- mean(d, na.rm = T)
  n_lambda[n] <- lambda$non_lagged
  
  e <- lapply(lambda$M_s_lagged[[1]], function(x) as.vector(x)) %>% bind_rows %>% t
  f <- corrr::correlate(e)
  g <- as.matrix(f[,-1])
  s_corr_hist[[n]] <- g
  s_corr_sum[n] <- mean(g, na.rm = T)
  s_lambda[n] <- lambda$lagged_s
  
  h <- lapply(lambda$M_g_lagged[[1]], function(x) as.vector(x)) %>% bind_rows %>% t
  i <- corrr::correlate(h)
  j <- as.matrix(i[,-1])
  g_corr_hist[[n]] <- j
  g_corr_sum[n] <- mean(j, na.rm = T)
  g_lambda[n] <- lambda$lagged_g
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
if(foption == "PF_neg_1yr") {
  lambda <- PF_lambdas(n_it = 10000, 
                       clim_sd = clim_sd[taskID], 
                       clim_corr = clim_corr[taskID], 
                       params_list = params_list, 
                       clim_params = list(s_temp = -1.233, 
                                          g_temp = -0.066,
                                          fpC_temp = -0.617,
                                          fnC_temp = -0.34))}

if(foption == "P_2yr") {
  lambda <- P_2yr(n_it = 10000, 
                  clim_sd = clim_sd[taskID], 
                  clim_corr = clim_corr[taskID])
} 

saveRDS(lambda, file = output)

Sys.time() - start
