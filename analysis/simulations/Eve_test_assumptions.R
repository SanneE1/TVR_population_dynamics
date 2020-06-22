suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ipmr))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))



start <- Sys.time()
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
  usage       = "Rscript %prog [options] output",
  option_list = Parsoptions,
  description = "\nan Run lagged ipmr functions",
  epilogue    = ""
)

cli <- parse_args(parser, positional_arguments = 1)

foption <- cli$options$foption
taskID <- as.integer(Sys.getenv("SGE_TASK_ID"))
output <- cli$args[1]

foption
taskID

source("/home/evers/lagged/ipmr_functions.R")

## load models -------------------------------------------------------------------------
survC <- readRDS("/data/gsclim/ipmr_lagged_test/HEQU_s_climate.rds")
growthC <- readRDS("/data/gsclim/ipmr_lagged_test/HEQU_g_climate.rds")

fp <- readRDS("/data/gsclim/ipmr_lagged_test/HEQU_fp.rds")
fn <- readRDS("/data/gsclim/ipmr_lagged_test/HEQU_fn.rds")
fpC <- readRDS("/data/gsclim/ipmr_lagged_test/HEQU_fp_climate.rds")
fnC <- readRDS("/data/gsclim/ipmr_lagged_test/HEQU_fn_climate.rds")



## create parameter list -------------------------------------------------------------------------

params_list <- list(
  s_int = fixef(survC)[1],
  s_slope = fixef(survC)[2],
  s_temp = fixef(survC)[5],
  g_int = fixef(growthC)[1],
  g_slope = fixef(growthC)[2],
  g_temp = fixef(growthC)[5],
  g_sd = sd(resid(growthC)),
  fp_int = fixef(fp)[1],
  fp_slope = fixef(fp)[2],
  fpC_int = fixef(fpC)[1],
  fpC_slope = fixef(fpC)[2],
  fpC_temp = fixef(fpC)[5],
  fn_int = fixef(fn)[1],
  fn_slope = fixef(fn)[2],
  fnC_int = fixef(fnC)[1],
  fnC_slope = fixef(fnC)[2],
  fnC_temp = fixef(fnC)[5],
  germ_mean = 0.1262968,
  germ_sd = 0.2725941,
  fd_mean = 1.178749,
  fd_sd = 0.8747764
)

clim_sd <- as.list(rep(seq(from = 0, to = 2, length.out = 16), 20))

if(foption == "P_1yr") {
  lambda <- P_1yr(clim_sd[[taskID]])
} 
if(foption == "P_2yr") {
  lambda <- P_2yr(clim_sd[[taskID]])
} 
if(foption == "P_neg_1yr"){
  lambda <- P_neg_1yr(clim_sd[[taskID]])
}
if(foption == "PF_1yr") {
  lambda <- PF_1yr(clim_sd[[taskID]])
}
if(foption == "PF_neg_1yr") {
  lambda <- PF_neg_1yr(clim_sd[[taskID]])
}

lambda

saveRDS(lambda, file = output)

Sys.time() - start
