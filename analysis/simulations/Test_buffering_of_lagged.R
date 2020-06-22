# Simulate the effect of lagged climate drivers on stochastic population growth rate
# This script if for climate independent F kernel, P kernel reacting to the same climate driver

library(ipmr)
library(lme4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)

source("ipmr_functions.R")

## load models -------------------------------------------------------------------------
survC <- readRDS("HEQU_s_climate.rds")
growthC <- readRDS("HEQU_g_climate.rds")

fp <- readRDS("HEQU_fp.rds")
fn <- readRDS("HEQU_fn.rds")
fpC <- readRDS("HEQU_fp_climate.rds")
fnC <- readRDS("HEQU_fn_climate.rds")



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


cores <- detectCores()
cl <- makeCluster(cores[1] - 2)
clusterExport(cl, c("P_neg_1yr", "P_1yr", "P_2yr", "PF_1yr", "params_list"))
clusterEvalQ(cl, c(library(ipmr), library(dplyr)))


clim_sd <- as.list(rep(seq(from = 0, to = 2, length.out = 16), 20))

lPneg <- parLapply(cl, clim_sd, P_neg_1yr) %>% bind_rows
lP1 <- parLapply(cl, clim_sd, P_1yr) %>% bind_rows
lP2 <- parLapply(cl, clim_sd, P_2yr) %>% bind_rows
lFP1 <- parLapply(cl, clim_sd, PF_1yr) %>% bind_rows

stopCluster(cl)


