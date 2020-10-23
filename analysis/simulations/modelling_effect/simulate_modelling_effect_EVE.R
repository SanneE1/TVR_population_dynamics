suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ipmr))
suppressPackageStartupMessages(library(rethinking))

start <- Sys.time()
start

parser <- OptionParser(
  usage       = "Rscript %prog output",
  option_list = Parsoptions,
  description = "\nan Run population estimation simulations",
  epilogue    = ""
)

cli <- parse_args(parser, positional_arguments = 1)
output <- cli$args[1]
taskID <- as.integer(Sys.getenv("SGE_TASK_ID"))

taskID


source("/home/evers/lagged_buffering/analysis/simulations/modelling_effect/population_simulation_functions.R")

clim_sd <- rep(c(1:2), 50)
clim_cor <- rep(seq(from = -0.9, to = 0.9, length.out = 5), each = 20)

clim_cor[taskID]
clim_sd[taskID]


df <- wrapper(init.pop.size = 1000,   
              n_yrs = 35,
              clim_corr = clim_cor[taskID],
              clim_sd = clim_sd[taskID])


saveRDS(df, file = output)
