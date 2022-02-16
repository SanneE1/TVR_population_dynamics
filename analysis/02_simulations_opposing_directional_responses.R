#
# Script to run the main simulations on effect of both lagged and recent responses to climate 
# in a matrix population model, where climate responses are in OPPOSING direction
#
# This script includes arguments so it can be submitted on the UFZ HPC "EVE"
# Also see "submit_simulation.sh" in the same folder for the job submission script
#

rm(list=ls())

set.seed(2)

library(dplyr)
library(tidyr)
library(popbio)
library(parallel)
library(faux)
library(boot)

# Get required arguments supplied during job submission
args = commandArgs(trailingOnly = T)

if(length(args)!=2) {  
stop("Provide (only) signal strength for the analysis", call.=F)
}

# Set the proportion of variance (p) explained by the climate driver in the temporal sequence.
# Here set to 1, 0.5, 0.25 or 0.05. See manuscript for more info
i = as.numeric(args[1])

# File location of the simulation script with SAME DIRECTIONAL responses to climate
output_dir <- args[2]

source_file <- "/gpfs1/data/lagged/analysis/01_simulations_same_directional_responses.R"

print(paste("signal strength =", i))
print(paste("output dir =", output_dir))
print(paste("source file =", source_file))


### Create line sourcing to
source_lines <- function(file, lines){
  source(textConnection(readLines(file)[lines]))
}

### Create mpm function with +P (s&g) and - F
mpm <- function(survival, growth, reproduction,
                clim_sd, sig.strength) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)

  # growth
  g_mean = 0.325
  g_sd = 0.118

  if(is.na(growth)) {
    mpm[2,1] <- g_mean
  } else {
    ## total deviation from mean = climate signal * signal strength & correction factor (partitioning at variance scale) + random noise * signal strength
    dev <- growth$clim * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + 
      rnorm(1,0, clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    ## Because of partitioning and correction factor above, the resulting distribution has a sd of clim_sd
    p <- pnorm(dev, mean = 0, sd = clim_sd)
    mpm[2,1] <- qbeta(p, (((g_mean*(1-g_mean))/(g_sd * clim_sd)^2) - 1) * g_mean,
                      (((g_mean*(1-g_mean))/(g_sd * clim_sd)^2) - 1) * (1 - g_mean))
  }

  # survival
  s_mean = 0.541
  s_sd = 0.135
  if(is.na(survival)) {
    mpm[2,2] <- s_mean
  } else {
    dev <- survival$clim * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1,0, clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd)
    mpm[2,2] <- qbeta(p, (((s_mean*(1-s_mean))/(s_sd * clim_sd)^2) - 1) * s_mean,
                      (((s_mean*(1-s_mean))/(s_sd * clim_sd)^2) - 1) * (1 - s_mean))
  }

  # reproduction
  f_mean = 0.788
  f_sd = 0.862
  if(is.na(reproduction)) {
    mpm[1,2] <- f_mean
  } else {
    dev <- (-reproduction$clim) * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1,0, clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd)
    mpm[1,2] <- qgamma(p, (f_mean^2)/(f_sd * clim_sd)^2, (f_mean)/(f_sd * clim_sd)^2)
  }

  return(mpm)
}

## Run simulations ----------------------------------------

# Get create_seq(), st.lamb(), n_it and set up parallel. Source lines from 01 script
print("retrieving functions and running simulations")
source_lines(source_file, c(39:68, 121:212))
