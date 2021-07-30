### Similar simulations as in 01 folder, but now vital rates respond in different direction to climate signals
rm(list=ls())
library(dplyr)
library(tidyr)
library(popbio)
library(parallel)
library(faux)
library(boot)

set.seed(2)

args = commandArgs(trailingOnly = T)

if(length(args)!=3) {  ## submit script provides two more arguments (output and source location)
stop("Provide (only) signal strength for the analysis", call.=F)
}

i = as.numeric(args[1])
source_file <- args[2]
output_dir <- args[3]

print(paste("signal strength =", i))
print(paste("source file =", source_file))
print(paste("output dir =", output_dir))

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
      growth$ran * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
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
    dev <- survival$clim * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + survival$ran * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
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
    dev <- (-reproduction$clim) * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + reproduction$ran * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd)
    mpm[1,2] <- qgamma(p, (f_mean^2)/(f_sd * clim_sd)^2, (f_mean)/(f_sd * clim_sd)^2)
  }

  return(mpm)
}

## Run simulations ----------------------------------------

# Get create_seq(), st.lamb(), n_it and set up parallel. Source lines from 01 folder
print("Getting create_seq(), st.lamb(), n_it and set up parallel")
source_lines(source_file, c(24:55, 109:147))

# Create climate sequences
print("creating temporal sequences")
source_lines(source_file, c(204:207))

# Run simulations for U and F matrix
print("running simulations")
source_lines(source_file, c(251:293))


