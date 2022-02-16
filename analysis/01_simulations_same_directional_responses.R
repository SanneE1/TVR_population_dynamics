#
# Script to run the main simulations on effect of both lagged and recent responses to climate 
# in a matrix population model, where climate responses are in the SAME direction
#
# This script includes arguments so it can be submitted on the UFZ HPC "EVE"
# Also see "submit_simulation.sh" in the same folder for the job submission script

rm(list=ls())

set.seed(2)

library(dplyr)
library(tidyr)
library(popbio)
library(parallel)
library(ggplot2)
library(faux)
library(boot)
library(purrr)

# Get required arguments supplied during job submission
args = commandArgs(trailingOnly = T)

if(length(args)!=2) {
  stop("Provide (only) signal strength for the analysis", call.=F)
}

# Set the proportion of variance (p) explained by the climate driver in the temporal sequence.
# Here set to 1, 0.5, 0.25 or 0.05. See manuscript for more info
i = as.numeric(args[1])

# Where to save the output
output_dir <- args[2]

# Print information for the log file
print(paste("sig.strength =", i))
print(paste("output dir =", output_dir))

# Number of iterations to run the simulations
n_it = 50000
print(paste("#iterations =", n_it))

# Create function that creates environmental sequence ------------------------------------------------------------------------------------

## Creates a sequence of climate anomalies (c in manuscript) with a specified standard deviation and autocorrelation.
## The function then creates another sequence of the same length and standard deviation, to be used as random noise.
## Next the two vectors are combined in two dataframes, one which is the original (and recent) dataframe and a second, 
## which is offset by a specified period to create a lagged sequence

create_seq <- function(n_it, clim_sd, clim_auto, lag) { 
  for(n in c(1:(n_it+lag))){ 
    if(n == 1) {
      seq <- rnorm(1)
    } else {
      seq[n] <- clim_auto * seq[n-1] + rnorm(1)
    }
  }
  seq <- scale(seq) * clim_sd

  lagged <- c(rep(NA, lag), seq)

  recent <- c(seq, rep(NA, lag))

  df <- data.frame(recent = recent,
                   lagged = lagged)
  return(df)
}

# Function to create MPM for each iteration. ------------------------------------------------------------------------------------

## This function creates a single matrix population model (MPM) based on the climate & random values given for each of the vital rates 
## Given the vital rate means and sds, as well as the climate sd, this function calculates the MPM cell values, using the method of moments 

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
    dev <- growth * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + 
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
    dev <- survival * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1,0, clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
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
    dev <- reproduction * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1,0, clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    mpm[1,2] <- qgamma(p, (f_mean^2)/(f_sd * clim_sd)^2, (f_mean)/(f_sd * clim_sd)^2)
  }
  
  return(mpm)  
}

# Function to create a MPM sequence from supplied climate and noise sequences and calculate stocastic lambda -----------------------------------

## This function will create a MPM for each iteration of the environmental sequence, and calculate the stochastic lambda of said MPM sequence 
st.lamb <- function(env_surv, env_growth, env_reproduction, 
                    clim_sd, clim_auto, sig.strength) {
  
  n_it = length(env_surv)
  
  env <- list(survival = env_surv,
              growth = env_growth,
              reproduction = env_reproduction,
              clim_sd = clim_sd,
              sig.strength = sig.strength)
  
  ### Get all mpm's
  mats <- purrr::pmap(env, mpm) %>% Filter(Negate(anyNA), .)
  mats <- mats %>% purrr::discard(function(x) all(x == 0))
  
  df <- data.frame(lambda = stoch.growth.rate(mats, maxt = n_it, verbose = F)$sim,
                   clim_sd = clim_sd,
                   clim_auto = clim_auto)
  
  return(list(df = df,
              mats = NA)) #sapply(mats, as.vector) %>% t %>% 
                #`colnames<-`(c("1,1", "2,1", "1,2", "2,2"))))
}

#---------------------------------------------------------------------------------------------------------
# Start of analyses
#---------------------------------------------------------------------------------------------------------

## set up sequences of climate standard deviation and autocorrelation values so that there are 30 duplicates 
## for each combination of sd & autocorrelation value
clim_sd <- rep(seq(from = 0.01, to = 2, length.out = 10), 90)
clim_auto <- rep(rep(c(-0.9,0,0.9), each = 10), 30)

# Set up parallel runs
cl <- makeForkCluster()
suppressForeignCheck(clusterEvalQ(cl, c(library(popbio), library(dplyr), library(purrr))))

# export objects to workers
suppressMessages(clusterExport(cl, c("i", "create_seq", "inv.logit", "mpm", "st.lamb", "clim_auto", "clim_sd", "n_it")))

#### Create main temporal sequences
lag_clim <- parLapplyLB(cl = cl,
                        as.list(c(1:900)), 
                        function(x) create_seq(n_it, clim_sd = clim_sd[x], clim_auto = clim_auto[x], lag = 1))
suppressMessages(clusterExport(cl, c("lag_clim")))

#### Lagged effect between U & F matrices
lag_p <- parLapplyLB(cl = cl,
                   lag_clim,
                   function(x) tryCatch(st.lamb(env_surv = x$lagged,
                                                env_growth = x$lagged,
                                                env_reproduction = x$recent,
                                                clim_sd = sd(x$recent, na.rm = T),
                                                clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                                sig.strength = i),
                                        error=function(err) NA)
)

lag_f <- parLapplyLB(cl = cl,
                   lag_clim,
                   function(x) tryCatch(st.lamb(env_surv = x$recent,
                                                env_growth = x$recent,
                                                env_reproduction = x$lagged,
                                                clim_sd = sd(x$recent, na.rm = T),
                                                clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                                sig.strength = i),
                                        error=function(err) NA)
)

lag_n2 <- parLapplyLB(cl = cl,
                    lag_clim,
                    function(x) tryCatch(st.lamb(env_surv = x$recent,
                                                 env_growth = x$recent,
                                                 env_reproduction = x$recent,
                                                 clim_sd = sd(x$recent, na.rm = T),
                                                 clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                                 sig.strength = i),
                                         error=function(err) NA)
)

lag_fp <- list("Umatrix" = lag_p, "Fmatrix" = lag_f, "none" = lag_n2)

print("done w/ UF sim")

saveRDS(lag_fp, file.path(output_dir, paste("mpm_", i, "_lagfp.RDS", sep = "")))

stopCluster(cl)


