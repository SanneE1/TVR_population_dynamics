rm(list=ls())
library(dplyr)
library(tidyr)
library(popbio)
library(parallel)
library(ggplot2)
library(faux)
library(boot)
set.seed(2)

args = commandArgs(trailingOnly = T)

if(length(args)!=2) {
  stop("Provide (only) signal strength for the analysis", call.=F)
}


i = as.numeric(args[1])
output_dir <- args[2]

print(paste("sig.strength =", i))
print(paste("output dir =", output_dir))

n_it = 50000
print(paste("#iterations =", n_it))

# Create function that creates environmental sequence ------------------------------------------------------------------------------------
### creates a sequence of climate anomalies whith a specified standard deviation and autocorrelation.
### the function then creates another sequence of the same length and standard deviation, to be used as random noise
### of each sequence two vectors are created, one which is the "original" and a second, which is offset by a specified period
### to create a lagged sequence
create_seq <- function(n_it, clim_sd, clim_auto, lag) { 
  for(n in c(1:(n_it+lag))){ 
    if(n == 1) {
      seq <- rnorm(1)
    } else {
      seq[n] <- clim_auto * seq[n-1] + rnorm(1)
    }
  }
  seq <- scale(seq) * clim_sd
  ran <- rnorm(n_it + lag, mean = 0, sd = clim_sd)
  
  lagged <- c(rep(NA, lag), seq)
  lagged_ran <- c(rep(NA, lag), ran)
  
  recent <- c(seq, rep(NA, lag))
  recent_ran <- c(ran, rep(NA, lag))
  
  df <- list(recent = data.frame(clim = recent,
                                 ran = recent_ran),
             lagged = data.frame(clim = lagged,
                                 ran = lagged_ran))
  return(df)
}

# Function to create MPM for each iteration. ------------------------------------------------------------------------------------
### This function takes SINGLE variables for the 3 climate drivers and 3 random noise values. 
### Using the method of moments, this function creates distribution of specified mean and sd values for each of the 3 vital rates. 
### The climate driver and the random noise for each VR are added to create a "deviation value". p calculates the probability of this "deviation value".
### We take that probability to calculate a value in the VR distributions

mpm <- function(survival, growth, reproduction, 
                clim_sd, sig.strength) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)
  
  # growth 
  g_mean = 0.325
  g_sd = 0.118
  
  if(is.na(growth$clim)) {
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
  if(is.na(survival$clim)) {
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
  if(is.na(reproduction$clim)) {
    mpm[1,2] <- f_mean
  } else {
    dev <- reproduction$clim * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + reproduction$ran * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    mpm[1,2] <- qgamma(p, (f_mean^2)/(f_sd * clim_sd)^2, (f_mean)/(f_sd * clim_sd)^2)
  }
  
  return(mpm)  
}

# Function to create a MPM sequence from supplied climate and noise sequences and calculate stocastic lambda -----------------------------------

st.lamb <- function(env_surv, env_growth, env_reproduction, 
                    clim_sd, clim_auto, sig.strength) {
  
  n_it = length(env_surv[,1])
  
  env <- list(survival = lapply(as.list(1:n_it), function(x) env_surv[x[1],]),
              growth = lapply(as.list(1:n_it), function(x) env_growth[x[1],]),
              reproduction = lapply(as.list(1:n_it), function(x) env_reproduction[x[1],]),
              clim_sd = clim_sd,
              sig.strength = sig.strength)
  
  ### Get all mpm's
  mats <- pmap(env, mpm) %>% Filter(Negate(anyNA), .)
  
  
  df <- data.frame(lambda = stoch.growth.rate(mats, maxt = n_it, verbose = F)$sim,
                   clim_sd = clim_sd,
                   clim_auto = clim_auto)
  
  return(list(df = df,
              mats = sapply(mats, as.vector) %>% t %>% 
                `colnames<-`(c("1,1", "2,1", "1,2", "2,2"))))
}

# Set up parallel runs
cl <- makeForkCluster(outfile = "")
## export libraries to workers
suppressMessages(clusterEvalQ(cl, c(library(popbio), library(dplyr), library(purrr))))


### autocorrelation in growth
clim_sd <- rep(seq(from = 0.01, to = 2, length.out = 10), 90)
clim_auto <- rep(rep(c(-0.9,0,0.9), each = 10), 30)

## export objects to workers
suppressMessages(clusterExport(cl, c("i", "create_seq", "inv.logit", "mpm", "st.lamb", "clim_auto", "clim_sd", "n_it")))


auto <- parLapply(cl = cl,
                  as.list(c(1:900)),
                  function(x) st.lamb(env_surv = create_seq(n_it = n_it, clim_sd[x], clim_auto[x], 0)[["recent"]],
                                      env_growth = create_seq(n_it = n_it, clim_sd[x], clim_auto[x], 0)[["recent"]],
                                      env_reproduction = data.frame(clim = rep(NA, n_it),
                                                                    ran = rep(NA, n_it)),
                                      clim_sd = clim_sd[x],
                                      clim_auto = clim_auto[x],
                                      sig.strength = i)
)

str(auto)

saveRDS(auto, file.path(output_dir, paste("mpm_", i, "_auto.RDS", sep = "")))



### Covariance all U vr
cov_clim <- function(x) {
  clim = rnorm_multi(n = n_it,
                     vars = 2,
                     mu = c(0,0),
                     sd = clim_sd[x],
                     r = clim_auto[x],
                     varnames = c("surv", "growth"))
  ran = rnorm(n = n_it, mean = 0, sd = clim_sd[x])
  
  df <- list(surv = data.frame(clim = clim$surv,
                               ran = ran),
             growth = data.frame(clim = clim$growth,
                                 ran = ran))
  return(df)
}

clim <- lapply(as.list(c(1:900)), cov_clim)

clusterExport(cl, c("clim"))

cov <- parLapply(cl = cl,
                 clim,
                 function(x) tryCatch(st.lamb(env_surv = x[["surv"]],
                                              env_growth = x[["growth"]],
                                              env_reproduction = data.frame(clim = rep(NA, n_it),
                                                                            ran = rep(NA, n_it)),
                                              clim_sd = sd(x[["surv"]]$clim),
                                              clim_auto = cor(x[["surv"]]$clim, x[["growth"]]$clim),
                                              sig.strength = i),
                                      error=function(err) NA)
)

str(cov)

saveRDS(cov, file.path(output_dir, paste("mpm_", i, "_cov.RDS", sep = "")))


#### Create main temporal sequences
lag_clim <- lapply(as.list(c(1:900)), function(x) create_seq(n_it, clim_sd = clim_sd[x], clim_auto = clim_auto[x], lag = 1))
suppressMessages(clusterExport(cl, c("lag_clim")))

#### Lagged effect within U matrix
lag_g <- parLapply(cl = cl,
                   lag_clim,
                   function(x) tryCatch(st.lamb(env_surv = x[["recent"]],
                                                env_growth = x[["lagged"]],
                                                env_reproduction = data.frame(clim = rep(NA, n_it),
                                                                              ran = rep(NA, n_it)),
                                                clim_sd = sd(x[["recent"]]$clim, na.rm = T),
                                                clim_auto = acf(x[["recent"]]$clim, plot = F, na.action = na.pass)$acf[2],
                                                sig.strength = i),
                                        error=function(err) NA)
)

lag_s <- parLapply(cl = cl,
                   lag_clim,
                   function(x) tryCatch(st.lamb(env_surv = x[["lagged"]],
                                                env_growth = x[["recent"]],
                                                env_reproduction = data.frame(clim = rep(NA, n_it),
                                                                              ran = rep(NA, n_it)),
                                                clim_sd = sd(x[["recent"]]$clim, na.rm = T),
                                                clim_auto = acf(x[["recent"]]$clim, plot = F, na.action = na.pass)$acf[2],
                                                sig.strength = i),
                                        error=function(err) NA)
)

lag_n <- parLapply(cl = cl,
                   lag_clim,
                   function(x) tryCatch(st.lamb(env_surv = x[["recent"]],
                                                env_growth = x[["recent"]],
                                                env_reproduction = data.frame(clim = rep(NA, n_it),
                                                                              ran = rep(NA, n_it)),
                                                clim_sd = sd(x[["recent"]]$clim, na.rm = T),
                                                clim_auto = acf(x[["recent"]]$clim, plot = F, na.action = na.pass)$acf[2],
                                                sig.strength = i),
                                        error=function(err) NA)
)

lag <- list("growth" = lag_g, "survival" = lag_s, "none" = lag_n)

str(lag)

saveRDS(lag, file.path(output_dir, paste("mpm_", i, "_lag.RDS", sep = "")))

#### Lagged effect between U & F matrices
lag_p <- parLapply(cl = cl,
                   lag_clim,
                   function(x) tryCatch(st.lamb(env_surv = x[["lagged"]],
                                                env_growth = x[["lagged"]],
                                                env_reproduction = x[["recent"]],
                                                clim_sd = sd(x[["recent"]]$clim, na.rm = T),
                                                clim_auto = acf(x[["recent"]]$clim, plot = F, na.action = na.pass)$acf[2],
                                                sig.strength = i),
                                        error=function(err) NA)
)

lag_f <- parLapply(cl = cl,
                   lag_clim,
                   function(x) tryCatch(st.lamb(env_surv = x[["recent"]],
                                                env_growth = x[["recent"]],
                                                env_reproduction = x[["lagged"]],
                                                clim_sd = sd(x[["recent"]]$clim, na.rm = T),
                                                clim_auto = acf(x[["recent"]]$clim, plot = F, na.action = na.pass)$acf[2],
                                                sig.strength = i),
                                        error=function(err) NA)
)

lag_n2 <- parLapply(cl = cl,
                    lag_clim,
                    function(x) tryCatch(st.lamb(env_surv = x[["recent"]],
                                                 env_growth = x[["recent"]],
                                                 env_reproduction = x[["recent"]],
                                                 clim_sd = sd(x[["recent"]]$clim, na.rm = T),
                                                 clim_auto = acf(x[["recent"]]$clim, plot = F, na.action = na.pass)$acf[2],
                                                 sig.strength = i),
                                         error=function(err) NA)
)

lag_fp <- list("Umatrix" = lag_p, "Fmatrix" = lag_f, "none" = lag_n2)

str(lag_fp)

saveRDS(lag_fp, file.path(output_dir, paste("mpm_", i, "_lagfp.RDS", sep = "")))

stopCluster(cl)


