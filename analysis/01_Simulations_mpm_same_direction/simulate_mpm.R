rm(list=ls())
library(dplyr)
library(tidyr)
library(popbio)
library(pbapply)
library(parallel)
library(ggplot2)
library(faux)
library(boot)
set.seed(2)

output_dir <- "results/01_Simulations_mpm_same_directions/rds/"

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

  lagged <- c(rep(NA, lag), seq)
  recent <- c(seq, rep(NA, lag))
  
  df <- data.frame(recent = recent,
                   lagged = lagged)
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
  
  ## survival Juveniles
  # sj_mean = 
  # sj_sd = 
  # if(is.na(survival)) {
  #   
  #   mpm[1,1] <- sj_mean
  #   
  # } else {
  #   
  #   dev <- survival .....  ## total deviation from mean = climate signal & random noise
  #   p <- pnorm(dev, mean = 0, sd = sqrt(clim_sd * 2)) ## Here I use sd = sqrt(2) because dev has variance var(growth) + var(noise). sqrt() to get stand. dev.
  #   q <- qnorm(p)
  #   mpm[1,1] <- qbeta(p, (((sj_mean*(1-sj_mean))/sj_sd^2) - 1) * sj_mean,
  #                     (((sj_mean*(1-sj_mean))/sj_sd^2) - 1) * (1 - sj_mean))
  #   
  # }
  
  # growth 
  g_mean = 0.325
  g_sd = 0.118
  
  if(is.na(growth)) {
    mpm[2,1] <- g_mean  
  } else {
    ## total deviation from mean = climate signal * signal strength & correction factor (partitioning at variance scale) + random noise * signal strength
    dev <- growth * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + 
      rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
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
    dev <- survival * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
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
    dev <- reproduction * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    mpm[1,2] <- qgamma(p, (f_mean^2)/(f_sd * clim_sd)^2, (f_mean)/(f_sd * clim_sd)^2)
  }
  
  return(mpm)  
}

# Function to create a MPM sequence from supplied climate and noise sequences and calculate stocastic lambda -----------------------------------

st.lamb <- function(env_surv, env_growth, env_reproduction, 
                    clim_sd, clim_auto, sig.strength) {
  
  n_it = length(env_surv)
  
  env <- list(survival = env_surv,
              growth = env_growth,
              reproduction = env_reproduction,
              clim_sd = clim_sd,
              sig.strength = sig.strength)
  
  ### Get all mpm's
  mats <- pmap(env, mpm) %>% Filter(Negate(anyNA), .)
  
  
  df <- data.frame(lambda = stoch.growth.rate(mats, maxt = n_it)$sim,
                   clim_sd = clim_sd,
                   clim_auto = clim_auto)
  
  return(list(df = df,
              mats = sapply(mats, as.vector) %>% t %>% 
                `colnames<-`(c("1,1", "2,1", "1,2", "2,2"))))
}

n_it = 5000

for(i in c(1, 0.5, 0.25, 0.05)) {   #### i means that I can run the same code for different climate effect strength without manual changes
  
  print(i)
  
  # Set up parallel runs
  no_cores <- detectCores()
  cl <- makeCluster(no_cores - 2)
  ## export libraries to workers
  clusterEvalQ(cl, c(library(popbio), library(dplyr), library(purrr)))
  
  
  ### autocorrelation in growth
  clim_sd <- rep(seq(from = 0.01, to = 2, length.out = 10), 90)
  clim_auto <- rep(rep(c(-0.9,0,0.9), each = 10), 30)
  
  ## export objects to workers
  clusterExport(cl, c("i", "create_seq", "inv.logit", "mpm", "st.lamb", "clim_auto", "clim_sd", "n_it"))
  
  
  auto <- pblapply(cl = cl,
                   as.list(c(1:900)),
                   function(x) tryCatch(st.lamb(env_surv = create_seq(n_it = n_it, clim_sd[x], clim_auto[x], 0)$recent,
                                                env_growth = create_seq(n_it = n_it, clim_sd[x], clim_auto[x], 0)$recent,
                                                env_reproduction = rep(NA,n_it),
                                                clim_sd = clim_sd[x],
                                                clim_auto = clim_auto[x],
                                                sig.strength = i),
                                        error=function(err) NA)
  )
  
  
  saveRDS(auto, file.path(output_dir, paste("mpm_", i, "_auto.RDS", sep = "")))
  
  
  
  ### Covariance all U vr
  clim <- lapply(as.list(c(1:900)), function(x)
    rnorm_multi(n = n_it,
                vars = 2,
                mu = c(0,0),
                sd = clim_sd[x],
                r = clim_auto[x],
                varnames = c("surv", "growth")) )
  
  clusterExport(cl, c("clim"))
  
  cov <- pblapply(cl = cl,
                  clim,
                  function(x) tryCatch(st.lamb(env_surv = x$surv,
                                               env_growth = x$growth,
                                               env_reproduction = rep(NA,length(x$surv)),
                                               clim_sd = sd(x$surv),
                                               clim_auto = cor(x$surv, x$growth),
                                               sig.strength = i),
                                       error=function(err) NA)
  )
  
  saveRDS(cov, file.path(output_dir, paste("mpm_", i, "_cov.RDS", sep = "")))
  
  
  #### Lagged effect within U matrix
  
  lag_clim <- lapply(as.list(c(1:900)), function(x) create_seq(n_it, clim_sd = clim_sd[x], clim_auto = clim_auto[x], lag = 1))
  
  clusterExport(cl, c("lag_clim"))
  
  lag_g <- pblapply(cl = cl,
                    lag_clim,
                    function(x) tryCatch(st.lamb(env_surv = x$recent,
                                                 env_growth = x$lagged,
                                                 env_reproduction = rep(NA,length(x$recent)),
                                                 clim_sd = sd(x$recent, na.rm = T),
                                                 clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                                 sig.strength = i),
                                         error=function(err) NA)
  )
  
  lag_s <- pblapply(cl = cl,
                    lag_clim,
                    function(x) tryCatch(st.lamb(env_surv = x$lagged,
                                                 env_growth = x$recent,
                                                 env_reproduction = rep(NA,length(x$recent)),
                                                 clim_sd = sd(x$recent, na.rm = T),
                                                 clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                                 sig.strength = i),
                                         error=function(err) NA)
  )
  
  lag_n <- pblapply(cl = cl,
                    lag_clim,
                    function(x) tryCatch(st.lamb(env_surv = x$recent,
                                                 env_growth = x$recent,
                                                 env_reproduction = rep(NA,length(x$recent)),
                                                 clim_sd = sd(x$recent, na.rm = T),
                                                 clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                                 sig.strength = i),
                                         error=function(err) NA)
  )
  
  lag <- list("growth" = lag_g, "survival" = lag_s, "none" = lag_n)
  
  saveRDS(lag, file.path(output_dir, paste("mpm_", i, "_lag.RDS", sep = "")))
  
  #### Lagged effect between U & F matrices
  
  lag_p <- pblapply(cl = cl,
                    lag_clim,
                    function(x) tryCatch(st.lamb(env_surv = x$lagged,
                                                 env_growth = x$lagged,
                                                 env_reproduction = x$recent,
                                                 clim_sd = sd(x$recent, na.rm = T),
                                                 clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                                 sig.strength = i),
                                         error=function(err) NA)
  )
  
  lag_f <- pblapply(cl = cl,
                    lag_clim,
                    function(x) tryCatch(st.lamb(env_surv = x$recent,
                                                 env_growth = x$recent,
                                                 env_reproduction = x$lagged,
                                                 clim_sd = sd(x$recent, na.rm = T),
                                                 clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                                 sig.strength = i),
                                         error=function(err) NA)
  )
  
  lag_n2 <- pblapply(cl = cl,
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
  
  saveRDS(lag_fp, file.path(output_dir, paste("mpm_", i, "_lagfp.RDS", sep = "")))
  
  stopCluster(cl)
  
}
