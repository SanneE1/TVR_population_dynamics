library(dplyr)
library(popbio)
library(pbapply)
library(parallel)
library(ggplot2)
library(faux)


inv_logit <- function(x) {
  return(
    1/(1 + exp(-(x)))
  )
}

### Create environmental sequence ----------------------------
create_seq <- function(n_it, clim_sd, clim_corr, lag) { 
  for(n in c(1:(n_it+lag))){ 
    if(n == 1) {
      seq <- rnorm(1)
    } else {
      seq[n] <- clim_corr * seq[n-1] + rnorm(1)
    }
  }
  seq <- scale(seq) * clim_sd
  lagged <- c(rep(NA, lag), seq)
  recent <- c(seq, rep(NA, lag))
  df <- data.frame(recent = recent,
                   lagged = lagged)
  return(df)
}


mpm <- function(survival, growth, reproduction) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)
  
  #growth
  mpm[2,1] <- inv_logit(growth)  ### inv_logit(0) = 0.5 (intercept)
  
  # survival
  mpm[2,2] <- inv_logit(survival)
  
  # reproduction 
  mpm[1,2] <- exp(1.2 + reproduction)
  
  return(mpm)  
}

# "Stocastic" mpm -----------------------------------

st.lamb <- function(env_surv, env_growth, env_reproduction) {
  
  n_it = length(env_surv)
  
  env <- data.frame(survival = env_surv,
                    growth = env_growth,
                    reproduction = env_reproduction)
  env <- env[complete.cases(env), ]
  
  env <- as.list(as.data.frame(t(env)))
  
  ### Get all mpm's
  mats <- lapply(env, function(x) mpm(x[1], x[2], x[3]))
  
  
  df = stoch.growth.rate(mats, maxt = 1000)$sim
  
  return(df)
}


for(i in c(-1,1)) {   #### i means that I can run the same code for different climate effects without manual changes
  
  print(i)
  
  # Set up parallel runs
  no_cores <- detectCores()
  cl <- makeCluster(no_cores - 2)
  ## export libraries to workers
  clusterEvalQ(cl, c(library(popbio), library(dplyr)))
  ## provide seed number to workers
  # clusterSetRNGStream(cl, iseed = 04103)
  
  ### autocorrelation in growth
  
  clim_sd <- rep(seq(from = 0, to = 2, length.out = 10), 90)
  clim_corr <- rep(rep(c(-0.9,0,0.9), each = 10), 30)
  
  ## export objects to workers
  clusterExport(cl, c("i", "create_seq", "inv_logit", "mpm", "st.lamb", "clim_corr", "clim_sd"))
  

  auto <- pblapply(cl = cl,
                   as.list(c(1:900)),
                   function(x) st.lamb(env_surv = rep(0, 50000) * i,
                                       env_growth = create_seq(n_it = 50000, clim_sd[x], clim_corr[x], 0) * i,
                                       env_reproduction = rep(0,50000) * i)
  )

  saveRDS(auto, paste("results/simulations/mpm/mpm_", i, "_auto.RDS", sep = ""))



  ### Covariance all vr
  clim <- lapply(as.list(c(1:900)), function(x)
    rnorm_multi(n = 50000,
                vars = 2,
                mu = c(0,0),
                sd = clim_sd[x],
                r = clim_corr[x],
                varnames = c("surv", "growth")) )

  clusterExport(cl, c("clim"))

  cov <- pblapply(cl = cl,
                  clim,
                  function(x) st.lamb(env_surv = x$surv * i,
                                      env_growth = x$growth * i,
                                      env_reproduction = rep(0,50000) * i)
  )

  saveRDS(cov, paste("results/simulations/mpm/mpm_", i, "_cov.RDS", sep = ""))

  #### Lagged effect within P "functions"
  
  lag_clim <- lapply(as.list(c(1:900)), function(x) create_seq(50000, clim_sd = clim_sd[x], clim_corr = clim_corr[x], lag = 1))
  
  clusterExport(cl, c("lag_clim"))
  
  lag_g <- pblapply(cl = cl,
                    lag_clim,
                    function(x) st.lamb(env_surv = x$recent * i,
                                        env_growth = x$lagged * i,
                                        env_reproduction = rep(0,length(x$recent)) * i)
  )

  lag_s <- pblapply(cl = cl,
                    lag_clim,
                    function(x) st.lamb(env_surv = x$lagged * i,
                                        env_growth = x$recent * i,
                                        env_reproduction = rep(0,length(x$recent)) * i)
  )

  lag_n <- pblapply(cl = cl,
                    lag_clim,
                    function(x) st.lamb(env_surv = x$recent * i,
                                        env_growth = x$recent * i,
                                        env_reproduction = rep(0,length(x$recent)) * i)
  )

  lag <- list("growth" = lag_g, "survival" = lag_s, "none" = lag_n)

  saveRDS(lag, paste("results/simulations/mpm/mpm_", i, "_lag.RDS", sep = ""))

  #### Lagged effect
  
  lag_p <- pblapply(cl = cl, 
                    lag_clim, 
                    function(x) st.lamb(env_surv = x$lagged * i,
                                        env_growth = x$lagged * i,
                                        env_reproduction = x$recent * i)
  )
  
  lag_f <- pblapply(cl = cl, 
                    lag_clim, 
                    function(x) st.lamb(env_surv = x$recent * i,
                                        env_growth = x$recent * i,
                                        env_reproduction = x$lagged * i)
  )
  
  lag_n <- pblapply(cl = cl, 
                    lag_clim, 
                    function(x) st.lamb(env_surv = x$recent * i,
                                        env_growth = x$recent * i,
                                        env_reproduction = x$recent * i)
  )
  
  lag_fp <- list("Pkernel" = lag_p, "Fkernel" = lag_f, "none" = lag_n)
  
  saveRDS(lag_fp, paste("results/simulations/mpm/mpm_", i, "_lagfp.RDS", sep = ""))
  
  stopCluster(cl)
}
