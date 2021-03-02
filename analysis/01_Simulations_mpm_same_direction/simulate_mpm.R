library(dplyr)
library(tidyr)
library(popbio)
library(pbapply)
library(parallel)
library(ggplot2)
library(faux)

output_dir <- "results/01_Simulations_mpm_same_directions/rds/"

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


mpm <- function(survival, growth, reproduction, clim_sd, sig.strength = 1) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)
  
  #growth                     Get the error term sd to reflect the sd of the environment sequence
  mpm[2,1] <- inv_logit(growth * sig.strength + ((1-sig.strength) * rnorm(1, 0, clim_sd)) ) ### inv_logit(0) = 0.5 (intercept)
  
  # survival
  mpm[2,2] <- inv_logit(survival + ((1-sig.strength) * rnorm(1, 0, clim_sd)) )
  
  # reproduction 
  mpm[1,2] <- exp(1.2 + reproduction + ((1-sig.strength) * rnorm(1, 0, clim_sd)) )
  
  return(mpm)  
}

# "Stocastic" mpm -----------------------------------

st.lamb <- function(env_surv, env_growth, env_reproduction, clim_sd, sig.strength) {
  
  n_it = length(env_surv)
  
  env <- data.frame(survival = env_surv,
                    growth = env_growth,
                    reproduction = env_reproduction)
  env <- env[complete.cases(env), ]
  
  env <- as.list(as.data.frame(t(env)))
  
  ### Get all mpm's
  mats <- lapply(env, function(x) mpm(x[1], x[2], x[3], clim_sd, sig.strength))
  
  
  df = stoch.growth.rate(mats, maxt = 1000)$sim
  
  return(list(df = df,
              mats = data.frame(t(lapply(mats, as.vector) %>% bind_rows)) %>% 
                `colnames<-`(c("1,1", "2,1", "1,2", "2,2"))))
}


for(i in c(1, 0.5, 0.25, 0.05)) {   #### i means that I can run the same code for different climate effect strength without manual changes
  
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
                   as.list(c(1:225)),
                   function(x) try(st.lamb(env_surv = rep(0, 50000),
                                       env_growth = create_seq(n_it = 50000, clim_sd[x], clim_corr[x], 0)$recent,
                                       env_reproduction = rep(0,50000),
                                       clim_sd = clim_sd[x],
                                       sig.strength = i))
  )
  
  auto_df <- lapply(auto, function(x) x$df) %>% unlist
  
  saveRDS(auto_df, file.path(output_dir, paste("mpm_", i, "_auto.RDS", sep = "")))



  ### Covariance all vr
  clim <- lapply(as.list(c(1:900)), function(x)
    rnorm_multi(n = 5000,
                vars = 2,
                mu = c(0,0),
                sd = clim_sd[x],
                r = clim_corr[x],
                varnames = c("surv", "growth")) )

  clusterExport(cl, c("clim"))

  cov <- pblapply(cl = cl,
                  clim,
                  function(x) st.lamb(env_surv = x$surv,
                                      env_growth = x$growth,
                                      env_reproduction = rep(0,5000),
                                      clim_sd = sd(x$surv),
                                      sig.strength = i)
  )
  cov_df <- lapply(cov, function(x) x$df) %>% unlist
  
  saveRDS(cov_df, paste("/mpm_", i, "_cov.RDS", sep = ""))

  #### Lagged effect within P "functions"
  
  lag_clim <- lapply(as.list(c(1:900)), function(x) create_seq(5000, clim_sd = clim_sd[x], clim_corr = clim_corr[x], lag = 1))
  
  clusterExport(cl, c("lag_clim"))
  
  lag_g <- pblapply(cl = cl,
                    lag_clim,
                    function(x) st.lamb(env_surv = x$recent,
                                        env_growth = x$lagged,
                                        env_reproduction = rep(0,length(x$recent)),
                                        clim_sd = sd(x$recent, na.rm = T),
                                        sig.strength = i)
  )

  lag_s <- pblapply(cl = cl,
                    lag_clim,
                    function(x) st.lamb(env_surv = x$lagged,
                                        env_growth = x$recent,
                                        env_reproduction = rep(0,length(x$recent)),
                                        clim_sd = sd(x$recent, na.rm = T),
                                        sig.strength = i)
  )

  lag_n <- pblapply(cl = cl,
                    lag_clim,
                    function(x) st.lamb(env_surv = x$recent,
                                        env_growth = x$recent,
                                        env_reproduction = rep(0,length(x$recent)),
                                        clim_sd = sd(x$recent, na.rm = T),
                                        sig.strength = i)
  )
  lagg_df <- lapply(lag_g, function(x) x$df) %>% unlist
  lags_df <- lapply(lag_s, function(x) x$df) %>% unlist
  lagn_df <- lapply(lag_n, function(x) x$df) %>% unlist
  
  lag <- list("growth" = lagg_df, "survival" = lags_df, "none" = lagn_df)

  saveRDS(lag, file.path(output_dir, paste("mpm_", i, "_lag.RDS", sep = "")))

  #### Lagged effect
  
  lag_p <- pblapply(cl = cl, 
                    lag_clim, 
                    function(x) st.lamb(env_surv = x$lagged,
                                        env_growth = x$lagged,
                                        env_reproduction = x$recent,
                                        clim_sd = sd(x$recent, na.rm = T),
                                        sig.strength = i)
  )
  
  lag_f <- pblapply(cl = cl, 
                    lag_clim, 
                    function(x) st.lamb(env_surv = x$recent,
                                        env_growth = x$recent,
                                        env_reproduction = x$lagged,
                                        clim_sd = sd(x$recent, na.rm = T),
                                        sig.strength = i)
  )
  
  lag_n2 <- pblapply(cl = cl, 
                    lag_clim, 
                    function(x) st.lamb(env_surv = x$recent,
                                        env_growth = x$recent,
                                        env_reproduction = x$recent,
                                        clim_sd = sd(x$recent, na.rm = T),
                                        sig.strength = i)
  )
  lagp_df <- lapply(lag_p, function(x) x$df) %>% unlist
  lagf_df <- lapply(lag_f, function(x) x$df) %>% unlist
  lagn2_df <- lapply(lag_n2, function(x) x$df) %>% unlist
  
  lag_fp <- list("Pkernel" = lagp_df, "Fkernel" = lagf_df, "none" = lagn2_df)
  
  saveRDS(lag_fp, file.path(output_dir, paste("mpm_", i, "_lagfp.RDS", sep = "")))
  
  stopCluster(cl)

}
