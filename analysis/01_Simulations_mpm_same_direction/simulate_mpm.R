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

## Alright kids, this is where it gets complicated: partitioning on variance scale and distributions on sd scale 
mpm <- function(survival, growth, reproduction, clim_sd, sig.strength) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)
  
  #growth                    
  if(is.na(growth)) {
    mpm[2,1] <- inv_logit(0) ### inv_logit(0) = 0.5 (intercept) 
  } else {
    mpm[2,1] <- inv_logit((growth * (sqrt(clim_sd^2 * sig.strength)/clim_sd)) + 
                            ((sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) * rnorm(1, 0, clim_sd)) 
    ) 
  }
  # survival
  if(is.na(survival)) {
    mpm[2,2] <- inv_logit(0)
  } else {
    mpm[2,2] <- inv_logit((survival * (sqrt(clim_sd^2 * sig.strength)/clim_sd)) + 
                            ((sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) * rnorm(1, 0, clim_sd)) 
    )
  }
  # reproduction 
  if(is.na(reproduction)) {
    mpm[1,2] <- exp(1.2)
  } else {
    mpm[1,2] <- exp(1.2 + 
                      (reproduction * (sqrt(clim_sd^2 * sig.strength)/clim_sd)) + 
                      ((sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) * rnorm(1, 0, clim_sd)) 
    )
  }
  return(mpm)  
}

# "Stocastic" mpm -----------------------------------

st.lamb <- function(env_surv, env_growth, env_reproduction, clim_sd, clim_auto, sig.strength) {
  
  n_it = length(env_surv)
  
  env <- data.frame(survival = env_surv,
                    growth = env_growth,
                    reproduction = env_reproduction)
  
  
  env <- as.list(as.data.frame(t(env)))
  
  ### Get all mpm's
  mats <- lapply(env, function(x) mpm(x[1], x[2], x[3], clim_sd, sig.strength))
  
  
  df <- data.frame(lambda = stoch.growth.rate(mats, maxt = n_it)$sim,
                  clim_sd = clim_sd,
                  clim_auto = clim_auto)
  
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
  
  clim_sd <- rep(seq(from = 0.01, to = 2, length.out = 10), 90)
  clim_auto <- rep(rep(c(-0.9,0,0.9), each = 10), 30)
  
  ## export objects to workers
  clusterExport(cl, c("i", "create_seq", "inv_logit", "mpm", "st.lamb", "clim_auto", "clim_sd"))


  auto <- pblapply(cl = cl,
                   as.list(c(1:900)),
                   function(x) try(st.lamb(env_surv = rep(NA, 5000),
                                       env_growth = create_seq(n_it = 5000, clim_sd[x], clim_auto[x], 0)$recent,
                                       env_reproduction = rep(NA,5000),
                                       clim_sd = clim_sd[x],
                                       clim_auto = clim_auto[x],
                                       sig.strength = i))
  )


  saveRDS(auto, file.path(output_dir, paste("mpm_", i, "_auto.RDS", sep = "")))



  ### Covariance all U vr
  clim <- lapply(as.list(c(1:900)), function(x)
    rnorm_multi(n = 5000,
                vars = 2,
                mu = c(0,0),
                sd = clim_sd[x],
                r = clim_auto[x],
                varnames = c("surv", "growth")) )

  clusterExport(cl, c("clim"))

  cov <- pblapply(cl = cl,
                  clim,
                  function(x) st.lamb(env_surv = x$surv,
                                      env_growth = x$growth,
                                      env_reproduction = rep(NA,5000),
                                      clim_sd = sd(x$surv),
                                      clim_auto = cor(x$surv, x$growth),
                                      sig.strength = i)
  )

  saveRDS(cov, file.path(output_dir, paste("mpm_", i, "_cov.RDS", sep = "")))
  

  #### Lagged effect within U matrix

  lag_clim <- lapply(as.list(c(1:900)), function(x) create_seq(5000, clim_sd = clim_sd[x], clim_auto = clim_auto[x], lag = 1))

  clusterExport(cl, c("lag_clim"))

  lag_g <- pblapply(cl = cl,
                    lag_clim,
                    function(x) st.lamb(env_surv = x$recent,
                                        env_growth = x$lagged,
                                        env_reproduction = rep(NA,length(x$recent)),
                                        clim_sd = sd(x$recent, na.rm = T),
                                        clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                        sig.strength = i)
  )

  lag_s <- pblapply(cl = cl,
                    lag_clim,
                    function(x) st.lamb(env_surv = x$lagged,
                                        env_growth = x$recent,
                                        env_reproduction = rep(NA,length(x$recent)),
                                        clim_sd = sd(x$recent, na.rm = T),
                                        clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                        sig.strength = i)
  )

  lag_n <- pblapply(cl = cl,
                    lag_clim,
                    function(x) st.lamb(env_surv = x$recent,
                                        env_growth = x$recent,
                                        env_reproduction = rep(NA,length(x$recent)),
                                        clim_sd = sd(x$recent, na.rm = T),
                                        clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                        sig.strength = i)
  )

  lag <- list("growth" = lag_g, "survival" = lag_s, "none" = lag_n)

  saveRDS(lag, file.path(output_dir, paste("mpm_", i, "_lag.RDS", sep = "")))

  #### Lagged effect between U & F matrices

  lag_p <- pblapply(cl = cl,
                    lag_clim,
                    function(x) st.lamb(env_surv = x$lagged,
                                        env_growth = x$lagged,
                                        env_reproduction = x$recent,
                                        clim_sd = sd(x$recent, na.rm = T),
                                        clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                        sig.strength = i)
  )

  lag_f <- pblapply(cl = cl,
                    lag_clim,
                    function(x) st.lamb(env_surv = x$recent,
                                        env_growth = x$recent,
                                        env_reproduction = x$lagged,
                                        clim_sd = sd(x$recent, na.rm = T),
                                        clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                        sig.strength = i)
  )

  lag_n2 <- pblapply(cl = cl,
                    lag_clim,
                    function(x) st.lamb(env_surv = x$recent,
                                        env_growth = x$recent,
                                        env_reproduction = x$recent,
                                        clim_sd = sd(x$recent, na.rm = T),
                                        clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                        sig.strength = i)
  )

  lag_fp <- list("Umatrix" = lag_p, "Fmatrix" = lag_f, "none" = lag_n2)

  saveRDS(lag_fp, file.path(output_dir, paste("mpm_", i, "_lagfp.RDS", sep = "")))

  stopCluster(cl)

}
