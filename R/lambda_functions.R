# Functions to simulate stochastic lambda for MPMs responding to climate (using the method of moments)

# Create function that creates environmental sequence ------------------------------------------------------------------------------------

## Creates a sequence of climate anomalies (c in manuscript) with a specified standard deviation and autocorrelation.
## The function then creates another sequence of the same length and standard deviation, to be used as random noise.
## Next the two vectors are combined in two dataframes, one which is the original ("recent") dataframe and a second, 
## which is offset by a specified period to create a "lagged" sequence

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

mpm <- function(mpm_df, survival, growth, reproduction, 
                clim_sd, sig.strength) {
  
  #correction factor (partitioning at variance scale) for climate (thetaP) and random noise (thetaP1 (theta[p-1]))
  thetaP <- sqrt(clim_sd^2 * sig.strength)/clim_sd
  thetaP1 <- sqrt(clim_sd^2 * (1-sig.strength))/clim_sd
  
  # survival juveniles
  sj_mean = mpm_df$Sj
  sj_sd = mpm_df$Sj_sd
  
  if(is.na(survival)) {
    Sj <- sj_mean
  } else {
    dev <- survival * thetaP + rnorm(1,0, clim_sd) * thetaP1 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    Sj <- qbeta(p, (((sj_mean*(1-sj_mean))/(sj_sd * clim_sd)^2) - 1) * sj_mean,
                (((sj_mean*(1-sj_mean))/(sj_sd * clim_sd)^2) - 1) * (1 - sj_mean))
  }
  
  # survival adults
  sa_mean = mpm_df$Sa
  sa_sd = mpm_df$Sa_sd
  
  if(is.na(survival)) {
    Sa <- sa_mean
  } else {
    dev <- survival * thetaP + rnorm(1,0, clim_sd) * thetaP1 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    Sa <- qbeta(p, (((sa_mean*(1-sa_mean))/(sa_sd * clim_sd)^2) - 1) * sa_mean,
                (((sa_mean*(1-sa_mean))/(sa_sd * clim_sd)^2) - 1) * (1 - sa_mean))
  }
  
  # reproduction 
  rho_mean = mpm_df$rho
  rho_sd = mpm_df$rho_sd
  
  if(is.na(reproduction)) {
    rho <- rho_mean
  } else {
    dev <- reproduction * thetaP + rnorm(1,0, clim_sd) * thetaP1 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    rho <- qgamma(p, (rho_mean^2)/(rho_sd * clim_sd)^2, (rho_mean)/(rho_sd * clim_sd)^2)
  }
  
  # creat matrix from vr
  mat <- matrix(0,2,2)
  mat[1,1] <- Sj*(1-mpm_df$gamma)
  mat[2,1] <- Sj*mpm_df$gamma
  mat[1,2] <- rho
  mat[2,2] <- Sa
  
  
  return(mat)  
}

## does the same as the mpm function above, but simulates opposing responses to climate in surv/growth and fecundity
mpm_o <- function(mpm_df, survival, growth, reproduction, 
                  clim_sd, sig.strength) {
  
  thetaP <- sqrt(clim_sd^2 * sig.strength)/clim_sd
  thetaP1 <- sqrt(clim_sd^2 * (1-sig.strength))/clim_sd
  
  # survival juveniles
  sj_mean = mpm_df$Sj
  sj_sd = mpm_df$Sj_sd
  
  if(is.na(survival)) {
    Sj <- sj_mean
  } else {
    dev <- survival * thetaP + rnorm(1,0, clim_sd) * thetaP1 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    Sj <- qbeta(p, (((sj_mean*(1-sj_mean))/(sj_sd * clim_sd)^2) - 1) * sj_mean,
                (((sj_mean*(1-sj_mean))/(sj_sd * clim_sd)^2) - 1) * (1 - sj_mean))
  }
  
  # survival adults
  sa_mean = mpm_df$Sa
  sa_sd = mpm_df$Sa_sd
  
  if(is.na(survival)) {
    Sa <- sa_mean
  } else {
    dev <- survival * thetaP + rnorm(1,0, clim_sd) * thetaP1 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    Sa <- qbeta(p, (((sa_mean*(1-sa_mean))/(sa_sd * clim_sd)^2) - 1) * sa_mean,
                (((sa_mean*(1-sa_mean))/(sa_sd * clim_sd)^2) - 1) * (1 - sa_mean))
  }
  
  # reproduction 
  rho_mean = mpm_df$rho
  rho_sd = mpm_df$rho_sd
  
  if(is.na(reproduction)) {
    rho <- rho_mean
  } else {
    dev <- (-reproduction) * thetaP + rnorm(1,0, clim_sd) * thetaP1 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    rho <- qgamma(p, (rho_mean^2)/(rho_sd * clim_sd)^2, (rho_mean)/(rho_sd * clim_sd)^2)
  }
  
  # creat matrix from vr
  mat <- matrix(0,2,2)
  mat[1,1] <- Sj*(1-mpm_df$gamma)
  mat[2,1] <- Sj*mpm_df$gamma
  mat[1,2] <- rho
  mat[2,2] <- Sa
  
  
  return(mat)  
}

## does the same as the other two, but for one vital rate, the signal strength is halved (decreasing the sensitivity of a vital rate to the climate)
mpm_weak <- function(mpm_df, survival, growth, reproduction, 
                     clim_sd, sig.strength, weak_vr) {
  
  #correction factor (partitioning at variance scale) for climate (thetaP) and random noise (thetaP1 (theta[p-1]))
  thetaP <- sqrt(clim_sd^2 * sig.strength)/clim_sd
  thetaP1 <- sqrt(clim_sd^2 * (1-sig.strength))/clim_sd
  
  weak_thetaP <- sqrt(clim_sd^2 * (0.5 * sig.strength))/clim_sd
  weak_thetaP1 <- sqrt(clim_sd^2 * (1-(0.5 * sig.strength)))/clim_sd
  
  # survival juveniles
  sj_mean = mpm_df$Sj
  sj_sd = mpm_df$Sj_sd
  
  if(is.na(survival)) {
    Sj <- sj_mean
  } else {
    if(weak_vr == "Sj"){
      dev <- survival * weak_thetaP + rnorm(1,0, clim_sd) * weak_thetaP1
    } else {
      dev <- survival * thetaP + rnorm(1,0, clim_sd) * thetaP1  
    } 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    Sj <- qbeta(p, (((sj_mean*(1-sj_mean))/(sj_sd * clim_sd)^2) - 1) * sj_mean,
                (((sj_mean*(1-sj_mean))/(sj_sd * clim_sd)^2) - 1) * (1 - sj_mean))
  }
  
  # survival adults
  sa_mean = mpm_df$Sa
  sa_sd = mpm_df$Sa_sd
  
  if(is.na(survival)) {
    Sa <- sa_mean
  } else {
    if(weak_vr == "Sa"){
      dev <- survival * weak_thetaP + rnorm(1,0, clim_sd) * weak_thetaP1
    } else {
      dev <- survival * thetaP + rnorm(1,0, clim_sd) * thetaP1  
    } 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    Sa <- qbeta(p, (((sa_mean*(1-sa_mean))/(sa_sd * clim_sd)^2) - 1) * sa_mean,
                (((sa_mean*(1-sa_mean))/(sa_sd * clim_sd)^2) - 1) * (1 - sa_mean))
  }
  
  # reproduction 
  rho_mean = mpm_df$rho
  rho_sd = mpm_df$rho_sd
  
  if(is.na(reproduction)) {
    rho <- rho_mean
  } else {
    if(weak_vr == "Sj"){
      dev <- reproduction * weak_thetaP + rnorm(1,0, clim_sd) * weak_thetaP1
    } else {
      dev <- reproduction * thetaP + rnorm(1,0, clim_sd) * thetaP1  
    } 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    rho <- qgamma(p, (rho_mean^2)/(rho_sd * clim_sd)^2, (rho_mean)/(rho_sd * clim_sd)^2)
  }
  
  # creat matrix from vr
  mat <- matrix(0,2,2)
  mat[1,1] <- Sj*(1-mpm_df$gamma)
  mat[2,1] <- Sj*mpm_df$gamma
  mat[1,2] <- rho
  mat[2,2] <- Sa
  
  
  return(mat)  
}

# Function to create a MPM sequence from supplied climate and noise sequences and calculate stocastic lambda -----------------------------------

## This function will create a MPM for each iteration of the environmental sequence, and calculate the stochastic lambda of said MPM sequence 
st.lamb <- function(mpm_df, env_surv, env_growth, env_reproduction, 
                    clim_sd, clim_auto, sig.strength, return.mpm = F) {
  
  n_it = length(env_surv)
  
  env <- list(
    mpm_df = list(mpm_df),
    survival = env_surv,
    growth = env_growth,
    reproduction = env_reproduction,
    clim_sd = clim_sd,
    sig.strength = sig.strength)
  
  ### Get all mpm's
  mats <- purrr::pmap(env, mpm) %>% Filter(Negate(anyNA), .)
  mats <- mats %>% purrr::discard(function(x) all(x == 0))
  
  lam_sim <- popbio::stoch.growth.rate(mats, maxt = length(mats), verbose = F)$sim
  
  df <- mpm_df %>% tibble::add_column(lambda = ifelse(!is.nan(lam_sim), lam_sim,
                                                      popbio::stoch.growth.rate(mats, maxt = length(mats), verbose = F)$approx),
                                      clim_sd = clim_sd,
                                      clim_auto = clim_auto,
                                      n_mats = nrow(mats),
                                      lam_calc = ifelse(!is.nan(lam_sim), "sim", "approx"))   ## Return for check that simulations produce right number of mats (i.e. not just mats with NaNs)
  
  if(return.mpm == F){ 
    return(df = df) 
  } else {
    return(list(df = df,
                mats = sapply(mats, as.vector) %>% t %>%
                  `colnames<-`(c("1,1", "2,1", "1,2", "2,2")))) 
  }
  
}

## does the same as the function above, only uses the mpm_o function rather than the mpm function
st.lamb_o <- function(mpm_df, env_surv, env_growth, env_reproduction, 
                      clim_sd, clim_auto, sig.strength, return.mpm = F) {
  
  n_it = length(env_surv)
  
  env <- list(
    mpm_df = list(mpm_df),
    survival = env_surv,
    growth = env_growth,
    reproduction = env_reproduction,
    clim_sd = clim_sd,
    sig.strength = sig.strength)
  
  ### Get all mpm's
  mats <- purrr::pmap(env, mpm_o) %>% Filter(Negate(anyNA), .)
  mats <- mats %>% purrr::discard(function(x) all(x == 0))
  
  df <- mpm_df %>% tibble::add_column(lambda = popbio::stoch.growth.rate(mats, maxt = n_it, verbose = F)$sim,
                                      clim_sd = clim_sd,
                                      clim_auto = clim_auto,
                                      n_mats = nrow(mats))
  
  if(return.mpm == F){ 
    return(df = df) 
  } else {
    return(list(df = df,
                mats = sapply(mats, as.vector) %>% t %>%
                  `colnames<-`(c("1,1", "2,1", "1,2", "2,2")))) 
  }
  
}

## for the mpm_weak function
st.lamb_weak <- function(mpm_df, env_surv, env_growth, env_reproduction, 
                         clim_sd, clim_auto, sig.strength, weak_vr, return.mpm = F) {
  
  n_it = length(env_surv)
  
  env <- list(
    mpm_df = list(mpm_df),
    survival = env_surv,
    growth = env_growth,
    reproduction = env_reproduction,
    clim_sd = clim_sd,
    sig.strength = sig.strength,
    weak_vr = weak_vr)
  
  ### Get all mpm's
  mats <- purrr::pmap(env, mpm_weak) %>% Filter(Negate(anyNA), .)
  mats <- mats %>% purrr::discard(function(x) all(x == 0))
  
  lam_sim <- popbio::stoch.growth.rate(mats, maxt = length(mats), verbose = F)$sim
  
  df <- mpm_df %>% tibble::add_column(lambda = ifelse(!is.nan(lam_sim), lam_sim,
                                                      popbio::stoch.growth.rate(mats, maxt = length(mats), verbose = F)$approx),
                                      clim_sd = clim_sd,
                                      clim_auto = clim_auto,
                                      n_mats = nrow(mats),
                                      lam_calc = ifelse(!is.nan(lam_sim), "sim", "approx"))   ## Return for check that simulations produce right number of mats (i.e. not just mats with NaNs)
  
  if(return.mpm == F){ 
    return(df = df) 
  } else {
    return(list(df = df,
                mats = sapply(mats, as.vector) %>% t %>%
                  `colnames<-`(c("1,1", "2,1", "1,2", "2,2")))) 
  }
  
}

