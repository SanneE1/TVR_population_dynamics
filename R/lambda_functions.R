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
  
  # growth 
  g_mean = mpm_df$gamma
  g_sd = sqrt(mpm_df$gamma_var)
  
  if(is.na(growth)) {
    gamma <- g_mean  
  } else {
    ## total deviation from mean = climate signal * signal strength & correction factor (partitioning at variance scale) + random noise * signal strength
    dev <- growth * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + 
      rnorm(1,0, clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    ## Because of partitioning and correction factor above, the resulting distribution has a sd of clim_sd
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    gamma <- qbeta(p, (((g_mean*(1-g_mean))/(g_sd * clim_sd)^2) - 1) * g_mean,
                   (((g_mean*(1-g_mean))/(g_sd * clim_sd)^2) - 1) * (1 - g_mean))
  }
  
  # survival juveniles
  sj_mean = mpm_df$Sj
  sj_sd = sqrt(mpm_df$Sj_var)
  
  if(is.na(survival)) {
    Sj <- sj_mean
  } else {
    dev <- survival * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1,0, clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    Sj <- qbeta(p, (((sj_mean*(1-sj_mean))/(sj_sd * clim_sd)^2) - 1) * sj_mean,
                (((sj_mean*(1-sj_mean))/(sj_sd * clim_sd)^2) - 1) * (1 - sj_mean))
  }
  
  # survival adults
  sa_mean = mpm_df$Sa
  sa_sd = sqrt(mpm_df$Sa_var)
  
  if(is.na(survival)) {
    Sa <- sa_mean
  } else {
    dev <- survival * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1,0, clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    Sa <- qbeta(p, (((sa_mean*(1-sa_mean))/(sa_sd * clim_sd)^2) - 1) * sa_mean,
                (((sa_mean*(1-sa_mean))/(sa_sd * clim_sd)^2) - 1) * (1 - sa_mean))
  }
  
  # reproduction 
  phi_mean = mpm_df$phi
  phi_sd = sqrt(mpm_df$phi_var)
  
  if(is.na(reproduction)) {
    phi <- phi_mean
  } else {
    dev <- reproduction * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1,0, clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    phi <- qgamma(p, (phi_mean^2)/(phi_sd * clim_sd)^2, (phi_mean)/(phi_sd * clim_sd)^2)
  }
  
  # creat matrix from vr
  mat <- matrix(0,2,2)
  mat[1,1] <- Sj*(1-gamma)
  mat[2,1] <- Sj*gamma
  mat[1,2] <- phi
  mat[2,2] <- Sa
  
  
  return(mat)  
}

## does the same as the mpm function above, but simulates opposing responses to climate in surv/growth and fecundity
mpm_o <- function(mpm_df, survival, growth, reproduction, 
                clim_sd, sig.strength) {
  
  # growth 
  g_mean = mpm_df$gamma
  g_sd = sqrt(mpm_df$gamma_var)
  
  if(is.na(growth)) {
    gamma <- g_mean  
  } else {
    ## total deviation from mean = climate signal * signal strength & correction factor (partitioning at variance scale) + random noise * signal strength
    dev <- growth * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + 
      rnorm(1,0, clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    ## Because of partitioning and correction factor above, the resulting distribution has a sd of clim_sd
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    gamma <- qbeta(p, (((g_mean*(1-g_mean))/(g_sd * clim_sd)^2) - 1) * g_mean,
                   (((g_mean*(1-g_mean))/(g_sd * clim_sd)^2) - 1) * (1 - g_mean))
  }
  
  # survival juveniles
  sj_mean = mpm_df$Sj
  sj_sd = sqrt(mpm_df$Sj_var)
  
  if(is.na(survival)) {
    Sj <- sj_mean
  } else {
    dev <- survival * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1,0, clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    Sj <- qbeta(p, (((sj_mean*(1-sj_mean))/(sj_sd * clim_sd)^2) - 1) * sj_mean,
                (((sj_mean*(1-sj_mean))/(sj_sd * clim_sd)^2) - 1) * (1 - sj_mean))
  }
  
  # survival adults
  sa_mean = mpm_df$Sa
  sa_sd = sqrt(mpm_df$Sa_var)
  
  if(is.na(survival)) {
    Sa <- sa_mean
  } else {
    dev <- survival * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1,0, clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    Sa <- qbeta(p, (((sa_mean*(1-sa_mean))/(sa_sd * clim_sd)^2) - 1) * sa_mean,
                (((sa_mean*(1-sa_mean))/(sa_sd * clim_sd)^2) - 1) * (1 - sa_mean))
  }
  
  # reproduction 
  phi_mean = mpm_df$phi
  phi_sd = sqrt(mpm_df$phi_var)
  
  if(is.na(reproduction)) {
    phi <- phi_mean
  } else {
    dev <- (-reproduction) * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1,0, clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    phi <- qgamma(p, (phi_mean^2)/(phi_sd * clim_sd)^2, (phi_mean)/(phi_sd * clim_sd)^2)
  }
  
  # creat matrix from vr
  mat <- matrix(0,2,2)
  mat[1,1] <- Sj*(1-gamma)
  mat[2,1] <- Sj*gamma
  mat[1,2] <- phi
  mat[2,2] <- Sa
  
  
  return(mat)  
}

# Function to create a MPM sequence from supplied climate and noise sequences and calculate stocastic lambda -----------------------------------

## This function will create a MPM for each iteration of the environmental sequence, and calculate the stochastic lambda of said MPM sequence 
st.lamb <- function(mpm_df, env_surv, env_growth, env_reproduction, 
                    clim_sd, clim_auto, sig.strength) {
  
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
  
  df <- mpm_df %>% tibble::add_column(lambda = popbio::stoch.growth.rate(mats, maxt = n_it, verbose = F)$sim,
                                      clim_sd = clim_sd,
                                      clim_auto = clim_auto)
  
  return(df = df) 
  # if you want to also return the matrices:
  # return(list(df = df,
  #        mats = sapply(mats, as.vector) %>% t %>% 
  # `colnames<-`(c("1,1", "2,1", "1,2", "2,2"))))
}

## does the same as the function above, only uses the mpm_o function rather than the mpm function
st.lamb_o <- function(mpm_df, env_surv, env_growth, env_reproduction, 
                    clim_sd, clim_auto, sig.strength) {
  
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
                                      clim_auto = clim_auto)
  
  return(df = df) 
  # if you want to also return the matrices:
  # return(list(df = df,
  #        mats = sapply(mats, as.vector) %>% t %>% 
  # `colnames<-`(c("1,1", "2,1", "1,2", "2,2"))))
}
