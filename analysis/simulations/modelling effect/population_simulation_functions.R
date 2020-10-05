### Functions for IBM simulations


## Vital rate functions

# survival function - logistic
s_z <- function(z, par, clim, yr) {
  
  x <- par$s_int + par$s_size * log(z) + par$s_temp * clim$recent[which(clim$yr == yr)]
  p <- 1/(1+exp(-(x)))
  
  return(p)
}

# growth function - linear

g_z1z <- function(z1, z, par, clim, yr) {
  
  mean <- exp(par$g_int + par$g_size * log(z) + par$g_temp * clim$lagged[which(clim$yr == yr)])  ## mean size next year
  p_grow <- dnorm(z1, mean = mean, sd = par$g_sd)   ### probability distribution that you grow to size z1 
  ### based on mean size and sd
  return(p_grow)
}

# Flower probability - logistic
fp_z <- function(z, par) {
  
  x <- par$fp_int + par$fp_size * log(z)
  p <- 1/(1+exp(-(x)))
  
  return(p)
}

# Seed production - exp
seed_z <- function(z, par) {
  
  n <- exp(par$seed_int + par$seed_size * log(z))
  
  return(n)
}

# recruit size
fd_z1 <- function(z1, par) {
  
  fd <- dnorm(z1, mean = par$fd_int, sd = par$fd_sd)    ## probability of size z1 recruits
  
  return(fd)
}



simulate <- function(init.pop.size, n_yrs, clim_sd, clim_corr) {
  
  clim <- data.frame(yr = c(-4:n_yrs),
                     recent = create_seq(n_it = n_yrs, clim_sd = clim_sd, clim_corr = clim_corr))
  a <- clim %>% mutate(yr = yr + 1) %>% rename(lagged = recent)
  clim <- left_join(clim, a)
  print(clim_sd)
  print(clim_corr)
  
  ## True parameters
  par.true <- data.frame(s_int = -0.229,
                         s_size = 1.077,
                         s_temp = 1.233,
                         g_int = 0.424,
                         g_size = 0.846,
                         g_temp = -0.066,
                         g_sd = 1.076,
                         fp_int = -3.970,
                         fp_size = 1.719,
                         fd_int = 1.178749,
                         fd_sd = 0.76,
                         seed_int = -0.6652477,
                         seed_size = 0.8747764,
                         germ_int = 0.1262968,
                         germ_sd = 0.2725941)  ### see if I can work germ_sd into the simulation
  
  
  ## Simulate a population
  
  # set initial population sizes
  z <- round(abs(rnorm(init.pop.size, mean = par.true$fd_int, sd = par.true$fd_sd)), digits = 0)
  
  # vectors to store pop size and mean size
  pop.size.t <- mean.z.t <- mean.z.repr.t <- n_recr.t <- n_seeds.t <- numeric(n_yrs)
  
  
  
  yr <- 1
  
  try(
  while (yr <= n_yrs && length(z) > 0) {
    
    # calculate current population size
    pop.size <- length(z)
    
    # calculate/generate survival chance, based on size z for each of the "individuals"
    # in population vector z
    surv <- rbinom(n = pop.size, prob = s_z(z, par.true, clim, yr), size = 1)
    
    # generate size of the surviving individuals in t+1
    df <- data.frame(z = z, surv = surv) %>% 
      mutate(z1 = ifelse(surv == 1, 
                         round(rnorm(n = pop.size, 
                                     mean = exp(par.true$g_int + par.true$g_size * log(z) + 
                                                  par.true$g_temp * clim$lagged[which(clim$yr == yr)]), 
                                     sd = par.true$g_sd), digits = 0), 
                         NA)) %>%
      mutate(z1 = replace(z1,
                          z1 <= 0,
                          1))
    
    # generate probability of flowering for surviving individuals
    df <- df %>%
      mutate(fp = ifelse(surv == 1, 
                         rbinom(n = pop.size,
                                prob = fp_z(z, par.true),
                                size = 1), 
                         NA))
    
    ## Calculate number of seeds produced
    df <- df %>%
      mutate(n_seeds = ifelse(fp == 1,
                              rpois(n = pop.size,
                                    seed_z(z, par.true)),
                              NA))
    
    ## Calculate number of recruits
    recr <- ifelse(sum(df$n_seeds, na.rm = T) == 0, 0, rbinom(1,
                                                              sum(df$n_seeds, na.rm = T),
                                                              par.true$germ_int))
    
    ## Calculate size distribution of recruits
    recz <- round(rnorm(recr, mean = par.true$fd_int, sd = par.true$fd_sd), digits = 0)
    recz[which(recz <= 0)] <- 1
    
    # store simulated data
    if (yr == 1) {
      sim.data <- data.frame(df, yr = 1, 
                             recent = clim$recent[clim$yr == 1], 
                             lagged = clim$lagged[clim$yr == 1])
    } else {
      sim.data <- rbind(sim.data, data.frame(df, yr = yr, 
                                             recent = clim$recent[which(clim$yr == yr)], 
                                             lagged = clim$lagged[which(clim$yr == yr)]))
    }
    
    z1 <- c(recz, df$z1[which(df$surv == 1)])
    
    ## store the population size and mean size
    pop.size.t[yr] <- length(df$z)
    mean.z.t[yr] <- mean(df$z)
    mean.z.repr.t[yr] <- mean(df$z[which(df$fp == 1)], na.rm = TRUE)
    n_seeds.t[yr] <- sum(df$n_seeds, na.rm = T)
    n_recr.t[yr] <- recr
    
    z <- z1
    yr <- yr + 1
    
    
  })
  
  model_lagged <- summary(glm(z1 ~ log(z) + lagged, family = "poisson", data = sim.data))  ## Right model
  model_recent <- summary(glm(z1 ~ log(z) + recent, family = "poisson", data = sim.data))  ## Wrong model
  
  # model_lagged$coefficients["lagged", c("Estimate", "Pr(>|z|)")]
  # model_recent$coefficients["recent", c("Estimate", "Pr(>|z|)")]
  # 
  df <- tibble(clim_corr = clim_corr,
               clim_sd = clim_sd,
               lagged_estimate = model_lagged$coefficients["lagged", c("Estimate")],
               lagged_p = model_lagged$coefficients["lagged", c("Pr(>|z|)")],
               recent_estimate = model_recent$coefficients["recent", c("Estimate")],
               recent_p = model_recent$coefficients["recent", c("Pr(>|z|)")],
               pop.sizes = list(pop.size.t))
  
  return(df)
}
