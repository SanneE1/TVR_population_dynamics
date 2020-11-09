## Wrapper to simulate data, calculate vital rates, and use these vital rates for lambda calculations

wrapper <- function(sample, clim_sd, clim_corr) {
  
  data <- simulate(sample, clim_sd, clim_corr)

  vr_lagged <- bayes_vr_lagged(data)
  vr_recent <- bayes_vr_recent(data)
  
  lagged_lambdas <- bayes_lambda(vr_lagged, clim_sd, clim_corr, n_it = 5000)
  recent_lambdas <- bayes_lambda(vr_recent, clim_sd, clim_corr, n_it = 5000)
  
  tibble(data$df, list(lagged_lambdas), list(recent_lambdas))
}


### Functions for IBM simulations

### Create climate sequence
create_seq <- function(n_it, clim_sd, clim_corr) { 
  for(n in c(1:(n_it+5))) 
    if(n == 1) {
      seq <- rnorm(1)
    } else {
      seq[n] <- clim_corr * seq[n-1] + rnorm(1)
    }
  seq <- scale(seq) * clim_sd
  return(seq)
}


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


## simulate dataset
simulate <- function(sample, clim_sd, clim_corr, init.pop.size = 100000, run_yrs = 1000) {
  
  clim <- data.frame(yr = c(-4:run_yrs),
                     recent = create_seq(n_it = run_yrs, clim_sd = clim_sd, clim_corr = clim_corr))
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
  pop.size.t <- mean.z.t <- mean.z.repr.t <- n_recr.t <- n_seeds.t <- numeric(run_yrs)
  
  
  
  yr <- 1
  
  try(
  while (yr <= run_yrs && length(z) > 0) {
    
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
      z_new <- recz
    } else {
      sim.data <- rbind(sim.data, data.frame(df, yr = yr, 
                                             recent = clim$recent[which(clim$yr == yr)], 
                                             lagged = clim$lagged[which(clim$yr == yr)]))
      z_new <- c(z_new, recz)
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
  
  full_data <- sim.data
  
  if( max(sim.data$yr) - round(max(sim.data$yr)*0.75) >= sample ) {
    start_sample <- round(max(sim.data$yr)*0.75)
    sample_range <- c(start_sample:(start_sample+sample))
    
    sim.data <- sim.data %>% filter(yr %in% sample_range)

  } else if (max(sim.data$yr) - round(max(sim.data$yr)*0.5) >= sample ) {
    start_sample <- round(max(sim.data$yr)*0.5)
    sample_range <- c(start_sample:(start_sample+sample))
    
    sim.data <- sim.data %>% filter(yr %in% sample_range)
  } else {
    start_sample <- round(max(sim.data$yr)-sample)
    sample_range <- c(start_sample:(start_sample+sample))
      
    sim.data <- sim.data %>% filter(yr %in% sample_range)
  }
  
  model_lagged <- summary(glm(z1 ~ log(z) + lagged, family = "poisson", data = sim.data))  ## Right model
  model_recent <- summary(glm(z1 ~ log(z) + recent, family = "poisson", data = sim.data))  ## Wrong model
  
  actual_lambda <- data.frame(N1 = c(pop.size.t, NA), N = c(NA, pop.size.t)) %>% 
                            mutate(lambda = N1/N) %>%
                            summarise(lambda = log(mean(lambda, na.rm = T)))
  
  df <- tibble(clim_corr = clim_corr,
               clim_sd = clim_sd,
               actual_lambda = actual_lambda$lambda,
               lagged_estimate = model_lagged$coefficients["lagged", c("Estimate")],
               lagged_p = model_lagged$coefficients["lagged", c("Pr(>|z|)")],
               recent_estimate = model_recent$coefficients["recent", c("Estimate")],
               recent_p = model_recent$coefficients["recent", c("Pr(>|z|)")],
               pop.sizes = list(pop.size.t[sample_range]),
               n_seeds.t = list(n_seeds.t[sample_range]),
               n_recr.t = list(n_recr.t[sample_range]),
               z_new = list(z_new[sample_range])
               )
  
  return(list(df = df, sim.data = sim.data, full.data = full_data))
  
}

## calculate vr models from simulated dataset
#using lagged climate for growth
bayes_vr_lagged <- function(data, n_sample = 100) {
  
  sim.data <- data$sim.data
  sim.data <- sim.data %>% filter(z != 0)
  
  data_mod <- list(
    yr_1 = sim.data$yr,
    z_1 = sim.data$z,
    surv = sim.data$surv,
    recent_1 = sim.data$recent,
    lagged_1 = sim.data$lagged,
    yr_2 = sim.data$yr[which(sim.data$surv == 1)],
    z_2 = sim.data$z[which(sim.data$surv == 1)],
    z1 = sim.data$z1[which(sim.data$surv == 1)],
    fp = sim.data$fp[which(sim.data$surv == 1)],
    recent_2 = sim.data$recent[which(sim.data$surv == 1)],
    lagged_2 = sim.data$lagged[which(sim.data$surv == 1)],
    yr_3 = sim.data$yr[which(sim.data$fp == 1)],
    z_3 = sim.data$z[which(sim.data$fp == 1)],
    n_seeds = sim.data$n_seeds[which(sim.data$fp == 1)],
    n_recr = data$df$n_recr.t[[1]], 
    n_seeds_t = data$df$n_seeds.t[[1]],
    z_new = data$df$z_new[[1]],
    N3 = length(sim.data$yr[which(sim.data$fp == 1)]),
    N2 = length(sim.data$yr[which(sim.data$surv == 1)]),
    N1 = length(sim.data$yr),
    N = length(data$df$z_new[[1]])
  )
  
  recent_mod <- cmdstan_model("/home/evers/lagged_buffering/analysis/simulations/modelling_effect/vr_model_recent.stan")
  recent_fit <- recent_mod$sample(data = data_mod,
                                  chains = 1,
                                  iter_sampling = 2500,
                                  iter_warmup = 2500)
  
  
  params <- posterior::as_draws_df(recent_fit$draws())
  params <- params[sample(2500, n_sample), ]
  
  return(params)
  
}


# using recent climate for growth
bayes_vr_recent <- function(data, n_sample = 100) {
  
  sim.data <- data$sim.data
  sim.data <- sim.data %>% filter(z != 0)
  
  
  data_mod <- list(
    yr_1 = sim.data$yr,
    z_1 = sim.data$z,
    surv = sim.data$surv,
    recent_1 = sim.data$recent,
    lagged_1 = sim.data$lagged,
    yr_2 = sim.data$yr[which(sim.data$surv == 1)],
    z_2 = sim.data$z[which(sim.data$surv == 1)],
    z1 = sim.data$z1[which(sim.data$surv == 1)],
    fp = sim.data$fp[which(sim.data$surv == 1)],
    recent_2 = sim.data$recent[which(sim.data$surv == 1)],
    lagged_2 = sim.data$lagged[which(sim.data$surv == 1)],
    yr_3 = sim.data$yr[which(sim.data$fp == 1)],
    z_3 = sim.data$z[which(sim.data$fp == 1)],
    n_seeds = sim.data$n_seeds[which(sim.data$fp == 1)],
    n_recr = data$df$n_recr.t[[1]], 
    n_seeds_t = data$df$n_seeds.t[[1]],
    z_new = data$df$z_new[[1]]
  )
  
  lagged_mod <- cmdstan_model("/home/evers/lagged_buffering/analysis/simulations/modelling_effect/vr_model_lagged.stan")
  lagged_fit <- lagged_mod$sample(data = data_mod,
                                  chains = 1,
                                  iter_sampling = 2500,
                                  iter_warmup = 2500)
  
  
  params <- posterior::as_draws_df(lagged_fit$draws())
  params <- params[sample(2500, n_sample), ]
  
  return(params)
  
}


# calculate lambda from calculated vitalrate models
bayes_lambda <- function(param.sample, clim_sd, clim_corr, type, n_mesh = 100, n_it = 10000) {
  
  init_pop_vec <- runif(n_mesh)
  environ_seq <- create_seq(n_it = n_it, clim_sd = clim_sd, clim_corr = clim_corr)
  
  ## Define environment -------------------------------------------------------------------------
  
  env_sampler <- function(environ_seq, iteration) {
    
    temp <- list("temp0" = environ_seq[iteration + 1],
                 "temp1" = environ_seq[iteration]
    )
    
    return(temp)
  }
  
  
  ## create custom functions -------------------------------------------------------------------------
  
  inv_logit <- function(x) {
    return(
      1/(1 + exp(-(x)))
    )
  }
  
  pois <- function(x) {
    return(
      exp(x)
    )
  }
  
  my_functions <- list(inv_logit = inv_logit,
                       pois = pois,
                       env_sampler = env_sampler)
  
  
  lambda <- c(1:length(param.sample$lp__))
  
  for(i in c(1:length(param.sample$lp__))) {
    svMisc::progress(i, progress.bar = T)
    
    params_list <- list(s_int = param.sample$s_int[i],
                        s_size = param.sample$s_size[i],
                        s_temp = param.sample$s_recent[i],
                        g_int = param.sample$g_int[i],
                        g_size = param.sample$g_size[i],
                        g_recent = param.sample$g_recent[i],
                        g_lagged = param.sample$g_lagged[i],
                        g_sd = param.sample$g_sd[i],
                        fp_int = param.sample$fp_int[i],
                        fp_size = param.sample$fp_size[i],
                        fd_int = param.sample$fd_int[i],
                        fd_sd = param.sample$fd_sd[i],
                        seed_int = param.sample$seed_int[i],
                        seed_size = param.sample$seed_size[i],
                        germ_int = param.sample$germ_int[i],
                        germ_sd = param.sample$germ_sd[i])
    
    if(is.null(params_list$g_lagged) == T) { params_list$g_lagged <- 0}
    if(is.null(params_list$g_recent) == T) { params_list$g_recent <- 0}
    
    ipm <-  init_ipm("simple_di_stoch_param") %>%
      define_kernel(
        name = "P",
        
        formula = s * g,
        family = "CC",
        
        s = inv_logit(s_int + s_slope * log(size_1) + s_recent * temp0),
        g = dnorm(size_2, mean = g_mean, sd = g_sd),
        g_mean = pois(g_int + g_slope * log(size_1)  + g_recent * temp0 + g_lagged * temp1),
        
        data_list = params_list,
        states = list(c('size')),
        
        has_hier_effs = FALSE,
        
        evict_cor = TRUE,
        evict_fun = truncated_distributions("norm", "g")
      ) %>%
      define_kernel(
        name = "F",
        
        formula = fp * fn * germ * fd,
        family = "CC",
        
        fp = inv_logit(fp_int + fp_slope * log(size_1)),
        fn = pois(seed_int + seed_slope * log(size_1)),
        germ = inv_logit(p_germ),
        fd = dnorm(size_2, mean = pois(fd_int), sd = fd_sd),
        
        data_list = params_list,
        states = list(c("size")),
        
        has_hier_effs = FALSE,
        
        evict_cor = TRUE,
        evict_fun = truncated_distributions("norm", "fd")
      ) %>%
      define_k(
        name = "K",
        family = "IPM",
        K = P + F,
        n_size_t_1 = K %*% n_size_t,
        data_list = list(),
        states = list(c("size")),
        has_hier_effs = FALSE,
        
        evict_cor = FALSE
      ) %>% 
      define_impl(
        make_impl_args_list(
          kernel_names = c("K", "P", "F"),
          int_rule = rep("midpoint", 3),
          dom_start = rep("size", 3),
          dom_end = rep("size", 3)
        )
      ) %>%
      define_domains( size = c(1, 115, n_mesh)
      ) %>% 
      define_env_state(
        env_values = env_sampler(environ_seq = environ_seq,
                                 iteration = t),
        data_list = list(
          environ_seq = environ_seq,
          env_sampler = env_sampler
        )
      ) %>%
      define_pop_state(
        pop_vectors = list(
          n_size_t = init_pop_vec
        )
      ) %>%
      make_ipm(usr_funs = my_functions,
               iterate = TRUE,
               iterations = n_it)
    
    lambda[i] <- lambda(ipm, "pop_size", "stochastic")
  }
  
  return(lambda)
}

