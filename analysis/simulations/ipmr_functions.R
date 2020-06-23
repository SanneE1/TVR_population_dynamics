
create_seq <- function(n_it, clim_sd, clim_corr) { 
  for(n in c(1:(n_it+5))) 
    if(n == 1) {
      seq <- rnorm(1,0, clim_sd)
    } else {
      seq[n] <- clim_corr * seq[n-1] + rnorm(1,0,clim_sd)
    }
  return(seq)
}


P_1yr <- function(n_it, clim_sd, clim_corr) {
  
  init_pop_vec <- runif(200)
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
  
  
  ### set up non-lagged ipm -------------------------------------------------------------------------
 
  clim_mod <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp0),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp0),
      
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
      fn = pois(fn_int + fn_slope * log(size_1)),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
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
    define_domains( size = c(1, 115, 200)
    ) %>% 
    define_env_state(
      env_params = env_sampler(environ_seq = environ_seq,
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
  
  
  lambdas <- tibble(clim_sd = clim_sd,
                    autocorrelation = clim_corr,
                    ### get lambda non-lagged ---------------------------------------------------------
                    non_lagged = lambda(clim_mod, "pop_size", "stochastic"), 
                    non_lagged_all = list(lambda(clim_mod, "pop_size", "all")))
  
  # remove clim_mod object to save memory
  rm(clim_mod)
  
  ## lagged ipm ----------------------------------------
  
  ipm_lagged <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp0),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp1),
      
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
      fn = pois(fn_int + fn_slope * log(size_1)),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
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
    define_domains( size = c(1, 115, 200))  %>%
    define_env_state(
      env_params = env_sampler(environ_seq = environ_seq,
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
             iterate = T, 
             iterations = n_it)

  
  lambdas$lagged <- lambda(ipm_lagged, "pop_size", "stochastic")
  lambdas$lagged_all <- list(lambda(ipm_lagged, "pop_size", "all"))
  
  
  return(lambdas)
}

P_neg_1yr <- function(n_it, clim_sd, clim_corr) {
  
  
  init_pop_vec <- runif(200)
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
  
  
  ### set up non-lagged ipm -------------------------------------------------------------------------
  
  clim_mod <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) - s_temp * temp0),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp0),
      
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
      fn = pois(fn_int + fn_slope * log(size_1)),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
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
    define_domains( size = c(1, 115, 200)
    ) %>% 
    define_env_state(
      env_params = env_sampler(environ_seq = environ_seq,
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
  
  
  ## lagged ipm ----------------------------------------
  
  ipm_lagged <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) - s_temp * temp0),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp1),
      
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
      fn = pois(fn_int + fn_slope * log(size_1)),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
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
    define_domains( size = c(1, 115, 200))  %>%
    define_env_state(
      env_params = env_sampler(environ_seq = environ_seq,
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
             iterate = T, 
             iterations = n_it)
  
  
  ### temperature sd of simulation -------------------------------------------------
  lambdas <- tibble(clim_sd = clim_sd,
                    autocorrelation = clim_corr,
                    ### get lambda non-lagged ---------------------------------------------------------
                    non_lagged = clim_mod$pop_state$lambda[-c(1:round(0.1*n_it))] %>% log %>% mean,   ## discard the first 10% as "burn-in" phase
                    non_lagged_all = list(lambda(clim_mod, "pop_size", "all")[-c(1:(n_it - 1000))]),
                    ### get lambda lagged  ---------------------------------------------------------
                    lagged = ipm_lagged$pop_state$lambda[-c(1:round(0.1*n_it))] %>% log %>% mean,   ## discard the first 10% as "burn-in" phase
                    lagged_all = list(lambda(ipm_lagged, "pop_size", "all")[-c(1:(n_it - 1000))])
  )
  
  return(lambdas)
}

P_2yr <- function(n_it, clim_sd, clim_corr) {
  
  
  init_pop_vec <- runif(200)
  environ_seq <- create_seq(n_it = n_it, clim_sd = clim_sd, clim_corr = clim_corr)
  
  ## Define environment -------------------------------------------------------------------------
  
  env_sampler <- function(environ_seq, iteration) {
    
    temp <- list("temp0" = environ_seq[iteration + 2],
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
  
  
  ### set up non-lagged ipm -------------------------------------------------------------------------
  
  clim_mod <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp0),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp0),
      
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
      fn = pois(fn_int + fn_slope * log(size_1)),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
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
    define_domains( size = c(1, 115, 200)
    ) %>% 
    define_env_state(
      env_params = env_sampler(environ_seq = environ_seq,
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
  
  
  ## lagged ipm ----------------------------------------
  
  ipm_lagged <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp0),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp1),
      
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
      fn = pois(fn_int + fn_slope * log(size_1)),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
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
    define_domains( size = c(1, 115, 200))  %>%
    define_env_state(
      env_params = env_sampler(environ_seq = environ_seq,
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
             iterate = T, 
             iterations = n_it)
  
  
  ### temperature sd of simulation -------------------------------------------------
  lambdas <- tibble(clim_sd = clim_sd,
                    autocorrelation = clim_corr,
                    ### get lambda non-lagged ---------------------------------------------------------
                    non_lagged = clim_mod$pop_state$lambda[-c(1:round(0.1*n_it))] %>% log %>% mean,   ## discard the first 10% as "burn-in" phase
                    non_lagged_all = list(lambda(clim_mod, "pop_size", "all")[-c(1:(n_it - 1000))]),
                    ### get lambda lagged  ---------------------------------------------------------
                    lagged = ipm_lagged$pop_state$lambda[-c(1:round(0.1*n_it))] %>% log %>% mean,   ## discard the first 10% as "burn-in" phase
                    lagged_all = list(lambda(ipm_lagged, "pop_size", "all")[-c(1:(n_it - 1000))])
  )
  
  return(lambdas)
}

PF_1yr <- function(n_it, clim_sd, clim_corr) {
  
  
  init_pop_vec <- runif(200)
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
  
  
  ### set up non-lagged ipm -------------------------------------------------------------------------
  
  clim_mod <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp0),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp0),
      
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
      
      fp = inv_logit(fpC_int + fpC_slope * log(size_1) + fpC_temp * temp0),
      fn = pois(fnC_int + fnC_slope * log(size_1) + fpC_temp * temp0),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
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
    define_domains( size = c(1, 115, 200)
    ) %>% 
    define_env_state(
      env_params = env_sampler(environ_seq = environ_seq,
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
  
  
  ## lagged ipm ----------------------------------------
  
  ipm_P_lagged <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp1),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp1),
      
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
      
      fp = inv_logit(fpC_int + fpC_slope * log(size_1) + fpC_temp * temp0),
      fn = pois(fnC_int + fnC_slope * log(size_1) + fpC_temp * temp0),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
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
    define_domains( size = c(1, 115, 200))  %>%
    define_env_state(
      env_params = env_sampler(environ_seq = environ_seq,
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
             iterate = T, 
             iterations = n_it)
  
  
  ipm_F_lagged <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp0),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp0),
      
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
      
      fp = inv_logit(fpC_int + fpC_slope * log(size_1) + fpC_temp * temp1),
      fn = pois(fnC_int + fnC_slope * log(size_1) + fpC_temp * temp1),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
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
    define_domains( size = c(1, 115, 200))  %>%
    define_env_state(
      env_params = env_sampler(environ_seq = environ_seq,
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
             iterate = T, 
             iterations = n_it)
  
  
  ### temperature sd of simulation -------------------------------------------------
  lambdas <- tibble(clim_sd = clim_sd,
                    autocorrelation = clim_corr,
                    ### get lambda non-lagged ---------------------------------------------------------
                    non_lagged = clim_mod$pop_state$lambda[-c(1:round(0.1*n_it))] %>% log %>% mean,   ## discard the first 10% as "burn-in" phase
                    non_lagged_all = list(lambda(clim_mod, "pop_size", "all")[-c(1:(n_it - 1000))]),
                    ### get lambda lagged  ---------------------------------------------------------
                    lagged_P = ipm_P_lagged$pop_state$lambda[-c(1:round(0.1*n_it))] %>% log %>% mean,   ## discard the first 10% as "burn-in" phase
                    lagged_P_all = list(lambda(ipm_P_lagged, "pop_size", "all")[-c(1:(n_it - 1000))]),
                    lagged_F = ipm_F_lagged$pop_state$lambda[-c(1:round(0.1*n_it))] %>% log %>% mean,   ## discard the first 10% as "burn-in" phase
                    lagged_F_all = list(lambda(ipm_F_lagged, "pop_size", "all")[-c(1:(n_it - 1000))])
  )
  
  return(lambdas)
}

PF_neg_1yr <- function(n_it, clim_sd, clim_corr) {
  
  
  init_pop_vec <- runif(200)
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
  
  
  ### set up non-lagged ipm -------------------------------------------------------------------------
  
  clim_mod <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp0),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp0),
      
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
      
      fp = inv_logit(fpC_int + fpC_slope * log(size_1) + fpC_temp * temp0),
      fn = pois(fnC_int + fnC_slope * log(size_1) + fpC_temp * temp0),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
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
    define_domains( size = c(1, 115, 200)
    ) %>% 
    define_env_state(
      env_params = env_sampler(environ_seq = environ_seq,
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
  
  
  ## lagged ipm ----------------------------------------
  
  ipm_P_lagged <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp1),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp1),
      
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
      
      fp = inv_logit(fpC_int + fpC_slope * log(size_1) + fpC_temp * temp0),
      fn = pois(fnC_int + fnC_slope * log(size_1) + fpC_temp * temp0),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
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
    define_domains( size = c(1, 115, 200))  %>%
    define_env_state(
      env_params = env_sampler(environ_seq = environ_seq,
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
             iterate = T, 
             iterations = n_it)
  
  
  ipm_F_lagged <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) - s_temp * temp0),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp0),
      
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
      
      fp = inv_logit(fpC_int + fpC_slope * log(size_1) + fpC_temp * temp1),
      fn = pois(fnC_int + fnC_slope * log(size_1) + fpC_temp * temp1),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
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
    define_domains( size = c(1, 115, 200))  %>%
    define_env_state(
      env_params = env_sampler(environ_seq = environ_seq,
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
             iterate = T, 
             iterations = n_it)
  
  
  ### temperature sd of simulation -------------------------------------------------
  lambdas <- tibble(clim_sd = clim_sd,
                    autocorrelation = clim_corr,
                    ### get lambda non-lagged ---------------------------------------------------------
                    non_lagged = clim_mod$pop_state$lambda[-c(1:round(0.1*n_it))] %>% log %>% mean,
                    non_lagged_all = list(lambda(clim_mod, "pop_size", "all")[-c(1:(n_it - 1000))]),
                    ### get lambda lagged  ---------------------------------------------------------
                    lagged_P = ipm_P_lagged$pop_state$lambda[-c(1:round(0.1*n_it))] %>% log %>% mean,
                    lagged_P_all = list(lambda(ipm_P_lagged, "pop_size", "all")[-c(1:(n_it - 1000))]),
                    lagged_F = ipm_F_lagged$pop_state$lambda[-c(1:round(0.1*n_it))] %>% log %>% mean,
                    lagged_F_all = list(lambda(ipm_F_lagged, "pop_size", "all")[-c(1:(n_it - 1000))])
  )
  
  return(lambdas)
}

P_neg_1yr_man <- function(n_it, clim_sd) {
  
  
  init_pop_vec <- runif(200)
  environ_seq <- rep(c(rep(0,99),2), round(n_it/100))
  
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
  
  
  ### set up non-lagged ipm -------------------------------------------------------------------------
  
  clim_mod <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) - s_temp * temp0),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp0),
      
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
      fn = pois(fn_int + fn_slope * log(size_1)),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
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
    define_domains( size = c(1, 115, 200)
    ) %>% 
    define_env_state(
      env_params = env_sampler(environ_seq = environ_seq,
                               iteration = t),
      data_list = list(
        environ_seq = environ_seq,
        env_sampler = env_sampler
      )
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_size_t = init_pop_vec,
        
      )
    ) %>%
    make_ipm(usr_funs = my_functions,
             iterate = TRUE,
             iterations = n_it,
             normalize_pop_size = FALSE)
  
  
  ## lagged ipm ----------------------------------------
  
  ipm_lagged <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) - s_temp * temp0),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp1),
      
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
      fn = pois(fn_int + fn_slope * log(size_1)),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
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
    define_domains( size = c(1, 115, 200))  %>%
    define_env_state(
      env_params = env_sampler(environ_seq = environ_seq,
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
             iterate = T, 
             iterations = n_it,
             normalize_pop_size = FALSE)
  
  
  ### temperature sd of simulation -------------------------------------------------
  lambdas <- tibble(clim_sd = clim_sd,
                    autocorrelation = clim_corr,
                    ### get lambda non-lagged ---------------------------------------------------------
                    non_lagged = clim_mod$pop_state$lambda[-c(1:round(0.1*n_it))] %>% log %>% mean,
                    non_lagged_all = list(lambda(clim_mod, "pop_size", "all")[-c(1:(n_it - 1000))]),
                    ### get lambda lagged  ---------------------------------------------------------
                    lagged = ipm_lagged$pop_state$lambda[-c(1:round(0.1*n_it))] %>% log %>% mean,
                    lagged_all = list(lambda(ipm_lagged, "pop_size", "all")[-c(1:(n_it - 1000))]),
                    ### population sizes
                    pop_non_lagged = list(clim_mod$pop_state$pop_state_size),
                    pop_lagged = list(ipm_lagged$pop_state$pop_state_size)
  )
  
  return(lambdas)
}


