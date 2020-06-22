

P_1yr <- function(temp_sd) {
  
  print(temp_sd)

  init_pop_vec <- runif(200)
  
  ## Define environment -------------------------------------------------------------------------
  
  sample_env <- function(env_params) {
    temp <- rnorm(1,
                  env_params$temp_mu,
                  env_params$temp_sd)
    out <- list(temp = temp)
    return(out)
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
                       sample_env = sample_env)
  
  
  ### Create function to run different climate sd's in parlapply ----------------------------------------------------------------
  
  
  env_params <- list(
    temp_mu = 0,
    temp_sd = temp_sd
  )
  
  ### set up non-lagged ipm -------------------------------------------------------------------------
  
  clim_mod <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp),
      
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
    define_domains( size = c(1, 115, 200)) 
  
  ### iterate non-lagged ipm --------------------------------------------------------
  
  ipm_non_lagged <- clim_mod %>%
    define_env_state(
      env_covs = sample_env(env_params),
      data_list = list(env_params = env_params)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_size_t = init_pop_vec
      )
    ) %>%
    make_ipm(usr_funs = my_functions,
             iterate = TRUE,
             iterations = 1000)
  
  
  
  #### lagged climate vector ------------------------------------------------------
  temp_lagged <- data.frame(n = c(1:1001), 
                            temp0 = rnorm(1001, sd = env_params$temp_sd))
  temp_lagged <- left_join(temp_lagged, 
                           temp_lagged %>%
                             mutate(n = n + 1) %>%
                             rename(temp1 = temp0)
  )
  
  print("temp_lagged")
  head(temp_lagged)

  params_list2 <- append(params_list,
                         list(temp0 = NA,
                              temp1 = NA))
  la <- NA
  pop_vect <- list(init_pop_vec)
  ### Iterate with lagged climate -----------------------------------------------------
  
  for(i in c(2:1001)){
    params_list2$temp0 <- temp_lagged$temp0[i]
    params_list2$temp1 <- temp_lagged$temp1[i]
    
    pop_vec <- pop_vect[[i-1]]
    ipm_lagged <- init_ipm("simple_di_stoch_param") %>%
      define_kernel(
        name = "P",
        
        formula = s * g,
        family = "CC",
        
        s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp0),
        g = dnorm(size_2, mean = g_mean, sd = g_sd),
        g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp1),
        
        data_list = params_list2,
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
        
        data_list = params_list2,
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
        env_covs = list(temp0 = params_list2$temp0,
                        temp1 = params_list2$temp1),
        data_list = list(params_list2 = params_list2)
      ) %>%
      define_pop_state(
        pop_vectors = list(
          n_size_t = pop_vec
        )
      ) %>%
      make_ipm(usr_funs = my_functions, 
               iterate = T, 
               iterations = 1)
    
    la <- append(la, ipm_lagged$pop_state$lambda)
    pop_vect <- append(pop_vect, list(ipm_lagged$pop_state$pop_state_size[,2]))
  }
  
  print("la")
  head(la)
  print("lambdas")
  head(ipm_non_lagged$pop_state$lambda)

  ### temperature sd of simulation -------------------------------------------------
  lambdas <- data.frame(temperature_sd = env_params$temp_sd,
                        ### get lambda non-lagged ---------------------------------------------------------
                        non_lagged = ipm_non_lagged$pop_state$lambda %>% log %>% mean,
                        non_lagged_sd = ipm_non_lagged$pop_state$lambda %>% log %>% sd,
                        ### get lambda lagged  ---------------------------------------------------------
                        lagged = la %>% log %>% mean(na.rm = T),
                        lagged_sd = la %>% log %>% sd(na.rm = T)
  )
  
  return(lambdas)
}

P_neg_1yr <- function(temp_sd) {
  
  init_pop_vec <- runif(200)
  
  ## Define environment -------------------------------------------------------------------------
  
  sample_env <- function(env_params) {
    temp <- rnorm(1,
                  env_params$temp_mu,
                  env_params$temp_sd)
    out <- list(temp = temp)
    return(out)
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
                       sample_env = sample_env)
  
  
  ### Create function to run different climate sd's in parlapply ----------------------------------------------------------------
  
  
  env_params <- list(
    temp_mu = 0,
    temp_sd = temp_sd
  )
  
  ### set up non-lagged ipm -------------------------------------------------------------------------
  
  clim_mod <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) - s_temp * temp),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp),
      
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
    define_domains( size = c(1, 115, 200)) 
  
  ### iterate non-lagged ipm --------------------------------------------------------
  
  ipm_non_lagged <- clim_mod %>%
    define_env_state(
      env_covs = sample_env(env_params),
      data_list = list(env_params = env_params)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_size_t = init_pop_vec
      )
    ) %>%
    make_ipm(usr_funs = my_functions,
             iterate = TRUE,
             iterations = 1000)
  
  
  
  #### lagged climate vector ------------------------------------------------------
  temp_lagged <- data.frame(n = c(1:1001), 
                            temp0 = rnorm(1001, sd = env_params$temp_sd))
  temp_lagged <- left_join(temp_lagged, 
                           temp_lagged %>%
                             mutate(n = n + 1) %>%
                             rename(temp1 = temp0)
  )
  
  params_list2 <- append(params_list,
                         list(temp0 = NA,
                              temp1 = NA))
  la <- NA
  pop_vect <- list(init_pop_vec)
  ### Iterate with lagged climate -----------------------------------------------------
  
  for(i in c(2:1001)){
    params_list2$temp0 <- temp_lagged$temp0[i]
    params_list2$temp1 <- temp_lagged$temp1[i]
    
    pop_vec <- pop_vect[[i-1]]
    ipm_lagged <- init_ipm("simple_di_stoch_param") %>%
      define_kernel(
        name = "P",
        
        formula = s * g,
        family = "CC",
        
        s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp0),
        g = dnorm(size_2, mean = g_mean, sd = g_sd),
        g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp1),
        
        data_list = params_list2,
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
        
        data_list = params_list2,
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
        env_covs = list(temp0 = params_list2$temp0,
                        temp1 = params_list2$temp1),
        data_list = list(params_list2 = params_list2)
      ) %>%
      define_pop_state(
        pop_vectors = list(
          n_size_t = pop_vec
        )
      ) %>%
      make_ipm(usr_funs = my_functions, 
               iterate = T, 
               iterations = 1)
    
    la <- append(la, ipm_lagged$pop_state$lambda)
    pop_vect <- append(pop_vect, list(ipm_lagged$pop_state$pop_state_size[,2]))
  }
  
  
  ### temperature sd of simulation -------------------------------------------------
  lambdas <- data.frame(temperature_sd = env_params$temp_sd,
                        ### get lambda non-lagged ---------------------------------------------------------
                        non_lagged = ipm_non_lagged$pop_state$lambda %>% log %>% mean,
                        non_lagged_sd = ipm_non_lagged$pop_state$lambda %>% log %>% sd,
                        ### get lambda lagged  ---------------------------------------------------------
                        lagged = la %>% log %>% mean(na.rm = T),
                        lagged_sd = la %>% log %>% sd(na.rm = T)
  )
  
  return(lambdas)
}

P_2yr <- function(temp_sd) {
  
  init_pop_vec <- runif(200)
  
  ## Define environment -------------------------------------------------------------------------
  
  sample_env <- function(env_params) {
    temp <- rnorm(1,
                  env_params$temp_mu,
                  env_params$temp_sd)
    out <- list(temp = temp)
    return(out)
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
                       sample_env = sample_env)
  
  
  ### Create function to run different climate sd's in parlapply ----------------------------------------------------------------
  
  
  env_params <- list(
    temp_mu = 0,
    temp_sd = temp_sd
  )
  
  ### set up non-lagged ipm -------------------------------------------------------------------------
  
  clim_mod <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp),
      
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
    define_domains( size = c(1, 115, 200)) 
  
  ### iterate non-lagged ipm --------------------------------------------------------
  
  ipm_non_lagged <- clim_mod %>%
    define_env_state(
      env_covs = sample_env(env_params),
      data_list = list(env_params = env_params)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_size_t = init_pop_vec
      )
    ) %>%
    make_ipm(usr_funs = my_functions,
             iterate = TRUE,
             iterations = 1000)
  
  
  
  #### lagged climate vector ------------------------------------------------------
  temp_lagged <- data.frame(n = c(1:1002), 
                            temp0 = rnorm(1002, sd = env_params$temp_sd))
  temp_lagged <- left_join(temp_lagged, 
                           temp_lagged %>%
                             mutate(n = n + 2) %>%
                             rename(temp1 = temp0)
  )
  
  params_list2 <- append(params_list,
                         list(temp0 = NA,
                              temp1 = NA))
  la <- NA
  pop_vect <- list(init_pop_vec,
                   init_pop_vec)
  ### Iterate with lagged climate -----------------------------------------------------
  
  for(i in c(3:1002)){
    params_list2$temp0 <- temp_lagged$temp0[i]
    params_list2$temp1 <- temp_lagged$temp1[i]
    
    pop_vec <- pop_vect[[i-1]]
    ipm_lagged <- init_ipm("simple_di_stoch_param") %>%
      define_kernel(
        name = "P",
        
        formula = s * g,
        family = "CC",
        
        s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp0),
        g = dnorm(size_2, mean = g_mean, sd = g_sd),
        g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp1),
        
        data_list = params_list2,
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
        
        data_list = params_list2,
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
        env_covs = list(temp0 = params_list2$temp0,
                        temp1 = params_list2$temp1),
        data_list = list(params_list2 = params_list2)
      ) %>%
      define_pop_state(
        pop_vectors = list(
          n_size_t = pop_vec
        )
      ) %>%
      make_ipm(usr_funs = my_functions, 
               iterate = T, 
               iterations = 1)
    
    la <- append(la, ipm_lagged$pop_state$lambda)
    pop_vect <- append(pop_vect, list(ipm_lagged$pop_state$pop_state_size[,2]))
  }
  
  
  ### temperature sd of simulation -------------------------------------------------
  lambdas <- data.frame(temperature_sd = env_params$temp_sd,
                        ### get lambda non-lagged ---------------------------------------------------------
                        non_lagged = ipm_non_lagged$pop_state$lambda %>% log %>% mean,
                        non_lagged_sd = ipm_non_lagged$pop_state$lambda %>% log %>% sd,
                        ### get lambda lagged  ---------------------------------------------------------
                        lagged = la %>% log %>% mean(na.rm = T),
                        lagged_sd = la %>% log %>% sd(na.rm = T)
  )
  
  return(lambdas)
}

PF_1yr <- function(temp_sd) {
  
  init_pop_vec <- runif(200)
  
  ## Define environment -------------------------------------------------------------------------
  
  sample_env <- function(env_params) {
    temp <- rnorm(1,
                  env_params$temp_mu,
                  env_params$temp_sd)
    out <- list(temp = temp)
    return(out)
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
                       sample_env = sample_env)
  
  
  ### Create function to run different climate sd's in parlapply ----------------------------------------------------------------
  
  
  env_params <- list(
    temp_mu = 0,
    temp_sd = temp_sd
  )
  
  ### set up non-lagged ipm -------------------------------------------------------------------------
  
  clim_mod <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp),
      
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
      
      fp = inv_logit(fpC_int + fpC_slope * log(size_1) + fpC_temp * temp),
      fn = pois(fnC_int + fnC_slope * log(size_1) + fnC_temp * temp),
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
    define_domains( size = c(1, 115, 200)) 
  
  ### iterate non-lagged ipm --------------------------------------------------------
  
  ipm_non_lagged <- clim_mod %>%
    define_env_state(
      env_covs = sample_env(env_params),
      data_list = list(env_params = env_params)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_size_t = init_pop_vec
      )
    ) %>%
    make_ipm(usr_funs = my_functions,
             iterate = TRUE,
             iterations = 1000)
  
  
  
  #### lagged climate vector ------------------------------------------------------
  temp_lagged <- data.frame(n = c(1:1001), 
                            temp0 = rnorm(1001, sd = env_params$temp_sd))
  temp_lagged <- left_join(temp_lagged, 
                           temp_lagged %>%
                             mutate(n = n + 1) %>%
                             rename(temp1 = temp0)
  )
  
  params_list2 <- append(params_list,
                         list(temp0 = NA,
                              temp1 = NA))
  la <- NA
  pop_vect_a <- list(init_pop_vec)
  
  ### Iterate with lagged climate -----------------------------------------------------
  
  ## P kernel lagged
  
  for(i in c(2:1001)){
    params_list2$temp0 <- temp_lagged$temp0[i]
    params_list2$temp1 <- temp_lagged$temp1[i]
    
    pop_vec <- pop_vect_a[[i-1]]
    ipm_lagged <- init_ipm("simple_di_stoch_param") %>%
      define_kernel(
        name = "P",
        
        formula = s * g,
        family = "CC",
        
        s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp1),
        g = dnorm(size_2, mean = g_mean, sd = g_sd),
        g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp1),
        
        data_list = params_list2,
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
        
        data_list = params_list2,
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
        env_covs = list(temp0 = params_list2$temp0,
                        temp1 = params_list2$temp1),
        data_list = list(params_list2 = params_list2)
      ) %>%
      define_pop_state(
        pop_vectors = list(
          n_size_t = pop_vec
        )
      ) %>%
      make_ipm(usr_funs = my_functions, 
               iterate = T, 
               iterations = 1)
    
    la <- append(la, ipm_lagged$pop_state$lambda)
    pop_vect_a <- append(pop_vect_a, list(ipm_lagged$pop_state$pop_state_size[,2]))
  }
  
  ## F kernel lagged
  
  lb <- NA
  pop_vect_b <- list(init_pop_vec)
  ### Iterate with lagged climate -----------------------------------------------------
  
  for(i in c(2:1001)){
    params_list2$temp0 <- temp_lagged$temp0[i]
    params_list2$temp1 <- temp_lagged$temp1[i]
    
    pop_vec <- pop_vect_b[[i-1]]
    ipm_lagged <- init_ipm("simple_di_stoch_param") %>%
      define_kernel(
        name = "P",
        
        formula = s * g,
        family = "CC",
        
        s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp0),
        g = dnorm(size_2, mean = g_mean, sd = g_sd),
        g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp0),
        
        data_list = params_list2,
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
        
        data_list = params_list2,
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
        env_covs = list(temp0 = params_list2$temp0,
                        temp1 = params_list2$temp1),
        data_list = list(params_list2 = params_list2)
      ) %>%
      define_pop_state(
        pop_vectors = list(
          n_size_t = pop_vec
        )
      ) %>%
      make_ipm(usr_funs = my_functions, 
               iterate = T, 
               iterations = 1)
    
    lb <- append(lb, ipm_lagged$pop_state$lambda)
    pop_vect_b <- append(pop_vect_b, list(ipm_lagged$pop_state$pop_state_size[,2]))
  }
  
  ### temperature sd of simulation -------------------------------------------------
  lambdas <- data.frame(temperature_sd = env_params$temp_sd,
                        ### get lambda non-lagged ---------------------------------------------------------
                        non_lagged = ipm_non_lagged$pop_state$lambda %>% log %>% mean,
                        non_lagged_sd = ipm_non_lagged$pop_state$lambda %>% log %>% sd,
                        ### get lambda lagged P kernel ---------------------------------------------------------
                        lagged_P = la %>% log %>% mean(na.rm = T),
                        lagged_P_sd = la %>% log %>% sd(na.rm = T),
                        ### get lambda lagged F kernel ---------------------------------------------------------
                        lagged_F = lb %>% log %>% mean(na.rm = T),
                        lagged_F_sd = lb %>% log %>% sd(na.rm = T)
  )
  
  return(lambdas)
}

PF_neg_1yr <- function(temp_sd) {
  
  init_pop_vec <- runif(200)
  
  ## Define environment -------------------------------------------------------------------------
  
  sample_env <- function(env_params) {
    temp <- rnorm(1,
                  env_params$temp_mu,
                  env_params$temp_sd)
    out <- list(temp = temp)
    return(out)
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
                       sample_env = sample_env)
  
  
  ### Create function to run different climate sd's in parlapply ----------------------------------------------------------------
  
  
  env_params <- list(
    temp_mu = 0,
    temp_sd = temp_sd
  )
  
  ### set up non-lagged ipm -------------------------------------------------------------------------
  
  clim_mod <- init_ipm("simple_di_stoch_param") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp),
      
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
      
      fp = inv_logit(fpC_int + fpC_slope * log(size_1) + fpC_temp * temp),
      fn = pois(fnC_int + fnC_slope * log(size_1) + fnC_temp * temp),
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
    define_domains( size = c(1, 115, 200)) 
  
  ### iterate non-lagged ipm --------------------------------------------------------
  
  ipm_non_lagged <- clim_mod %>%
    define_env_state(
      env_covs = sample_env(env_params),
      data_list = list(env_params = env_params)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_size_t = init_pop_vec
      )
    ) %>%
    make_ipm(usr_funs = my_functions,
             iterate = TRUE,
             iterations = 1000)
  
  
  
  #### lagged climate vector ------------------------------------------------------
  temp_lagged <- data.frame(n = c(1:1001), 
                            temp0 = rnorm(1001, sd = env_params$temp_sd))
  temp_lagged <- left_join(temp_lagged, 
                           temp_lagged %>%
                             mutate(n = n + 1) %>%
                             rename(temp1 = temp0)
  )
  
  params_list2 <- append(params_list,
                         list(temp0 = NA,
                              temp1 = NA))
  la <- NA
  pop_vect_a <- list(init_pop_vec)
  
  ### Iterate with lagged climate -----------------------------------------------------
  
  ## P kernel lagged
  
  for(i in c(2:1001)){
    params_list2$temp0 <- temp_lagged$temp0[i]
    params_list2$temp1 <- temp_lagged$temp1[i]
    
    pop_vec <- pop_vect_a[[i-1]]
    ipm_lagged <- init_ipm("simple_di_stoch_param") %>%
      define_kernel(
        name = "P",
        
        formula = s * g,
        family = "CC",
        
        s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp1),
        g = dnorm(size_2, mean = g_mean, sd = g_sd),
        g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp1),
        
        data_list = params_list2,
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
        
        data_list = params_list2,
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
        env_covs = list(temp0 = params_list2$temp0,
                        temp1 = params_list2$temp1),
        data_list = list(params_list2 = params_list2)
      ) %>%
      define_pop_state(
        pop_vectors = list(
          n_size_t = pop_vec
        )
      ) %>%
      make_ipm(usr_funs = my_functions, 
               iterate = T, 
               iterations = 1)
    
    la <- append(la, ipm_lagged$pop_state$lambda)
    pop_vect_a <- append(pop_vect_a, list(ipm_lagged$pop_state$pop_state_size[,2]))
  }
  
  ## F kernel lagged
  
  lb <- NA
  pop_vect_b <- list(init_pop_vec)
  ### Iterate with lagged climate -----------------------------------------------------
  
  for(i in c(2:1001)){
    params_list2$temp0 <- temp_lagged$temp0[i]
    params_list2$temp1 <- temp_lagged$temp1[i]
    
    pop_vec <- pop_vect_b[[i-1]]
    ipm_lagged <- init_ipm("simple_di_stoch_param") %>%
      define_kernel(
        name = "P",
        
        formula = s * g,
        family = "CC",
        
        s = inv_logit(s_int + s_slope * log(size_1) - s_temp * temp0),
        g = dnorm(size_2, mean = g_mean, sd = g_sd),
        g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp0),
        
        data_list = params_list2,
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
        
        data_list = params_list2,
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
        env_covs = list(temp0 = params_list2$temp0,
                        temp1 = params_list2$temp1),
        data_list = list(params_list2 = params_list2)
      ) %>%
      define_pop_state(
        pop_vectors = list(
          n_size_t = pop_vec
        )
      ) %>%
      make_ipm(usr_funs = my_functions, 
               iterate = T, 
               iterations = 1)
    
    lb <- append(lb, ipm_lagged$pop_state$lambda)
    pop_vect_b <- append(pop_vect_b, list(ipm_lagged$pop_state$pop_state_size[,2]))
  }
  
  ### temperature sd of simulation -------------------------------------------------
  lambdas <- data.frame(temperature_sd = env_params$temp_sd,
                        ### get lambda non-lagged ---------------------------------------------------------
                        non_lagged = ipm_non_lagged$pop_state$lambda %>% log %>% mean,
                        non_lagged_sd = ipm_non_lagged$pop_state$lambda %>% log %>% sd,
                        ### get lambda lagged P kernel ---------------------------------------------------------
                        lagged_P = la %>% log %>% mean(na.rm = T),
                        lagged_P_sd = la %>% log %>% sd(na.rm = T),
                        ### get lambda lagged F kernel ---------------------------------------------------------
                        lagged_F = lb %>% log %>% mean(na.rm = T),
                        lagged_F_sd = lb %>% log %>% sd(na.rm = T)
  )
  
  return(lambdas)
}

P_neg_1yr_man <- function(temp_sd) {
  
  init_pop_vec <- runif(200) * 5
  

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
                       pois = pois)
  
  #### (lagged) climate vector ------------------------------------------------------
  temp_lagged <- data.frame(n = c(1:150), 
                            temp0 = c(rep(0,100), temp_sd, rep(0,49)))
  temp_lagged <- left_join(temp_lagged, 
                           temp_lagged %>%
                             mutate(n = n + 1) %>%
                             rename(temp1 = temp0)
  )
  
  params_list2 <- append(params_list,
                         list(temp0 = NA,
                              temp1 = NA))
  la <- NA
  pop_vect <- list(init_pop_vec)
  type <- NA  
  temp0_val <- NA
  temp1_val <- NA
  pop_n <- NA
  
  ### Iterate with lagged climate -----------------------------------------------------
  for(j in c("lagged", "non-lagged")) {

  for(i in c(2:151)){
    
    if(j == "lagged"){
      params_list2$temp0 <- temp_lagged$temp0[i]
      params_list2$temp1 <- temp_lagged$temp1[i]
    }
    if(j == "non-lagged"){
      params_list2$temp0 <- temp_lagged$temp0[i]
      params_list2$temp1 <- temp_lagged$temp0[i]
    }    
    
    pop_vec <- pop_vect[[i-1]]
    ipm_lagged <- init_ipm("simple_di_stoch_param") %>%
      define_kernel(
        name = "P",
        
        formula = s * g,
        family = "CC",
        
        s = inv_logit(s_int + s_slope * log(size_1) - s_temp * temp0),
        g = dnorm(size_2, mean = g_mean, sd = g_sd),
        g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp1),
        
        data_list = params_list2,
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
        
        data_list = params_list2,
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
        env_covs = list(temp0 = params_list2$temp0,
                        temp1 = params_list2$temp1),
        data_list = list(params_list2 = params_list2)
      ) %>%
      define_pop_state(
        pop_vectors = list(
          n_size_t = pop_vec
        )
      ) %>%
      make_ipm(usr_funs = my_functions, 
               iterate = T, 
               iterations = 1,
               normalize_pop_size = FALSE)
    
    la <- append(la, ipm_lagged$pop_state$lambda)
    pop_vect <- append(pop_vect, list(ipm_lagged$pop_state$pop_state_size[,2]))
    pop_n <- append(pop_n, sum(ipm_lagged$pop_state$pop_state_size[,2]))
    type <- append(type, j)
    temp0_val <- append(temp0_val, params_list2$temp0)
    temp1_val <- append(temp1_val, params_list2$temp1)
  
  } 
  }
  
  ### temperature sd of simulation -------------------------------------------------
  lambdas <- list(lambdas = la,
                  pop_n = pop_n,
                  pop_vect = pop_vect,
                  type = type,
                  value_temp0 = temp0_val,
                  value_temp1 = temp1_val)
  
  return(lambdas)
}

