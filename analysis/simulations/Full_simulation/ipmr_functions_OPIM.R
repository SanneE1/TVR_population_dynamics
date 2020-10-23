
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


P_lambdas <- function(n_it, clim_sd, clim_corr, params_list, clim_params, n_mesh = 200, save_K = FALSE, n_save_K = 0.1) {
  
  init_pop_vec <- runif(n_mesh)
  environ_seq <- create_seq(n_it = n_it, clim_sd = clim_sd, clim_corr = clim_corr)
  
  params_list <- append(params_list, clim_params)
  
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
      fn = rnorm(1, mean = fn_int, sd = fn_slope),
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
             iterations = n_it,
             report_progress = T)
  
  message("ipm 1 done")
  lambdas <- tibble(clim_sd = clim_sd,
                    autocorrelation = clim_corr,
                    s_temp = clim_params$s_temp,
                    g_temp = clim_params$g_temp,
                    ### get lambda non-lagged ---------------------------------------------------------
                    non_lagged = lambda(clim_mod, "pop_size", "stochastic"), 
                    non_lagged_all = list(lambda(clim_mod, "pop_size", "all")),
                    n_env_seq = list(clim_mod$env_seq))
  
  message("trying to assign K matrices to tibble")
  
  if(save_K == T) {
    lambdas$M_non_lagged <- list(clim_mod$iterators[c((n_it - (n_it * n_save_K)):n_it)])
  }
  message("done")
  
  # remove clim_mod object to save memory
  rm(clim_mod)
  
  ## lagged s ipm ----------------------------------------
  message("starting 2nd ipm")  
  ipm_s_lagged <- init_ipm("simple_di_stoch_param") %>%
    
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp1),
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
      fn = rnorm(1, mean = fn_int, sd = fn_slope),
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
    define_domains( size = c(1, 115, n_mesh))  %>%
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
  
  message("done")
  lambdas$lagged_s <- lambda(ipm_s_lagged, "pop_size", "stochastic")
  lambdas$lagged_s_all <- list(lambda(ipm_s_lagged, "pop_size", "all"))
  lambdas$s_env_seq <- list(ipm_s_lagged$env_seq)
  
  message("assigning 2nd set K kernels")
  if(save_K == T){
    lambdas$M_s_lagged <- list(ipm_s_lagged$iterators[c((n_it - (n_it * n_save_K)):n_it)])
  }
  message("done")
  # remove clim_mod object to save memory
  rm(ipm_s_lagged)
  
  ## lagged g ipm ----------------------------------------
  message("starting 3rd ipm")
  ipm_g_lagged <- init_ipm("simple_di_stoch_param") %>%
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
      fn = rnorm(1, mean = fn_int, sd = fn_slope),
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
    define_domains( size = c(1, 115, n_mesh))  %>%
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
  
  message("done")
  lambdas$lagged_g <- lambda(ipm_g_lagged, "pop_size", "stochastic")
  lambdas$lagged_g_all <- list(lambda(ipm_g_lagged, "pop_size", "all"))
  lambdas$g_env_seq <- list(ipm_g_lagged$env_seq)
  
  if(save_K == T) {
    lambdas$M_g_lagged <- list(ipm_g_lagged$iterators[c((n_it - (n_it * n_save_K)):n_it)])
  }
  message("done assigning 3rd K kernels")
  
  return(lambdas)
}

# PF_lambdas <- function(n_it, clim_sd, clim_corr, params_list, clim_params) {
#   
#   init_pop_vec <- runif(200)
#   environ_seq <- create_seq(n_it = n_it, clim_sd = clim_sd, clim_corr = clim_corr)
#   
#   params_list <- append(params_list, clim_params)
#   
#   ## Define environment -------------------------------------------------------------------------
#   
#   env_sampler <- function(environ_seq, iteration) {
#     
#     temp <- list("temp0" = environ_seq[iteration + 1],
#                  "temp1" = environ_seq[iteration]
#     )
#     
#     return(temp)
#   }
#   
#   
#   ## create custom functions -------------------------------------------------------------------------
#   
#   inv_logit <- function(x) {
#     return(
#       1/(1 + exp(-(x)))
#     )
#   }
#   
#   pois <- function(x) {
#     return(
#       exp(x)
#     )
#   }
#   
#   my_functions <- list(inv_logit = inv_logit,
#                        pois = pois,
#                        env_sampler = env_sampler)
#   
#   
#   ### set up non-lagged ipm -------------------------------------------------------------------------
#   
#   clim_mod <- init_ipm("simple_di_stoch_param") %>%
#     define_kernel(
#       name = "P",
#       
#       formula = s * g,
#       family = "CC",
#       
#       s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp0),
#       g = dnorm(size_2, mean = g_mean, sd = g_sd),
#       g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp0),
#       
#       data_list = params_list,
#       states = list(c('size')),
#       
#       has_hier_effs = FALSE,
#       
#       evict_cor = TRUE,
#       evict_fun = truncated_distributions("norm", "g")
#     ) %>%
#     define_kernel(
#       name = "F",
#       
#       formula = fp * fn * germ * fd,
#       family = "CC",
#       
#       fp = inv_logit(fpC_int + fpC_slope * log(size_1) + fpC_temp * temp0),
#       fn = pois(fnC_int + fnC_slope * log(size_1) + fnC_temp * temp0),
#       germ = germ_mean,
#       fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
#       
#       data_list = params_list,
#       states = list(c("size")),
#       
#       has_hier_effs = FALSE,
#       
#       evict_cor = TRUE,
#       evict_fun = truncated_distributions("norm", "fd")
#     ) %>%
#     define_k(
#       name = "K",
#       family = "IPM",
#       K = P + F,
#       n_size_t_1 = K %*% n_size_t,
#       data_list = list(),
#       states = list(c("size")),
#       has_hier_effs = FALSE,
#       
#       evict_cor = FALSE
#     ) %>% 
#     define_impl(
#       make_impl_args_list(
#         kernel_names = c("K", "P", "F"),
#         int_rule = rep("midpoint", 3),
#         dom_start = rep("size", 3),
#         dom_end = rep("size", 3)
#       )
#     ) %>%
#     define_domains( size = c(1, 115, 200)
#     ) %>% 
#     define_env_state(
#       env_params = env_sampler(environ_seq = environ_seq,
#                                iteration = t),
#       data_list = list(
#         environ_seq = environ_seq,
#         env_sampler = env_sampler
#       )
#     ) %>%
#     define_pop_state(
#       pop_vectors = list(
#         n_size_t = init_pop_vec
#       )
#     ) %>%
#     make_ipm(usr_funs = my_functions,
#              iterate = TRUE,
#              iterations = n_it)
#   
#   
#   lambdas <- tibble(clim_sd = clim_sd,
#                     autocorrelation = clim_corr,
#                     s_temp = clim_params$s_temp,
#                     g_temp = clim_params$g_temp,
#                     fpC_temp  = clim_params$fpC_temp,
#                     fnC_temp = clim_params$fnC_temp,
#                     ### get lambda non-lagged ---------------------------------------------------------
#                     non_lagged = lambda(clim_mod, "pop_size", "stochastic"), 
#                     non_lagged_all = list(lambda(clim_mod, "pop_size", "all")))
#   
#   # remove clim_mod object to save memory
#   rm(clim_mod)
#   
#   ## lagged P ipm ----------------------------------------
#   
#   ipm_P_lagged <- init_ipm("simple_di_stoch_param") %>%
#     define_kernel(
#       name = "P",
#       
#       formula = s * g,
#       family = "CC",
#       
#       s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp1),
#       g = dnorm(size_2, mean = g_mean, sd = g_sd),
#       g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp1),
#       
#       data_list = params_list,
#       states = list(c('size')),
#       
#       has_hier_effs = FALSE,
#       
#       evict_cor = TRUE,
#       evict_fun = truncated_distributions("norm", "g")
#     ) %>%
#     define_kernel(
#       name = "F",
#       
#       formula = fp * fn * germ * fd,
#       family = "CC",
#       
#       fp = inv_logit(fpC_int + fpC_slope * log(size_1) + fpC_temp * temp0),
#       fn = pois(fnC_int + fnC_slope * log(size_1) + fnC_temp * temp0),
#       germ = germ_mean,
#       fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
#       
#       data_list = params_list,
#       states = list(c("size")),
#       
#       has_hier_effs = FALSE,
#       
#       evict_cor = TRUE,
#       evict_fun = truncated_distributions("norm", "fd")
#     ) %>%
#     define_k(
#       name = "K",
#       family = "IPM",
#       K = P + F,
#       n_size_t_1 = K %*% n_size_t,
#       data_list = list(),
#       states = list(c("size")),
#       has_hier_effs = FALSE,
#       
#       evict_cor = FALSE
#     ) %>% 
#     define_impl(
#       make_impl_args_list(
#         kernel_names = c("K", "P", "F"),
#         int_rule = rep("midpoint", 3),
#         dom_start = rep("size", 3),
#         dom_end = rep("size", 3)
#       )
#     ) %>%
#     define_domains( size = c(1, 115, 200))  %>%
#     define_env_state(
#       env_params = env_sampler(environ_seq = environ_seq,
#                                iteration = t),
#       data_list = list(
#         environ_seq = environ_seq,
#         env_sampler = env_sampler
#       )
#     ) %>%
#     define_pop_state(
#       pop_vectors = list(
#         n_size_t = init_pop_vec
#       )
#     ) %>%
#     make_ipm(usr_funs = my_functions, 
#              iterate = T, 
#              iterations = n_it)
#   
#   
#   lambdas$lagged_P <- lambda(ipm_P_lagged, "pop_size", "stochastic")
#   lambdas$lagged_P_all <- list(lambda(ipm_P_lagged, "pop_size", "all"))
#   
#   # remove clim_mod object to save memory
#   rm(ipm_P_lagged)
#   
#   ## lagged F ipm ----------------------------------------
#   
#   ipm_F_lagged <- init_ipm("simple_di_stoch_param") %>%
#     define_kernel(
#       name = "P",
#       
#       formula = s * g,
#       family = "CC",
#       
#       s = inv_logit(s_int + s_slope * log(size_1) + s_temp * temp0),
#       g = dnorm(size_2, mean = g_mean, sd = g_sd),
#       g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * temp0),
#       
#       data_list = params_list,
#       states = list(c('size')),
#       
#       has_hier_effs = FALSE,
#       
#       evict_cor = TRUE,
#       evict_fun = truncated_distributions("norm", "g")
#     ) %>%
#     define_kernel(
#       name = "F",
#       
#       formula = fp * fn * germ * fd,
#       family = "CC",
#       
#       fp = inv_logit(fpC_int + fpC_slope * log(size_1) + fpC_temp * temp1),
#       fn = pois(fnC_int + fnC_slope * log(size_1) + fnC_temp * temp1),
#       germ = germ_mean,
#       fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
#       
#       data_list = params_list,
#       states = list(c("size")),
#       
#       has_hier_effs = FALSE,
#       
#       evict_cor = TRUE,
#       evict_fun = truncated_distributions("norm", "fd")
#     ) %>%
#     define_k(
#       name = "K",
#       family = "IPM",
#       K = P + F,
#       n_size_t_1 = K %*% n_size_t,
#       data_list = list(),
#       states = list(c("size")),
#       has_hier_effs = FALSE,
#       
#       evict_cor = FALSE
#     ) %>% 
#     define_impl(
#       make_impl_args_list(
#         kernel_names = c("K", "P", "F"),
#         int_rule = rep("midpoint", 3),
#         dom_start = rep("size", 3),
#         dom_end = rep("size", 3)
#       )
#     ) %>%
#     define_domains( size = c(1, 115, 200))  %>%
#     define_env_state(
#       env_params = env_sampler(environ_seq = environ_seq,
#                                iteration = t),
#       data_list = list(
#         environ_seq = environ_seq,
#         env_sampler = env_sampler
#       )
#     ) %>%
#     define_pop_state(
#       pop_vectors = list(
#         n_size_t = init_pop_vec
#       )
#     ) %>%
#     make_ipm(usr_funs = my_functions, 
#              iterate = T, 
#              iterations = n_it)
#   
#   
#   lambdas$lagged_F <- lambda(ipm_F_lagged, "pop_size", "stochastic")
#   lambdas$lagged_F_all <- list(lambda(ipm_F_lagged, "pop_size", "all"))
#   
#   return(lambdas)
# }
# 
