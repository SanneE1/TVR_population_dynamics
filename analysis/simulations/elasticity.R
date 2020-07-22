### Run elasticity analysis

library(ipmr)
library(dplyr)
library(ggplot2)


params_list <- list(
  s_int = -0.229,
  s_slope = 1.077,
  s_temp = 1.233,
  g_int = 0.424,
  g_slope = 0.846,
  g_temp = -0.066,
  g_sd = 1.076,
  fp_int = -3.970,
  fp_slope = 1.719,
  fpC_int = -3.859385,
  fpC_slope = 1.719768,
  fpC_temp = -0.6169492,
  fn_int = -0.6652477,
  fn_slope = 0.7048809,
  fnC_int = -0.5661762,
  fnC_slope = 0.7048782,
  fnC_temp = -0.3398345,
  germ_mean = 0.1262968,
  germ_sd = 0.2725941,
  fd_mean = 1.178749,
  fd_sd = 0.8747764
)


# ipm set up

init_pop_vec <- runif(200)

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


## create dataframes for ipm run

pertubation <- 0.1 
levels <- c(-1, 0, 1) ## climate anomaly levels at which to test sensitivity/elasticity


lambda <- expand.grid(list(param = c("As", "Ag", "Afp", "Afn"), 
                           levels = levels,
                           pertubation = c(pertubation, 0, -pertubation))) 
lambda$lambda <- NA

params_list <- append(params_list, list(
  As = 0,
  Ag = 0,
  Afp = 0,
  Afn = 0
))

for (j in c(1:length(lambda[,1]))) {
  
  pertub_params <- append(params_list, list(levels = lambda$levels[j]))
  pertub_params[as.character(lambda$param[j])] <- lambda$pertubation[j]
  
  ipm <- init_ipm("simple_di_det") %>%
    define_kernel(
      name = "P",
      
      formula = s * g,
      family = "CC",
      
      s = inv_logit(s_int + s_slope * log(size_1) + s_temp * levels + As),
      g = dnorm(size_2, mean = g_mean, sd = g_sd),
      g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * levels + Ag),
      
      data_list = pertub_params,
      states = list(c('size')),
      
      has_hier_effs = FALSE,
      
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "g")
    ) %>%
    define_kernel(
      name = "F",
      
      formula = fp * fn * germ * fd,
      family = "CC",
      
      fp = inv_logit(fp_int + fp_slope * log(size_1) + Afp),
      fn = pois(fn_int + fn_slope * log(size_1) + Afn),
      germ = germ_mean,
      fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
      
      data_list = pertub_params,
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
    make_ipm(usr_funs = my_functions)
  
  
  lambda$lambda[j] <- ipmr::lambda(ipm = ipm, comp_method = "eigen")
  
} 

sensitivity <- lambda %>% group_by(param, levels) %>% 
  pivot_wider(names_from = pertubation, values_from = lambda) %>%
  mutate(sens = (`0.1` - `-0.1`) / 2*pertubation)



ggplot(sensitivity, aes(x = param, y = sens, group = levels)) + 
  geom_bar(aes(fill = as.factor(levels)), stat = "identity", position = position_dodge()) + 
  scale_fill_brewer(palette = "Greys") + 
  theme_classic() + 
  labs(fill = "Temperature \nanomaly")
