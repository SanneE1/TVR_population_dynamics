### Run elasticity analysis for the ipm (first part) and mpm (second part) with only P kernel sensitive to climate

library(ipmr)
library(dplyr)
library(ggplot2)
library(tidyr)


#-------------------------------------------------
# IPM elasticity
#-------------------------------------------------
# 
# params_list <- list(
#   s_int = -0.229,
#   s_slope = 1.077,
#   s_temp = 1.233,
#   g_int = 0.424,
#   g_slope = 0.846,
#   g_temp = -0.066,
#   g_sd = 1.076,
#   fp_int = -3.970,
#   fp_slope = 1.719,
#   # fpC_int = -3.859385,
#   # fpC_slope = 1.719768,
#   # fpC_temp = -0.6169492,
#   fn_int = -0.6652477,
#   fn_slope = 0.7048809,
#   # fnC_int = -0.5661762,
#   # fnC_slope = 0.7048782,
#   # fnC_temp = -0.3398345,
#   germ_mean = 0.1262968,
#   germ_sd = 0.2725941,
#   fd_mean = 1.178749,
#   fd_sd = 0.8747764
# )
# 
# 
# # ipm set up
# 
# init_pop_vec <- runif(200)
# 
# ## create custom functions -------------------------------------------------------------------------
# 
# inv_logit <- function(x) {
#   return(
#     1/(1 + exp(-(x)))
#   )
# }
# 
# pois <- function(x) {
#   return(
#     exp(x)
#   )
# }
# 
# my_functions <- list(inv_logit = inv_logit,
#                      pois = pois)
# 
# 
# ## create dataframes for ipm run
# 
# pertubation <- 0.1 
# levels <- c(-1, 0, 1) ## climate anomaly levels at which to test sensitivity/elasticity
# 
# 
# lambda <- expand.grid(list(param = names(params_list), 
#                            levels = levels,
#                            pertubation = c(pertubation, 0, -pertubation))) 
# lambda$lambda <- NA
# 
# 
# for (j in c(1:length(lambda$param))) {
#   
#   pertub_params <- params_list
#   pertub_params[[lambda$param[j]]] <- pertub_params[[lambda$param[j]]] + lambda$pertubation[j]
#   pertub_params <- append(pertub_params, list(levels = lambda$levels[j]))
#   
#   ipm <- init_ipm("simple_di_det") %>%
#     define_kernel(
#       name = "P",
#       
#       formula = s * g,
#       family = "CC",
#       
#       s = inv_logit(s_int + s_slope * log(size_1) + s_temp * levels),
#       g = dnorm(size_2, mean = g_mean, sd = g_sd),
#       g_mean = pois(g_int + g_slope * log(size_1)  + g_temp * levels),
#       
#       data_list = pertub_params,
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
#       fp = inv_logit(fp_int + fp_slope * log(size_1)),
#       fn = pois(fn_int + fn_slope * log(size_1)),
#       germ = germ_mean,
#       fd = dnorm(size_2, mean = fd_mean, sd = fd_sd),
#       
#       data_list = pertub_params,
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
#     define_domains( size = c(1, 115, 50)
#     ) %>%
#     make_ipm(usr_funs = my_functions)
#   
#   
#   lambda$lambda[j] <- ipmr::lambda(ipm = ipm, comp_method = "eigen")
#   
# } 

sensitivity <- lambda %>% group_by(param, levels) %>% 
  pivot_wider(names_from = pertubation, values_from = lambda) %>%
  mutate(sens = (`0.1` - `-0.1`) / 2*pertubation)

sensitivity <- left_join(sensitivity, 
                         tibble::rownames_to_column(as.data.frame(
                           params_list %>% bind_rows() %>% t)) %>% rename(param = rowname,
                                                                          value = V1) ) %>%
  mutate(elasticity = sens * abs(value) / pertubation)
 

sensitivity$vr <- factor(rep(c("surv", "surv", "surv",
                    "growth", "growth", "growth", "growth",
                    "fp", "fp",
                    "fn", "fn",
                    "germ", "germ",
                    "fd", "fd"), 3))

ggplot(sensitivity, aes(x = param, y = elasticity, group = levels)) + 
  geom_bar(aes(fill = as.factor(levels)), stat = "identity", position = position_dodge()) + 
  scale_fill_brewer(palette = "Greens") + 
  scale_x_discrete(limits = c("s_int", "s_slope", "s_temp", "g_int", "g_slope", "g_temp", "g_sd",
                              "fp_int", "fp_slope", "fn_int", "fn_slope", "germ_mean", "germ_sd",
                              "fd_mean", "fd_sd")) +
  theme_classic() + 
  labs(fill = "Temperature \nanomaly")




#-------------------------------------------------
# MPM elasticity (climate insensitive)
#-------------------------------------------------


# Transforms all values below/above limits in min/max size

# Survival function
s_z <- function(z, pars) {
  x <- pars$s_int + pars$s_size * log(z)
  p <- 1/(1+exp(-(x)))

  return(p)
}

# Growth function
g_z1z <- function(z1, z, pars) {


  mean <- exp(pars$g_int + pars$g_size * log(z))  ## mean size next year
  p_grow <- dnorm(z1, mean = mean, sd = pars$g_sd)   ### probability distribution that you grow to size z1
  ### based on mean size and sd
  return(p_grow)
}

# Survival * Growth
P_z1z <- function(z1, z, pars) {
  p <- s_z(z, pars) * g_z1z(z1, z, pars)
  return(p)
}

# Flower probability
fp_z <- function(z, pars) {

  x <- pars$fp_int + pars$fp_size * log(z)
  p <- 1/(1+exp(-(x)))

  return(p)
}
# Seed production
seed_z <- function(z, pars) {

  n <- exp(pars$seed_int + pars$seed_size * log(z))

  return(n)
}
# recruit size
fd_z1 <- function(z1, pars) {

  fd <- dnorm(z1, mean = pars$fd_int, sd = pars$fd_sd)    ## probability of size z1 recruits

  return(fd)
}

# Fecundity function
fec <- function(z1, z, pars) {
  fec <- fp_z(z, pars) * seed_z(z, pars) * pars$germ_int * fd_z1(z1, pars)
  return(fec)
}


pars <- data.frame(s_int = -0.229,
                   s_size = 1.077,
                   g_int = 0.424,
                   g_size = 0.846,
                   g_sd = 1.076,
                   fp_int = -3.970,
                   fp_size = 1.719,
                   fd_int = 1.178749,
                   fd_sd = 0.76,
                   seed_int = -0.6652477,
                   seed_size = 0.8747764,
                   germ_int = 0.1262968,
                   germ_sd = 0.2725941,

                   mat_siz = 100,
                   L = 1,
                   U = 115

)


kernel <- function(pars){

  n   <- pars$mat_siz
  L   <- pars$L
  U   <- pars$U
  #these are the upper and lower integration limits
  h   <- (U-L)/n                   #Bin size
  b   <- L+c(0:n)*h                #Lower boundaries of bins
  y   <- 0.5*(b[1:n]+b[2:(n+1)])   #Bins' midpoints


  # Fertility matrix
  Fmat <- matrix(0,(n),(n))
  Fmat <- outer(y, y, fec, pars) * h


  # Survival matrix
  smat <- matrix(0,(n), (n))
  smat <- s_z(y, pars)

  # Growth matrix
  gmat <- matrix(0, (n), (n))
  gmat <- outer(y,y, g_z1z, pars) * h

  # Growth/survival matrix
  Tmat <- gmat # place holder

  # fix eviction of offspring
  for(i in 1:(pars$mat_siz/2)) {
    gmat[1,i]<- gmat[1,i]+1-sum(gmat[,i])
    Tmat[,i]<-gmat[,i]*smat[i]
  }
  # fix eviction of large adults
  for(i in (pars$mat_siz/2+1):pars$mat_siz) {
    gmat[n,i]<-gmat[n,i]+1-sum(gmat[,i])
    Tmat[,i]<-gmat[,i]*smat[i]
  }
  
  # Full kernel
  k_yx <- Fmat + Tmat

  return(list(k_yx = k_yx,
              smat = smat,
              gmat = gmat,
              Tmat = Tmat,
              Fmat = Fmat))

}

k <- kernel(pars)
lambda <- Re(eigen(k$k_yx)$values[1])
elast <- matrix(0, pars$mat_siz, pars$mat_siz)
pert <- 0.05

for (r in c(1:pars$mat_siz)) {
  for (c in c(1:pars$mat_siz)) {
    a <- k$k_yx[r,c]
    da <- a + pert
    dmat <- k$k_yx
    dmat[r,c] <- da
    dlambda <- Re(eigen(dmat)$values[1])
    elast[r,c] <- (a/lambda) * ((dlambda - lambda) / pert)
  }

}

save(k, elast, file = "results/simulations/elasticity_mpm.RData")

load("results/simulations/elasticity_mpm.RData")

plot <- elast %>% 
  as_tibble %>% 
  tibble::rowid_to_column(var = "X") %>% 
  gather(key = "Y", value = "elasticity", -1) %>%
  mutate(Y = as.numeric(gsub("V", "", Y))) %>%
  ggplot(aes(X, Y, fill = elasticity)) + geom_tile() + 
  theme_minimal() + coord_flip() + scale_x_reverse() + 
  scale_fill_gradient('elasticity', limits=c(min(elast), max(elast)), 
                      breaks = round(seq(0,0.07, length.out = 7), digits = 3),  high = "blue", low = "white")

ggplot2::ggsave(filename = "results/mpm_elasticity.png", plot)

