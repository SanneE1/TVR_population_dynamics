rm(list=ls())
library(patchwork)

#-------------------------------------------
# Source functions from main simulations
#-------------------------------------------

#Set signal strength to 0.5
i = 0.5

### Create line sourcing to 
source_lines <- function(file, lines){
  source(textConnection(readLines(file)[lines]))
}

### Get necesary functions. Source lines from 01 folder
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(2:36, 107:134))

### Change output directory to current folder
output_dir <- "results/05_Simulations_mpm_non_linear_averaging/"


#-----------------------------------
# High survival rate
#-----------------------------------

### Create mpm function with higher end survival/growth probablility (~0.7)
mpm <- function(survival, growth, reproduction, 
                clim_sd, sig.strength) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)
  
  ## survival Juveniles
  # sj_mean = 
  # sj_sd = 
  # if(is.na(survival)) {
  #   
  #   mpm[1,1] <- sj_mean
  #   
  # } else {
  #   
  #   dev <- survival .....  ## total deviation from mean = climate signal & random noise
  #   p <- pnorm(dev, mean = 0, sd = sqrt(clim_sd * 2)) ## Here I use sd = sqrt(2) because dev has variance var(growth) + var(noise). sqrt() to get stand. dev.
  #   q <- qnorm(p)
  #   mpm[1,1] <- qbeta(p, (((sj_mean*(1-sj_mean))/sj_sd^2) - 1) * sj_mean,
  #                     (((sj_mean*(1-sj_mean))/sj_sd^2) - 1) * (1 - sj_mean))
  #   
  # }
  
  # growth 
  g_mean = 0.7
  g_sd = 0.118
  
  if(is.na(growth)) {
    mpm[2,1] <- g_mean  
  } else {
    ## total deviation from mean = climate signal * signal strength & correction factor (partitioning at variance scale) + random noise * signal strength
    dev <- growth * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + 
      rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    ## Because of partitioning and correction factor above, the resulting distribution has a sd of clim_sd
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    mpm[2,1] <- qbeta(p, (((g_mean*(1-g_mean))/(g_sd * clim_sd)^2) - 1) * g_mean,
                      (((g_mean*(1-g_mean))/(g_sd * clim_sd)^2) - 1) * (1 - g_mean))
  }
  
  # survival
  s_mean = 0.7
  s_sd = 0.135
  if(is.na(survival)) {
    mpm[2,2] <- s_mean
  } else {
    dev <- survival * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    mpm[2,2] <- qbeta(p, (((s_mean*(1-s_mean))/(s_sd * clim_sd)^2) - 1) * s_mean,
                      (((s_mean*(1-s_mean))/(s_sd * clim_sd)^2) - 1) * (1 - s_mean))
  }
  
  # reproduction 
  f_mean = 0.788
  f_sd = 0.862
  if(is.na(reproduction)) {
    mpm[1,2] <- f_mean
  } else {
    dev <- reproduction * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    mpm[1,2] <- qgamma(p, (f_mean^2)/(f_sd * clim_sd)^2, (f_mean)/(f_sd * clim_sd)^2)
  }
  
  return(mpm)  
}


# Set up parallel & climate sequences
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(139:152, 196:200))

# Run simulation within U and between U&F matrices
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(200:235, 238:274))

saveRDS(lag, file.path(output_dir, "high_survival_lag.RDS"))
saveRDS(lag_fp, file.path(output_dir, "high_survival_lagfp.RDS"))

stopCluster(cl)


#-----------------------------------
# Low survival rate
#-----------------------------------

### Create mpm function with lower end survival probablility (~0.23)
mpm <- function(survival, growth, reproduction, 
                clim_sd, sig.strength) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)
  
  ## survival Juveniles
  # sj_mean = 
  # sj_sd = 
  # if(is.na(survival)) {
  #   
  #   mpm[1,1] <- sj_mean
  #   
  # } else {
  #   
  #   dev <- survival .....  ## total deviation from mean = climate signal & random noise
  #   p <- pnorm(dev, mean = 0, sd = sqrt(clim_sd * 2)) ## Here I use sd = sqrt(2) because dev has variance var(growth) + var(noise). sqrt() to get stand. dev.
  #   q <- qnorm(p)
  #   mpm[1,1] <- qbeta(p, (((sj_mean*(1-sj_mean))/sj_sd^2) - 1) * sj_mean,
  #                     (((sj_mean*(1-sj_mean))/sj_sd^2) - 1) * (1 - sj_mean))
  #   
  # }
  
  # growth 
  g_mean = 0.2
  g_sd = 0.118
  
  if(is.na(growth)) {
    mpm[2,1] <- g_mean  
  } else {
    ## total deviation from mean = climate signal * signal strength & correction factor (partitioning at variance scale) + random noise * signal strength
    dev <- growth * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + 
      rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    ## Because of partitioning and correction factor above, the resulting distribution has a sd of clim_sd
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    mpm[2,1] <- qbeta(p, (((g_mean*(1-g_mean))/(g_sd * clim_sd)^2) - 1) * g_mean,
                      (((g_mean*(1-g_mean))/(g_sd * clim_sd)^2) - 1) * (1 - g_mean))
  }
  
  # survival
  s_mean = 0.2
  s_sd = 0.135
  if(is.na(survival)) {
    mpm[2,2] <- s_mean
  } else {
    dev <- survival * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    mpm[2,2] <- qbeta(p, (((s_mean*(1-s_mean))/(s_sd * clim_sd)^2) - 1) * s_mean,
                      (((s_mean*(1-s_mean))/(s_sd * clim_sd)^2) - 1) * (1 - s_mean))
  }
  
  # reproduction 
  f_mean = 0.788
  f_sd = 0.862
  if(is.na(reproduction)) {
    mpm[1,2] <- f_mean
  } else {
    dev <- reproduction * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    mpm[1,2] <- qgamma(p, (f_mean^2)/(f_sd * clim_sd)^2, (f_mean)/(f_sd * clim_sd)^2)
  }
  
  return(mpm)  
}

# Set up parallel & climate sequences
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(139:152, 196:200))

# Run simulation within U and between U&F matrices
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(200:235, 238:274))

saveRDS(lag, file.path(output_dir, "low_survival_lag.RDS"))
saveRDS(lag_fp, file.path(output_dir, "low_survival_lagfp.RDS"))

stopCluster(cl)


#-------------------------------------------
# Plot results
#-------------------------------------------

### load original simulation (with survival mean at 0.5) 
medium_lag <- lapply(readRDS(file.path("results/01_Simulations_mpm_same_directions/", "rds", "mpm_0.5_lag.RDS")),
                     function(x) lapply(x, function(y) y$df) %>% bind_rows %>%
                       mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
  bind_rows(., .id = "type") %>%
  ggplot(.) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+
  labs(colour = "Lag type", title = "survival&growth mean set to 0.5") +
  facet_grid(cols = vars(auto_cat))

medium_fp <- lapply(readRDS(file.path("results/01_Simulations_mpm_same_directions/", "rds", paste0("mpm_", i, "_lagfp.RDS"))),
                    function(x) lapply(x, function(y) y$df) %>% bind_rows %>% 
                      mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
  bind_rows(., .id = "type") %>%
  ggplot(.) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  labs(colour = "Lag type", title = "survival&growth mean set to 0.5") +
  facet_grid(cols = vars(auto_cat))

high_lag <- lapply(readRDS(file.path(output_dir, "high_survival_lag.RDS")),
                   function(x) lapply(x, function(y) y$df) %>% bind_rows %>%
                     mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
  bind_rows(., .id = "type") %>%
  ggplot(.) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+
  labs(colour = "Lag type", title = "survival&growth mean set to 0.7") +
  facet_grid(cols = vars(auto_cat))

high_fp <- lapply(readRDS(file.path(output_dir, "high_survival_lagfp.RDS")),
                  function(x) lapply(x, function(y) y$df) %>% bind_rows %>% 
                    mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
  bind_rows(., .id = "type") %>%
  ggplot(.) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  labs(colour = "Lag type", title = "survival&growth mean set to 0.7") +
  facet_grid(cols = vars(auto_cat))

low_lag <- lapply(readRDS(file.path(output_dir, "low_survival_lag.RDS")),
                  function(x) lapply(x, function(y) y$df) %>% bind_rows %>%
                    mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
  bind_rows(., .id = "type") %>%
  ggplot(.) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+
  labs(colour = "Lag type", title = "survival&growth mean set to 0.27") +
  facet_grid(cols = vars(auto_cat))

low_fp <- lapply(readRDS(file.path(output_dir, "low_survival_lagfp.RDS")),
                 function(x) lapply(x, function(y) y$df) %>% bind_rows %>% 
                   mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
  bind_rows(., .id = "type") %>%
  ggplot(.) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  labs(colour = "Lag type", title = "survival&growth mean set to 0.27") +
  facet_grid(cols = vars(auto_cat))


lag <- high_lag + medium_lag + low_lag + plot_layout(nrow = 3, guides = "collect") + plot_annotation(title = "Lagged climate in growth or survival at different mean survival/growth pr", 
                                                                                              subtitle = "survival climate effect = pos \ngrowth climate effect= pos")

fp <- high_fp + medium_fp + low_fp + plot_layout(nrow = 3, guides = "collect") + plot_annotation(title = "Lagged climate in U or F matrix at different mean survival/growth pr", 
                                                                                                    subtitle = "survival climate effect = pos \ngrowth climate effect= pos")

ggsave(filename = file.path(output_dir, "lambda_comparison_non_lin_avg_sg.tiff"), plot = lag)
ggsave(filename = file.path(output_dir, "lambda_comparison_non_lin_avg_fp.tiff"), plot = fp)

