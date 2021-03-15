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
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(1:33, 51:76))

### Change output directory to current folder
output_dir <- "results/05_Simulations_mpm_non_linear_averaging/"


#-----------------------------------
# High survival rate
#-----------------------------------

### Create mpm function with higher end survival probablility (~0.27)
mpm <- function(survival, growth, reproduction, clim_sd, sig.strength = 1) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)
  
  #growth                     Get the error term sd to reflect the sd of the environment sequence
  mpm[2,1] <- inv_logit(growth * sig.strength + ((1-sig.strength) * rnorm(1, 0, clim_sd)) ) ### inv_logit(0) = 0.5 (intercept)
  
  # survival/stasis
  mpm[2,2] <- inv_logit(1 + survival * sig.strength + ((1-sig.strength) * rnorm(1, 0, clim_sd)) )
  
  # reproduction 
  mpm[1,2] <- exp(1.2 + reproduction * sig.strength + ((1-sig.strength) * rnorm(1, 0, clim_sd)) )
  
  return(mpm)  
}


# Set up parallel
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(81:97))

# Run simulation within P Kernel and between P&F Kernels
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(138:175, 178:211))

saveRDS(lag, file.path(output_dir, "high_survival_lag.RDS"))
saveRDS(lag_fp, file.path(output_dir, "high_survival_lagfp.RDS"))

stopCluster(cl)


#-----------------------------------
# Low survival rate
#-----------------------------------

### Create mpm function with lower end survival probablility (~0.23)
mpm <- function(survival, growth, reproduction, clim_sd, sig.strength = 1) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)
  
  #growth                     Get the error term sd to reflect the sd of the environment sequence
  mpm[2,1] <- inv_logit(growth * sig.strength + ((1-sig.strength) * rnorm(1, 0, clim_sd)) ) ### inv_logit(0) = 0.5 (intercept)
  
  # survival/stasis
  mpm[2,2] <- inv_logit(-1 + survival * sig.strength + ((1-sig.strength) * rnorm(1, 0, clim_sd)) )
  
  # reproduction 
  mpm[1,2] <- exp(1.2 + reproduction  * sig.strength + ((1-sig.strength) * rnorm(1, 0, clim_sd)) )
  
  return(mpm)  
}

# Set up parallel
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(81:97))

# Run simulation within P Kernel and between P&F Kernels
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(138:175, 178:211))

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
  labs(colour = "Lag type", title = "survival mean set to 0.5") +
  facet_grid(cols = vars(auto_cat))  

medium_fp <- lapply(readRDS(file.path("results/01_Simulations_mpm_same_directions/", "rds", "mpm_0.5_lagfp.RDS")),
                    function(x) lapply(x, function(y) y$df) %>% bind_rows %>% 
                      mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
  bind_rows(., .id = "type") %>%
  ggplot(.) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  labs(colour = "Lag type", title = "survival mean set to 0.5") +
  facet_grid(cols = vars(auto_cat))

high_lag <- lapply(readRDS(file.path("results/05_Simulations_mpm_non_linear_averaging/high_survival_lag.RDS")),
                   function(x) lapply(x, function(y) y$df) %>% bind_rows %>% 
                     mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
  bind_rows(., .id = "type") %>%
  ggplot(.) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  labs(colour = "Lag type", title = "survival mean set to 0.73") +
  facet_grid(cols = vars(auto_cat))  
high_fp <- lapply(readRDS(file.path("results/05_Simulations_mpm_non_linear_averaging/high_survival_lagfp.RDS")),
                  function(x) lapply(x, function(y) y$df) %>% bind_rows %>% 
                    mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
  bind_rows(., .id = "type") %>%
  ggplot(.) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  labs(colour = "Lag type", title = "survival mean set to 0.5") +
  facet_grid(cols = vars(auto_cat))

low_lag <- lapply(readRDS(file.path("results/05_Simulations_mpm_non_linear_averaging/low_survival_lag.RDS")),
                  function(x) lapply(x, function(y) y$df) %>% bind_rows %>% 
                    mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
  bind_rows(., .id = "type") %>%
  ggplot(.) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  labs(colour = "Lag type", title = "survival mean set to 0.27") +
  facet_grid(cols = vars(auto_cat))  
low_fp <- lapply(readRDS(file.path("results/05_Simulations_mpm_non_linear_averaging/low_survival_lagfp.RDS")),
                 function(x) lapply(x, function(y) y$df) %>% bind_rows %>% 
                   mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
  bind_rows(., .id = "type") %>%
  ggplot(.) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  labs(colour = "Lag type", title = "survival mean set to 0.5") +
  facet_grid(cols = vars(auto_cat))


lag <- high_lag + medium_lag + low_lag + plot_layout(nrow = 3, guides = "collect") + plot_annotation(title = "Lagged climate in growth or survival at different mean survival pr", 
                                                                                              subtitle = "survival climate effect = pos \ngrowth climate effect= neg")

fp <- high_fp + medium_fp + low_fp + plot_layout(nrow = 3, guides = "collect") + plot_annotation(title = "Lagged climate in P or F kernel at different mean survival pr", 
                                                                                                    subtitle = "survival climate effect = pos \ngrowth climate effect= neg")

ggsave(filename = file.path(output_dir, "lambda_comparison_non_lin_avg_sg.tiff"), plot = lag)
ggsave(filename = file.path(output_dir, "lambda_comparison_non_lin_avg_fp.tiff"), plot = fp)

