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
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(1:33, 67:92))

### Change output directory to current folder
output_dir <- "results/05_Simulations_mpm_non_linear_averaging/"


#-----------------------------------
# High survival rate
#-----------------------------------

### Create mpm function with higher end survival probablility (~0.27)
mpm <- function(survival, growth, reproduction, clim_sd, sig.strength) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)
  
  #growth                    
  if(is.na(growth)) {
    mpm[2,1] <- inv_logit(0) ### inv_logit(0) = 0.5 (intercept) 
  } else {
    mpm[2,1] <- inv_logit((growth * (sqrt(clim_sd^2 * sig.strength)/clim_sd)) + 
                            ((sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) * rnorm(1, 0, clim_sd)) 
    ) 
  }
  # survival
  if(is.na(survival)) {
    mpm[2,2] <- inv_logit(0)
  } else {
    mpm[2,2] <- inv_logit((1 + 
                             survival * (sqrt(clim_sd^2 * sig.strength)/clim_sd)) + 
                            ((sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) * rnorm(1, 0, clim_sd)) 
    )
  }
  # reproduction 
  if(is.na(reproduction)) {
    mpm[1,2] <- exp(1.2)
  } else {
    mpm[1,2] <- exp(1.2 + 
                      (reproduction * (sqrt(clim_sd^2 * sig.strength)/clim_sd)) + 
                      ((sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) * rnorm(1, 0, clim_sd)) 
    )
  }
  return(mpm)  
}


# Set up parallel
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(98:113))

# Run simulation within U matrix and between U&F matrices
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(154:191, 194:227))

saveRDS(lag, file.path(output_dir, "high_survival_lag.RDS"))
saveRDS(lag_fp, file.path(output_dir, "high_survival_lagfp.RDS"))

stopCluster(cl)


#-----------------------------------
# Low survival rate
#-----------------------------------

### Create mpm function with lower end survival probablility (~0.23)
mpm <- function(survival, growth, reproduction, clim_sd, sig.strength) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)
  
  #growth                    
  if(is.na(growth)) {
    mpm[2,1] <- inv_logit(0) ### inv_logit(0) = 0.5 (intercept) 
  } else {
    mpm[2,1] <- inv_logit((growth * (sqrt(clim_sd^2 * sig.strength)/clim_sd)) + 
                            ((sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) * rnorm(1, 0, clim_sd)) 
    ) 
  }
  # survival
  if(is.na(survival)) {
    mpm[2,2] <- inv_logit(0)
  } else {
    mpm[2,2] <- inv_logit((-1 +
                             survival * (sqrt(clim_sd^2 * sig.strength)/clim_sd)) + 
                            ((sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) * rnorm(1, 0, clim_sd)) 
    )
  }
  # reproduction 
  if(is.na(reproduction)) {
    mpm[1,2] <- exp(1.2)
  } else {
    mpm[1,2] <- exp(1.2 + 
                      (reproduction * (sqrt(clim_sd^2 * sig.strength)/clim_sd)) + 
                      ((sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) * rnorm(1, 0, clim_sd)) 
    )
  }
  return(mpm)  
}

# Set up parallel
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(98:113))

# Run simulation within U matrix and between U&F matrices
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(154:191, 194:227))

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

medium_fp <- lapply(readRDS(file.path("results/01_Simulations_mpm_same_directions/", "rds", paste0("mpm_", i, "_lagfp.RDS"))),
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

