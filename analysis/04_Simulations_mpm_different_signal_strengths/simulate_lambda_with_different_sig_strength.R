### simulate the effect of recent and lagged climate drivers when signal strength of different VR depends on the 
### sensitivity of the VR

### Create line sourcing to 
source_lines <- function(file, lines){
  source(textConnection(readLines(file)[lines]))
}

## source necessary function from main file
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(1:33, 51:76))

## overwrite output directory to right folder
output_dir <- "results/04_simulations_mpm_different_signal_strengths/"

### Create MPM function that deals with different signal strenghts
mpm <- function(survival, growth, reproduction, clim_sd, sig.strength = 1) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)
  
  #growth                     Get the error term sd to reflect the sd of the environment sequence
  mpm[2,1] <- inv_logit((growth * sig.strength[1]) + ((1-sig.strength[1]) * rnorm(1, 0, clim_sd)) ) ### inv_logit(0) = 0.5 (intercept)
  
  # survival
  mpm[2,2] <- inv_logit((survival * sig.strength[3]) + ((1-sig.strength[3]) * rnorm(1, 0, clim_sd)) )
  
  # reproduction 
  mpm[1,2] <- exp(1.2 + (reproduction * sig.strength[2]) + ((1-sig.strength[2]) * rnorm(1, 0, clim_sd)) )
  
  return(mpm)  
}


### calculate elasticities
A <- matrix2(c(0,exp(1.2),0.5,0.5), stages = c("Juvenile", "Adult"))
E <- elasticity(A)


### Elasticities of VR need to be negatively correlated with the climate variance in VR (and thus sig.strength)
sig.strength <- 1-as.vector(E)[-1]
sig.strength <- sig.strength/sum(sig.strength)


### Simulate Lambda under differnt SD and autocorrelation
i = sig.strength ## sig.strength needs to be put to i as i was used to loop through sig.strenghts in the original file

source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(82:97, 138:175, 178:211))

### Save output
saveRDS(lag, file.path(output_dir, "mpm_diff_ss_lag.RDS"))
saveRDS(lag_fp, file.path(output_dir, "mpm_diff_ss_lagfp.RDS"))

### Plot lambda's

lag_df <- lapply(lag, function(x) lapply(x, function(y) y$df) %>%
                                           bind_rows) %>%
  bind_rows(., .id = "type") %>%
  mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))

lag_p <- ggplot(lag_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  labs(colour = "Lag type", title = "Lagged climate in growth or survival", 
       subtitle = paste("sig.strength surv:", round(sig.strength[3], digits = 2), "growth:",round(sig.strength[1], digits = 2),
                        "fec:", round(sig.strength[2], digits = 2) )) +
  facet_grid(cols = vars(auto_cat))

lagpf_df <- lapply(lag_fp, function(x) lapply(x, function(y) y$df) %>%
                     bind_rows) %>%
  bind_rows(., .id = "type") %>%
  mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))

lagpf_p <- ggplot(lagpf_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  labs(colour = "Lag type", title = "Lagged climate in P or F", 
       subtitle = paste("sig.strength surv:", round(sig.strength[3], digits = 2), "growth:",round(sig.strength[1], digits = 2),
                        "fec:", round(sig.strength[2], digits = 2) )) +
  facet_grid(cols = vars(auto_cat)) + scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))

ggsave(file.path(output_dir, "lambda_plot_lag.tiff"), plot = lag_p)
ggsave(file.path(output_dir, "lambda_plot_lagfp.tiff"), plot = lagpf_p)
