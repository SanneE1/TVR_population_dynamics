### Similar simulations as in 01 folder, but now vital rates respond in different direction to climate signals

library(dplyr)
library(tidyr)
library(popbio)
library(pbapply)
library(parallel)
library(ggplot2)
library(faux)
library(patchwork)


output_dir <- "results/02_Simulations_mpm_opposing_directions/pos_s_neg_g/"

#Set signal strength to 0.5
i = 0.5

### Create line sourcing to 
source_lines <- function(file, lines){
  source(textConnection(readLines(file)[lines]))
}

### Get necesary functions. Source lines from 01 folder
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(10:33, 51:76))

### Create mpm function with +s and -g

mpm <- function(survival, growth, reproduction, clim_sd, sig.strength = 1) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)
  
  #growth                     Get the error term sd to reflect the sd of the environment sequence
  mpm[2,1] <- inv_logit(- (growth * sig.strength) + ((1-sig.strength) * rnorm(1, 0, clim_sd)) ) ### inv_logit(0) = 0.5 (intercept)
  
  # survival/stasis
  mpm[2,2] <- inv_logit( (survival * sig.strength) + ((1-sig.strength) * rnorm(1, 0, clim_sd)) )
  
  # reproduction 
  mpm[1,2] <- exp(1.2 + (reproduction * sig.strength) + ((1-sig.strength) * rnorm(1, 0, clim_sd)) )
  
  return(mpm)  
}


## Run simulations ----------------------------------------

# Set up parallel
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(81:97))

# Run covariance simulation and within F Kernel
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(114:177))


# Plot results neg/pos

cov <- lapply(cov, function(x) x$df) %>% bind_rows
lag <- lapply(lag, function(x) lapply(x, function(y) y$df) %>% bind_rows)

cov_df <- cov %>% mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))
cov_p <- ggplot(cov_df) + geom_smooth(aes(x = clim_sd, y = lambda, linetype = as.factor(auto_cat)), colour = "grey30") + 
  labs(linetype = "Climate \nautocorrelation", title = "Covariance in survival & growth", subtitle = "survival climate effect = pos \ngrowth climate effect= neg")

lag_df <- lapply(lag, function(x) x %>% mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
  bind_rows(., .id = "type")

lag_p <- ggplot(lag_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  labs(colour = "Lag type", title = "Lagged climate in growth or survival", subtitle = "survival climate effect = pos \ngrowth climate effect= neg") +
  facet_grid(cols = vars(auto_cat))  


## Load results and plot original pos/pos simulations

orig_dir <- "results/01_Simulations_mpm_same_directions/"

### Load results
pospos <- file.path(orig_dir, "rds")

cov_pos <- lapply(list.files(pospos, pattern = "cov.RDS", full.names = T), readRDS)
lag_pos <- lapply(list.files(pospos, pattern = "lag.RDS", full.names = T), readRDS)

sig.strength <- regmatches(list.files(pospos, pattern = "cov.RDS"), 
                           regexec("mpm\\_(.+)\\_", list.files(pospos, pattern = "auto.RDS"))) %>% 
  lapply(., function(x) x[2])

 cov_pos <- lapply(cov_pos[[which(sig.strength == i)]], function(x) x$df) %>% bind_rows %>% 
   mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9"))) %>%
   ggplot(.) + 
   geom_smooth(aes(x = clim_sd, y = lambda, linetype = as.factor(auto_cat)), colour = "grey30") + 
   labs(linetype = "Climate \nautocorrelation", title = "Covariance in survival & growth", 
        subtitle = "survival climate effect = pos \ngrowth climate effect= pos")
 
 
lag_pos <- lapply(lag_pos[[which(sig.strength == i)]], function(x) lapply(x, function(y) y$df) %>% bind_rows %>% 
               mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
   bind_rows(., .id = "type") %>%
  ggplot(.) +
  geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  labs(colour = "Lag type", title = "Lagged climate in growth or survival",
       subtitle = "survival climate effect = pos \ngrowth climate effect= pos") +
  facet_grid(cols = vars(auto_cat)) 



ggsave(filename = file.path(output_dir, "cov_plot_sg.png"), 
       cov_p + cov_pos + plot_layout(guides = "collect"))
ggsave(filename = file.path(output_dir, "lag_plot_sg.png"), 
       lag_p / lag_pos + plot_layout(guides = "collect"))
