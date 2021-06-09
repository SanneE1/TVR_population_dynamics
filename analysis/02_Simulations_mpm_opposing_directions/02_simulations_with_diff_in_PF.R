### Similar simulations as in 01 folder, but now vital rates respond in different direction to climate signals
rm(list=ls())
library(dplyr)
library(tidyr)
library(popbio)
library(pbapply)
library(parallel)
library(ggplot2)
library(faux)
library(boot)
library(patchwork)

set.seed(2)

output_dir <- "results/02_Simulations_mpm_opposing_directions/pos_P_neg_F/"

#Set signal strength to 0.5
i = 0.5

### Create line sourcing to 
source_lines <- function(file, lines){
  source(textConnection(readLines(file)[lines]))
}

### Get necesary functions. Source lines from 01 folder
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(14:35, 107:134))

### Create mpm function with +P (s&g) and - F
mpm <- function(survival, growth, reproduction, 
                clim_sd, sig.strength) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)
  
  # growth 
  g_mean = 0.325
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
  s_mean = 0.541
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
    dev <- (-reproduction) * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    mpm[1,2] <- qgamma(p, (f_mean^2)/(f_sd * clim_sd)^2, (f_mean)/(f_sd * clim_sd)^2)
  }
  
  return(mpm)  
}

## Run simulations ----------------------------------------

# Set up parallel
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(139:152))

# Create climate sequences
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(196:200))

# Run simulations for U and F matrix
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(238:277))


# Plot results ----------------------------------------
# lag_fp <- readRDS(file.path(output_dir, paste0("mpm_", i, "_lagfp.RDS")))

lag_pf <- lapply(lag_fp, function(x) 
  lapply(x, function(y) y$df) %>% bind_rows  %>% 
    mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
  bind_rows(., .id = "type")


lagpf_p <- ggplot(lag_pf %>% filter(type != "Fmatrix")) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  labs(colour = "Lag type", title = "Lagged climate in U or F matrix", 
       subtitle = "U climate effect = pos \nF climate effect = neg") +
  facet_grid(cols = vars(auto_cat)) + scale_colour_manual(values = c("#E7B800", "#FC4E07"))



## Get plot of original (both positive) simulation
orig_dir <- "results/01_Simulations_mpm_same_directions/"

### Load results
pospos <- file.path(orig_dir, "rds")

lagpf_pos <- lapply(list.files(pospos, pattern = "lagfp.RDS", full.names = T), readRDS)
sig.strength <- regmatches(list.files(pospos, pattern = "auto.RDS"), 
                           regexec("mpm\\_(.+)\\_", list.files(pospos, pattern = "auto.RDS"))) %>% 
  lapply(., function(x) x[2])
lagpf_pos <- lagpf_pos[[which(sig.strength == i)]]

lagpf_pos <- lapply(lagpf_pos, function(x) lapply(x, function(y) y$df) %>% bind_rows %>% 
                    mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
  bind_rows(., .id = "type") %>% filter(type != "Fmatrix") %>%
  ggplot(.) +
  geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  labs(colour = "Lag type", title = "Lagged climate in U or F matrix",
       subtitle = "U climate effect = pos \nF climate effect= pos") +
  facet_grid(cols = vars(auto_cat)) + scale_colour_manual(values = c("#E7B800", "#FC4E07")) 

ggsave(filename = file.path(output_dir, "lag_plot_PF.png"), 
       (lagpf_p / lagpf_pos) * coord_cartesian(ylim = c(-0.55, -0.15)) + plot_layout(guides = "collect"))
