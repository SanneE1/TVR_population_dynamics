### Script to simulate a population based on "true" parameters and climate effects

library(dplyr)
library(tidyr)
library(ggplot2)


source("analysis/simulations/ipmr_functions.R")
source("analysis/simulations/modelling effect/population_simulation_functions.R")

clim_sd <- rep(seq(from = 0.1, to = 2, length.out = 5), 50)
clim_cor <- rep(seq(from = -0.9, to = 0.9, length.out = 5), each = 50)

df <- lapply(as.list(c(1:250)), function(n) tryCatch(simulate(init.pop.size = 1000,   ## tryCatch allows lapply to keep going in case of an error (population/model crashes)
                                                    n_yrs = 35,
                                                    clim_corr = clim_cor[n],
                                                    clim_sd = clim_sd[n]), error = function(e) NULL) 
)

df <- df %>% bind_rows %>%
  pivot_longer(cols = contains("estimate"), names_to = "type", values_to = "estimate")

ggplot(df) + 
  geom_smooth(aes(x = clim_sd, y = estimate, colour = type)) + 
  facet_grid(cols = vars(clim_corr)) + 
  geom_hline(aes(yintercept = -0.066))

