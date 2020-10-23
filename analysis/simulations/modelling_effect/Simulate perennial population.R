set.seed(13)
### Script to simulate a population based on "true" parameters and climate effects

library(dplyr)
library(tidyr)
library(ggplot2)
library(ipmr)
library(rethinking)
library(pbapply)


source("analysis/simulations/modelling_effect/population_simulation_functions.R")

clim_sd <- rep(c(0.1, 0.5, 1, 1.5, 2), 5)
clim_cor <- rep(seq(from = -0.9, to = 0.9, length.out = 5), each = 5)


cores <- detectCores()
cl <- makeCluster(cores-2)
clusterExport(cl = cl, ls())
clusterEvalQ(cl, c(library(rethinking), library(ipmr), library(dplyr), library(tidyr)))

start <- Sys.time()
start

df <- pblapply(cl = cl, as.list(c(1:5)), function(n) wrapper(init.pop.size = 1000,   ## tryCatch allows lapply to keep going in case of an error (population/model crashes)
                                                    n_yrs = 35,
                                                    clim_corr = clim_cor[n],
                                                    clim_sd = clim_sd[n]) 
)

stopCluster(cl)

saveRDS("results/simulations/Simulated_pop.rds")
df <- readRDS("results/simulations/Simulated_pop.rds")

# end <- Sys.time()
end


df1 <- df %>% bind_rows %>%
  unnest(c(`list(lagged_lambdas)`, `list(recent_lambdas)`)) %>%
  rename(lagged_lambdas = `list(lagged_lambdas)`,
         recent_lambdas = `list(recent_lambdas)`) %>% 
  select(clim_corr, clim_sd, lagged_lambdas, recent_lambdas) %>%
  pivot_longer(contains("lambdas"), names_to = "type", values_to = "estimated_lambda")

df2 <- df %>% bind_rows %>%
  pivot_longer(contains("estimate"), names_to = "type", values_to = "estimate_coefficient")

ggplot(df1) + geom_smooth(aes(x = clim_corr, y = estimated_lambda, colour = type)) 

ggplot(df2) +
  geom_smooth(aes(x = clim_sd, y = estimate_coefficient, colour = type)) +
  facet_grid(cols = vars(clim_corr)) +
  geom_hline(aes(yintercept = -0.066))

