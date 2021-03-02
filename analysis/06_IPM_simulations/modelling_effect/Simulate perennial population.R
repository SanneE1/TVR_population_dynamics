set.seed(13)
### Script to simulate a population based on "true" parameters and climate effects

library(dplyr)
library(tidyr)
library(ggplot2)
library(ipmr)
library(pbapply)
library(cmdstanr)
library(parallel)

source("analysis/simulations/modelling_effect/population_simulation_functions.R")

clim_sd <- rep(c(0.1, 0.75, 1.5), 5)
clim_cor <- rep(seq(from = -0.9, to = 0.9, length.out = 5), each = 3)


cores <- detectCores()
cl <- makeCluster(cores-2)
clusterExport(cl = cl, ls())
clusterEvalQ(cl, c(library(cmdstanr), library(ipmr), library(dplyr), library(tidyr)))

start <- Sys.time()
start

df <- pblapply(cl = cl, as.list(c(1:15)), function(n) wrapper(sample = 35,   ## Number of consecutive years to follow
                                                              clim_corr = clim_cor[n],
                                                              clim_sd = clim_sd[n]) 
)

end <- Sys.time()
end - start
stopCluster(cl)

saveRDS(df, "results/simulations/Simulated_pop_3.rds")
df <- readRDS("results/simulations/Simulated_pop_1.rds") %>% rbind(.,readRDS("results/simulations/Simulated_pop_2.rds") )

end <- Sys.time()
end


df1 <- df %>% #bind_rows %>%
  unnest(c(`list(lagged_lambdas)`, `list(recent_lambdas)`)) %>%
  rename(lagged_lambdas = `list(lagged_lambdas)`,
         recent_lambdas = `list(recent_lambdas)`) %>% 
  select(-c(pop.sizes, n_seeds.t, n_recr.t, z_new)) %>%
  pivot_longer(contains("lambdas"), names_to = "type", values_to = "estimated_lambda") %>%
  mutate(diff_lambda = estimated_lambda - actual_lambda)
  

df2 <- df %>% bind_rows %>%
  pivot_longer(contains("model"), names_to = "type", values_to = "estimate_coefficient")


ggplot(df1) + geom_boxplot(aes(x = as.factor(clim_sd), y = estimated_lambda, colour = type), position = position_dodge(width = 1)) +
  ylab("Estimated Log Lambda") + xlab("Climate autocorrelation") + facet_grid(cols = vars(clim_corr))

ggplot(df1) + geom_boxplot(aes(x = as.factor(clim_sd), y = diff_lambda, colour = type), position = position_dodge(width = 1)) +
  ylab("Difference in lambdas \nestimated and moddeled") + xlab("Climate autocorrelation") + facet_grid(cols = vars(clim_corr))



ggplot(df2) +
  geom_smooth(aes(x = clim_sd, y = estimate_coefficient, colour = type)) +
  facet_grid(cols = vars(clim_corr)) +
  geom_hline(aes(yintercept = -0.066)) +
  ylab("Estimated climate coefficient") + xlab("Climate standard deviation")

