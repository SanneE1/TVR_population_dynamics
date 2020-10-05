
library(ipmr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)


df <- readRDS("results/simulations/P_1yr_HEQU_original_mesh_100.RDS") %>%
  rename(clim_corr = `clim_corr[taskID]`,
         clim_sd = `clim_sd[taskID]`)

a <- df %>% pivot_longer(cols = contains("lambda"), names_to = "type", values_to = "lambda") %>%
  select(clim_corr, clim_sd, type, lambda)
a$type <- gsub("_lambda", "", a$type)

b <- df %>% pivot_longer(cols = contains("corr_sum"), names_to = "type", values_to = "cell_corr_sum") %>%
  select(clim_corr, clim_sd, type, cell_corr_sum)
b$type <- gsub("_corr_sum", "", b$type)

c <- df %>% pivot_longer(cols = contains("corr_sd"), names_to = "type", values_to = "cell_corr_sd") %>%
  select(clim_corr, clim_sd, type, cell_corr_sd)
c$type <- gsub("_corr_sd", "", c$type)

data <- cbind(a, b$cell_corr_sum, c$cell_corr_sd) %>%
  rename(cell_corr_sum = `b$cell_corr_sum`,
         cell_corr_sd = `c$cell_corr_sd`)
data$type <- factor(data$type, levels = c("g", "s", "n"))



l <- ggplot(data)+ 
  geom_smooth(aes(x = clim_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
 facet_grid(rows = vars(clim_corr)) +
  ylim(c(-0.3,0.02))

cor <- ggplot(data) +
  geom_smooth(aes(x = clim_sd,
                  y = cell_corr_sum,
                  color = type),
              size = 1) +
  facet_grid(rows = vars(clim_corr)) +
  ylim(c(0.05, 0.4))

sd <- ggplot(data) +
  geom_smooth(aes(x = clim_sd,
                  y = cell_corr_sd,
                  color = type),
              size = 1) +
  facet_grid(rows = vars(clim_corr)) + 
  ylim(c(0.25,0.5))



####  Growth climate effect 2x as big as original


df1 <- readRDS("results/simulations/P_1yr_HEQU_growth2_mesh_100.RDS") %>%
  rename(clim_corr = `clim_corr[taskID]`,
         clim_sd = `clim_sd[taskID]`)

a1 <- df1 %>% pivot_longer(cols = contains("lambda"), names_to = "type", values_to = "lambda") %>%
  select(clim_corr, clim_sd, type, lambda)
a1$type <- gsub("_lambda", "", a1$type)

b1 <- df1 %>% pivot_longer(cols = contains("corr_sum"), names_to = "type", values_to = "cell_corr_sum") %>%
  select(clim_corr, clim_sd, type, cell_corr_sum)
b1$type <- gsub("_corr_sum", "", b1$type)

c1 <- df1 %>% pivot_longer(cols = contains("corr_sd"), names_to = "type", values_to = "cell_corr_sd") %>%
  select(clim_corr, clim_sd, type, cell_corr_sd)
c1$type <- gsub("_corr_sd", "", c1$type)

data1 <- cbind(a1, b1$cell_corr_sum, c1$cell_corr_sd) %>%
  rename(cell_corr_sum = `b1$cell_corr_sum`,
         cell_corr_sd = `c1$cell_corr_sd`)
data1$type <- factor(data1$type, levels = c("g", "s", "n"))



l1 <- ggplot(data1)+ 
  geom_smooth(aes(x = clim_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
  facet_grid(rows = vars(clim_corr)) +
  ylim(c(-0.3,0.02))

cor1 <- ggplot(data1) +
  geom_smooth(aes(x = clim_sd,
                  y = cell_corr_sum,
                  color = type),
              size = 1) +
  facet_grid(rows = vars(clim_corr)) +
  ylim(c(0.05, 0.4))

sd1 <- ggplot(data1) +
  geom_smooth(aes(x = clim_sd,
                  y = cell_corr_sd,
                  color = type),
              size = 1) +
  facet_grid(rows = vars(clim_corr)) + 
  ylim(c(0.25,0.5))

(l + cor + sd) / (l1 + cor1 + sd1)



##################

diff <- df %>%
  mutate(lam_s = s_lambda - n_lambda,
         lam_g = g_lambda - n_lambda,
         corr_sum_s = s_corr_sum - n_corr_sum,
         corr_sum_g = g_corr_sum - n_corr_sum,
         corr_sd_s = s_corr_sd - n_corr_sd,
         corr_sd_g = g_corr_sd - n_corr_sd) %>%
  select(clim_corr, clim_sd,
         lam_s, lam_g,
         corr_sum_s, corr_sum_g,
         corr_sd_s, corr_sd_g) 

diff_lambda <- diff %>%
  pivot_longer(cols = contains("lam_"), names_to = "type", values_to = "lambda_diff")
diff_lambda$type <- gsub("lam_", "", diff_lambda$type)

diff_sum <- diff %>%
  pivot_longer(cols = contains("corr_sum_"), names_to = "type", values_to = "corr_sum_diff")
diff_sum$type <- gsub("corr_sum_", "", diff_sum$type)

diff_sd <- diff %>%
  pivot_longer(cols = contains("corr_sd_"), names_to = "type", values_to = "corr_sd_diff")
diff_sd$type <- gsub("corr_sd_", "", diff_sd$type)


diff_l <- ggplot(diff_lambda)+ 
  geom_smooth(aes(x = clim_sd,
                  y = lambda_diff, 
                  color = type), 
              size = 1) + 
  facet_grid(rows = vars(clim_corr)) +
  ggtitle("difference with non_lagged lambda")

diff_csum <- ggplot(diff_sum)+ 
  geom_smooth(aes(x = clim_sd,
                  y = corr_sum_diff, 
                  color = type), 
              size = 1) + 
  facet_grid(rows = vars(clim_corr)) +
  ggtitle("difference with non_lagged cell correlation sum")

diff_csd <- ggplot(diff_sd)+ 
  geom_smooth(aes(x = clim_sd,
                  y = corr_sd_diff, 
                  color = type), 
              size = 1) + 
  facet_grid(rows = vars(clim_corr))  +
  ggtitle("difference with non_lagged cell correlation sd")

l + diff_l + diff_csum + diff_csd
