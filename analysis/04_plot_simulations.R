# Script to produce figures 1 & 2 of the manuscript

rm(list=ls())

set.seed(2)

library(tidyverse)
# library(tidyr)
library(popbio)
# library(parallel)
# library(faux)
# library(boot)
library(patchwork)
library(ggforce)

# --------------------------------------------
# Plot simulation lambda's
# --------------------------------------------

df <- read.csv(file.path("results/lambdas_life_histories.csv")) %>% 
  mutate(auto_cat = as.factor(auto_cat)) %>%
  pivot_longer(c("Umatrix", "none"), names_to = "lag_type", values_to = "lambda")


for(j in c("positive", "negative")){
  
plot_all <- ggplot(df %>% filter(vr_cov == j)) +
  geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(lag_type), linetype = auto_cat), se = T) +
  geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(lag_type), shape = auto_cat),
             position = position_dodge(width = 0.1)) +
  scale_colour_manual(name = "Simulation type",
                      values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      labels = c("none" = "control", "Umatrix" = "MCD")) +
  scale_shape(name = "Autocorrelation") +
  scale_linetype(name = "Autocorrelation") +
  ylab("stochastic log lambda") + xlab(~ paste(sigma[c], " of environmental sequence")) +   
  theme(legend.position = "bottom")
  

for(i in c(1:49)){
    
    ggsave(filename = file.path("results", "lambda_climSD_plots" , 
                                j, paste0("lambda_climsd_all_ID_sigstrength_", i, ".tiff")), 
                                plot = plot_all +
                                  facet_wrap_paginate(~ sig.strength + lh_id, 
                                                      labeller = labeller(sig.strength = ~ paste("sig.strength:", .),
                                                                          lh_id = ~ paste("life history id: ", .)),
                                                      nrow = 4, ncol = 3, scales = "free", page = i), 
                                width = 190, height = 277, units = "mm" 
    )
    
  }
}

plot_one <- df %>% filter(lh_id == 60) %>% 
  mutate(vr_cov = factor(vr_cov, levels = c("positive", "negative"), labels = c("Positive covariation", "Negative covariation"))) %>%
  ggplot(.) +
  geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(lag_type), shape = auto_cat),
             position = position_dodge(width = 0.15), alpha = 0.1) +
  geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(lag_type), linetype = auto_cat), se = T) +
  scale_colour_manual(name = "Simulation type",
                      values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      labels = c("none" = "control", "Umatrix" = "MCD")) +
  scale_shape(name = "Autocorrelation") +
  scale_linetype(name = "Autocorrelation") +
  ylab("stochastic log lambda") + xlab(~ paste(sigma[c], " of environmental sequence")) +   
  theme(legend.position = "bottom") + 
  facet_row(vars(vr_cov), labeller = )





##  -------------------------------------
##  Plot a lambda time sequence of simulation with i = 0.5 (50% of temporal variance explained by climate driver)
##  -------------------------------------
set.seed(2)

source("R/lambda_functions.R")

i = 0.5

life_history <- 60
lh_df <- read.csv("results/life_histories_df.csv")
lh_df <- lh_df[life_history,]

seq0 <- create_seq(n_it = 10000, clim_sd = 1, clim_auto = 0, lag = 1)
seqneg <- create_seq(n_it = 10000, clim_sd = 1, clim_auto = -0.9, lag = 1)
seqpos <- create_seq(n_it = 10000, clim_sd = 1, clim_auto = 0.9, lag = 1)

U0 <- st.lamb(mpm_df = lh_df,
              env_surv = seq0$lagged,
              env_growth = seq0$lagged,
              env_reproduction = seq0$recent,
              clim_sd = sd(seq0$recent, na.rm = T),
              clim_auto = acf(seq0$recent, plot = F, na.action = na.pass)$acf[2],
              sig.strength = i, 
              return.mpm = T)$mats
c0 <- st.lamb(mpm_df = lh_df,
              env_surv = seq0$recent,
              env_growth = seq0$recent,
              env_reproduction = seq0$recent,
              clim_sd = sd(seq0$recent, na.rm = T),
              clim_auto = acf(seq0$recent, plot = F, na.action = na.pass)$acf[2],
              sig.strength = i, 
              return.mpm = T)$mats

Uneg <- st.lamb(mpm_df = lh_df,
                env_surv = seqneg$lagged,
                env_growth = seqneg$lagged,
                env_reproduction = seqneg$recent,
                clim_sd = sd(seqneg$recent, na.rm = T),
                clim_auto = acf(seqneg$recent, plot = F, na.action = na.pass)$acf[2],
                sig.strength = i, 
                return.mpm = T)$mats
cneg <- st.lamb(mpm_df = lh_df,
                env_surv = seqneg$recent,
                env_growth = seqneg$recent,
                env_reproduction = seqneg$recent,
                clim_sd = sd(seqneg$recent, na.rm = T),
                clim_auto = acf(seqneg$recent, plot = F, na.action = na.pass)$acf[2],
                sig.strength = i, 
                return.mpm = T)$mats

Upos <- st.lamb(mpm_df = lh_df,
                env_surv = seqpos$lagged,
                env_growth = seqpos$lagged,
                env_reproduction = seqpos$recent,
                clim_sd = sd(seqpos$recent, na.rm = T),
                clim_auto = acf(seqpos$recent, plot = F, na.action = na.pass)$acf[2],
                sig.strength = i, 
                return.mpm = T)$mats
cpos <- st.lamb(mpm_df = lh_df,
                env_surv = seqpos$recent,
                env_growth = seqpos$recent,
                env_reproduction = seqpos$recent,
                clim_sd = sd(seqpos$recent, na.rm = T),
                clim_auto = acf(seqpos$recent, plot = F, na.action = na.pass)$acf[2],
                sig.strength = i, 
                return.mpm = T)$mats

l_seq0 <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(U0))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                      none = lapply(as.list(as.data.frame(t(c0))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                      time = c(1:10002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

l_seqneg <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(Uneg))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        none = lapply(as.list(as.data.frame(t(cneg))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        time = c(1:10002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

l_seqpos <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(Upos))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        none = lapply(as.list(as.data.frame(t(cpos))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        time = c(1:10002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

### zoom in on 35 years
time_series0 <- ggplot(l_seq0) +
  geom_line(aes(y = lambda, x = time, colour = type), show.legend = F) +
  xlim(95,130) + ylim(c(-3,1.5)) + ylab("log lambda") +
  labs(subtitle = "0 autocorrelation") +
  scale_colour_manual(values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      labels = c("none" = "control", "Umatrix" = "MCD"))


time_seriesneg <- ggplot(l_seqneg) +
  geom_line(aes(y = lambda, x = time, colour = type), show.legend = F) +
  xlim(95,130) + ylim(c(-3,1.5)) + ylab("log lambda") +
  labs(subtitle = "-0.9 autocorrelation") +
  scale_colour_manual(values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      label = c("Umatrix" = "MCD", "none" = "control"))

time_seriespos <- ggplot(l_seqpos) +
  geom_line(aes(y = lambda, x = time, colour = type), show.legend = F) +
  xlim(95,130) + ylim(c(-3,1.5)) + ylab("log lambda") +
  labs(subtitle = "0.9 autocorrelation") +
  scale_colour_manual(values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      label = c("none" = "control", "Umatrix" = "MCD"))

int_ann_0 <- l_seq0 %>%
  group_by(type) %>%
  arrange(time) %>%
  mutate(diff = (lambda %>% exp) - (lag(lambda, default = first(lambda)) %>% exp),
         type = factor(type, levels = c("Umatrix", "none"))) %>%
  ggplot(.) +
  geom_density(aes(x = diff, fill = type, colour = type), alpha = 0.5) +
  scale_fill_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "MCD"),
                    values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  scale_colour_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "MCD"),
                      values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  xlab("interannual difference") + coord_cartesian(xlim = c(-1, 1))

int_ann_neg <- l_seqneg %>%
  group_by(type) %>%
  arrange(time) %>%
  mutate(diff = (lambda %>% exp) - (lag(lambda, default = first(lambda)) %>% exp),
         type = factor(type, levels = c("Umatrix", "none"))) %>%
  ggplot(.) +
  geom_density(aes(x = diff, fill = type, colour = type), alpha = 0.5, show.legend = F)  +
  scale_fill_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "MCD"),
                    values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  scale_colour_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "MCD"),
                      values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  xlab("interannual difference") + coord_cartesian(c(-1, 1))

int_ann_pos <- l_seqpos %>%
  group_by(type) %>%
  arrange(time) %>%
  mutate(diff = (lambda %>% exp) - (lag(lambda, default = first(lambda)) %>% exp),
         type = factor(type, levels = c("Umatrix", "none"))) %>%
  ggplot(.) +
  geom_density(aes(x = diff, fill = type, colour = type), alpha = 0.5, show.legend = F) +
  scale_fill_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "MCD"),
                    values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  scale_colour_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "MCD"),
                      values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  xlab("interannual difference") + coord_cartesian(c(-1, 1))



a <- wrap_plots((time_seriespos / time_series0 / time_seriesneg), (int_ann_pos / int_ann_0 / int_ann_neg),
                byrow = F,
                ncol = 2,
                widths = c(3,1),
                guides = "collect") & theme(legend.position = "bottom") & theme_minimal()


a[[1]] <- a[[1]] + plot_layout(tag_level = "new")
a[[2]] <- a[[2]] + plot_layout(tag_level = "new")


time_diff <- a + plot_annotation(tag_levels = c('A', 'i'), tag_sep = '.') & theme(legend.position = "bottom")


png(filename = file.path("results", "timeseries&diff.png"), type = "cairo",
    width = 9.22, height = 6.5, units = "in", res = 600)
print(time_diff)
dev.off()
