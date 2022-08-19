# Script to produce figures 1 & 2 of the manuscript

rm(list=ls())

set.seed(2)

library(dplyr)
library(tidyr)
library(popbio)
library(parallel)
library(faux)
library(boot)
library(ggplot2)
library(patchwork)

### Create line sourcing to
source_lines <- function(file, lines){
  source(textConnection(readLines(file)[lines]))
}

# --------------------------------------------
# Plot simulation lambda's
# --------------------------------------------

label_auto <- c(
  "0.9" = "0.9 autocorrelation",
  "0" = "0 autocorrelation",
  "-0.9" = "-0.9 autocorrelation"
)

for(i in c(1, 0.5, 0.25, 0.05)){
  
  lag_same <- readRDS(file.path("results/01_Simulations_mpm_same_directions/", 
                                paste("mpm_", i, "_lagfp.RDS", sep = ""))) %>%
    lapply(., function(x) lapply(x, function(y) y$df) %>% bind_rows %>% 
             mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))
    ) %>%
    bind_rows(., .id = "type") %>% mutate(type = factor(type, levels = c("none", "Umatrix", "Fmatrix")))
  
  lag_opp <- readRDS(file.path("results/02_Simulations_mpm_opposing_directions/", 
                               paste("mpm_", i, "_lagfp.RDS", sep = ""))) %>%
    lapply(., function(x) lapply(x, function(y) y$df) %>% bind_rows %>% 
             mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))
    ) %>%
    bind_rows(., .id = "type") %>% mutate(type = factor(type, levels = c("none", "Umatrix", "Fmatrix")))
  
  
  plot_same <- lag_same %>% filter(type != "Fmatrix") %>%
    ggplot(.) + 
    geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)), se = F)+ 
    geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(type)), position = position_dodge(width = 0.1))+ 
    labs(title = "Same directional response of submatrices to climate driver") +
    ylab("stochastic log lambda") + xlab(~ paste(sigma[c], " of environmental sequence")) + 
    ylim(-0.6,-0.1)+
    facet_grid(cols = vars(auto_cat), labeller = labeller(auto_cat = label_auto)) + 
    scale_colour_manual(name = "Simulation type",
                        values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                        labels = c("none" = "control", "Umatrix" = "MCD")) + theme_minimal() +
    theme(legend.position = "bottom")
  
  plot_opp <- lag_opp %>% filter(type != "Fmatrix") %>%
    ggplot(.) + 
    geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)), se = F)+ 
    geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(type)), position = position_dodge(width = 0.1))+ 
    labs(title = "Opposing directional response of submatrices to climate driver") +
    ylab("stochastic log lambda") + xlab(~ paste(sigma[c], " of environmental sequence")) + 
    ylim(-0.6,-0.1)+
    facet_grid(cols = vars(auto_cat), labeller = labeller(auto_cat = label_auto)) + 
    scale_colour_manual(name = "Simulation type",
                        values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                        labels = c("none" = "control", "Umatrix" = "MCD")) + theme_minimal() +
    theme(legend.position = "bottom")
  
  plot <- (plot_same / plot_opp ) + 
    plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & 
    theme(legend.position = "bottom")
  
    png(filename = file.path("results", paste0("Comparison_plot_", i, ".png")),
       width = 6, height = 7, type = "cairo", units = "in", res = 600)
    print(plot)
    dev.off()
  
  
}

rm(lag_same, lag_opp, plot_same, plot_opp, plot)

plot_UFc <- readRDS("results/01_Simulations_mpm_same_directions/mpm_0.5_lagfp.RDS") %>%
  lapply(., function(x) lapply(x, function(y) y$df) %>% bind_rows %>% 
           mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))
  ) %>%
  bind_rows(., .id = "type") %>% mutate(type = factor(type, levels = c("none", "Umatrix", "Fmatrix"))) %>% 
  ggplot(.) + 
  geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)), se = F)+ 
  geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(type)), position = position_dodge(width = 0.1))+ 
  labs(title = "Same directional response of submatrices to climate driver") +
  ylab("stochastic log lambda") + xlab(~ paste(sigma[c], " of environmental sequence")) + 
  ylim(-0.6,-0.1)+
  facet_grid(cols = vars(auto_cat), labeller = labeller(auto_cat = label_auto)) + 
  scale_colour_manual(name = "Simulation type",
                      values = c("none" = "#0072B2", "Umatrix" = "#E69F00", "Fmatrix" = "#CC79A7"),
                      labels = c("none" = "control", "Umatrix" = "U matrix lagged", "Fmatrix" = "F matrix lagged")) + 
  theme_minimal() +
  theme(legend.position = "bottom")

png(filename = file.path("results", "01_Simulations_mpm_same_directions","UFc_Comparison_plot_0.5.png"),
    width = 6, height = 4, type = "cairo", units = "in", res = 600)
print(plot_UFc)
dev.off()


##  -------------------------------------
##  Plot a lambda time sequence of simulation with i = 0.5 (50% of temporal variance explained by climate driver)
##  -------------------------------------
source_lines("analysis/01_simulations_same_directional_responses.R", c(43:120))

st.lamb <- function(env_surv, env_growth, env_reproduction, 
                    clim_sd = 2, clim_auto, sig.strength = 0.5) {
  
  n_it = length(env_surv)
  
  env <- list(survival = env_surv,
              growth = env_growth,
              reproduction = env_reproduction,
              clim_sd = clim_sd,
              sig.strength = sig.strength)
  
  ### Get all mpm's
  mats <- purrr::pmap(env, mpm) %>% Filter(Negate(anyNA), .)
  mats <- mats %>% purrr::discard(function(x) all(x == 0))
  
  df <- data.frame(lambda = stoch.growth.rate(mats, maxt = n_it, verbose = F)$sim,
                   clim_sd = clim_sd,
                   clim_auto = clim_auto)
  
  return(list(df = df,
              mats = sapply(mats, as.vector) %>% t %>% 
  `colnames<-`(c("1,1", "2,1", "1,2", "2,2"))))
}

seq0 <- create_seq(n_it = 50000, clim_sd = 2, clim_auto = 0, lag = 1)
seqneg <- create_seq(n_it = 50000, clim_sd = 2, clim_auto = -0.9, lag = 1)
seqpos <- create_seq(n_it = 50000, clim_sd = 2, clim_auto = 0.9, lag = 1)

U0 <- st.lamb(env_surv = seq0$lagged,
              env_growth = seq0$lagged,
              env_reproduction = seq0$recent,
              clim_sd = 2, clim_auto = 0, sig.strength = 0.5)$mats
c0 <- st.lamb(env_surv = seq0$recent,
              env_growth = seq0$recent,
              env_reproduction = seq0$recent,
              clim_sd = 2, clim_auto = 0, sig.strength = 0.5)$mats

Uneg <- st.lamb(env_surv = seqneg$lagged,
              env_growth = seqneg$lagged,
              env_reproduction = seqneg$recent,
              clim_sd = 2, clim_auto = -0.9, sig.strength = 0.5)$mats
cneg <- st.lamb(env_surv = seqneg$recent,
              env_growth = seqneg$recent,
              env_reproduction = seqneg$recent,
              clim_sd = 2, clim_auto = -0.9, sig.strength = 0.5)$mats

Upos <- st.lamb(env_surv = seqpos$lagged,
              env_growth = seqpos$lagged,
              env_reproduction = seqpos$recent,
              clim_sd = 2, clim_auto = 0.9, sig.strength = 0.5)$mats
cpos <- st.lamb(env_surv = seqpos$recent,
              env_growth = seqpos$recent,
              env_reproduction = seqpos$recent,
              clim_sd = 2, clim_auto = 0.9, sig.strength = 0.5)$mats

l_seq0 <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(U0))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                      none = lapply(as.list(as.data.frame(t(c0))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                      time = c(1:50002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

l_seqneg <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(Uneg))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        none = lapply(as.list(as.data.frame(t(cneg))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        time = c(1:50002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

l_seqpos <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(Upos))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        none = lapply(as.list(as.data.frame(t(cpos))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        time = c(1:50002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

### zoom in on 35 years
time_series0 <- ggplot(l_seq0) + 
  geom_line(aes(y = lambda, x = time, colour = type), show.legend = F) + 
  xlim(1095,1130) + ylim(c(-3,1.5)) + ylab("log lambda") + 
  labs(subtitle = "0 autocorrelation") + 
  scale_colour_manual(values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      labels = c("none" = "control", "Umatrix" = "MCD"))


time_seriesneg <- ggplot(l_seqneg) + 
  geom_line(aes(y = lambda, x = time, colour = type), show.legend = F) + 
  xlim(1095,1130) + ylim(c(-3,1.5)) + ylab("log lambda") + 
  labs(subtitle = "-0.9 autocorrelation") + 
  scale_colour_manual(values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      label = c("Umatrix" = "MCD", "none" = "control"))

time_seriespos <- ggplot(l_seqpos) + 
  geom_line(aes(y = lambda, x = time, colour = type), show.legend = F) + 
  xlim(1095,1130) + ylim(c(-3,1.5)) + ylab("log lambda") + 
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
  xlab("interannual difference") + coord_cartesian(xlim = c(-3, 3))

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
  xlab("interannual difference") + coord_cartesian(c(-3, 3))

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
  xlab("interannual difference") + coord_cartesian(c(-3, 3))



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
