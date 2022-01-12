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
    ylab("stochastic log lambda") + xlab("SD of environmental sequence") + 
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
    ylab("stochastic log lambda") + xlab("SD of environmental sequence") + 
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


##  -------------------------------------
##  Plot a lambda time sequence of simulation with i = 0.5 (50% of temporal variance explained by climate driver)
##  -------------------------------------
lagpf_0.5 <- readRDS(list.files("results/01_Simulations_mpm_same_directions/", pattern = "0.5_lagfp.RDS", full.names = T))


l_seq0 <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(lagpf_0.5$Umatrix[[20]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                      none = lapply(as.list(as.data.frame(t(lagpf_0.5$none[[20]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                      time = c(1:50002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

l_seqneg <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(lagpf_0.5$Umatrix[[10]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        none = lapply(as.list(as.data.frame(t(lagpf_0.5$none[[10]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        time = c(1:50002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

l_seqpos <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(lagpf_0.5$Umatrix[[30]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        none = lapply(as.list(as.data.frame(t(lagpf_0.5$none[[30]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        time = c(1:50002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

### zoom in on 35 years
time_series0 <- ggplot(l_seq0) + 
  geom_line(aes(y = lambda, x = time, colour = type), show.legend = F) + 
  xlim(1095,1130) + ylim(c(-7.5,2)) + ylab("log lambda") + 
  labs(subtitle = "0 autocorrelation") + 
  scale_colour_manual(values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      labels = c("none" = "control", "Umatrix" = "MCD"))


time_seriesneg <- ggplot(l_seqneg) + 
  geom_line(aes(y = lambda, x = time, colour = type), show.legend = F) + 
  xlim(1095,1130) + ylim(c(-7.5,2)) + ylab("log lambda") + 
  labs(subtitle = "-0.9 autocorrelation") + 
  scale_colour_manual(values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      label = c("Umatrix" = "MCD", "none" = "control"))

time_seriespos <- ggplot(l_seqpos) + 
  geom_line(aes(y = lambda, x = time, colour = type), show.legend = F) + 
  xlim(1095,1130) + ylim(c(-7.5,2)) + ylab("log lambda") + 
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
