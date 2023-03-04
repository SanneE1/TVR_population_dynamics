# Script to produce figures 1 & 2 of the manuscript

rm(list=ls())

library(tidyverse)
library(popbio)
library(patchwork)
library(ggforce)
library(ggrepel)
library(grid)

source("R/plot_functions.R")
source("R/lambda_functions.R")

##  -------------------------------------
##  Plot a lambda time sequence of simulation with i = 0.5 (50% of temporal variance explained by climate driver)
##  Figure 2
##  -------------------------------------
set.seed(2)

i = 0.5

life_history <- 60
lh_df <- read.csv("results/life_histories_df.csv")
lh_df <- lh_df[life_history,]

seq0 <- create_seq(n_it = 10000, clim_sd = 0.5, clim_auto = 0, lag = 1)
seqneg <- create_seq(n_it = 10000, clim_sd = 1, clim_auto = -0.6, lag = 1)
seqpos <- create_seq(n_it = 10000, clim_sd = 1, clim_auto = 0.6, lag = 1)

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
                      labels = c("none" = "control", "Umatrix" = "TVR"))


time_seriesneg <- ggplot(l_seqneg) +
  geom_line(aes(y = lambda, x = time, colour = type), show.legend = F) +
  xlim(95,130) + ylim(c(-3,1.5)) + ylab("log lambda") +
  labs(subtitle = "-0.6 autocorrelation") +
  scale_colour_manual(values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      label = c("Umatrix" = "TVR", "none" = "control"))

time_seriespos <- ggplot(l_seqpos) +
  geom_line(aes(y = lambda, x = time, colour = type), show.legend = F) +
  xlim(95,130) + ylim(c(-3,1.5)) + ylab("log lambda") +
  labs(subtitle = "0.6 autocorrelation") +
  scale_colour_manual(values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      label = c("none" = "control", "Umatrix" = "TVR"))

int_ann_0 <- l_seq0 %>%
  group_by(type) %>%
  arrange(time) %>%
  mutate(diff = (lambda %>% exp) - (lag(lambda, default = first(lambda)) %>% exp),
         type = factor(type, levels = c("Umatrix", "none"))) %>%
  ggplot(.) +
  geom_density(aes(x = diff, fill = type, colour = type), alpha = 0.5) +
  scale_fill_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "TVR"),
                    values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  scale_colour_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "TVR"),
                      values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  xlab("interannual difference") + coord_cartesian(xlim = c(-2, 2))

int_ann_neg <- l_seqneg %>%
  group_by(type) %>%
  arrange(time) %>%
  mutate(diff = (lambda %>% exp) - (lag(lambda, default = first(lambda)) %>% exp),
         type = factor(type, levels = c("Umatrix", "none"))) %>%
  ggplot(.) +
  geom_density(aes(x = diff, fill = type, colour = type), alpha = 0.5, show.legend = F)  +
  scale_fill_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "TVR"),
                    values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  scale_colour_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "TVR"),
                      values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  xlab("interannual difference") + coord_cartesian(c(-2, 2))

int_ann_pos <- l_seqpos %>%
  group_by(type) %>%
  arrange(time) %>%
  mutate(diff = (lambda %>% exp) - (lag(lambda, default = first(lambda)) %>% exp),
         type = factor(type, levels = c("Umatrix", "none"))) %>%
  ggplot(.) +
  geom_density(aes(x = diff, fill = type, colour = type), alpha = 0.5, show.legend = F) +
  scale_fill_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "TVR"),
                    values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  scale_colour_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "TVR"),
                      values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  xlab("interannual difference") + coord_cartesian(c(-2, 2))



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





# --------------------------------------------
# Plot simulation lambda's (Figure 3 and appendix)
# --------------------------------------------

df <- read.csv("results/lambdas_life_histories.csv") %>%
  left_join(., 
            read.csv("results/life_histories_df.csv")) %>%
  rename(life.expect_o = life.expect,
         iteroparity_o = iteroparity) %>%
  mutate(auto_cat = cut(clim_auto, breaks = 7, labels = c(-0.6, -0.3, -0.1, 0, 0.1, 0.3, 0.6)),
         life.expect = scale(life.expect_o),
         iteroparity = scale(iteroparity_o)) 

for(j in c("positive", "negative")){
  
  plot_all <- ggplot(df %>% filter(vr_cov == j)) +
    geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(lag_type), linetype = auto_cat), se = T) +
    geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(lag_type), shape = auto_cat),
               position = position_dodge(width = 0.1)) +
    scale_colour_manual(name = "Simulation type",
                        values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                        labels = c("none" = "control", "Umatrix" = "TVR")) +
    scale_shape(name = "Autocorrelation") +
    scale_linetype(name = "Autocorrelation") +
    ylab("stochastic log lambda") + xlab(~ paste(sigma[V], " of environmental sequence")) +
    theme(legend.position = "bottom")
  
  
  for(i in c(1:49)){
    
    ggsave(filename = file.path("results", "lambda_climSD_plots" ,
                                j, paste0("lambda_climsd_all_ID_sigstrength_", i, ".png")),
           plot = plot_all +
             facet_wrap_paginate(~ sig.strength + lh_id,
                                 labeller = labeller(sig.strength = ~ paste("sig.strength:", .),
                                                     lh_id = ~ paste("life history id: ", .)),
                                 nrow = 4, ncol = 3, scales = "free", page = i),
           width = 190, height = 277, units = "mm", dpi = 200,
    )
    
  }
}

plot_one <- plot_w_vrcov_function(id = 60)

ggsave(plot = plot_one, filename = file.path("results", "lh_id_60_lambda_vs_climsd.tiff"),
       width = 7, height = 6, dpi = 400)

## --------------------------------------------
## Plot simulation lambda's for different ends of the life history traits 
## (Supplement)
## --------------------------------------------

lh_df <- read.csv("results/life_histories_df.csv")

### Graphical representation of the lifehistory traits chosen for the next graphs
lh_plot <- ggplot(lh_df, aes(x = log(life.expect), y = iteroparity, label = lh_id)) + 
  geom_point() + 
  geom_text_repel(data = subset(lh_df, lh_id %in% c(60,
                                                    144, 92)),
                  nudge_y = 4.3,
                  box.padding = 1.5,
                  point.padding = 0,
                  direction = "x") +
  geom_text_repel(data = subset(lh_df, lh_id %in% c(60, 44, 111)),
                  nudge_x = -1,
                  box.padding = 0.15,
                  point.padding = 0,
                  direction = "y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, .12))) +
  scale_x_continuous(expand = expansion(mult = c(0.05, .12))) +
  viridis::scale_colour_viridis() + 
  theme_minimal() + theme(legend.position = "bottom")

ggsave("results/life_history_trait_distribution.png", lh_plot,
       width = 5, height = 3)

## graph to compare simulation results for different trait "extremes"
L_low <- plot_one_function(id = 92) + ylab("stochastic log lambda")+ xlab(~ paste(sigma[c], " of environmental sequence"))
L_high <- plot_one_function(id = 144) + ylab("stochastic log lambda")

S_low <- plot_one_function(id = 111) + xlab(~ paste(sigma[c], " of environmental sequence")) 
S_high <- plot_one_function(id = 44)


lifehist_plot <- plot_spacer() + wrap_elements(textGrob("Life Expectancy", gp = gpar(fontsize = 16))) + wrap_elements(textGrob("Degree of iteroparity", gp = gpar(fontsize = 16))) +
  wrap_elements(textGrob("High", rot = 90, gp = gpar(fontsize = 16))) + L_high + S_high +  
  wrap_elements(textGrob("low", rot = 90, gp = gpar(fontsize = 16))) + L_low + S_low +  
  plot_layout(nrow = 3, guides = "collect", 
              widths = c(1,9,9),
              heights = c(2,9,9)) & theme( legend.direction = "vertical")


ggsave(lifehist_plot, filename = "results/lifehistory_comparison_plot.png",
       width = 8, height = 6)


### ---------------------------------------------------------------------------------------------
### Heat map predictions life history space
### Figure 5
### ---------------------------------------------------------------------------------------------

rel_dec <- readRDS("results/models_rel_decr.rds")
rel_dec_pos <- rel_dec[[1]]
rel_dec_neg <- rel_dec[[2]]


lh_df <- read.csv("results/life_histories_df.csv") %>% 
  mutate(life.expect_s = scale(log(life.expect)),
         iteroparity_s = scale(iteroparity))


new_df <- expand.grid(life.expect_real = seq(from = min(log(lh_df$life.expect)), to = max(log(lh_df$life.expect)), length.out = 20),
                      iteroparity_real = seq(from = min(lh_df$iteroparity), to = max(lh_df$iteroparity), length.out = 20),
                      auto_cat = c(-0.6,0,0.6),
                      sig.strength = 0.5,
                      lh_id = 1) %>%
  mutate(life.expect = (life.expect_real - attr(lh_df$life.expect_s, 'scaled:center'))/ attr(lh_df$life.expect_s, 'scaled:scale'),
         iteroparity = (iteroparity_real - attr(lh_df$iteroparity_s, 'scaled:center'))/ attr(lh_df$iteroparity_s, 'scaled:scale'))


### 
new_df$Positive <- predict(rel_dec_pos, newdata = new_df, re.form = NA)
new_df$Negative <- predict(rel_dec_neg, newdata = new_df, re.form = NA)

new_df <- pivot_longer(new_df, cols = c("Positive", "Negative"), names_to = "vr_cor", values_to = "rel_diff")
new_df$vr_cor <- factor(new_df$vr_cor, levels = c("Positive", "Negative"))

heat <- ggplot(new_df) + 
  geom_tile(aes(x = life.expect_real, y = iteroparity_real, fill = rel_diff)) +
  geom_contour(aes(x = life.expect_real, y = iteroparity_real, z = rel_diff),
               binwidth = 0.5, colour = "white") +
  metR::geom_text_contour(aes(x = life.expect_real, y = iteroparity_real, z = rel_diff),
                          binwidth = 0.5, rotate = F) +
  geom_point(data = lh_df, 
             aes(x = log(life.expect), y = iteroparity), alpha = 0.2) +
  facet_grid(rows = vars(vr_cor), cols = vars(auto_cat)) +
  scale_fill_gradient2(low = "#0072B2",
                       mid = "#999999",
                       high = "#E69F00",
                       midpoint = 0, name = "Relative \ndifference",
                       limits = c(-2, 2),
                       oob = scales::squish) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Vital rate correlation", breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Climate autocorrelation levels", breaks = NULL, labels = NULL)) +
  ylab("Degree of iteroparity") + xlab("Ln(Life expectancy)") +
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(colour = "grey20"),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12))

heat

ggsave(plot = heat, filename = "results/lh_heatmap.png",
       width = 8, height = 6, dpi = 400)


### ---------------------------------------------------------------------------------------------
### Plot model coefficients (Figures 4 and 6)
### ---------------------------------------------------------------------------------------------

tb1 <- readRDS("results/models_lambda.rds")
tb2 <- readRDS("results/models_rel_decr.rds")

# a <- confint(tb1[[1]])
# b <- confint(tb1[[2]])
# c <- confint(tb2[[1]])
# d <- confint(tb2[[2]])
# 
# saveRDS(list(a = a, b = b, c = c, d = d), "results/model_confint.rds")

confint <- readRDS("results/model_confint.rds")
a <- confint$a
b <- confint$b
c <- confint$c
d <- confint$d


df_pos <- as.data.frame(a) %>% rownames_to_column(var = "cf")
df_pos$type <- "positive"
df_pos$mean <- (df_pos$`2.5 %` + df_pos$`97.5 %`)/2

df_neg <- as.data.frame(b) %>% rownames_to_column(var = "cf")
df_neg$type <- "negative"
df_neg$mean <- (df_neg$`2.5 %` + df_neg$`97.5 %`)/2

df_coef <- rbind(df_pos, df_neg) %>% filter(!(cf %in% c(".sig01", ".sigma"))) %>%
  mutate(cf_f = factor(cf, levels = c("(Intercept)", "clim_sd", "I(clim_sd^2)", "clim_auto",
                                      "sig.strength", "lag_typeUmatrix", "clim_sd:lag_typeUmatrix", 
                                      "clim_auto:lag_typeUmatrix", "sig.strength:lag_typeUmatrix")
  ),
  type = factor(type, levels = c("positive", "negative")))

full_labels <- rev(expression("Intercept", "Climate std"~(sigma[c]), "Climate std"~(sigma[c]^2),
                          "Climate autocorrelation"~(r[K]), "Climate signal strength (p)",
                          "TVR simulation type", "TVR - "~sigma[c]~" interaction",
                          "TVR - "~r[K]~" interaction", "TVR - p interaction"))
select_labels <- rev(expression("Climate autocorrelation"~(r[K]), "Climate signal strength (p)",
                            "TVR simulation type", "TVR - "~sigma[c]~" interaction",
                            "TVR - "~r[K]~" interaction", "TVR - p interaction"))



coef_plot_full <- ggplot(df_coef) +
  geom_vline(aes(xintercept = 0)) +
  geom_pointrange(aes(x = mean, xmin = `2.5 %`, xmax = `97.5 %`, y = cf_f, colour = type), 
                  position = ggstance::position_dodgev(height = 0.5)) +
  scale_y_discrete(labels = full_labels,
                   limits = rev) +
  scale_colour_manual(values = c("positive" = "#000000","negative" = "#D55E00"),
                      name = "Vital rate \ncorrelation") +
  xlab(expression("Estimated effect on"~ln(lambda[s]))) +
  theme_minimal() +
  theme(axis.title.y = element_blank())


coef_plot <- ggplot(df_coef %>% filter(!(cf %in% c("(Intercept)", "clim_sd", "I(clim_sd^2)")))) +
  geom_vline(aes(xintercept = 0)) +
  geom_pointrange(aes(x = mean, xmin = `2.5 %`, xmax = `97.5 %`, y = cf_f, colour = type), 
                  position = ggstance::position_dodgev(height = 0.5)) +
  scale_y_discrete(labels = select_labels,
                   limits = rev) +
  scale_colour_manual(values = c("positive" = "#000000","negative" = "#D55E00"),
                      name = "Vital rate \ncorrelation") +
  xlab(expression("Estimated effect on"~ln(lambda[s]))) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(hjust = 0.95))

ggsave("results/model_lambda_plot.png", plot = coef_plot_full, width = 7, height = 4)
ggsave("results/model_lambda_plot_selected.png", plot = coef_plot, width = 7, height = 3.5, )




df_pos <- as.data.frame(c) %>% rownames_to_column(var = "cf")
df_pos$type <- "positive"
df_pos$mean <- (df_pos$`2.5 %` + df_pos$`97.5 %`)/2

df_neg <- as.data.frame(d) %>% rownames_to_column(var = "cf")
df_neg$type <- "negative"
df_neg$mean <- (df_neg$`2.5 %` + df_neg$`97.5 %`)/2

df_coef2 <- rbind(df_pos, df_neg) %>% filter(!(cf %in% c(".sig01", ".sigma"))) %>%
  mutate(cf_f = factor(cf, levels = c("(Intercept)", "life.expect", "iteroparity",
                                      "auto_cat", "sig.strength", "life.expect:auto_cat", "iteroparity:auto_cat")
  ),
  type = factor(type, levels = c("positive", "negative")))

full_labels2 <- rev(expression("Intercept", "ln(life expectancy)", "Degree of iteroparity",
                              "Climate autocorrelation"~(r[K]), "Climate signal strength (p)",
                              r[K]~": ln(life expectancy) interaction", r[K]~": Iteroparity interaction"
                              ))


coef_plot_diff <- ggplot(df_coef2) +
  geom_vline(aes(xintercept = 0)) +
  geom_pointrange(aes(x = mean, xmin = `2.5 %`, xmax = `97.5 %`, y = cf_f, colour = type), 
                  position = ggstance::position_dodgev(height = 0.5)) +
  scale_y_discrete(labels = full_labels2,
                   limits = rev) +
  scale_colour_manual(values = c("positive" = "#000000","negative" = "#D55E00"),
                      name = "Vital rate \ncorrelation") +
  xlab(expression("Estimated effect on ln(relative difference in"~lambda[s]~")" )) +
  theme_minimal() +
  theme(axis.title.y = element_blank())


ggsave("results/model_rel_diff.png", plot = coef_plot_diff, width = 7, height = 4)
