library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(patchwork)
library(lme4)
library(ipmr)

source("ipmr_functions.R")

P1 <- as.list(file.path("P-1yr", list.files("P-1yr/")))
P1neg <- as.list(file.path("P-neg-1yr", list.files("P-neg-1yr/")))
P2 <- as.list(file.path("P-2yr", list.files("P-2yr/")))
PF <- as.list(file.path("PF-1yr", list.files("PF-1yr/")))
PF_neg <- as.list(file.path("PF-neg-1yr", list.files("PF-neg-1yr/")))

P1 <- lapply(P1, readRDS) %>% bind_rows
P1neg <- lapply(P1neg, readRDS) %>% bind_rows
P2 <- lapply(P2, readRDS) %>% bind_rows
PF <- lapply(PF, readRDS) %>% bind_rows
PF_neg <- lapply(PF_neg, readRDS) %>% bind_rows


FP <- PF %>% select(temperature_sd,non_lagged, non_lagged_sd, lagged_P.1, lagged_P_sd.1) %>%
  rename(lagged = lagged_P.1,
         lagged_sd = lagged_P_sd.1)

PF <- PF %>% select(temperature_sd, non_lagged, non_lagged_sd, lagged_P, lagged_P_sd) %>%
  rename(lagged = lagged_P,
         lagged_sd = lagged_P_sd)


FP_neg <- PF_neg %>% select(temperature_sd,non_lagged, non_lagged_sd, lagged_F, lagged_F_sd) %>%
  rename(lagged = lagged_F,
         lagged_sd = lagged_F_sd)

PF_neg <- PF_neg %>% select(temperature_sd, non_lagged, non_lagged_sd, lagged_P, lagged_P_sd) %>%
  rename(lagged = lagged_P,
         lagged_sd = lagged_P_sd)


lams <- lapply(list(P1 = P1, 
                    P1neg = P1neg, 
                    P2 = P2, 
                    PF = PF, 
                    FP = FP, 
                    PF_neg = PF_neg, 
                    FP_neg = FP_neg), 
               function(df) df %>% 
                 pivot_longer(cols = c("non_lagged", "lagged"), 
                              names_to = "type", 
                              values_to = "lambda") %>%
                 rename(lagged = lagged_sd,
                        non_lagged = non_lagged_sd) %>%
                 group_by(type) %>%
                 pivot_longer(cols = c("non_lagged", "lagged"),
                              names_to = "type_sd",
                              values_to = "sd") %>%
                 subset(type == type_sd)
               
)


plotP1 <- ggplot(lams$P1) + 
  geom_ribbon(data = lams$P1 %>% group_by(temperature_sd, type) %>% mutate(lambda = mean(lambda),
                                                                       sd = mean(sd)),
              aes(x = temperature_sd,
                  ymin = lambda - sd,
                  ymax = lambda + sd,
                  fill = type),
              alpha = 0.3) + 
  geom_smooth(aes(x = temperature_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
  ylim(c(-1.25, 0.75)) + 
  annotation_custom(tableGrob(data.frame(vitalrate = c("survival", "growth", "fp", "fn"), 
                                         "parameter value" = c(1.233434, -0.06588578, NA, NA))),
                              xmin = 0, xmax = 1, ymin = -1.25, ymax = -0.75) +
  ggtitle("P kernel lagged 1yr, F kernel insensitive")


plotP1neg <- ggplot(lams$P1neg) + 
  geom_ribbon(data = lams$P1neg %>% group_by(temperature_sd, type) %>% mutate(lambda = mean(lambda),
                                                                       sd = mean(sd)),
              aes(x = temperature_sd,
                  ymin = lambda - sd,
                  ymax = lambda + sd,
                  fill = type),
              alpha = 0.3) + 
  geom_smooth(aes(x = temperature_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
  ylim(c(-1.25, 0.75)) + 
  annotation_custom(tableGrob(data.frame(vitalrate = c("survival", "growth", "fp", "fn"), 
                                         "parameter value" = c(-1.233434, -0.06588578, NA, NA))),
                    xmin = 0, xmax = 1, ymin = -1.25, ymax = -0.75)+
  ggtitle("P kernel lagged 1yr, F kernel insensitive")


plotP2 <- ggplot(lams$P2) + 
  geom_ribbon(data = lams$P2 %>% group_by(temperature_sd, type) %>% mutate(lambda = mean(lambda),
                                                                       sd = mean(sd)),
              aes(x = temperature_sd,
                  ymin = lambda - sd,
                  ymax = lambda + sd,
                  fill = type),
              alpha = 0.3) + 
  geom_smooth(aes(x = temperature_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
  ylim(c(-1.25, 0.75)) + 
  annotation_custom(tableGrob(data.frame(vitalrate = c("survival", "growth", "fp", "fn"), 
                                         "parameter value" = c(1.233434, -0.06588578, NA, NA))),
                    xmin = 0, xmax = 1, ymin = -1.25, ymax = -0.75) +
  ggtitle("P kernel lagged 2yr, F kernel insensitive")


plotPF <- ggplot(lams$PF) + 
  geom_ribbon(data = lams$PF %>% group_by(temperature_sd, type) %>% mutate(lambda = mean(lambda),
                                                                       sd = mean(sd)),
              aes(x = temperature_sd,
                  ymin = lambda - sd,
                  ymax = lambda + sd,
                  fill = type),
              alpha = 0.3) + 
  geom_smooth(aes(x = temperature_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
  ylim(c(-1.25, 0.75)) + 
  annotation_custom(tableGrob(data.frame(vitalrate = c("survival", "growth", "fp", "fn"), 
                                         "parameter value" = c(1.233434, -0.06588578, -0.6169492, -0.3398345))),
                    xmin = 0, xmax = 1, ymin = -1.25, ymax = -0.75) +
  ggtitle("Lagged P kernel, F climate sensitive")


plotFP <- ggplot(lams$FP) + 
  geom_ribbon(data = lams$FP %>% group_by(temperature_sd, type) %>% mutate(lambda = mean(lambda),
                                                                       sd = mean(sd)),
              aes(x = temperature_sd,
                  ymin = lambda - sd,
                  ymax = lambda + sd,
                  fill = type),
              alpha = 0.3) + 
  geom_smooth(aes(x = temperature_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
  ylim(c(-1.25, 0.75)) +
  annotation_custom(tableGrob(data.frame(vitalrate = c("survival", "growth", "fp", "fn"), 
                                         "parameter value" = c(1.233434, -0.06588578, -0.6169492, -0.3398345))),
                    xmin = 0, xmax = 1, ymin = -1.25, ymax = -0.75) +
  ggtitle("Lagged F kernel, P climate sensitive")

plotPF_neg <- ggplot(lams$PF_neg) + 
  geom_ribbon(data = lams$PF_neg %>% group_by(temperature_sd, type) %>% mutate(lambda = mean(lambda),
                                                                           sd = mean(sd)),
              aes(x = temperature_sd,
                  ymin = lambda - sd,
                  ymax = lambda + sd,
                  fill = type),
              alpha = 0.3) + 
  geom_smooth(aes(x = temperature_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
  ylim(c(-1.25, 0.75)) + 
  annotation_custom(tableGrob(data.frame(vitalrate = c("survival", "growth", "fp", "fn"), 
                                         "parameter value" = c(-1.233434, -0.06588578, -0.6169492, -0.3398345))),
                    xmin = 0, xmax = 1, ymin = -1.25, ymax = -0.75) +
  ggtitle("Lagged P kernel, F climate sensitive")


plotFP_neg <- ggplot(lams$FP_neg) + 
  geom_ribbon(data = lams$FP_neg %>% group_by(temperature_sd, type) %>% mutate(lambda = mean(lambda),
                                                                           sd = mean(sd)),
              aes(x = temperature_sd,
                  ymin = lambda - sd,
                  ymax = lambda + sd,
                  fill = type),
              alpha = 0.3) + 
  geom_smooth(aes(x = temperature_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
  ylim(c(-1.25, 0.75)) +
  annotation_custom(tableGrob(data.frame(vitalrate = c("survival", "growth", "fp", "fn"), 
                                         "parameter value" = c(-1.233434, -0.06588578, -0.6169492, -0.3398345))),
                    xmin = 0, xmax = 1, ymin = -1.25, ymax = -0.75) +
  ggtitle("Lagged F kernel, P climate sensitive")


(plotP1 + plotP1neg + plotP2) / (plotPF + plotFP + plot_spacer()) / (plotPF_neg + plotFP_neg + plot_spacer())



#### follow population numbers

survC <- readRDS("HEQU_s_climate.rds")
growthC <- readRDS("HEQU_g_climate.rds")

fp <- readRDS("HEQU_fp.rds")
fn <- readRDS("HEQU_fn.rds")
fpC <- readRDS("HEQU_fp_climate.rds")
fnC <- readRDS("HEQU_fn_climate.rds")



## create parameter list -------------------------------------------------------------------------

params_list <- list(
  s_int = fixef(survC)[1],
  s_slope = fixef(survC)[2],
  s_temp = fixef(survC)[5],
  g_int = fixef(growthC)[1],
  g_slope = fixef(growthC)[2],
  g_temp = fixef(growthC)[5],
  g_sd = sd(resid(growthC)),
  fp_int = fixef(fp)[1],
  fp_slope = fixef(fp)[2],
  fpC_int = fixef(fpC)[1],
  fpC_slope = fixef(fpC)[2],
  fpC_temp = fixef(fpC)[5],
  fn_int = fixef(fn)[1],
  fn_slope = fixef(fn)[2],
  fnC_int = fixef(fnC)[1],
  fnC_slope = fixef(fnC)[2],
  fnC_temp = fixef(fnC)[5],
  germ_mean = 0.1262968,
  germ_sd = 0.2725941,
  fd_mean = 1.178749,
  fd_sd = 0.8747764
)



manual_P2 <- P_neg_1yr_man(2) %>% bind_cols
manual_P_2 <- P_neg_1yr_man(-2) %>% bind_cols


Man_pos <- ggplot(manual_P2[-1,]) +
  geom_line(aes(x = rep(1:150,2), y = pop_n, colour = type, linetype = type), size = 2) +
  geom_vline(xintercept = 100)+
  scale_x_continuous(breaks = c(98, 100, 102))+
  coord_cartesian(xlim = c(97,103), ylim = c(50, 150)) +
  ggtitle("Climate anomaly at yr 100 = 2")

Man_neg <- ggplot(manual_P_2[-1,]) +
  geom_line(aes(x = rep(1:150,2), y = pop_n, colour = type, linetype = type), size = 2) +
  geom_vline(xintercept = 100)+
  scale_x_continuous(breaks = c(98, 100, 102)) +
  coord_cartesian(xlim = c(97,103), ylim = c(50,150))  +
  ggtitle("Climate anomaly at yr 100 = -2")

Man_pos / Man_neg

