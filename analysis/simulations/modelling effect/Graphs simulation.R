library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(patchwork)
library(lme4)
library(ipmr)

source("analysis/simulations/ipmr_functions.R")

P1 <- readRDS("results/simulations/P_1yr_Merged.RDS") 
# P1neg <- readRDS("results/simulations/P_neg_Merged.RDS") %>% select(clim_sd, autocorrelation, non_lagged, lagged)
# P2 <- readRDS("results/simulations/P_2yr_Merged.RDS") %>% select(clim_sd, autocorrelation, non_lagged, lagged)
# PF <- readRDS("results/simulations/PF_1yr_Merged.RDS")
# PF_neg <- readRDS("results/simulations/PF_neg_Merged.RDS")
# 
# FP <- PF %>% select(clim_sd, autocorrelation, non_lagged, lagged_F) 
# 
# PF <- PF %>% select(clim_sd, autocorrelation, non_lagged, lagged_P) 
# 
# FP_neg <- PF_neg %>% select(clim_sd, autocorrelation, non_lagged, lagged_F) 
# 
# PF_neg <- PF_neg %>% select(clim_sd, autocorrelation, non_lagged, lagged_P) 

lams <- lapply(list(P1 = P1,
                    P1neg = P1neg, 
                    P2 = P2,
                    PF = PF,
                    FP = FP,
                    PF_neg = PF_neg,
                    FP_neg = FP_neg),
               function(df) df %>% 
                 pivot_longer(cols = -c("clim_sd", "autocorrelation"), 
                              names_to = "type", 
                              values_to = "lambda")
)

plotP1 <- ggplot(lams$P1) + 
  geom_smooth(aes(x = clim_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
  facet_grid(cols = vars(autocorrelation)) + 
  ylim(c(-0.6, 0.1)) +
    ggtitle("P kernel lagged 1yr, F kernel insensitive") +
  theme(legend.position = "bottom")

P1table <- tableGrob(data.frame(vitalrate = c("survival", "growth", "fp", "fn"), 
                                "parameter value" = c(1.233434, -0.06588578, NA, NA),
                                lagged = c("no", "yes", NA, NA))) 


plotP1neg <- ggplot(lams$P1neg) + 
  geom_smooth(aes(x = clim_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
  facet_grid(cols = vars(autocorrelation)) + ylim(c(-0.6, 0.1)) +
  ggtitle("P kernel lagged 1yr, F kernel insensitive") +
  theme(legend.position = "none")

P1negtable <- tableGrob(data.frame(vitalrate = c("survival", "growth", "fp", "fn"), 
                                   "parameter value" = c(-1.233434, -0.06588578, NA, NA),
                                   lagged = c("no", "yes", NA, NA)))


plotFP <- ggplot(lams$FP)+ 
  geom_smooth(aes(x = clim_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
  facet_grid(cols = vars(autocorrelation)) + ylim(c(-0.6, 0.1)) +
  ggtitle("P climate sensitive, F kernel lagged") +
  theme(legend.position = "none")

FPtable <- tableGrob(data.frame(vitalrate = c("survival", "growth", "fp", "fn"), 
                                "parameter value" = c(1.233434, -0.06588578, -0.6169492, -0.3398345),
                                lagged = c("no", "no", "yes", "yes")))


plotPF <- ggplot(lams$FP)+ 
  geom_smooth(aes(x = clim_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
  facet_grid(cols = vars(autocorrelation)) + ylim(c(-0.6, 0.1)) +
  ggtitle("P kernel lagged, F climate sensitive") +
  theme(legend.position = "none")

PFtable <- tableGrob(data.frame(vitalrate = c("survival", "growth", "fp", "fn"), 
                                "parameter value" = c(1.233434, -0.06588578, -0.6169492, -0.3398345),
                                lagged = c("yes", "yes", "no", "no")))


plotFPneg <- ggplot(lams$FP_neg)+ 
  geom_smooth(aes(x = clim_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
  facet_grid(cols = vars(autocorrelation)) + ylim(c(-0.6, 0.1)) +
  ggtitle("P climate sensitive, F kernel lagged") +
  theme(legend.position = "none")

FPnegtable <- tableGrob(data.frame(vitalrate = c("survival", "growth", "fp", "fn"), 
                                "parameter value" = c(-1.233434, -0.06588578, -0.6169492, -0.3398345),
                                lagged = c("no", "no", "yes", "yes")))


plotPFneg <- ggplot(lams$FP_neg)+ 
  geom_smooth(aes(x = clim_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
  facet_grid(cols = vars(autocorrelation)) + ylim(c(-0.6, 0.1)) +
  ggtitle("P kernel lagged, F climate sensitive") +
  theme(legend.position = "none")

PFnegtable <- tableGrob(data.frame(vitalrate = c("survival", "growth", "fp", "fn"), 
                                "parameter value" = c(-1.233434, -0.06588578,  -0.6169492, -0.3398345),
                                lagged = c("yes", "yes", "no", "no")))




Plot <- plotP1 + P1table + plotP1neg + P1negtable +
  plotPF + PFtable + plotPFneg + PFnegtable + 
  plotFP + FPtable + plotFPneg + FPnegtable +
  plot_layout(ncol = 4, widths = c(2,1,2,1), guides = "collect") + 
  plot_annotation(theme = theme(legend.position = "bottom")) & theme(text = element_text(size = 18))


ggsave(Plot, filename = "results/simulations/Summary_plot.png", width = 20, height = 9)
