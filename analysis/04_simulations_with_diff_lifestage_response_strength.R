rm(list = ls())
set.seed(2)

library(dplyr)
library(ggplot2)

# Number of iterations to run the simulations
n_it = 10000

# Set the proportion of variance (p) explained by the climate driver in the temporal sequence.
i = 0.5

# Use lifehistory id = 60 which is also presented in the main manuscript
lh_df <- read.csv("results/life_histories_df.csv")
lh_df <- lh_df[60,]

# Source script with R functions required
source("R/lambda_functions.R")


#---------------------------------------------------------------------------------------------------------
# Start of analyses
#---------------------------------------------------------------------------------------------------------

## set up sequences of climate standard deviation and autocorrelation values so that there are 30 duplicates 
## for each combination of sd & autocorrelation value

df_clim <- expand.grid(clim_sd = seq(from = 0.01, to = 1, length.out = 5),
                       clim_auto =  c(-0.6,0,0.6),
                       rep = c(1:30))

#### Create main temporal sequences
lag_clim <- lapply(as.list(c(1:nrow(df_clim))), 
                   function(x) create_seq(n_it, 
                                          clim_sd = df_clim$clim_sd[x], 
                                          clim_auto = df_clim$clim_auto[x], 
                                          lag = 1))

#### Lagged effect between U & F matrices
lag_p_Sj <- lapply(lag_clim,
                function(x) st.lamb_weak(mpm_df = lh_df,
                                    env_surv = x$lagged,
                                    env_growth = x$lagged,
                                    env_reproduction = x$recent,
                                    clim_sd = sd(x$recent, na.rm = T),
                                    clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                    sig.strength = i, weak_vr = "Sj")
) %>% bind_rows() %>%
  tibble::add_column(lag_type = "Umatrix",
                     weak_vr = "Sj")

lag_n2_Sj <- lapply(lag_clim,
                 function(x) st.lamb_weak(mpm_df = lh_df,
                                     env_surv = x$recent,
                                     env_growth = x$recent,
                                     env_reproduction = x$recent,
                                     clim_sd = sd(x$recent, na.rm = T),
                                     clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                     sig.strength = i, weak_vr = "Sj")
) %>% bind_rows() %>%
  tibble::add_column(lag_type = "none",
                     weak_vr = "Sj")

#### Lagged effect between U & F matrices
lag_p_Sa <- lapply(lag_clim,
                   function(x) st.lamb_weak(mpm_df = lh_df,
                                            env_surv = x$lagged,
                                            env_growth = x$lagged,
                                            env_reproduction = x$recent,
                                            clim_sd = sd(x$recent, na.rm = T),
                                            clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                            sig.strength = i, weak_vr = "Sa")
) %>% bind_rows() %>%
  tibble::add_column(lag_type = "Umatrix",
                     weak_vr = "Sa")

lag_n2_Sa <- lapply(lag_clim,
                    function(x) st.lamb_weak(mpm_df = lh_df,
                                             env_surv = x$recent,
                                             env_growth = x$recent,
                                             env_reproduction = x$recent,
                                             clim_sd = sd(x$recent, na.rm = T),
                                             clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                             sig.strength = i, weak_vr = "Sa")
) %>% bind_rows() %>%
  tibble::add_column(lag_type = "none",
                     weak_vr = "Sa")

#### Lagged effect between U & F matrices
lag_p_rho <- lapply(lag_clim,
                   function(x) st.lamb_weak(mpm_df = lh_df,
                                            env_surv = x$lagged,
                                            env_growth = x$lagged,
                                            env_reproduction = x$recent,
                                            clim_sd = sd(x$recent, na.rm = T),
                                            clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                            sig.strength = i, weak_vr = "rho")
) %>% bind_rows() %>%
  tibble::add_column(lag_type = "Umatrix",
                     weak_vr = "rho")

lag_n2_rho <- lapply(lag_clim,
                    function(x) st.lamb_weak(mpm_df = lh_df,
                                             env_surv = x$recent,
                                             env_growth = x$recent,
                                             env_reproduction = x$recent,
                                             clim_sd = sd(x$recent, na.rm = T),
                                             clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                             sig.strength = i, weak_vr = "rho")
) %>% bind_rows() %>%
  tibble::add_column(lag_type = "none",
                     weak_vr = "rho")


#### Lagged effect between U & F matrices
lag_p_all <- lapply(lag_clim,
                   function(x) st.lamb_weak(mpm_df = lh_df,
                                            env_surv = x$lagged,
                                            env_growth = x$lagged,
                                            env_reproduction = x$recent,
                                            clim_sd = sd(x$recent, na.rm = T),
                                            clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                            sig.strength = i, weak_vr = "all")
) %>% bind_rows() %>%
  tibble::add_column(lag_type = "Umatrix",
                     weak_vr = "all equal")

lag_n2_all <- lapply(lag_clim,
                    function(x) st.lamb_weak(mpm_df = lh_df,
                                             env_surv = x$recent,
                                             env_growth = x$recent,
                                             env_reproduction = x$recent,
                                             clim_sd = sd(x$recent, na.rm = T),
                                             clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                             sig.strength = i, weak_vr = "all")
) %>% bind_rows() %>%
  tibble::add_column(lag_type = "none",
                     weak_vr = "all equal")



lag_df <- rbind(lag_p_Sj, lag_n2_Sj, lag_p_Sa, lag_n2_Sa, lag_p_rho, lag_n2_rho, lag_p_all, lag_n2_all)

write.csv(lag_df, file.path("results", "lhID_60_diff_lifestage_resp_lambda_vs_climsd.csv"), 
          row.names = F)

# lag_df <- read.csv(file.path("results", "lhID_60_diff_windows_lambda_vs_climsd.csv"))

df <- lag_df %>% 
  mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c(-0.6, 0, 0.6))) 

plot <- ggplot(df) +
  geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(lag_type), linetype = auto_cat), se = F) +
  scale_colour_manual(name = "Simulation type",
                      values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      labels = c("none" = "control", "Umatrix" = "MCD")) +
  scale_linetype(name = "Autocorrelation") +
  ylab("stochastic log lambda") + xlab(~ paste(sigma[c], " of environmental sequence")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin=margin(),
        text = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) +
  facet_wrap(vars(weak_vr)) 

ggsave(plot = plot, filename = file.path("results", "lhID_60_diff_lifestage_resp_lambda_vs_climsd.png"),
       width = 6, height = 6, dpi = 400)


