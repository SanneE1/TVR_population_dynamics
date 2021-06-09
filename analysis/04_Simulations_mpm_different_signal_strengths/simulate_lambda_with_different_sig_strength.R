### simulate the effect of recent and lagged climate drivers when signal strength of different VR depends on the 
### sensitivity of the VR
rm(list=ls())
library(dplyr)
library(tidyr)
library(popbio)
library(pbapply)
library(parallel)
library(ggplot2)
library(faux)
library(boot)
set.seed(2)

n_it = 5000

### Create line sourcing to 
source_lines <- function(file, lines){
  source(textConnection(readLines(file)[lines]))
}

## output directory to right folder
output_dir <- "results/04_simulations_mpm_different_signal_strengths/"

# Create function that creates environmental sequence ------------------------------------------------------------------------------------
create_seq <- function(n_it, clim_sd, clim_auto, lag) { 
  for(n in c(1:(n_it+lag))){ 
    if(n == 1) {
      seq <- rnorm(1)
    } else {
      seq[n] <- clim_auto * seq[n-1] + rnorm(1)
    }
  }
  seq <- scale(seq) * clim_sd
  
  lagged <- c(rep(NA, lag), seq)
  recent <- c(seq, rep(NA, lag))
  
  df <- data.frame(recent = recent,
                   lagged = lagged)
  return(df)
}

# Create MPM function that deals with different signal strengths
mpm <- function(survival, growth, reproduction, 
                clim_sd, sig.strength) {
  ## Basic mpm
  mpm <- matrix(0, nrow = 2, ncol = 2)
  
  ## survival Juveniles
  # sj_mean = 
  # sj_sd = 
  # if(is.na(survival)) {
  #   
  #   mpm[1,1] <- sj_mean
  #   
  # } else {
  #   
  #   dev <- survival .....  ## total deviation from mean = climate signal & random noise
  #   p <- pnorm(dev, mean = 0, sd = sqrt(clim_sd * 2)) ## Here I use sd = sqrt(2) because dev has variance var(growth) + var(noise). sqrt() to get stand. dev.
  #   q <- qnorm(p)
  #   mpm[1,1] <- qbeta(p, (((sj_mean*(1-sj_mean))/sj_sd^2) - 1) * sj_mean,
  #                     (((sj_mean*(1-sj_mean))/sj_sd^2) - 1) * (1 - sj_mean))
  #   
  # }
  
  # growth 
  g_mean = 0.325
  g_sd = 0.118
  
  if(is.na(growth)) {
    mpm[2,1] <- g_mean  
  } else {
    ## total deviation from mean = climate signal * signal strength & correction factor (partitioning at variance scale) + random noise * signal strength
    dev <- growth * (sqrt(clim_sd^2 * sig.strength[1])/clim_sd) + 
      rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength[1]))/clim_sd) 
    ## Because of partitioning and correction factor above, the resulting distribution has a sd of clim_sd
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    mpm[2,1] <- qbeta(p, (((g_mean*(1-g_mean))/(g_sd * clim_sd)^2) - 1) * g_mean,
                      (((g_mean*(1-g_mean))/(g_sd * clim_sd)^2) - 1) * (1 - g_mean))
  }
  
  # survival
  s_mean = 0.541
  s_sd = 0.135
  if(is.na(survival)) {
    mpm[2,2] <- s_mean
  } else {
    dev <- survival * (sqrt(clim_sd^2 * sig.strength[3])/clim_sd) + rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength[3]))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    mpm[2,2] <- qbeta(p, (((s_mean*(1-s_mean))/(s_sd * clim_sd)^2) - 1) * s_mean,
                      (((s_mean*(1-s_mean))/(s_sd * clim_sd)^2) - 1) * (1 - s_mean))
  }
  
  # reproduction 
  f_mean = 0.788
  f_sd = 0.862
  if(is.na(reproduction)) {
    mpm[1,2] <- f_mean
  } else {
    dev <- reproduction * (sqrt(clim_sd^2 * sig.strength[2])/clim_sd) + rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength[2]))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    mpm[1,2] <- qgamma(p, (f_mean^2)/(f_sd * clim_sd)^2, (f_mean)/(f_sd * clim_sd)^2)
  }
  
  return(mpm)  
}

# Function to create a MPM sequence from supplied climate and noise sequences and calculate stocastic lambda -----------------------------------

st.lamb <- function(env_surv, env_growth, env_reproduction, 
                    clim_sd, clim_auto, sig.strength) {
  
  n_it = length(env_surv)
  
  env <- list(survival = env_surv,
              growth = env_growth,
              reproduction = env_reproduction,
              clim_sd = clim_sd,
              sig.strength = list(sig.strength))
  
  ### Get all mpm's
  mats <- pmap(env, mpm) %>% Filter(Negate(anyNA), .)
  
  
  df <- data.frame(lambda = stoch.growth.rate(mats, maxt = n_it)$sim,
                   clim_sd = clim_sd,
                   clim_auto = clim_auto)
  
  return(list(df = df,
              mats = sapply(mats, as.vector) %>% t %>% 
                `colnames<-`(c("1,1", "2,1", "1,2", "2,2"))))
}


### calculate elasticities
A <- matrix2(c(0,0.788,0.325,0.541), stages = c("Juvenile", "Adult"))
E <- elasticity(A)


### Elasticities of VR need to be negatively correlated with the climate variance in VR (and thus sig.strength)
sig.strength <- 1-as.vector(E)[-1]
sig.strength <- sig.strength/sum(sig.strength)


### Simulate Lambda under differnt SD and autocorrelation
i = sig.strength ## sig.strength needs to be put to i as i was used to loop through sig.strenghts in the original file

# Set up parallel and run U&F simulations
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", 
             c(139:152, 196:200, 238:274))
stopCluster(cl)

### Save output
saveRDS(lag_fp, file.path(output_dir, "mpm_diff_ss_laguf.RDS"))

### Plot lambda's

lag_df <- lapply(lag_fp, function(x) lapply(x, function(y) y$df) %>%
                                           bind_rows) %>%
  bind_rows(., .id = "type") %>%
  mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))

lag_p <- ggplot(lag_df %>% filter(type != "Fmatrix")) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(type)), position = position_dodge(width = 0.1))+ 
  labs(colour = "Lag type", title = "Lagged climate in growth or survival", 
       subtitle = paste("sig.strength surv:", round(sig.strength[3], digits = 2), "growth:",round(sig.strength[1], digits = 2),
                        "fec:", round(sig.strength[2], digits = 2) )) +
  facet_grid(cols = vars(auto_cat))

lagpf_df <- lapply(lag_fp, function(x) lapply(x, function(y) y$df) %>%
                     bind_rows) %>%
  bind_rows(., .id = "type") %>%
  mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))

lagpf_p <- ggplot(lagpf_df %>% filter(type != "Fmatrix")) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
  geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(type)), position = position_dodge(width = 0.1))+ 
  labs(colour = "Lag type", title = "Lagged climate in U or F matrix", 
       subtitle = paste("sig.strength surv:", round(sig.strength[3], digits = 2), "growth:",round(sig.strength[1], digits = 2),
                        "fec:", round(sig.strength[2], digits = 2) )) +
  facet_grid(cols = vars(auto_cat)) + scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))

ggsave(file.path(output_dir, "lambda_plot_lag.tiff"), plot = lag_p)
ggsave(file.path(output_dir, "lambda_plot_lagfp.tiff"), plot = lagpf_p)
