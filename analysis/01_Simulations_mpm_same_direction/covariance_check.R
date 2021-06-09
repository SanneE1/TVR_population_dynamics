rm(list=ls())
library(dplyr)
library(tidyr)
library(popbio)
library(pbapply)
library(parallel)
library(ggplot2)
library(faux)
library(boot)
library(patchwork)

set.seed(2)

output_dir <- "results/01_Simulations_mpm_same_directions/"


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
    dev <- growth * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + 
      rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
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
    dev <- survival * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
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
    dev <- reproduction * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + rnorm(1, mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd) 
    p <- pnorm(dev, mean = 0, sd = clim_sd) 
    mpm[1,2] <- qgamma(p, (f_mean^2)/(f_sd * clim_sd)^2, (f_mean)/(f_sd * clim_sd)^2)
  }
  
  return(mpm)  
}


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

#---------------------------------------------------------------------
# Effect of lagged climate drivers on interannual variance in lambda
#---------------------------------------------------------------------

a <- purrr::pmap(list(survival = c(0,0,2,0,0), 
                      growth = c(0,0,2,0,0), 
                      reproduction = c(0,0,0,2,0),
                      clim_sd = 2,
                      sig.strength = 1), 
                 mpm)
la <- lapply(a, lambda) %>% unlist


b <- purrr::pmap(list(survival = c(0,0,2,0,0), 
                      growth = c(0,0,2,0,0), 
                      reproduction = c(0,0,2,0,0),
                      clim_sd = 2,
                      sig.strength = 1), 
                 mpm)
lb <- lapply(b, lambda) %>% unlist

df <- data.frame(climate = c(0,0,2,0,0),
                 recent = lb,
                 lagged = la) %>%
  mutate(diff_r = recent - lag(recent, default = first(recent)),
         diff_l = lagged - lag(lagged, default = first(lagged)))

ggplot(df) + 
  geom_point(aes(x = as.integer(row.names(df)), y=climate), size = 2) +
  geom_line(aes(x = as.numeric(row.names(df)), y = recent), colour = "#E7B800", size = 1) + 
  geom_line(aes(x = as.numeric(row.names(df)), y = lagged), colour = "#FC4E07", size = 1) +
  xlab("timestep") + ylab("climate anomaly & lambda") + theme_minimal()


sum(abs(df$diff_r))
sum(abs(df$diff_l))


#---------------------------------------------------------------------
# Difference in VR covariance under recent or mixed climate drivers
#---------------------------------------------------------------------
output_dir <- "results/01_Simulations_mpm_same_directions/"
location <- file.path(output_dir, "rds")

sig.strength <- regmatches(list.files(location, pattern = "auto.RDS"), 
                           regexec("mpm\\_(.+)\\_", list.files(location, pattern = "auto.RDS"))) %>% 
  lapply(., function(x) x[2])

lagpf_f <- lapply(list.files(location, pattern = "lagfp.RDS", full.names = T), readRDS)
U <- lapply(lagpf_f[[which(sig.strength == 0.5)]]$Umatrix, function(x) x$mats)
C <- lapply(lagpf_f[[which(sig.strength == 0.5)]]$none, function(x) x$mats)

U1 <- tibble(covariance = lapply(U, function(x) data.frame(cov = cov(x) %>% as.vector,
                                      from = rep(c("JuvSurv", "Progression", "Fecundity", "Survival"), each = 4),
                                      to = rep(c("JuvSurv", "Progression", "Fecundity", "Survival"), 4))),
            type = "Umatrix lagged",
            clim_sd = rep(seq(from = 0.01, to = 2, length.out = 10), 90),
            clim_auto = rep(rep(c(-0.9,0,0.9), each = 10), 30)
) %>% unnest(cols = covariance)

C1 <- tibble(covariance = lapply(C, function(x) data.frame(cov = cov(x) %>% as.vector,
                                                           from = rep(c("JuvSurv", "Progression", "Fecundity", "Survival"), each = 4),
                                                           to = rep(c("JuvSurv", "Progression", "Fecundity", "Survival"), 4))),
             type = "control",
             clim_sd = rep(seq(from = 0.01, to = 2, length.out = 10), 90),
             clim_auto = rep(rep(c(-0.9,0,0.9), each = 10), 30)
) %>% unnest(cols = covariance)

df1 <- rbind(U1, C1) %>% filter(clim_auto == 0 & clim_sd == 2)


PS <- ggplot(df1 %>% filter(from == "Progression" & to == "Survival")) + 
  geom_density(aes(x = cov, colour = type, fill = type), alpha = 0.5) +
  ggtitle("progression & survival")

PF <- ggplot(df1 %>% filter(from == "Progression" & to == "Fecundity")) + 
  geom_density(aes(x = cov, colour = type, fill = type), alpha = 0.5) +
  ggtitle("progression & fecundity")

FS <- ggplot(df1 %>% filter(from == "Fecundity" & to == "Survival")) + 
  geom_density(aes(x = cov, colour = type, fill = type), alpha = 0.5) +
  ggtitle("fecundity & survival")

plot <- PF + FS + PS + guide_area() + 
  plot_layout(nrow = 2, guides = "collect") + 
  plot_annotation(title = "Covariance between:")

ggsave(filename = file.path(output_dir, "Covariance_VR_auto0_sd2.tiff"), plot = plot,
       width = 5, height = 5)
