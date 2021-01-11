library(ggplot2)
library(dplyr)
library(patchwork)


sum_plot <- function(auto, cov, lag, lag_pf){
  ### Climate variables
  clim_sd <- rep(seq(from = 0, to = 2, length.out = 10), 90)
  clim_corr <- rep(rep(c(-0.9,0,0.9), each = 10), 30)
  ## Plot
  auto_df <- data.frame(clim_sd = clim_sd,
                        clim_corr = clim_corr,
                        lambda = unlist(auto))
    auto_p <- ggplot(auto_df) + geom_smooth(aes(x = clim_sd, y = lambda, linetype = as.factor(clim_corr)), colour = "grey30") + 
      labs(linetype = "Climate autocorrelation", title = "Autocorrelation \nin growth")
  
  cov_df <- data.frame(clim_sd = clim_sd,
                       clim_corr = clim_corr,
                       lambda = unlist(cov))
  cov_p <- ggplot(cov_df) + geom_smooth(aes(x = clim_sd, y = lambda, linetype = as.factor(clim_corr)), colour = "grey30")+ 
    labs(linetype = "Climate autocorrelation", title = "Covariance in \nsurvival & growth")
  
  lag_df <- data.frame(clim_sd = clim_sd,
                       clim_corr = clim_corr,
                       lambda = unlist(lag),
                       type = rep(names(lag), each = length(lag[[1]])))
    lag_p <- ggplot(lag_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
      labs(colour = "Lag type", title = "Lagged climate in growth \nor survival") +
      facet_grid(cols = vars(clim_corr))
  
   lagpf_df <- data.frame(clim_sd = clim_sd,
                        clim_corr = clim_corr,
                        lambda = unlist(lag_pf),
                        type = rep(names(lag_pf), each = length(lag_pf[[1]])))
   lagpf_p <- ggplot(lagpf_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
     labs(colour = "Lag type", title = "Lagged climate in P \nor F") +
     facet_grid(cols = vars(clim_corr))
    
  return(list(auto = auto_p,
              cov = cov_p, 
              lag = lag_p,
              lagpf = lagpf_p))
}

### Load results

# positive -> s = inv_logit(climate), g = inv_logit(climate), f = exp(1.2)
auto <- readRDS("results/simulations/mpm/mpm_1_auto.RDS")
cov <- readRDS("results/simulations/mpm/mpm_1_cov.RDS")
lag <- readRDS("results/simulations/mpm/mpm_1_lag.RDS")
lagpf <- readRDS("results/simulations/mpm/mpm_1_lagfp.RDS")

plots1 <- sum_plot(auto, cov, lag, lagpf)

mat1 <- matrix(nrow = 2, ncol = 2, dimnames = list(c("seedlings", "reproductive"), c("seedlings", "reproductive")))
mat1[2,1] <- "inv_logit(0 + climate)"
mat1[2,2] <- "inv_logit(0 + climate)"
mat1[1,2] <- "exp(1.2 + climate)"

## Autocorrelation and covariance effects
auto1 <- plots1$auto + plots1$cov + plot_layout(guides = "collect") + plot_annotation(title = "Positive climate effect")

## Lagged effects
lag1 <- plots1$lag / plots1$lagpf + 
  plot_layout(guides = "collect") + plot_annotation(title = "Positive climate effect")


# negative -> s = inv_logit(-climate), g = inv_logit(-climate), f = exp(1.2)
auto_n <- readRDS("results/simulations/mpm/mpm_-1_auto.RDS")
cov_n <- readRDS("results/simulations/mpm/mpm_-1_cov.RDS")
lag_n <- readRDS("results/simulations/mpm/mpm_-1_lag.RDS")
lagpf_n <- readRDS("results/simulations/mpm/mpm_-1_lagfp.RDS")

mat2 <- matrix(nrow = 2, ncol = 2, dimnames = list(c("seedlings", "reproductive"), c("seedlings", "reproductive")))
mat2[2,1] <- "inv_logit(0 - climate)"
mat2[2,2] <- "inv_logit(0 - climate)"
mat2[1,2] <- "exp(1.2 - climate)"

plots2 <- sum_plot(auto_n, cov_n, lag_n, lagpf_n)

## Autocorrelation and covariance effects
auto2 <- plots2$auto + plots2$cov

## Lagged effects
lag2 <- plots2$lag / plots2lagpf + 
  plot_layout(guides = "collect") + plot_annotation(title = "Negative climate effects")


# half positive effect -> s = inv_logit(-climate * 0.5), g = inv_logit(-climate * 0.5), f = exp(1.2)
auto_s <- readRDS("results/simulations/mpm/mpm_0.5_auto.RDS")
cov_s <- readRDS("results/simulations/mpm/mpm_0.5_cov.RDS")
lag_s <- readRDS("results/simulations/mpm/mpm_0.5_lag.RDS")
lagpf_s <- readRDS("results/simulations/mpm/mpm_0.5_lagfp.RDS")


mat3 <- matrix(nrow = 2, ncol = 2, dimnames = list(c("seedlings", "reproductive"), c("seedlings", "reproductive")))
mat3[2,1] <- "inv_logit(0 + climate * 0.5)"
mat3[2,2] <- "inv_logit(0 + climate * 0.5)"
mat3[1,2] <- "exp(1.2 + climate * 0.5)"


plots3<- sum_plot(auto_s, cov_s, lag_s, lagpf_s)

## Autocorrelation and covariance effects
auto3 <- plots3$auto+ plots3$cov

## Lagged effects
lag3 <- plots3$lag / plots3$lagpf + 
  plot_layout(guides = "collect") + plot_annotation(title = "Negative climate effects")
