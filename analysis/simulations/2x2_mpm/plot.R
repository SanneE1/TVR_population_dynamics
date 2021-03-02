library(ggplot2)
library(dplyr)
library(patchwork)



## lambda summary plot -------------------------------------

sum_plot <- function(auto, cov, lag, lag_pf, sig.strength, save = T){
  ### Climate variables
  clim_sd <- rep(seq(from = 0, to = 2, length.out = 10), 90)
  clim_corr <- rep(rep(c(-0.9,0,0.9), each = 10), 30)
  ## Plot
  auto_df <- data.frame(clim_sd = clim_sd,
                        clim_corr = clim_corr,
                        lambda = unlist(auto))
    auto_p <- ggplot(auto_df) + geom_smooth(aes(x = clim_sd, y = lambda, linetype = as.factor(clim_corr)), colour = "grey30") + 
      labs(linetype = "Climate \nautocorrelation", title = "Autocorrelation \nin growth") 
  
  cov_df <- data.frame(clim_sd = clim_sd,
                       clim_corr = clim_corr,
                       lambda = unlist(cov))
  cov_p <- ggplot(cov_df) + geom_smooth(aes(x = clim_sd, y = lambda, linetype = as.factor(clim_corr)), colour = "grey30")+ 
    labs(linetype = "Climate \nautocorrelation", title = "Covariance in \nsurvival & growth")
  
  lag_df <- data.frame(clim_sd = clim_sd,
                       clim_corr = clim_corr,
                       lambda = unlist(lag),
                       type = rep(names(lag), each = length(lag[[1]])))
    lag_p <- ggplot(lag_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
      labs(colour = "Lag type", title = "Lagged climate in growth or survival") +
      facet_grid(cols = vars(clim_corr))  
  
   lagpf_df <- data.frame(clim_sd = clim_sd,
                        clim_corr = clim_corr,
                        lambda = unlist(lag_pf),
                        type = rep(names(lag_pf), each = length(lag_pf[[1]])))
   lagpf_p <- ggplot(lagpf_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
     labs(colour = "Lag type", title = "Lagged climate in P or F") +
     facet_grid(cols = vars(clim_corr)) + scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
  
  plot_AC <- auto_p + cov_p + plot_layout(guides = "collect") + plot_annotation(title = paste("Climate signal strenth =", sig.strength))
  plot_lag <-  (lag_p / lagpf_p) + plot_layout(guides = "collect") + plot_annotation(title = paste("Climate signal strenth =", sig.strength))
  
  if(save) {
  ggsave(filename = file.path("results/simulations/mpm/", paste0("auto_cov_plot_signal_strength_", sig.strength, ".png")), 
         plot_AC)
    ggsave(filename = file.path("results/simulations/mpm/", paste0("lag_plot_signal_strength_", sig.strength, ".png")), 
           plot_lag)
  }
  
  return(list(plot_AC = plot_AC, 
              plot_lag = plot_lag))
}

## function to create dataset formate needed and plot cell values --------------------------

plot_cell <- function(mats, title) {
  values <- rbind(low = mats[[2]]$mats, medium = mats[[6]]$mats, high = mats[[10]]$mats) %>% 
    tibble::rownames_to_column(., "SD") %>% 
    mutate(SD = regmatches(SD, regexpr("^[[:lower:]]{3}", SD))) %>% 
    pivot_longer(cols = -SD, names_to = "cell", values_to = "value")
  
  plot <- ggplot(values, aes(value, fill = SD, colour = SD)) + 
    geom_histogram(position = "dodge") + facet_wrap(vars(cell), scales = "free") +
    scale_y_log10() + ggtitle(title)
  
  return(plot)
}


### Load results
auto <- lapply(list.files("results/simulations/mpm/", pattern = "auto.RDS", full.names = T), readRDS)
cov <- lapply(list.files("results/simulations/mpm/", pattern = "cov.RDS", full.names = T), readRDS)
lag <- lapply(list.files("results/simulations/mpm/", pattern = "lag.RDS", full.names = T), readRDS)
lagpf <- lapply(list.files("results/simulations/mpm/", pattern = "lagfp.RDS", full.names = T), readRDS)
sig.strength <- regmatches(list.files("results/simulations/mpm/", pattern = "auto.RDS"), regexec("mpm\\_(.+)\\_", list.files("results/simulations/mpm/", pattern = "auto.RDS"))) %>% lapply(., function(x) x[2])
m <- pmap(list(auto, cov, lag, lagpf, sig.strength), sum_plot)



## Save histograms of cell values
pdf(paste("cell_value_plots_sigstrength", i, ".pdf"))

plot_cell(auto, "autocorrelation")
plot_cell(cov, "covariation")

plot_cell(lag_s, "survival lagged")
plot_cell(lag_g, "growth lagged")
plot_cell(lag_n, "neither vital rate lagged")

plot_cell(lag_p, "P lagged")
plot_cell(lag_f, "F lagged")
plot_cell(lag_n2, "Neither submatrices lagged")

dev.off()

}