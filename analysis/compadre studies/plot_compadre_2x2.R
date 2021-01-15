library(dplyr)
library(ggplot2)
library(patchwork)
#-----------------------------------------------------------
# plot results
#-----------------------------------------------------------

# Summary plot function for easy looping
sum_plot <- function(lag, lag_pf){
  ### Climate variables
  clim_sd <- rep(seq(from = 0, to = 2, length.out = 10), 90)
  clim_corr <- rep(rep(c(-0.9,0,0.9), each = 10), 30)
  
  lag_df <- data.frame(clim_sd = clim_sd,
                       clim_corr = clim_corr,
                       lambda = unlist(lag),
                       type = rep(names(lag), each = length(lag[[1]])))
  lag_p <- ggplot(lag_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    labs(colour = "Lag type", title = "Lagged climate in growth or survival") +
    facet_grid(cols = vars(clim_corr)) + 
    scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + theme(legend.position = "bottom",
                                                                             plot.title = element_text(hjust = 0.5))
  
  lagpf_df <- data.frame(clim_sd = clim_sd,
                         clim_corr = clim_corr,
                         lambda = unlist(lag_pf),
                         type = rep(names(lag_pf), each = length(lag_pf[[1]])))
  lagpf_p <- ggplot(lagpf_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    labs(colour = "Lag type", title = "Lagged climate in P or F") +
    facet_grid(cols = vars(clim_corr)) + 
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))
  
  return(list(lag = lag_p,
              lagpf = lagpf_p))
}


### Load species
species <- stringr::str_match(list.files("results/compadre_2x2/", pattern = "_lag.RDS"), "mpm_(.*?)_lag.RDS")[,2]


for(i in species) {
  
  lag <- readRDS(paste0("results/compadre_2x2/mpm_", i, "_lag.RDS"))
  lagpf <- readRDS(paste0("results/compadre_2x2/mpm_", i, "_lagfp.RDS"))
  
plots <- sum_plot(lag, lagpf)

plots <- plots$lag / plots$lagpf / guide_area() +
  plot_layout(guides = "collect", heights = c(3,3,1)) +
  plot_annotation(title = i) 

ggsave(plots, filename = paste0("results/compadre_2x2/plots_", i, ".tiff"))

}


