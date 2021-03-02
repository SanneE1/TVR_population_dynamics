library(tidyverse)
library(patchwork)

files <- list.files(path = "results/compadre/", pattern = "mpm_", full.names = T)

for(i in c(1:length(files))) {
  
  df <- readRDS(files[i]) 
  df <- lapply( df, function(x) x %>% mutate(clim_sd = rep(seq(from = 0.01, to = 2, length.out = 10), 90),
                                      clim_auto = rep(rep(c(-0.9,0,0.9), each = 10), 30)))
  
  species <- stringr::str_match(files[i], "mpm_(.*?)_laguf.RDS")[,2]
  
  df1 <- df %>% bind_rows()
  
  plot1 <- ggplot(df1) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    labs(colour = "Lag type", title = "Lagged climate in growth or survival") +
    facet_grid(cols = vars(clim_auto)) + xlim(c(0,2)) +
    scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + theme(legend.position = "bottom",
                                                                             plot.title = element_text(hjust = 0.5)) +
    plot_annotation(title = species)
  
  df2 <- lapply(df, function(x) x %>% pivot_wider(names_from = type, values_from = lambda, values_fn = list))
  df3 <- lapply(df2, function(x) unnest(x, cols = last_col()) %>% pivot_longer(cols = last_col() ,names_to = "type", values_to = "lambda"))
    bind_rows() %>%
    mutate(diff_U = (Umatrix - None)/abs(None))


  plot2 <- ggplot(df3) + geom_smooth(aes(x = clim_sd, y = diff_U, colour = as.factor(clim_auto))) +
    xlim(c(0,2)) + labs(colour = "autocorrelation",
                        title = paste("Proportional differences in \nlambda for", species, "\nwhen Umatrix responds to lagged climate"),
                        x = "Climate SD",
                        y = "Relative difference in lambda") +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))
  
  name1 <- paste0("lambdas_", species, ".tiff")
  name2 <- paste0("lambdas_", species, "_reldiff.tiff")
  
  direct <-"results/compadre/plots/"
  
  ggsave(plot1, filename = file.path(direct, name1))
  # ggsave(plot2, filename = file.path(direct, name2))
  
  
}




