library(tidyverse)
library(patchwork)

output_dir <- "results/06_COMPADRE_studies/"
files <- list.files(path = file.path(output_dir, "rds"), pattern = "mpm_", full.names = T)

for(i in c(1:length(files))) {
  
  results <- readRDS(files[i]) 
  species <- stringr::str_match(files[i], "mpm_(.*?)_laguf.RDS")[,2]

  df <- lapply( results, function(x) x %>% bind_rows) %>% 
    bind_rows(., .id = "type") %>% 
    mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))
  
  plot1 <- ggplot(df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    labs(colour = "Lag type", title = "Lagged climate in growth or survival") +
    facet_grid(cols = vars(auto_cat)) + xlim(c(0,2)) +
    scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + theme(legend.position = "bottom",
                                                                             plot.title = element_text(hjust = 0.5)) +
    plot_annotation(title = species)
  
  df2 <- df %>% 
    pivot_wider(names_from = type, values_from = "lambda", values_fn = list) %>%
    unnest(cols = c(Umatrix, Fmatrix, None)) %>% 
    mutate(rel.diff.U = (Umatrix-None)/abs(None))


  plot2 <- ggplot(df2) + geom_smooth(aes(x = clim_sd, y = rel.diff.U, colour = as.factor(auto_cat))) +
    xlim(c(0,2)) + labs(colour = "autocorrelation",
                        title = paste("Proportional differences in \nlambda for", species, "\nwhen Umatrix responds to lagged climate"),
                        x = "Climate SD",
                        y = "Relative difference in lambda") +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))
  
  name1 <- paste0("lambdas_", species, ".tiff")
  name2 <- paste0("lambdas_", species, "_reldiff.tiff")
  
  ggsave(plot1, filename = file.path(output_dir, "lambda_plots", name1))
  ggsave(plot2, filename = file.path(output_dir, "difference_plots", name2))
  
  
}




