library(tidyverse)
library(patchwork)

output_dir <- "results/06_COMPADRE_studies/"
files <- list.files(path = file.path(output_dir, "rds"), pattern = "mpm_", full.names = T)

plot_cell <- function(mats, title) {
  
  df <- data.frame(clim_sd = rep(seq(from = 0.01, to = 2, length.out = 10), 90),
                   clim_auto = rep(rep(c(-0.9,0,0.9), each = 10), 30)
  )
  
  # Randomly select a sequence with a low and high sd
  low_sd <- sample(which(round(df$clim_sd, digits = 3) == 0.231), 1)
  high_sd <- sample(which(round(df$clim_sd, digits = 3) == 2), 1)
  
  
  values <- rbind(low = mats[[low_sd]]$mats, high = mats[[high_sd]]$mats) %>% 
    tibble::rownames_to_column(., "SD") %>% 
    mutate(SD = regmatches(SD, regexpr("^[[:lower:]]{3}", SD))) %>% 
    pivot_longer(cols = -SD, names_to = "cell", values_to = "value")
  
  ggplot(values, aes(value, fill = SD, colour = SD)) + 
    geom_histogram(position = "dodge") + facet_wrap(vars(cell), scales = "free") +
    scale_y_log10() + ggtitle(title)
}

for(i in c(1:length(files))) {
  
  results <- readRDS(files[i]) 
  df <- lapply( results, function(x) x %>% mutate(clim_sd = rep(seq(from = 0.01, to = 2, length.out = 10), 90),
                                      clim_auto = rep(rep(c(-0.9,0,0.9), each = 10), 30)))
  
  species <- stringr::str_match(files[i], "mpm_(.*?)_laguf.RDS")[,2]
  
  df1 <- df %>% bind_rows() 
  
  plot1 <- ggplot(df1) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    labs(colour = "Lag type", title = "Lagged climate in growth or survival") +
    facet_grid(cols = vars(clim_auto)) + xlim(c(0,2)) +
    scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + theme(legend.position = "bottom",
                                                                             plot.title = element_text(hjust = 0.5)) +
    plot_annotation(title = species)
  
  df2 <- df1 %>% 
    pivot_wider(names_from = type, values_from = "lambda", values_fn = list) %>%
    unnest(cols = c(Umatrix, Fmatrix, None)) %>% 
    mutate(rel.diff.U = (Umatrix-None)/abs(None))


  plot2 <- ggplot(df2) + geom_smooth(aes(x = clim_sd, y = rel.diff.U, colour = as.factor(clim_auto))) +
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
  
  ## Save histograms of cell values 
  ## functions to create dataset formate needed, plot cell values and save to pdf --------------------------
  
  results
  
  
}




