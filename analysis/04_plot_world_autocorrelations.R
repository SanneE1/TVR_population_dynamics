# Script to plot the world maps of precipitation and temperature autocorrelation

rm(ls())

library(tidyverse)
library(raster)
library(rasterVis)
library(patchwork)

set.seed(2)

temp <- raster("data/tmean_cor.tif")
precip <- raster("data/prec_cor.tif")
world <- map_data("world")
locations <- read.csv("data/species_authors.csv")


temp_plot <- gplot(temp) + 
  geom_tile(aes(fill = value)) +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "grey40", fill = "transparent", size = 0.3
  ) +
  geom_point(data = locations, aes(x = Lon, y = Lat), shape = 17, size = 2) +
  scale_fill_gradient2(low="blue", high="red", mid = "white", 
                       limits = c(-0.94, 0.94), na.value="transparent", 
                       name = "Autocorrelation value") + 
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(xlim = c(-180, 170), ylim = c(-90, 80)) +
  ylab("Latitude") + xlab("Longitude")  + guides(legend = "Autocorrelation value") +
  theme_bw() + 
  theme(text = element_text(size = 18), 
        legend.position = "bottom",
        panel.background = element_rect(fill = "gray95"))


precip_plot <- gplot(precip) + 
  geom_tile(aes(fill = value), na.rm = T) +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "grey40", fill = "transparent", size = 0.3
  ) +
  geom_point(data = locations, aes(x = Lon, y = Lat), shape = 17, size = 2) +
  scale_fill_gradient2(low="blue", high="red", mid = "white", 
                       limits = c(-0.94, 0.94), na.value="transparent", 
                       name = "Autocorrelation value") + 
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(xlim = c(-180, 170), ylim = c(-90, 80)) +
  ylab("Latitude") + xlab("Longitude")  + guides(legend = "Autocorrelation value") +
  theme_bw() + 
  theme(text = element_text(size = 18), 
        legend.position = "bottom",
        panel.background = element_rect(fill = "gray95"))

plot2 <- temp_plot / precip_plot + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

ggsave(plot2, filename = "results/autocorrelation_temp_precip.png",
       width = 9, height = 10)





