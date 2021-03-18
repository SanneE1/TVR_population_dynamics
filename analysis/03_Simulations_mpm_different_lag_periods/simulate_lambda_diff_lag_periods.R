### Run MPM simulations for different lag periods
# other packages are sourced below
library(patchwork)

### Create line sourcing function 
source_lines <- function(file, lines){
  source(textConnection(readLines(file)[lines]))
}

## Get all the packages & functions for simulations
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(1:92))

# Overwrite output directory to current section
output_dir <- "results/03_Simulations_mpm_different_lag_periods/"

## Run lag-analysis on signal strength 0.5
i = 0.5

# Set up parallel and get clim_sd&clim_cor sequences
source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(98:113))

# Run simulations for different lag periods 1-4 yrs
for(l in c(1:4)) {
  
  ### Create climate sequences with different lag periods
  lag_clim <- lapply(as.list(c(1:900)), function(x) create_seq(5000, clim_sd = clim_sd[x], clim_auto = clim_auto[x], lag = l))
  
  ## Run simulations with lag in U/F matrix
  source_lines("analysis/01_Simulations_mpm_same_direction/simulate_mpm.R", c(194:226))
  
  ## Save under different names
  saveRDS(lag_fp, file.path(output_dir, paste("mpm_", l, "_lagfp.RDS", sep = "")))
  
  
}

# stop clusters
stopCluster(cl)


#### Read RDS files of simulations and plot results
format_df <- function(df) {
  lapply(df, function(x) lapply(x, function(y) y$df) %>% bind_rows %>% 
           mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
    bind_rows(., .id = "type")
}


lag_fp1 <- readRDS(file.path(output_dir, "mpm_1_lagfp.RDS")) %>% format_df %>% mutate(lag = 1)
lag_fp2 <- readRDS(file.path(output_dir, "mpm_2_lagfp.RDS")) %>% format_df %>% mutate(lag = 2)
lag_fp3 <- readRDS(file.path(output_dir, "mpm_3_lagfp.RDS")) %>% format_df %>% mutate(lag = 3)
lag_fp4 <- readRDS(file.path(output_dir, "mpm_4_lagfp.RDS")) %>% format_df %>% mutate(lag = 4)


plot_line <- function(df) {
   ggplot(df) +
    geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    labs(colour = "Lag type") +
    facet_grid(cols = vars(auto_cat)) + scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
}

l1 <- lag_fp1 %>% plot_line
l2 <- lag_fp2 %>% plot_line
l3 <- lag_fp3 %>% plot_line
l4 <- lag_fp4 %>% plot_line


line_plot <- l1 + l2 + l3 + l4 + plot_layout(nrow = 4, guides = "collect") + 
  plot_annotation(tag_levels = '1', tag_suffix = "yr lag",
                  title = "Lagged climate in P or F kernel", 
                  subtitle = "survival climate effect = pos \ngrowth climate effect= pos"
                  )

ggsave(plot = line_plot, file.path(output_dir, "comparative_line_plot.tiff"))


rel_line_plot1 <- rbind(lag_fp1, lag_fp2, lag_fp3, lag_fp4) %>% 
  pivot_wider(values_from = lambda, names_from = type) %>%
  mutate(rel.F.diff = (Fmatrix - none)/abs(none)) %>%
  ggplot(.) +
  geom_smooth(aes(x = clim_sd, y = rel.F.diff, colour = as.factor(lag)))+ 
  labs(colour = "# of yrs lagged") +
  facet_grid(cols = vars(auto_cat)) + 
  theme(legend.position = "bottom")

rel_line_plot2 <- rbind(lag_fp1, lag_fp2, lag_fp3, lag_fp4) %>% 
  pivot_wider(values_from = lambda, names_from = type) %>%
  mutate(rel.F.diff = (Fmatrix - none)/abs(none)) %>%
  ggplot(.) +
  geom_smooth(aes(x = clim_sd, y = rel.F.diff, colour = auto_cat))+ 
  labs(colour = "Autocorrelation") +
  facet_grid(cols = vars(lag)) + 
  theme(legend.position = "bottom")

ggsave(plot = rel_line_plot1, file.path(output_dir, "rel_line_plot1.tiff"))
ggsave(plot = rel_line_plot2, file.path(output_dir, "rel_line_plot2.tiff"))

box_plot <- rbind(lag_fp1, lag_fp2, lag_fp3, lag_fp4) %>% 
  pivot_wider(values_from = lambda, names_from = type) %>%
  mutate(rel.F.diff = (Fmatrix - none)/abs(none)) %>% 
  filter(round(clim_sd, digits = 1) == 1.1) %>%
  ggplot(.) +
  geom_boxplot(aes(x = auto_cat, y = rel.F.diff, fill = as.factor(lag))) +
  xlab("Climate autocorrelation") + ylab("Relative difference with lagged climate driver") +
  labs(fill = "# of yrs lagged", title = "Simulated stochastic lambda at 1.12 climate SD") + 
  theme_minimal() + geom_hline(yintercept = 0, size = 1)

ggsave(plot = box_plot, file.path(output_dir, "comparative_rel_box_plot.tiff"))


