rm(list=ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(popbio)
library(purrr)

output_dir <- "/gpfs1/data/lagged/results/01_Simulations_mpm_same_directions/"


##  -------------------------------------
##  lambda summary plot
##  -------------------------------------
label_auto <- c(
  "0.9" = "0.9 autocorrelation",
  "0" = "0 autocorrelation",
  "-0.9" = "-0.9 autocorrelation"
)

sum_plot <- function(file, ss){

  
  a = file[[1]]
  c = file[[2]]
  p = file[[3]]
  f = file[[4]]
  
  ### Climate variables
  clim_sd <- rep(seq(from = 0.01, to = 2, length.out = 10), 90)
  clim_corr <- rep(rep(c(-0.9,0,0.9), each = 10), 30)
  ## Plot
  auto_p <- ggplot(a) + geom_smooth(aes(x = clim_sd, y = lambda, linetype = as.factor(clim_auto)), colour = "grey30") + 
    labs(linetype = "Climate \nautocorrelation", title = "Autocorrelation in growth") +
    xlab("Climate standard deviation") + ylab("log lambda") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  cov_df <- c %>% mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))
  cov_p <- ggplot(cov_df) + geom_smooth(aes(x = clim_sd, y = lambda, linetype = as.factor(auto_cat)), colour = "grey30") + 
    labs(linetype = "Climate \ncovariance", title = "Covariance in \nsurvival & growth") +
    xlab("Climate standard deviation") + ylab("log lambda") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  lag_df <- lapply(p, function(x) x %>% mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
    bind_rows(., .id = "type")
  
  lag_p <- ggplot(lag_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)), se = F) + 
    geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(type)), position = position_dodge(width = 0.1)) + 
    labs(colour = "Lag type", title = "Lagged climate in growth or survival") +
    xlab("Climate standard deviation") + ylab("log lambda") +
    facet_grid(cols = vars(auto_cat), labeller = as_labeller(c("-0.9" = "Autocorrelation = -0.9",
                                                               "0" = "Autocorrelation = 0",
                                                               "0.9" = "Autocorrelation = 0.9"))) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  lagpf_df <- lapply(f, function(x) x %>% mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
    bind_rows(., .id = "type") %>% mutate(type = factor(type, levels = c("none", "Umatrix", "Fmatrix")))
  
  lagpf_p <- ggplot(lagpf_df) + 
    geom_smooth(aes(x = clim_sd, y = lambda, colour = type), se = F)+ 
    geom_point(aes(x = clim_sd, y = lambda, colour = type), position = position_dodge(width = 0.1))+ 
    labs(subtitle = "Same directional response of submatrices to climate driver") +
    ylab("stochastic log lambda") + xlab("SD of environmental sequence") + 
    facet_grid(cols = vars(auto_cat), labeller = labeller(auto_cat = label_auto)) + 
    scale_colour_manual(name = "Simulation type",
                        values = c("Umatrix" = "#E69F00", "none" = "#0072B2", "Fmatrix" = "#CC79A7"),
                        labels = c("none" = "control", "Umatrix" = "U_MCD", "Fmatrix" = "F_MCD")) + theme_minimal() +
    theme(legend.position = "bottom")
  
  plot_AC <- auto_p + cov_p + plot_annotation(title = paste("Climate signal strenth =", ss))
  plot_lag <-  (lag_p / lagpf_p) + plot_annotation(title = paste("Climate signal strenth =", ss))
  
  
  lagf <- lagpf_df %>% filter(type != "Fmatrix") %>%
    ggplot(.) + 
    geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)), se = F)+ 
    geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(type)), position = position_dodge(width = 0.1))+ 
    labs(subtitle = "Same directional response of submatrices to climate driver") +
    ylab("stochastic log lambda") + xlab("SD of environmental sequence") + 
    facet_grid(cols = vars(auto_cat), labeller = labeller(auto_cat = label_auto)) + 
    scale_colour_manual(name = "Simulation type",
                        values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                        labels = c("none" = "control", "Umatrix" = "MCD")) + theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(filename = file.path(output_dir, "lambda_plots", paste0("auto_cov_plot_", ss, ".png")), 
           plot_AC, width = 8.22, height = 5.53, units = "in")
  ggsave(filename = file.path(output_dir, "lambda_plots", paste0("lag_plot_", ss, ".png")), 
           plot_lag, width = 8.22, height = 5.53, units = "in")
  ggsave(filename = file.path(output_dir, "lambda_plots", paste0("UF_lag_plot_", ss, ".png")), 
         lagpf_p, width = 8.22, height = 4, units = "in")
  ggsave(filename = file.path(output_dir, "lambda_plots", paste0("U_lag_plot_", ss, ".png")),
         lagf, width = 8.22, height = 4)
  
  
  diff <- lagpf_df %>% pivot_wider(names_from = "type", values_from = lambda) %>%
    mutate(diff_U = (Umatrix - none)/ abs(none)) %>%
    ggplot(.) + geom_boxplot(aes(x = factor(round(clim_sd, digits = 1)), 
                                 y = diff_U, fill = auto_cat), 
                           position = "dodge") +
      geom_hline(aes(yintercept = 0)) + 
    ylab("Proportional difference in log lambda") +
    xlab("Climate standard deviation") + 
    labs(fill = "Climate \nautocorrelation", title = "Relative difference", 
         subtitle = "between control and simulations with lag in U-matrix") + 
    theme_minimal() +
    theme(legend.position = "bottom")
    
  
  ggsave(filename = file.path(output_dir, "lambda_plots", paste0("UF_relative_diff_", ss, ".png")),
         diff, width = 8.22, height = 5.53, units = "in")
  

}

### Load results
location <- file.path(output_dir, "rds")

print("loading RDSs done")
files_1 <- lapply(list.files(location, pattern = "1", full.names = T), readRDS)

print("formating autocorrelation df")
files_1[[1]] <- lapply(files_1[[1]], function(x) x$df) %>% bind_rows

print("format covariance df")
files_1[[2]] <- lapply(files_1[[2]], function(x) x$df) %>% bind_rows

print("format sg simulation df")
files_1[[3]] <- lapply(files_1[[3]], function(x) lapply(x, function(y) y$df) %>% bind_rows)

print("format fp simulation df")
files_1[[4]] <- lapply(files_1[[4]], function(x) lapply(x, function(y) y$df) %>% Filter(Negate(anyNA), .) %>% bind_rows)

sum_plot(files_1, ss = 1)
rm(files_1)

files_0.5 <- lapply(list.files(location, pattern = "0.5", full.names = T), readRDS)

print("formating autocorrelation df")
files_0.5[[1]] <- lapply(files_0.5[[1]], function(x) x$df) %>% bind_rows

print("format covariance df")
files_0.5[[2]] <- lapply(files_0.5[[2]], function(x) x$df) %>% bind_rows

print("format sg simulation df")
files_0.5[[3]] <- lapply(files_0.5[[3]], function(x) lapply(x, function(y) y$df) %>% bind_rows)

print("format fp simulation df")
files_0.5[[4]] <- lapply(files_0.5[[4]], function(x) lapply(x, function(y) y$df) %>% Filter(Negate(anyNA), .) %>% bind_rows)

sum_plot(file = files_0.5, ss = "0.5")
rm(files_0.5)

files_0.25 <- lapply(list.files(location, pattern = "0.25", full.names = T), readRDS)

print("formating autocorrelation df")
files_0.25[[1]] <- lapply(files_0.25[[1]], function(x) x$df) %>% bind_rows

print("format covariance df")
files_0.25[[2]] <- lapply(files_0.25[[2]], function(x) x$df) %>% bind_rows

print("format sg simulation df")
files_0.25[[3]] <- lapply(files_0.25[[3]], function(x) lapply(x, function(y) y$df) %>% bind_rows)

print("format fp simulation df")
files_0.25[[4]] <- lapply(files_0.25[[4]], function(x) lapply(x, function(y) y$df) %>% Filter(Negate(anyNA), .) %>% bind_rows)

sum_plot(files_0.25, ss = "0.25")
rm(files_0.25)


files_0.05 <- lapply(list.files(location, pattern = "0.05", full.names = T), readRDS)

print("formating autocorrelation df")
files_0.05[[1]] <- lapply(files_0.05[[1]], function(x) x$df) %>% bind_rows

print("format covariance df")
files_0.05[[2]] <- lapply(files_0.05[[2]], function(x) x$df) %>% bind_rows

print("format sg simulation df")
files_0.05[[3]] <- lapply(files_0.05[[3]], function(x) lapply(x, function(y) y$df) %>% bind_rows)

print("format fp simulation df")
files_0.05[[4]] <- lapply(files_0.05[[4]], function(x) lapply(x, function(y) y$df) %>% Filter(Negate(anyNA), .) %>% bind_rows)

sum_plot(files_0.05, ss = "0.05")
rm(files_0.05)

##  -------------------------------------
##  Plot a lambda time sequence of simulation with i = 0.5 (50% of temporal variance explained by climate driver)
##  -------------------------------------
lagpf_0.5 <- readRDS(list.files(location, pattern = "0.5_lagfp.RDS", full.names = T))


l_seq0 <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(lagpf_0.5$Umatrix[[20]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                     none = lapply(as.list(as.data.frame(t(lagpf_0.5$none[[20]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                     time = c(1:50002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

l_seqneg <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(lagpf_0.5$Umatrix[[10]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                      none = lapply(as.list(as.data.frame(t(lagpf_0.5$none[[10]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                      time = c(1:50002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

l_seqpos <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(lagpf_0.5$Umatrix[[30]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        none = lapply(as.list(as.data.frame(t(lagpf_0.5$none[[30]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        time = c(1:50002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

### zoom in on 35 years
time_series0 <- ggplot(l_seq0) + 
  geom_line(aes(y = lambda, x = time, colour = type), show.legend = F) + 
  xlim(1095,1130) + ylim(c(-7.5,2)) + ylab("log lambda") + 
  labs(subtitle = "0 autocorrelation") + 
  scale_colour_manual(values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      labels = c("none" = "control", "Umatrix" = "MCD"))


time_seriesneg <- ggplot(l_seqneg) + 
  geom_line(aes(y = lambda, x = time, colour = type), show.legend = F) + 
  xlim(1095,1130) + ylim(c(-7.5,2)) + ylab("log lambda") + 
  labs(subtitle = "-0.9 autocorrelation") + 
  scale_colour_manual(values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      label = c("Umatrix" = "MCD", "none" = "control"))

time_seriespos <- ggplot(l_seqpos) + 
  geom_line(aes(y = lambda, x = time, colour = type), show.legend = F) + 
  xlim(1095,1130) + ylim(c(-7.5,2)) + ylab("log lambda") + 
  labs(subtitle = "0.9 autocorrelation") + 
  scale_colour_manual(values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                      label = c("none" = "control", "Umatrix" = "MCD"))

int_ann_0 <- l_seq0 %>% 
  group_by(type) %>% 
  arrange(time) %>%
  mutate(diff = (lambda %>% exp) - (lag(lambda, default = first(lambda)) %>% exp),
         type = factor(type, levels = c("Umatrix", "none"))) %>%
  ggplot(.) +
  geom_density(aes(x = diff, fill = type, colour = type), alpha = 0.5) +
  scale_fill_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "MCD"), 
                    values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +  
  scale_colour_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "MCD"), 
                      values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  xlab("interannual difference") + coord_cartesian(c(-3, 3))

int_ann_neg <- l_seqneg %>% 
  group_by(type) %>% 
  arrange(time) %>% 
  mutate(diff = (lambda %>% exp) - (lag(lambda, default = first(lambda)) %>% exp),
         type = factor(type, levels = c("Umatrix", "none"))) %>%
  ggplot(.) +
  geom_density(aes(x = diff, fill = type, colour = type), alpha = 0.5, show.legend = F)  +
  scale_fill_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "MCD"), 
                    values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +  
  scale_colour_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "MCD"), 
                      values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  xlab("interannual difference") + coord_cartesian(c(-3, 3))

int_ann_pos <- l_seqpos %>% 
  group_by(type) %>% 
  arrange(time) %>% 
  mutate(diff = (lambda %>% exp) - (lag(lambda, default = first(lambda)) %>% exp),
         type = factor(type, levels = c("Umatrix", "none"))) %>%
  ggplot(.) +
  geom_density(aes(x = diff, fill = type, colour = type), alpha = 0.5, show.legend = F) +
  scale_fill_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "MCD"), 
                    values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +  
  scale_colour_manual(name = "Simulation", label = c("none" = "control", "Umatrix" = "MCD"), 
                      values = c("Umatrix" = "#E69F00", "none" = "#0072B2")) +
  xlab("interannual difference") + coord_cartesian(c(-3, 3))



a <- wrap_plots((time_seriespos / time_series0 / time_seriesneg), (int_ann_pos / int_ann_0 / int_ann_neg), 
           byrow = F, 
           ncol = 2, 
           widths = c(3,1), 
           guides = "collect") & theme(legend.position = "bottom") & theme_minimal() 


a[[1]] <- a[[1]] + plot_layout(tag_level = "new")
a[[2]] <- a[[2]] + plot_layout(tag_level = "new")


time_diff <- a + plot_annotation(tag_levels = c('A', 'i'), tag_sep = '.') & theme(legend.position = "bottom")


ggsave(filename = file.path(output_dir, "lambda_timeseries", "timeseries&diff.tiff"), 
       plot = time_diff, width = 9.22, height = 6.5)


rm(lagpf_0.5)

##  -------------------------------------
##  Save histograms of cell values 
##  -------------------------------------
## functions to create dataset formate needed, plot cell values and save to pdf --------------------------

plot_cell <- function(mats, title) {
  
  df <- data.frame(clim_sd = rep(seq(from = 0.01, to = 2, length.out = 10), 90),
                   clim_auto = rep(rep(c(-0.9,0,0.9), each = 10), 30)
  )
  
  # Randomly select a sequence with a low and high sd
  low_sd <- sample(which(round(df$clim_sd, digits = 3) == 0.231 & df$clim_auto == 0), 1)
  mid_sd <- sample(which(round(df$clim_sd, digits = 3) == 1.116 & df$clim_auto == 0), 1)
  high_sd <- sample(which(round(df$clim_sd, digits = 3) == 2 & df$clim_auto == 0), 1)
  
  
  values <- rbind(`0.2` = mats[[low_sd]]$mats %>% as.data.frame,
                  `1.1` = mats[[mid_sd]]$mats %>% as.data.frame,
                  `2.0` = mats[[high_sd]]$mats %>% as.data.frame) %>% 
    tibble::rownames_to_column(., "SD") %>% 
    mutate(SD = regmatches(SD, regexpr("^.{3}", SD))) %>% 
    pivot_longer(cols = -SD, names_to = "cell", values_to = "value")
  
  ggplot(values, aes(value, fill = SD, colour = SD)) + 
    geom_histogram(position = "dodge") + facet_wrap(vars(cell), scales = "free") +
    scale_y_log10() + ggtitle(title)
}


## print cellvalues for all signalstrengths -------------------- FOR LOOP OR FUNCTION CREATES A PDF THAT ISN'T READABLE, BUT MANUALLY WORKS JUST FINE?? ------------------
pdf(file.path(output_dir, "cell_value_plots", paste0("sigstrength_1.pdf")))

plot_cell(readRDS(file.path(output_dir, "rds", "mpm_1_auto.RDS")), "autocorrelation")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_1_cov.RDS")), "covariation")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_1_lag.RDS"))$survival, "survival lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_1_lag.RDS"))$growth, "growth lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_1_lag.RDS"))$none, "neither vital rate lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_1_lagfp.RDS"))$Umatrix, "U lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_1_lagfp.RDS"))$Fmatrix, "F lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_1_lagfp.RDS"))$none, "Neither submatrices lagged")

dev.off()


pdf(file.path(output_dir, "cell_value_plots", paste0("sigstrength_0.5.pdf")))

plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.5_auto.RDS")), "autocorrelation")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.5_cov.RDS")), "covariation")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.5_lag.RDS"))$survival, "survival lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.5_lag.RDS"))$growth, "growth lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.5_lag.RDS"))$none, "neither vital rate lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.5_lagfp.RDS"))$Umatrix, "U lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.5_lagfp.RDS"))$Fmatrix, "F lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.5_lagfp.RDS"))$none, "Neither submatrices lagged")

dev.off()

pdf(file.path(output_dir, "cell_value_plots", paste0("sigstrength_0.25.pdf")))

plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.25_auto.RDS")), "autocorrelation")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.25_cov.RDS")), "covariation")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.25_lag.RDS"))$survival, "survival lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.25_lag.RDS"))$growth, "growth lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.25_lag.RDS"))$none, "neither vital rate lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.25_lagfp.RDS"))$Umatrix, "U lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.25_lagfp.RDS"))$Fmatrix, "F lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.25_lagfp.RDS"))$none, "Neither submatrices lagged")

dev.off()


pdf(file.path(output_dir, "cell_value_plots", paste0("sigstrength_0.05.pdf")))

plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.05_auto.RDS")), "autocorrelation")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.05_cov.RDS")), "covariation")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.05_lag.RDS"))$survival, "survival lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.05_lag.RDS"))$growth, "growth lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.05_lag.RDS"))$none, "neither vital rate lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.05_lagfp.RDS"))$Umatrix, "U lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.05_lagfp.RDS"))$Fmatrix, "F lagged")
plot_cell(readRDS(file.path(output_dir, "rds", "mpm_0.05_lagfp.RDS"))$none, "Neither submatrices lagged")

dev.off()
