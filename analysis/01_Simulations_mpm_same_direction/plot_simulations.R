rm(list=ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(popbio)

output_dir <- "results/01_Simulations_mpm_same_directions/"

### Load results
location <- file.path(output_dir, "rds")

auto_f <- lapply(list.files(location, pattern = "auto.RDS", full.names = T), readRDS)
cov_f <- lapply(list.files(location, pattern = "cov.RDS", full.names = T), readRDS)
lag_f <- lapply(list.files(location, pattern = "lag.RDS", full.names = T), readRDS)
lagpf_f <- lapply(list.files(location, pattern = "lagfp.RDS", full.names = T), readRDS)
sig.strength <- regmatches(list.files(location, pattern = "auto.RDS"), 
                           regexec("mpm\\_(.+)\\_", list.files(location, pattern = "auto.RDS"))) %>% 
  lapply(., function(x) x[2])

##  -------------------------------------
##  lambda summary plot
##  -------------------------------------

sum_plot <- function(n, save = T){
  
  ## I tried using pmap and the object lists, but that gave:
  # "no applicable method for 'mutate_' applied to an object of class "list" " 
  # couldn't figured it out, and this way works too, but at some point I should get more familiar with pmap....
  a <- auto[[n]]
  c <- cov[[n]]
  p <- lag[[n]]
  f <- lagpf[[n]]
  ss <- sig.strength[[n]]
  
  ### Climate variables
  clim_sd <- rep(seq(from = 0.01, to = 2, length.out = 10), 90)
  clim_corr <- rep(rep(c(-0.9,0,0.9), each = 10), 30)
  ## Plot
  auto_p <- ggplot(a) + geom_smooth(aes(x = clim_sd, y = lambda, linetype = as.factor(clim_auto)), colour = "grey30") + 
    labs(linetype = "Climate \nautocorrelation", title = "Autocorrelation \nin growth") +
  xlab("Climate standard deviation") + ylab("log lambda")
  
  cov_df <- c %>% mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))
  cov_p <- ggplot(cov_df) + geom_smooth(aes(x = clim_sd, y = lambda, linetype = as.factor(auto_cat)), colour = "grey30") + 
    labs(linetype = "Climate \ncovariance", title = "Covariance in \nsurvival & growth") +
  xlab("Climate standard deviation") + ylab("log lambda")
  
  lag_df <- lapply(p, function(x) x %>% mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
    bind_rows(., .id = "type")
  
  lag_p <- ggplot(lag_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type))) + 
    geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(type)), position = position_dodge(width = 0.1)) + 
    labs(colour = "Lag type", title = "Lagged climate in growth or survival") +
    xlab("Climate standard deviation") + ylab("log lambda") +
    facet_grid(cols = vars(auto_cat), labeller = as_labeller(c("-0.9" = "Autocorrelation = -0.9",
                                                               "0" = "Autocorrelation = 0",
                                                               "0.9" = "Autocorrelation = 0.9")))
  
  lagpf_df <- lapply(f, function(x) x %>% mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
    bind_rows(., .id = "type")
  
  lagpf_p <- ggplot(lagpf_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(type)), position = position_dodge(width = 0.1))+ 
    labs(colour = "Lag type", title = "Lagged climate in U or F matrix") + 
    xlab("Climate standard deviation") + ylab("log lambda") +
    facet_grid(cols = vars(auto_cat), labeller = as_labeller(c("-0.9" = "Autocorrelation = -0.9",
                                                               "0" = "Autocorrelation = 0",
                                                               "0.9" = "Autocorrelation = 0.9"))) + 
    scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) 
  
  plot_AC <- auto_p + cov_p + plot_layout(guides = "collect") + plot_annotation(title = paste("Climate signal strenth =", ss))
  plot_lag <-  (lag_p / lagpf_p) + plot_layout(guides = "collect") + plot_annotation(title = paste("Climate signal strenth =", ss))
  
  
  lagf <- lagpf_df %>% filter(type != "Fmatrix") %>%
    ggplot(.) + 
    geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(type)), position = position_dodge(width = 0.1))+ 
    labs(colour = "Lag type", title = "Lagged climate in U") + 
    xlab("Climate standard deviation") + ylab("log lambda") +
    facet_grid(cols = vars(auto_cat), labeller = as_labeller(c("-0.9" = "Autocorrelation = -0.9",
                                                               "0" = "Autocorrelation = 0",
                                                               "0.9" = "Autocorrelation = 0.9"))) + 
    scale_colour_manual(values = c("#E7B800", "#FC4E07")) 
  
  
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
    labs(fill = "Climate \nautocorrelation", title = "Relative difference", subtitle = "between control and simulations with lag in U-matrix")
  
  ggsave(filename = file.path(output_dir, "lambda_plots", paste0("UF_relative_diff_", ss, ".png")),
         diff, width = 8.22, height = 5.53, units = "in")
  

}


auto <- lapply(auto_f, function(x) lapply(x, function(y) y$df) %>% bind_rows)  
cov <- lapply(cov_f, function(x) lapply(x, function(y) y$df) %>% bind_rows)
lag <- lapply(lag_f, function(x) lapply(x, function(y) lapply(y, function(z) z$df) %>% bind_rows))
lagpf <- lapply(lagpf_f, function(x) lapply(x, function(y) lapply(y, function(z) z[[1]]) %>% Filter(Negate(anyNA), .) %>% bind_rows))

m <- lapply(as.list(c(1:4)), sum_plot)



##  -------------------------------------
##  Plot a lambda time sequence 
##  -------------------------------------
l_seq0 <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(lagpf_f[[3]]$Umatrix[[20]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                     none = lapply(as.list(as.data.frame(t(lagpf_f[[3]]$none[[20]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                     time = c(1:5002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

l_seqneg <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(lagpf_f[[3]]$Umatrix[[10]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                      none = lapply(as.list(as.data.frame(t(lagpf_f[[3]]$none[[10]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                      time = c(1:5002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

l_seqpos <- data.frame( Umatrix = lapply(as.list(as.data.frame(t(lagpf_f[[3]]$Umatrix[[30]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        none = lapply(as.list(as.data.frame(t(lagpf_f[[3]]$none[[30]]$mats))), function(x) log(lambda(matrix2(x, byrow = F)))) %>% bind_rows() %>% t,
                        time = c(1:5002)
) %>% pivot_longer(cols = -time, names_to = "type", values_to = "lambda")

### zoom in on 35 years
time_series0 <- ggplot(l_seq0) + 
  geom_line(aes(y = lambda, x = time, colour = type)) + 
  coord_cartesian(xlim = c(1095,1130)) + ylab("log lambda") + 
  labs(subtitle = "0 autocorrelation, standard deviation = 2")

time_seriesneg <- ggplot(l_seqneg) + 
  geom_line(aes(y = lambda, x = time, colour = type)) + 
  coord_cartesian(xlim = c(1095,1130)) + ylab("log lambda") + 
  labs(subtitle = "-0.9 autocorrelation, standard deviation = 2")

time_seriespos <- ggplot(l_seqpos) + 
  geom_line(aes(y = lambda, x = time, colour = type)) + 
  coord_cartesian(xlim = c(1095,1130)) + ylab("log lambda") + 
  labs(subtitle = "0.9 autocorrelation, standard deviation = 2")

ggsave(filename = file.path(output_dir, "lambda_timeseries", "timeseries_0.tiff"), time_series0,width = 9.22, height = 6.5)
ggsave(filename = file.path(output_dir, "lambda_timeseries", "timeseries_neg.tiff"), time_seriesneg, width = 9.22, height = 6.5)
ggsave(filename = file.path(output_dir, "lambda_timeseries", "timeseries_pos.tiff"), time_seriespos, width = 9.22, height = 6.5)

time_series <- time_seriespos + time_series0 + time_seriesneg + plot_layout(nrow = 3, guides = "collect") + 
  plot_annotation(title = "lambda time series")

ggsave(filename = file.path(output_dir, "lambda_timeseries", "summary_timeseries.tiff"), 
       plot = time_series, width = 9.22, height = 6.5)

int_ann_0 <- l_seq0 %>% 
  group_by(type) %>% 
  arrange(time) %>% 
  mutate(diff = (lambda %>% exp) - (lag(lambda, default = first(lambda)) %>% exp)) %>%
  ggplot(.) +
  geom_density(aes(x = diff, fill = type, colour = type), alpha = 0.5) +
  xlab("interannual difference") + coord_cartesian(c(-3, 3))

int_ann_neg <- l_seqneg %>% 
  group_by(type) %>% 
  arrange(time) %>% 
  mutate(diff = (lambda %>% exp) - (lag(lambda, default = first(lambda)) %>% exp)) %>%
  ggplot(.) +
  geom_density(aes(x = diff, fill = type, colour = type), alpha = 0.5)  +
  xlab("interannual difference") + coord_cartesian(c(-3, 3))

int_ann_pos <- l_seqpos %>% 
  group_by(type) %>% 
  arrange(time) %>% 
  mutate(diff = (lambda %>% exp) - (lag(lambda, default = first(lambda)) %>% exp)) %>%
  ggplot(.) +
  geom_density(aes(x = diff, fill = type, colour = type), alpha = 0.5) +
  xlab("interannual difference") + coord_cartesian(c(-3, 3))

time_diff <-  (time_seriespos + time_series0 + time_seriesneg) * ylim(c(-7.5,2)) + 
  int_ann_pos + int_ann_0 + int_ann_neg + 
  plot_layout(nrow = 3, guides = "collect", byrow = F, widths = c(4,1)) + 
  plot_annotation(title = "lambda over time") 

ggsave(filename = file.path(output_dir, "lambda_timeseries", "timeseries&diff.tiff"), 
       plot = time_diff, width = 9.22, height = 6.5)



##  -------------------------------------
##  Save histograms of cell values 
##  -------------------------------------
## functions to create dataset formate needed, plot cell values and save to pdf --------------------------

plot_cell <- function(mats, title) {
  
  df <- data.frame(clim_sd = rep(seq(from = 0.01, to = 2, length.out = 10), 90),
                   clim_auto = rep(rep(c(-0.9,0,0.9), each = 10), 30)
  )
  
  # Randomly select a sequence with a low and high sd
  # low_sd <- sample(which(round(df$clim_sd, digits = 3) == 0.231 & df$clim_auto == 0), 1)
  mid_sd <- sample(which(round(df$clim_sd, digits = 3) == 1.116 & df$clim_auto == 0), 1)
  high_sd <- sample(which(round(df$clim_sd, digits = 3) == 2 & df$clim_auto == 0), 1)
  
  
  values <- rbind(#`0.2` = mats[[low_sd]]$mats %>% as.data.frame,
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
i = 1
pdf(file.path(output_dir, "cell_value_plots", paste0("sigstrength_", sig.strength[[i]], ".pdf")))

plot_cell(auto_f[[i]], "autocorrelation")
plot_cell(cov_f[[i]], "covariation")
plot_cell(lag_f[[i]]$survival, "survival lagged")
plot_cell(lag_f[[i]]$growth, "growth lagged")
plot_cell(lag_f[[i]]$none, "neither vital rate lagged")
plot_cell(lagpf_f[[i]]$Umatrix, "U lagged")
plot_cell(lagpf_f[[i]]$Fmatrix, "F lagged")
plot_cell(lagpf_f[[i]]$none, "Neither submatrices lagged")

dev.off()


i=2
pdf(file.path(output_dir, "cell_value_plots", paste0("sigstrength_", sig.strength[[i]], ".pdf")))

plot_cell(auto_f[[i]], "autocorrelation")
plot_cell(cov_f[[i]], "covariation")
plot_cell(lag_f[[i]]$survival, "survival lagged")
plot_cell(lag_f[[i]]$growth, "growth lagged")
plot_cell(lag_f[[i]]$none, "neither vital rate lagged")
plot_cell(lagpf_f[[i]]$Umatrix, "U lagged")
plot_cell(lagpf_f[[i]]$Fmatrix, "F lagged")
plot_cell(lagpf_f[[i]]$none, "Neither submatrices lagged")

dev.off()


i=3
pdf(file.path(output_dir, "cell_value_plots", paste0("sigstrength_", sig.strength[[i]], ".pdf")))

plot_cell(auto_f[[i]], "autocorrelation")
plot_cell(cov_f[[i]], "covariation")
plot_cell(lag_f[[i]]$survival, "survival lagged")
plot_cell(lag_f[[i]]$growth, "growth lagged")
plot_cell(lag_f[[i]]$none, "neither vital rate lagged")
plot_cell(lagpf_f[[i]]$Umatrix, "U lagged")
plot_cell(lagpf_f[[i]]$Fmatrix, "F lagged")
plot_cell(lagpf_f[[i]]$none, "Neither submatrices lagged")

dev.off()


i=4
pdf(file.path(output_dir, "cell_value_plots", paste0("sigstrength_", sig.strength[[i]], ".pdf")))

plot_cell(auto_f[[i]], "autocorrelation")
plot_cell(cov_f[[i]], "covariation")
plot_cell(lag_f[[i]]$survival, "survival lagged")
plot_cell(lag_f[[i]]$growth, "growth lagged")
plot_cell(lag_f[[i]]$none, "neither vital rate lagged")
plot_cell(lagpf_f[[i]]$Umatrix, "U lagged")
plot_cell(lagpf_f[[i]]$Fmatrix, "F lagged")
plot_cell(lagpf_f[[i]]$none, "Neither submatrices lagged")

dev.off()
