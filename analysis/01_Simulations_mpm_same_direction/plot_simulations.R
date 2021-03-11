library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

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


## lambda summary plot -------------------------------------

sum_plot <- function(n, save = T){
  
  ## I tried using pmap and the object lists, but that gave:
  # "no applicable method for 'mutate_' applied to an object of class "list" " 
  # couldn't figured it out, and this way works too, but at some point I should get more familiar with pmap....
  a <- auto[[n]]
  c <- cov[[n]]
  p <- lag[[n]]
  f <- lagpf[[n]]
  sig.strength <- sig.strength[[n]]
  
  ### Climate variables
  clim_sd <- rep(seq(from = 0.01, to = 2, length.out = 10), 90)
  clim_corr <- rep(rep(c(-0.9,0,0.9), each = 10), 30)
  ## Plot
  auto_p <- ggplot(a) + geom_smooth(aes(x = clim_sd, y = lambda, linetype = as.factor(clim_auto)), colour = "grey30") + 
    labs(linetype = "Climate \nautocorrelation", title = "Autocorrelation \nin growth") 
  
  cov_df <- c %>% mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))
  cov_p <- ggplot(cov_df) + geom_smooth(aes(x = clim_sd, y = lambda, linetype = as.factor(auto_cat)), colour = "grey30") + 
    labs(linetype = "Climate \nautocorrelation", title = "Covariance in \nsurvival & growth")
  
  lag_df <- lapply(p, function(x) x %>% mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
    bind_rows(., .id = "type")
  
  lag_p <- ggplot(lag_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    labs(colour = "Lag type", title = "Lagged climate in growth or survival") +
    facet_grid(cols = vars(auto_cat))  
  
  lagpf_df <- lapply(f, function(x) x %>% mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))) %>%
    bind_rows(., .id = "type")
  
  lagpf_p <- ggplot(lagpf_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    labs(colour = "Lag type", title = "Lagged climate in P or F") +
    facet_grid(cols = vars(auto_cat)) + scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
  
  plot_AC <- auto_p + cov_p + plot_layout(guides = "collect") + plot_annotation(title = paste("Climate signal strenth =", sig.strength))
  plot_lag <-  (lag_p / lagpf_p) + plot_layout(guides = "collect") + plot_annotation(title = paste("Climate signal strenth =", sig.strength))
  
  if(save) {
    ggsave(filename = file.path(output_dir, "lambda_plots", paste0("auto_cov_plot_signal_strength_", sig.strength, ".png")), 
           plot_AC)
    ggsave(filename = file.path(output_dir, "lambda_plots", paste0("lag_plot_signal_strength_", sig.strength, ".png")), 
           plot_lag)
  }
  
  return(list(plot_AC = plot_AC, 
              plot_lag = plot_lag))
}


auto <- lapply(auto_f, function(x) lapply(x, function(y) y$df) %>% bind_rows)  
cov <- lapply(cov_f, function(x) lapply(x, function(y) y$df) %>% bind_rows)
lag <- lapply(lag_f, function(x) lapply(x, function(y) lapply(y, function(z) z$df) %>% bind_rows))
lagpf <- lapply(lagpf_f, function(x) lapply(x, function(y) lapply(y, function(z) z$df) %>% bind_rows))

# m <- purrr::pmap(list(auto, cov, lag, lagpf, sig.strength), sum_plot)  ----- see commented section in function. Get a mutate error, even though manual run through works just fine

m <- lapply(as.list(c(1:4)), sum_plot)

## Save histograms of cell values 
## functions to create dataset formate needed, plot cell values and save to pdf --------------------------

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


# ## print for all signalstrengths -------------------- FOR LOOP CREATES A PDF THAT ISN'T READABLE, BUT MANUALLY WORKS JUST FINE?? ------------------
# for(i in c(1:length(sig.strength))) {
#   pdf(file.path(output_dir, "cell_value_plots", paste0("sigstrength_", sig.strength[[i]], ".pdf")))
#   
#   plot_cell(auto_f[[i]], "autocorrelation")
#   plot_cell(cov_f[[i]], "covariation")
#   
#   plot_cell(lag_f[[i]]$survival, "survival lagged")
#   plot_cell(lag_f[[i]]$growth, "growth lagged")
#   plot_cell(lag_f[[i]]$none, "neither vital rate lagged")
#   
#   plot_cell(lagpf_f[[i]]$Pkernel, "P lagged")
#   plot_cell(lagpf_f[[i]]$Fkernel, "F lagged")
#   plot_cell(lagpf_f[[i]]$none, "Neither submatrices lagged")
#   
#   dev.off()
# 
# }

i = 1
pdf(file.path(output_dir, "cell_value_plots", paste0("sigstrength_", sig.strength[[i]], ".pdf")))

plot_cell(auto_f[[i]], "autocorrelation")
plot_cell(cov_f[[i]], "covariation")
plot_cell(lag_f[[i]]$survival, "survival lagged")
plot_cell(lag_f[[i]]$growth, "growth lagged")
plot_cell(lag_f[[i]]$none, "neither vital rate lagged")
plot_cell(lagpf_f[[i]]$Pkernel, "P lagged")
plot_cell(lagpf_f[[i]]$Fkernel, "F lagged")
plot_cell(lagpf_f[[i]]$none, "Neither submatrices lagged")

dev.off()


i=2
pdf(file.path(output_dir, "cell_value_plots", paste0("sigstrength_", sig.strength[[i]], ".pdf")))

plot_cell(auto_f[[i]], "autocorrelation")
plot_cell(cov_f[[i]], "covariation")
plot_cell(lag_f[[i]]$survival, "survival lagged")
plot_cell(lag_f[[i]]$growth, "growth lagged")
plot_cell(lag_f[[i]]$none, "neither vital rate lagged")
plot_cell(lagpf_f[[i]]$Pkernel, "P lagged")
plot_cell(lagpf_f[[i]]$Fkernel, "F lagged")
plot_cell(lagpf_f[[i]]$none, "Neither submatrices lagged")

dev.off()


i=3
pdf(file.path(output_dir, "cell_value_plots", paste0("sigstrength_", sig.strength[[i]], ".pdf")))

plot_cell(auto_f[[i]], "autocorrelation")
plot_cell(cov_f[[i]], "covariation")
plot_cell(lag_f[[i]]$survival, "survival lagged")
plot_cell(lag_f[[i]]$growth, "growth lagged")
plot_cell(lag_f[[i]]$none, "neither vital rate lagged")
plot_cell(lagpf_f[[i]]$Pkernel, "P lagged")
plot_cell(lagpf_f[[i]]$Fkernel, "F lagged")
plot_cell(lagpf_f[[i]]$none, "Neither submatrices lagged")

dev.off()


i=4
pdf(file.path(output_dir, "cell_value_plots", paste0("sigstrength_", sig.strength[[i]], ".pdf")))

plot_cell(auto_f[[i]], "autocorrelation")
plot_cell(cov_f[[i]], "covariation")
plot_cell(lag_f[[i]]$survival, "survival lagged")
plot_cell(lag_f[[i]]$growth, "growth lagged")
plot_cell(lag_f[[i]]$none, "neither vital rate lagged")
plot_cell(lagpf_f[[i]]$Pkernel, "P lagged")
plot_cell(lagpf_f[[i]]$Fkernel, "F lagged")
plot_cell(lagpf_f[[i]]$none, "Neither submatrices lagged")

dev.off()
