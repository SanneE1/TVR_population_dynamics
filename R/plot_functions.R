### Plot functions

plot_w_vrcov_function <- function(id) {
  df %>%
    filter(lh_id == id & auto_cat %in% c("-0.6", "0", "0.6") & sig.strength == 0.5) %>%
    mutate(vr_cov = factor(vr_cov, levels = c("positive", "negative"), labels = c("A: Positive covariation", "B: Negative covariation"))) %>%
    ggplot(.) +
    geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(lag_type), linetype = auto_cat), se = F) +
    scale_colour_manual(name = "Simulation type",
                        values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                        labels = c("none" = "control", "Umatrix" = "TVR")) +
    scale_linetype(name = "Autocorrelation") +
    ylab("stochastic log lambda") + xlab(~ paste(sigma[c], " of environmental sequence")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          legend.margin=margin(),
          text = element_text(size = 14),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12)) +
    facet_row(vars(vr_cov))
}


plot_one_function <- function(id) {
  df %>%
    filter(lh_id == id & auto_cat %in% c("-0.6", "0", "0.6") & sig.strength == 0.5 & vr_cov == "positive") %>%
    ggplot(.) +
    geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(lag_type), linetype = auto_cat), se = F) +
    scale_colour_manual(name = "Simulation type",
                        values = c("Umatrix" = "#E69F00", "none" = "#0072B2"),
                        labels = c("none" = "control", "Umatrix" = "TVR")) +
    scale_linetype(name = "Autocorrelation") +
    xlab("") + ylab("") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          legend.margin=margin(),
          text = element_text(size = 14),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12)) +
    theme(legend.position = "bottom")
  
}
