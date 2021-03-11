library(tidyverse)
library(ggplot2)
library(pbapply)
library(popbio)

output_dir <- "results/06_COMPADRE_studies/actual_sd_auto/"

source_lines <- function(file, lines){
  source(textConnection(readLines(file)[lines]))
}

### Get the climate sd and autocorrelation for the populations with Lat/Lon
climate <- read.csv("data/All_populations.csv") %>% 
  pivot_wider(values_from = value, names_from = variable) %>%
  mutate(tmean = (tmin + tmax)/2) %>%
  group_by(SpeciesAuthor, MatrixPopulation, month) %>%
  mutate(prec_a = scale(prec),
         tmean_a = scale(tmean)) %>%
  ungroup(month) %>%
  group_by(year, .add=T) %>%
  mutate(an_prec = mean(prec_a),
         an_tmean = mean(tmean_a)) %>%
  ungroup(year) %>%
  summarise(prec_sd = sd(an_prec),
            prec_auto = acf(an_prec, plot=F, na.action = na.pass)$acf[2],
            tmean_sd = sd(an_tmean),
            tmean_auto = acf(an_tmean, plot=F, na.action = na.pass)$acf[2])

### Load selected species
species <- read.csv("data/species_authors.csv") %>% filter(!is.na(Lat))

### create main file/dataframe
df <- left_join(species, climate)

df$prcp_diffU_0.1 <- NA
df$prcp_diffU_0.5 <- NA
df$prcp_diffU_1 <- NA

df$temp_diffU_0.1 <- NA
df$temp_diffU_0.5 <- NA
df$temp_diffU_1 <- NA


### Run lambda simulations using matrices, sd and auto of specific populations
### 30 repetitions for 3 different climate signal strengths

load("data/COMPADRE_v.X.X.X.4.RData")

for(sp in c(1:nrow(df))) {
  
  i = df$SpeciesAuthor[sp]
  j = df$MatrixPopulation[sp]
  
  source_lines("analysis/06_COMPADRE_studies/05_simulate_lambda.R", c(30:137, 147:169))
  
  
  ## Create climate sequences with population's sd and autocorrelation 
  lag_prcp <- lapply(as.list(rep(c(1,0.5,0.1), 30)), function(x) 
    create_seq(clim_sd = df$prec_sd[sp], clim_corr = df$prec_auto[sp], sig.strength = x) %>% filter(complete.cases(.)))
  
  lag_temp <- lapply(as.list(rep(c(1,0.5,0.1), 30)), function(x) 
    create_seq(clim_sd = df$tmean_sd[sp], clim_corr = df$tmean_auto[sp], sig.strength = x) %>% filter(complete.cases(.)))
  
  # Run simulations for precipitation
  lag_clim <- lag_prcp
  
  ## Only run U and none lagged simulation
  source_lines("analysis/06_COMPADRE_studies/05_simulate_lambda.R", c(171:180, 189:196))
  
  df1 <- left_join(
    lag_u %>% group_by(sig.strength) %>%
      summarise(Ulambda = mean(lambda)),
    lag_n %>% group_by(sig.strength) %>%
      summarise(Nlambda = mean(lambda))
  ) %>%
    mutate(diffU = (Ulambda - Nlambda)/abs(Nlambda))
  
  plot_df <- rbind(
    lag_u %>% select(lambda, sig.strength) %>% mutate(type = "Umatrix"),
    lag_n %>% select(lambda, sig.strength) %>% mutate(type = "None")
  ) %>% mutate(climate = "Precipitation")
  
  # Run simulation for temperature
  lag_clim <- lag_temp
  
  ## Only run U and none lagged simulation
  source_lines("analysis/06_COMPADRE_studies/05_simulate_lambda.R", c(171:180, 189:196))
  
  df2 <- left_join(
    lag_u %>% group_by(sig.strength) %>%
      summarise(Ulambda = mean(lambda)),
    lag_n %>% group_by(sig.strength) %>%
      summarise(Nlambda = mean(lambda))
  ) %>%
    mutate(diffU = (Ulambda - Nlambda)/abs(Nlambda))
  
  df3 <- cbind(pivot_wider(df1, diffU, names_from = sig.strength, names_prefix = "prcp_diffU_", values_from = diffU),
               pivot_wider(df2, diffU, names_from = sig.strength, names_prefix = "temp_diffU_", values_from = diffU))
  
  df[which(df$SpeciesAuthor == i & df$MatrixPopulation == j), c(15:20)] <- df3
  
  plot_df <- rbind(plot_df,
                   rbind(lag_u %>% select(lambda, sig.strength) %>% mutate(type = "Umatrix"),
                         lag_n %>% select(lambda, sig.strength) %>% mutate(type = "None")
                   ) %>% mutate(climate = "Temperature"))
  
  
  facet_text <- data.frame(
    label = c(paste("sd =", round(df$prec_sd[sp],3), "autocor =", round(df$prec_auto[sp],3)), 
              paste("sd =", round(df$tmean_sd[sp],3), "autocor =", round(df$tmean_auto[sp],3))),
    climate = c("Precipitation", "Temperature"),
    x = c(2,2)
  )
  
  plot <- ggplot(plot_df) + geom_boxplot(aes(x = as.character(sig.strength), y = lambda, fill = type)) + 
    facet_grid(cols = vars(climate)) + 
    geom_text(data = facet_text,
              aes(x = x, y = Inf, label = label), vjust = 1.2) +
    xlab("climate signal strength") + ylab("log lambda")
  
  ggsave(file.path(output_dir, paste0("lambda_comparison_", i, "_", j, ".tiff")), plot = plot)
  
}








