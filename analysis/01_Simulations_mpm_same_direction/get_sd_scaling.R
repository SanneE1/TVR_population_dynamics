library(tidyverse)
library(ggplot2)
library(pbapply)
library(popbio)
library(boot)
library(parallel)

# ---------------------------------------------------------------
# load COMPADRE
# ---------------------------------------------------------------

load("data/COMPADRE_v.6.21.1.0.RData")

# ---------------------------------------------------------------
# get sApropos collapsed matrices
# ---------------------------------------------------------------

col_mat <- read.csv("data/vitalRatesCollapsed.csv") %>% 
  subset( !( survJuv > 1 | survAdu > 1 |
               retrAdu > 1 | progJuv > 1 | 
               cloJuv  > 1 | cloJuvJuv > 1 |
               cloAdu  > 1 | cloAduAdu > 1 ) ) %>%
  mutate(U11 = survJuv * (1-progJuv),
         U12 = retrAdu * survAdu,
         U21 = progJuv * survJuv,
         U22 = survAdu * (1-retrAdu),
         F11 = fecJuv * fecJuvJuv,   # fecJuvJuv = proportion of fecJuv that goes into F[1,1]
         F12 = fecAdu * (1-fecAduAdu),
         F21 = fecJuv * (1-fecJuvJuv),
         F22 = fecAdu * fecAduAdu 
           ) %>%
  # mutate(survJuv = logit(survJuv),
  #        progJuv = logit(progJuv),
  #        survAdu = logit(survAdu),
  #        fecAdu = log(fecAdu)) %>%
  group_by(SpeciesAuthor, MatrixPopulation) %>%
  summarise(survJuv_mean = mean(survJuv, na.rm = T),
            survJuv_sd = sd(survJuv, na.rm = T),
            progJuv_mean = mean(progJuv, na.rm = T),
            progJuv_sd = sd(progJuv, na.rm = T),
            survAdu_mean = mean(survAdu, na.rm = T),
            survAdu_sd = sd(survAdu, na.rm = T),
            fecAdu_mean = mean(fecAdu, na.rm = T),
            fecAdu_sd = sd(fecAdu, na.rm = T),
            U11_mean = mean(U11, na.rm = T),
            U11_sd = sd(U11, na.rm = T),
            U12_mean = mean(U12, na.rm = T),
            U12_sd = sd(U12, na.rm = T),
            U21_mean = mean(U21, na.rm = T),
            U21_sd = sd(U21, na.rm = T),
            U22_mean = mean(U22, na.rm = T),
            U22_sd = sd(U22, na.rm = T),
            F11_mean = mean(F11, na.rm = T),
            F11_sd = sd(F11, na.rm = T),
            F12_mean = mean(F12, na.rm = T),
            F12_sd = sd(F12, na.rm = T),
            F21_mean = mean(F21, na.rm = T),
            F21_sd = sd(F21, na.rm = T),
            F22_mean = mean(F22, na.rm = T),
            F22_sd = sd(F22, na.rm = T)) %>%
  filter(across(c(3:last_col()), ~ !is.na(.x)))
  

ggplot(col_mat) + geom_point(aes(x = survAdu_sd, y = fecAdu_sd, colour = progJuv_sd))
ggplot(col_mat) + geom_point(aes(x = survAdu_mean, y = fecAdu_mean, colour = progJuv_mean))



ggplot(col_mat) + geom_point(aes(x = survJuv_mean, y = survJuv_sd))
ggplot(col_mat) + geom_point(aes(x = survAdu_mean, y = survAdu_sd))
ggplot(col_mat) + geom_point(aes(x = progJuv_mean, y = progJuv_sd))
ggplot(col_mat) + geom_point(aes(x = fecAdu_mean, y = fecAdu_sd)) + coord_cartesian(xlim = c(0,25), ylim = c(0,25))

summary(lm(survJuv_sd ~ survJuv_mean, data = col_mat))
summary(lm(survAdu_sd ~ survAdu_mean, data = col_mat))
summary(lm(fecAdu_sd ~ fecAdu_mean, data = col_mat))


summary(lm(fecAdu_sd ~ survAdu_sd + survJuv_sd, data = col_mat))
summary(lm(survJuv_sd ~ fecAdu_sd, data = col_mat))

summary(lm(survAdu_sd ~ fecAdu_sd + survJuv_sd, data = col_mat))
summary(lm(survAdu_sd ~ survAdu_mean, data = col_mat))



#-----------------------------------------------------------
# simulate lambda with collapsed mpm for selected species
#-----------------------------------------------------------
# for comparison of COMPADRE results
n_it = 5000

# get cell values
df <- col_mat %>% filter(SpeciesAuthor == "Lesquerella_ovalifolia"|
                     SpeciesAuthor =="Cirsium_pitcheri_8"|
                     (SpeciesAuthor =="Dicerandra_frutescens" & MatrixPopulation == "Oak-hickory scrub 3 (Pop 19)")
) %>% select(SpeciesAuthor, MatrixPopulation, contains("U", ignore.case = F), contains("F", ignore.case = F))

cell_values <- pivot_longer(df, cols = -c("SpeciesAuthor", "MatrixPopulation"),
                             names_pattern = "^(.{3})_(.{2,4})", names_to = c("cell", "type")) %>%
  pivot_wider(names_from = "type", values_from = "value") %>% split(., f = .$SpeciesAuthor)

# set up climate sequences
create_seq <- function(n_it, clim_sd, clim_auto, lag) { 
  for(n in c(1:(n_it+lag))){ 
    if(n == 1) {
      seq <- rnorm(1)
    } else {
      seq[n] <- clim_auto * seq[n-1] + rnorm(1)
    }
  }
  seq <- scale(seq) * clim_sd
  lagged <- c(rep(NA, lag), seq)
  recent <- c(seq, rep(NA, lag))
  df <- data.frame(recent = recent,
                   lagged = lagged)
  return(df)
}

clim_sd <- rep(seq(from = 0.01, to = 2, length.out = 10), 90)
clim_auto <- rep(rep(c(-0.9,0,0.9), each = 10), 30)

lag_clim <- lapply(as.list(c(1:900)), function(x) 
  create_seq(n_it, clim_sd = clim_sd[x], clim_auto = clim_auto[x], lag = 1))


# Stocastic mpm -----------------------------------

st.lamb <- function(env_U, env_F, 
                    clim_sd, clim_auto, sig.strength) {
  
  n_it = length(env_U)
  
  env <- list(U_clim = env_U,
              F_clim = env_F,
              clim_sd = rep(clim_sd, n_it),
              sig.strength = rep(sig.strength, n_it))
  
  ### Get all mpm's
  mats <- pmap(env, mpm) 
  
  
  df <- data.frame(lambda = stoch.growth.rate(mats, maxt = n_it, verbose = F)$sim,
                   clim_sd = clim_sd,
                   clim_auto = clim_auto,
                   sig.strength = sig.strength)
  
  return(df)
}

#-----------------------------------------------------------
# Simulate lambda
#-----------------------------------------------------------

for(sp in c(1:length(cell_values))) {
  
  Ucell_values <- cell_values[[sp]] %>% filter(., grepl("^U", cell))
  Fcell_values <- cell_values[[sp]] %>% filter(., grepl("^F", cell))
  
  # create species specific mpm
  mpm <- function(U_clim, F_clim, sig.strength, clim_sd) {
    
    devU <- U_clim * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + 
      rnorm(length(Ucell_values$mean), mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd)
    pU <- pnorm(devU, mean = 0, sd = clim_sd) 
    Umat <- qbeta(pU, 
                  (((Ucell_values$mean*(1-Ucell_values$mean))/(Ucell_values$sd * clim_sd)^2) - 1) * Ucell_values$mean,
                  (((Ucell_values$mean*(1-Ucell_values$mean))/(Ucell_values$sd * clim_sd)^2) - 1) * (1 - Ucell_values$mean)) %>% 
      replace_na(., 0)
    
    devF <- F_clim * (sqrt(clim_sd^2 * sig.strength)/clim_sd) + 
      rnorm(length(Fcell_values$mean), mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd)
    pF <- pnorm(devF, mean = 0, sd = clim_sd) 
    Fmat <- qgamma(pF, 
                   (Fcell_values$mean^2)/(Fcell_values$sd * clim_sd)^2, 
                   (Fcell_values$mean)/(Fcell_values$sd * clim_sd)^2) %>% 
      replace_na(., 0)
    
    Amat <- Umat + Fmat
    mpm <- matrix2(Amat)
    
    return(mpm)  
  }
  
  # Set up parallel runs
  no_cores <- detectCores()
  cl <- makeCluster(no_cores - 2)
  ## export libraries to workers
  clusterEvalQ(cl, c(library(popbio), library(tidyverse), library(purrr), library(boot)))
  clusterExport(cl, c("st.lamb", "mpm", "lag_clim",
                      "Ucell_values", "Fcell_values"))
  
  
  plag_u <- pblapply(cl = cl,
                     lag_clim, 
                     function(x) st.lamb(env_U = x$lagged,
                                         env_F = x$recent,
                                         clim_sd = sd(x$recent, na.rm = T),
                                         clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                         sig.strength = 1)
  ) 
  
  plag_n <- pblapply(cl = cl,
                     lag_clim, 
                     function(x) st.lamb(env_U = x$recent,
                                         env_F = x$recent,
                                         clim_sd = sd(x$recent, na.rm = T),
                                         clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                         sig.strength = 1)
  ) 
  stopCluster(cl)

  data <- rbind(plag_u %>% bind_rows %>% mutate(type = "Umatrix"),
                plag_n %>% bind_rows %>% mutate(type = "None")) %>% 
    mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))
  
  
  origin <- readRDS(list.files(path = file.path("results", "06_COMPADRE_studies", "rds"), 
                     pattern = paste("mpm", df$SpeciesAuthor[sp], df$MatrixPopulation[sp], sep = "_"), 
                     full.names = T)) %>%
    lapply(., function(x) lapply(x, function(y) y %>% 
                                                mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")))
                                              ) %>% bind_rows) %>%
    bind_rows(., .id = "type")
  
  plot_df <- rbind(collapsed = data, 
        original = origin)
  
  lagpf_p <- ggplot(data) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    geom_point(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    labs(colour = "Lag type", title = "Lagged climate in U or F matrix") + 
    xlab("Climate standard deviation") + ylab("log lambda") +
    facet_grid(cols = vars(auto_cat), labeller = as_labeller(c("-0.9" = "Autocorrelation = -0.9",
                                                               "0" = "Autocorrelation = 0",
                                                               "0.9" = "Autocorrelation = 0.9"))) + 
    scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
    ggtitle(names(cell_values)[sp]) + theme(legend.position = "bottom")
  
  ggsave(file.path("results", paste0(names(cell_values)[sp], "_collapsed_mpm_simulations.tiff")), lagpf_p,
         width = 7, height = 4)
  
  }



