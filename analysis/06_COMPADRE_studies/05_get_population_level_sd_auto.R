library(tidyverse)
library(ggplot2)
library(pbapply)
library(popbio)
library(scales)

output_dir <- "results/06_COMPADRE_studies/actual_sd_auto/"

#-----------------------------------------------------------
# Create functions
#-----------------------------------------------------------

logit <- function(x) log(x/(1-x))

inv_logit <- function(x) {
  return(
    1/(1 + exp(-(x)))
  )
}

### Create environmental sequence ----------------------------
create_seq <- function(clim_sd, clim_corr, sig.strength, lag = 1, n_it = 5000) { 
  for(n in c(1:(n_it+lag))){ 
    if(n == 1) {
      seq <- rnorm(1)
    } else {
      seq[n] <- clim_corr * seq[n-1] + rnorm(1)
    }
  }
  seq <- scale(seq) * clim_sd
  lagged <- c(rep(NA, lag), seq)
  recent <- c(seq, rep(NA, lag))
  df <- data.frame(recent = recent,
                   lagged = lagged,
                   sig.strength = sig.strength)
  return(df)
}

# "Stocastic" mpm -----------------------------------

st.lamb <- function(growth, reproduction, clim_auto, sig.strength){
  
  n_it = length(growth)
  
  env <- data.frame(growth = growth,
                    reproduction = reproduction)
  env <- split(env, seq(nrow(env)))
  
  ### Get all mpm's
  mats <- lapply(env, function(x) mpm(U_clim = x$growth, F_clim = x$reproduction, sig.strength = sig.strength))
  
  df <- data.frame(lambda = stoch.growth.rate(mats, maxt = 5000, verbose = F)$sim,
                   clim_auto = clim_auto,
                   sig.strength = sig.strength)
  
  return(df)
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

df$prcp_diffU_0.05 <- NA
df$prcp_diffU_0.25 <- NA
df$prcp_diffU_0.5 <- NA
df$prcp_diffU_1 <- NA

df$temp_diffU_0.05 <- NA
df$temp_diffU_0.25 <- NA
df$temp_diffU_0.5 <- NA
df$temp_diffU_1 <- NA


### Run lambda simulations using matrices, sd and auto of specific populations
### 30 repetitions for 3 different climate signal strengths

load("data/COMPADRE_v.X.X.X.4.RData")

for(sp in c(1:nrow(df))) {
  
  i = df$SpeciesAuthor[sp]
  j = df$MatrixPopulation[sp]
  
  if(!is.na(j)){
    id <- which(compadre$metadata$SpeciesAuthor == i & 
                  compadre$metadata$MatrixPopulation == j &
                  compadre$metadata$MatrixComposite == "Individual")
  } else {
    id <- which(compadre$metadata$SpeciesAuthor == i & 
                  compadre$metadata$MatrixComposite == "Individual")
    
  }
  
  #-----------------------------------------------------------
  # get mpm's 
  #-----------------------------------------------------------
  
  # Retrieve all full matrices
  Amats <- lapply(as.list(id), function(x) as.vector(compadre$mat[x][[1]]$matA)) %>% bind_cols
  ## matrix cell position in vector -> c([1,1], [2,1], etc., [1,2], [2,2], etc.)
  Umats <- lapply(as.list(id), function(x) as.vector(compadre$mat[x][[1]]$matU)) %>% bind_cols 
  Fmats <- lapply(as.list(id), function(x) as.vector(compadre$mat[x][[1]]$matF)) %>% bind_cols
  
  # get dimension of the current current study/species_author. All matrices should be of the same size!!
  dim <- unique(compadre$metadata$MatrixDimension[id])
  
  if(length(dim) != 1) {
    stop("different sized matrices in same study")
  }
  
  print(paste("dim =", dim))
  
  
  ## Get mean and standard deviation for each cell in the matrices
  Acell_values <- Amats %>%
    summarise(mean = apply(., 1, mean),
              sd = apply(., 1, sd))%>%
    mutate(across(everything(), ~ replace(., is.na(.), 0)))
  
  Ucell_values <- Umats %>% 
    summarise(mean = apply(., 1, mean),
              sd = apply(., 1, sd))%>%
    mutate(across(everything(), ~ replace(., is.na(.), 0)))
  
  Fcell_values <- Fmats %>% 
    summarise(mean = apply(., 1, mean),
              sd = apply(., 1, sd))%>%
    mutate(across(everything(), ~ replace(., is.na(.), 0)))
  
  
  #-----------------------------------------------------------
  # create species specific 
  #-----------------------------------------------------------
  mpm <- function(U_clim, F_clim, sig.strength) {
    
    Umat <- Ucell_values$mean + (Ucell_values$sd * U_clim) * (sqrt(Ucell_values$sd^2 * sig.strength)/Ucell_values$sd) + 
      rnorm(n = length(Ucell_values$sd),
            mean = 0, sd = (Ucell_values$sd * (sqrt(Ucell_values$sd^2 * (1-sig.strength))/Ucell_values$sd))) 
    Fmat <- Fcell_values$mean + (Fcell_values$sd * F_clim) * (sqrt(Fcell_values$sd^2 * sig.strength)/Fcell_values$sd) + 
      rnorm(n = length(Fcell_values$sd),
            mean = 0, sd = (Fcell_values$sd * (sqrt(Fcell_values$sd^2 * (1-sig.strength))/Fcell_values$sd))) 
    
    ### nan's get produced in the above calculations, as some cells don't have values and thus there's a 
    ### division with 0 in the correction factors. The first case_when statement corrects this back to 0
    
    Umat <- case_when(is.nan(Umat) ~ 0,
                      Umat < 0 ~ 0,
                      Umat > 1 ~ 1,
                      between(Umat, 0, 1) ~ Umat)
    
    Fmat <- case_when(is.nan(Fmat) ~ 0,
                      Fmat <= 0 ~ 0,
                      Fmat > 0 ~ Fmat)
    Amat <- Umat + Fmat
    mpm <- matrix(Amat, nrow = dim)
    
    return(mpm)  
  }
  
  ## Create climate sequences with population's sd and autocorrelation 
  lag_prcp <- lapply(as.list(rep(c(1,0.5,0.25, 0.05), 30)), function(x) 
    create_seq(clim_sd = 1, clim_corr = df$prec_auto[sp], sig.strength = x) %>% filter(complete.cases(.)))
  
  lag_temp <- lapply(as.list(rep(c(1,0.5,0.25, 0.05), 30)), function(x) 
    create_seq(clim_sd = 1, clim_corr = df$tmean_auto[sp], sig.strength = x) %>% filter(complete.cases(.)))
  
  
  #-----------------------------------------------------------
  # Run simulations for precipitation
  #-----------------------------------------------------------
  #### Lagged effect in U or F matrix
  
  plag_u <- pblapply(lag_prcp, 
                     function(x) st.lamb(growth = x$lagged,
                                         reproduction = x$recent,
                                         clim_auto = df$prec_auto[sp],
                                         sig.strength = unique(x$sig.strength))
  ) 
  
  # plag_f <- pblapply(lag_prcp, 
  #                   function(x) st.lamb(growth = x$recent,
  #                                       reproduction = x$lagged,
  #                                       clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
  #                                       sig.strength = unique(x$sig.strength))
  # ) 
  
  plag_n <- pblapply(lag_prcp, 
                     function(x) st.lamb(growth = x$recent,
                                         reproduction = x$recent,
                                         clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                         sig.strength = unique(x$sig.strength))
  ) 
  
  #-----------------------------------------------------------
  # Run simulations for Temperature
  #-----------------------------------------------------------
  #### Lagged effect in U or F matrix
  
  
  tlag_u <- pblapply(lag_temp, 
                     function(x) st.lamb(growth = x$lagged,
                                         reproduction = x$recent,
                                         clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                         sig.strength = unique(x$sig.strength))
  ) 
  
  # tlag_f <- pblapply(lag_temp, 
  #                    function(x) st.lamb(growth = x$recent,
  #                                        reproduction = x$lagged,
  #                                        clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
  #                                        sig.strength = unique(x$sig.strength))
  # ) 
  
  tlag_n <- pblapply(lag_temp, 
                     function(x) st.lamb(growth = x$recent,
                                         reproduction = x$recent,
                                         clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                         sig.strength = unique(x$sig.strength))
  ) 
  
  
  ### Get differences in lambda between Ulagged and none for precipitation
  df1 <- left_join(
    plag_u %>% bind_rows() %>% 
      group_by(sig.strength) %>%
      summarise(Ulambda = mean(lambda, na.rm = T)),
    plag_n %>% bind_rows() %>% 
      group_by(sig.strength) %>%
      summarise(Nlambda = mean(lambda, na.rm = T))
  ) %>%
    mutate(diffU = (Ulambda - Nlambda)/abs(Nlambda))
  
  ### Get differences in lambda between Ulagged and none for temperature
  df2 <- left_join(
    tlag_u %>% bind_rows() %>% 
      group_by(sig.strength) %>%
      summarise(Ulambda = mean(lambda, na.rm = T)),
    tlag_n %>% bind_rows() %>% 
      group_by(sig.strength) %>%
      summarise(Nlambda = mean(lambda, na.rm = T))
  ) %>%
    mutate(diffU = (Ulambda - Nlambda)/abs(Nlambda))
  
  ### Combine datasets
  df3 <- cbind(pivot_wider(df1, diffU, names_from = sig.strength, names_prefix = "prcp_diffU_", values_from = diffU),
               pivot_wider(df2, diffU, names_from = sig.strength, names_prefix = "temp_diffU_", values_from = diffU))
  
  ### Assign differences to main data.frame
  df[which(df$SpeciesAuthor == i & df$MatrixPopulation == j), c(15:22)] <- df3
  
  
  ### Set up dataframe for lambda plots
  plot_df <- rbind(
    plag_u %>% bind_rows %>% select(lambda, sig.strength) %>% mutate(type = "Umatrix"),
    plag_n %>% bind_rows %>% select(lambda, sig.strength) %>% mutate(type = "None")
  ) %>% mutate(climate = "Precipitation")
  
  plot_df <- rbind(plot_df,
                   rbind(tlag_u %>% bind_rows %>% select(lambda, sig.strength) %>% mutate(type = "Umatrix"),
                         tlag_n %>% bind_rows %>% select(lambda, sig.strength) %>% mutate(type = "None")
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
  
  
  
  
  ## Run mpm sequences again
  # Randomly select climate sequences to check the cell values with sig.strength=1
  n = sample(c(1:30),1)
  
  growth <- lag_prcp[[n]]$recent
  reproduction <- lag_prcp[[n]]$lagged
  
  n_it = length(growth)
  env <- data.frame(growth = growth,
                    reproduction = reproduction)
  env <- env[complete.cases(env), ]
  env <- split(env, seq(nrow(env)))
  ### Get all Umat and Fmat values
  matsU <- function(climate, sig.strength) {
    Umat <- Ucell_values$mean + (Ucell_values$sd * climate$growth) * (sqrt(Ucell_values$sd^2 * sig.strength)/Ucell_values$sd) + 
      rnorm(n = length(Ucell_values$sd),
            mean = 0, sd = (Ucell_values$sd * (sqrt(Ucell_values$sd^2 * (1-sig.strength))/Ucell_values$sd))) 
    Umat <- case_when(is.nan(Umat) ~ 0,
                      Umat < 0 ~ 0,
                      Umat > 1 ~ 1,
                      between(Umat, 0, 1) ~ Umat)
    return(Umat)
  }
  
  matsF <- function(climate, sig.strength) {
    Fmat <- Fcell_values$mean + (Fcell_values$sd * climate$reproduction) * (sqrt(Fcell_values$sd^2 * sig.strength)/Fcell_values$sd) + 
      rnorm(n = length(Fcell_values$sd),
            mean = 0, sd = (Fcell_values$sd * (sqrt(Fcell_values$sd^2 * (1-sig.strength))/Fcell_values$sd))) 
    Fmat <- case_when(is.nan(Fmat) ~ 0,
                      Fmat <= 0 ~ 0,
                      Fmat > 0 ~ Fmat)
    return(Fmat)
  }
  
  
  ## Set up dataframes to compare simulated cell values to observed
  U1 <- lapply(env, function(x) matsU(x, sig.strength=1))
  U1 <- data.frame(t(lapply(U1, as.vector) %>% bind_rows)) %>% 
    pivot_longer(cols = everything(), names_to = "cell", values_to = "value") %>%
    mutate(cell = as.integer(gsub("X", "", cell)),
           matrix = "U",
           sig.strength = 1)
  U0.5 <- lapply(env, function(x) matsU(x, sig.strength=0.5))
  U0.5 <- data.frame(t(lapply(U0.5, as.vector) %>% bind_rows)) %>% 
    pivot_longer(cols = everything(), names_to = "cell", values_to = "value") %>%
    mutate(cell = as.integer(gsub("X", "", cell)),
           matrix = "U",
           sig.strength = 0.5) 
  
  F1 <- lapply(env, function(x) matsF(x, sig.strength = 1))
  F1 <- data.frame(t(lapply(F1, as.vector) %>% bind_rows)) %>% 
    pivot_longer(cols = everything(), names_to = "cell", values_to = "value") %>%
    mutate(cell = as.integer(gsub("X", "", cell)),
           matrix = "F",
           sig.strength = 1)                
  F0.5 <- lapply(env, function(x) matsF(x, sig.strength = 0.5))
  F0.5 <- data.frame(t(lapply(F0.5, as.vector) %>% bind_rows)) %>% 
    pivot_longer(cols = everything(), names_to = "cell", values_to = "value") %>%
    mutate(cell = as.integer(gsub("X", "", cell)),
           matrix = "F",
           sig.strength = 0.5)              
  
  actual_sd <- rbind(Ucell_values %>% select(sd) %>% mutate(matrix = "U") %>% rename(actual_sd = sd) %>%
                       tibble::rownames_to_column(var = "cell") %>% mutate(cell = as.integer(cell)),
                     Fcell_values %>% select(sd) %>% mutate(matrix = "F") %>% rename(actual_sd = sd) %>%
                       tibble::rownames_to_column(var = "cell") %>% mutate(cell = as.integer(cell)))
  
  U <- Umats %>% t %>% as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "cell", values_to = "value") %>%
    mutate(cell = as.integer(gsub("V", "", cell)))
  
  Fobs <- Fmats  %>% t %>% as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "cell", values_to = "value") %>%
    mutate(cell = as.integer(gsub("V", "", cell)))
  
  # plot sd of simulations and observed
  sds <- rbind(U1,U0.5,F1,F0.5) %>%
    group_by(matrix, sig.strength, cell) %>%
    summarise(sd = sd(value, na.rm = T)) %>%
    left_join(., actual_sd) %>% ungroup() %>%
    pivot_longer(cols = c("sd", "actual_sd"), names_to = "sd_type", values_to = "sd") %>% as.data.frame %>%
    ggplot(.) + 
    geom_bar(aes(x = matrix, y = sd, fill = factor(sd_type)), stat = "identity", position = "dodge") +
    facet_wrap(vars(cell), scales = "free", dir = "v") +
    theme(legend.position = "bottom") + scale_fill_discrete(name = "", labels = c("Observed sd", "Simulated sd"))
  
  
  ## sd distributions with different sig.strengths
  valuesUsig <- ggplot() + 
    geom_histogram(data = U1, aes(x = value, y = stat(width*density), fill = factor(sig.strength))) + 
    geom_histogram(data = U0.5, aes(x = value, y = -stat(width*density), fill = factor(sig.strength))) +
    facet_wrap(vars(cell), scales = "free", dir = "v") +
    scale_y_continuous(labels = percent_format()) + labs(title = paste(i, j), fill = "Climate signal strength") +
    ylab("in % of matrices") + xlab("cell value") + theme(legend.position = "bottom")
  
  valuesFsig <- ggplot() + 
    geom_histogram(data = F1, aes(x = value, y = stat(width*density), fill = factor(sig.strength))) + 
    geom_histogram(data = F0.5, aes(x = value, y = -stat(width*density), fill = factor(sig.strength))) +
    facet_wrap(vars(cell), scales = "free", dir = "v") +
    scale_y_continuous(labels = percent_format()) + labs(title = paste(i, j), fill = "Climate signal strength") +
    ylab("in % of matrices") + xlab("cell value") + theme(legend.position = "bottom")
  
  ## Compare simulated distribution with observed
  valuesUobs <- ggplot() + 
    geom_histogram(data = U, aes(x = value, y = stat(width*density), fill = "Observed")) + 
    geom_histogram(data = U0.5, aes(x = value, y = -stat(width*density), fill = "Simulated")) +
    facet_wrap(vars(cell), scales = "free", dir = "v") +
    scale_y_continuous(labels = percent_format()) + labs(title = "Observed vs. simulated matrix cell values" , 
                                                         subtitle = paste(i, j), fill = "U matrix cell values") +
    ylab("") + xlab("cell value") + theme(legend.position = "bottom")
  
  valuesFobs <- ggplot() + 
    geom_histogram(data = Fobs, aes(x = value, y = stat(width*density), fill = "Observed")) + 
    geom_histogram(data = F0.5, aes(x = value, y = -stat(width*density), fill = "Simulated")) +
    facet_wrap(vars(cell), scales = "free", dir = "v") +
    scale_y_continuous(labels = percent_format()) + labs(title = "Observed vs. simulated F matrix cell values" , 
                                                         subtitle = paste(i, j), fill = "F matrix cell values") +
    ylab("") + xlab("cell value") + theme(legend.position = "bottom")
  
  
  if(is.na(j)) {
    filename = file.path(output_dir, paste0("cell_values_", i, "_plots.pdf"))
  } else{
    filename = file.path(output_dir, paste0("cell_values_", i, "_", j, "_plots.pdf"))
  }
  
  
  pdf(filename)
  
  valuesUobs
  valuesFobs
  sds
  valuesUsig
  valuesFsig
  
  dev.off()
  
}








