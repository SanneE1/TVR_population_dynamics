rm(list = ls())
library(tidyverse)
library(ggplot2)
library(pbapply)
library(popbio)
library(scales)
library(boot)
library(parallel)
library(Rage)
library(rstatix)

n_it = 5000

output_dir <- "results/06_COMPADRE_studies/actual_sd_auto/"

#-----------------------------------------------------------
# Create functions
#-----------------------------------------------------------

### Create environmental sequence ----------------------------
### creates a sequence of climate anomalies whith a specified standard deviation and autocorrelation.
### the function then creates another sequence of the same length and standard deviation, to be used as random noise
### of each sequence two vectors are created, one which is the "original" and a second, which is offset by a specified period
### to create a lagged sequence
create_seq <- function(n_it, clim_sd, clim_auto, sig.strength, lag = 1) { 
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
                   lagged = lagged,
                   sig.strength = sig.strength)
  return(df)
}

# "Stocastic" mpm -----------------------------------

st.lamb <- function(env_U, env_F, 
                    clim_sd, clim_auto, sig.strength) {
  
  n_it = length(env_U)
  
  env <- list(U_clim = env_U,
              F_clim = env_F,
              clim_sd = rep(clim_sd, n_it),
              sig.strength = rep(sig.strength, n_it))
  
  ### Get all mpm's
  mats <- pmap(env, mpm) #%>% Filter(Negate(anyNA), .)
  
  
  df <- data.frame(lambda = stoch.growth.rate(mats, maxt = n_it, verbose = F)$sim,
                   clim_sd = clim_sd,
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

df$Pmean_c_0.05 <- NA
df$Pmean_U_0.05 <- NA
df$Sig_p_0.05 <- NA

df$Pmean_c_0.25 <- NA
df$Pmean_U_0.25 <- NA
df$Sig_p_0.25 <- NA

df$Pmean_c_0.50 <- NA
df$Pmean_U_0.50 <- NA
df$Sig_p_0.50 <- NA

df$Pmean_c_1 <- NA
df$Pmean_U_1 <- NA
df$Sig_p_1 <- NA

df$Tmean_c_0.05 <- NA
df$Tmean_U_0.05 <- NA
df$Sig_t_0.05 <- NA

df$Tmean_c_0.25 <- NA
df$Tmean_U_0.25 <- NA
df$Sig_t_0.25 <- NA

df$Tmean_c_0.50 <- NA
df$Tmean_U_0.50 <- NA
df$Sig_t_0.50 <- NA

df$Tmean_c_1 <- NA
df$Tmean_U_1 <- NA
df$Sig_t_1 <- NA

df$cov_sg <- NA
df$cov_sf <- NA
df$cov_gf <- NA


### Run lambda simulations using matrices, sd and auto of specific populations
### 30 repetitions for 3 different climate signal strengths

load("data/COMPADRE_v.6.21.1.0.RData")

for(sp in c(1:nrow(df))) {
  
  if(is.na(df$prec_auto[sp])) next   ## if prec_auto is na, prec_sd, tmean_auto & tmean_sd are also na
  
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
    mpm <- matrix(Amat, nrow = dim)

    return(mpm)
  }

  ## Create climate sequences with population's sd and autocorrelation
  lag_prcp <- lapply(as.list(rep(c(1,0.5,0.25, 0.05), 30)), function(x)
    create_seq(n_it = n_it, clim_sd = 1, clim_auto = df$prec_auto[sp], sig.strength = x) %>%
      filter(complete.cases(.)))

  lag_temp <- lapply(as.list(rep(c(1,0.5,0.25, 0.05), 30)), function(x)
    create_seq(n_it = n_it, clim_sd = 1, clim_auto = df$tmean_auto[sp], sig.strength = x) %>%
      filter(complete.cases(.)))


  # #-----------------------------------------------------------
  # # Run simulations for precipitation
  # #-----------------------------------------------------------
  # # Set up parallel runs
  # no_cores <- detectCores()
  # cl <- makeCluster(no_cores - 2)
  # ## export libraries to workers
  # clusterEvalQ(cl, c(library(popbio), library(tidyverse), library(purrr), library(boot)))
  # clusterExport(cl, c("st.lamb", "mpm", "lag_prcp", "lag_temp", "dim",
  #                     "Ucell_values", "Fcell_values", "df", "sp"))
  # 
  # 
  # #### Lagged effect in U or F matrix
  # 
  # plag_u <- pblapply(cl = cl,
  #                    lag_prcp, 
  #                    function(x) st.lamb(env_U = x$lagged,
  #                                        env_F = x$recent,
  #                                        clim_sd = 1,
  #                                        clim_auto = df$prec_auto[sp],
  #                                        sig.strength = unique(x$sig.strength))
  # ) 
  # 
  # plag_n <- pblapply(cl = cl,
  #                    lag_prcp, 
  #                    function(x) st.lamb(env_U = x$recent,
  #                                        env_F = x$recent,
  #                                        clim_sd = 1,
  #                                        clim_auto = df$prec_auto[sp],
  #                                        sig.strength = unique(x$sig.strength))
  # ) 
  # 
  # #-----------------------------------------------------------
  # # Run simulations for Temperature
  # #-----------------------------------------------------------
  # #### Lagged effect in U or F matrix
  # 
  # 
  # tlag_u <- pblapply(cl = cl, 
  #                    lag_temp, 
  #                    function(x) st.lamb(env_U = x$lagged,
  #                                        env_F = x$recent,
  #                                        clim_sd = 1,
  #                                        clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
  #                                        sig.strength = unique(x$sig.strength))
  # ) 
  # 
  # 
  # tlag_n <- pblapply(cl = cl,
  #                    lag_temp, 
  #                    function(x) st.lamb(env_U = x$recent,
  #                                        env_F = x$recent,
  #                                        clim_sd = 1,
  #                                        clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
  #                                        sig.strength = unique(x$sig.strength))
  # ) 
  # 
  # stopCluster(cl)
  # 
  # save(plag_u, plag_n, tlag_u, tlag_n, file = file.path(output_dir, "rds", paste0("simulations_", i,"_", j, ".RData"))) 
  # 
  
  load(file.path(output_dir, "rds", paste0("simulations_", i,"_", j, ".RData")))
  
  
  ### Get differences in lambda between Ulagged and none for precipitation
  df1 <- rbind(
    plag_u %>% bind_rows() %>% mutate(type = "Umatrix"), 
    plag_n %>% bind_rows() %>% mutate(type = "control"))
  
  df1 %>% group_by(sig.strength, type) %>% shapiro_test(lambda)
  test1 <- df1 %>% group_by(sig.strength) %>% wilcox_test(lambda ~ type) %>% add_significance() 
  test1b <- df1 %>% group_by(sig.strength, type) %>% get_summary_stats(lambda)
  df[which(df$SpeciesAuthor == i & df$MatrixPopulation == j), c("Sig_p_0.05", "Sig_p_0.25", "Sig_p_0.50", "Sig_p_1")] <- round(test1$p, digits = 3)
  df[which(df$SpeciesAuthor == i & df$MatrixPopulation == j), 
     c("Pmean_c_0.05", "Pmean_U_0.05", "Pmean_c_0.25", "Pmean_U_0.25", "Pmean_c_0.50", "Pmean_U_0.50", "Pmean_c_1", "Pmean_U_1")] <- round(test1b$mean, digits = 3)
  
  
  ### Get differences in lambda between Ulagged and none for temperature
  df2 <- rbind(
    tlag_u %>% bind_rows() %>% mutate(type = "Umatrix"), 
    tlag_n %>% bind_rows() %>% mutate(type = "control"))
  
  df2 %>% group_by(sig.strength, type) %>% shapiro_test(lambda)
  test2 <- df2 %>% group_by(sig.strength) %>% wilcox_test(lambda ~ type) %>% add_significance()
  test2b <- df2 %>% group_by(sig.strength, type) %>% get_summary_stats(lambda)
  df[which(df$SpeciesAuthor == i & df$MatrixPopulation == j), c("Sig_t_0.05", "Sig_t_0.25", "Sig_t_0.50", "Sig_t_1")] <- round(test2$p, digits = 3)
  df[which(df$SpeciesAuthor == i & df$MatrixPopulation == j), 
     c("Tmean_c_0.05", "Tmean_U_0.05", "Tmean_c_0.25", "Tmean_U_0.25", "Tmean_c_0.50", "Tmean_U_0.50", "Tmean_c_1", "Tmean_U_1")] <- round(test2b$mean, digits = 3)
  
  
  ### Get survival, growth and fecundity to calculate covariance of the vital rates
  u_mat <- lapply(as.list(id), function(x) compadre$mat[x][[1]]$matU %>% `row.names<-`(colnames(compadre$mat[x][[1]]$matU)))
  f_mat <- lapply(as.list(id), function(x) compadre$mat[id][[1]]$matF %>%
                    `row.names<-`(colnames(compadre$mat[id][[1]]$matU)) %>%
                    `colnames<-`(colnames(compadre$mat[id][[1]]$matU)) )

  gr <- map(u_mat, vr_growth) %>% unlist()
  su <- map(u_mat, vr_survival) %>% unlist()
  fe <- map2(u_mat, f_mat, vr_fecundity) %>% unlist()

  df$cov_sg[which(df$SpeciesAuthor == i & df$MatrixPopulation == j)] <- cov(gr, su, use = "pairwise.complete.obs")
  df$cov_sf[which(df$SpeciesAuthor == i & df$MatrixPopulation == j)] <- cov(gr, fe, use = "pairwise.complete.obs")
  df$cov_gf[which(df$SpeciesAuthor == i & df$MatrixPopulation == j)] <- cov(su, fe, use = "pairwise.complete.obs")


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
    label = c(paste("SD =", round(df$prec_sd[sp],3), "autocor =", round(df$prec_auto[sp],3)),
              paste("SD =", round(df$tmean_sd[sp],3), "autocor =", round(df$tmean_auto[sp],3))),
    climate = c("Precipitation", "Temperature"),
    x = c(2,2)
  )

  plot <- ggplot(plot_df) + geom_boxplot(aes(x = as.character(sig.strength), y = lambda, fill = type)) +
    facet_grid(cols = vars(climate)) +
    geom_text(data = facet_text,
              aes(x = x, y = Inf, label = label), vjust = 1.2, hjust = 0.35, size = 4.5) +
    xlab("climate signal strength") + ylab("log lambda") + 
    scale_fill_manual(name = "Simulation", label = c("None" = "control", "Umatrix" = "MCD"), 
                      values = c("Umatrix" = "#E69F00", "None" = "#0072B2")) + 
    theme_minimal() + theme(legend.position = "bottom",
                            strip.background = element_rect(fill = "grey80"),
                            text = element_text(size = 16))
    

  ggsave(file.path(output_dir, paste0("lambda_comparison_", i, "_", j, ".tiff")), plot = plot,
         width = 5.5, height = 6.5)




  # ## Run mpm sequences again
  # # Randomly select climate sequences to check the cell values with sig.strength=1
  # n1 = sample(seq(1,length(lag_prcp), by = 4),1)
  # n0.5 = sample(seq(2,length(lag_prcp), by = 4),1)
  # 
  # ### Get all Umat and Fmat values
  # matsU <- function(U_clim, clim_sd, sig.strength) {
  #   devU <- U_clim * (sqrt(clim_sd^2 * sig.strength)/clim_sd) +
  #     rnorm(length(Ucell_values$mean), mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd)
  #   pU <- pnorm(devU, mean = 0, sd = clim_sd)
  #   Umat <- qbeta(pU,
  #                 (((Ucell_values$mean*(1-Ucell_values$mean))/(Ucell_values$sd * clim_sd)^2) - 1) * Ucell_values$mean,
  #                 (((Ucell_values$mean*(1-Ucell_values$mean))/(Ucell_values$sd * clim_sd)^2) - 1) * (1 - Ucell_values$mean))
  #   Umat <- as.data.frame(t(Umat))
  #   colnames(Umat) <- c(1:dim^2)
  # 
  #   return(Umat)
  # }
  # 
  # matsF <- function(F_clim, clim_sd, sig.strength) {
  #   devF <- F_clim * (sqrt(clim_sd^2 * sig.strength)/clim_sd) +
  #     rnorm(length(Fcell_values$mean), mean = 0, sd = clim_sd) * (sqrt(clim_sd^2 * (1-sig.strength))/clim_sd)
  #   pF <- pnorm(devF, mean = 0, sd = clim_sd)
  #   Fmat <- qgamma(pF,
  #                  (Fcell_values$mean^2)/(Fcell_values$sd * clim_sd)^2,
  #                  (Fcell_values$mean)/(Fcell_values$sd * clim_sd)^2) %>% replace_na(., 0)
  #   Fmat <- as.data.frame(t(Fmat))
  #   colnames(Fmat) <- c(1:dim^2)
  # 
  #   return(Fmat)
  # }
  # 
  # 
  # ## Set up dataframes to compare simulated cell values to observed
  # U1 <- pmap(list(U_clim = lag_prcp[[n1]]$recent,
  #                 clim_sd = 1,
  #                 sig.strength = lag_prcp[[n1]]$sig.strength),
  #            matsU) %>% bind_rows %>%
  #   pivot_longer(cols = everything(), names_to = "cell", values_to = "value") %>%
  #   mutate(cell = as.integer(gsub("X", "", cell)),
  #          matrix = "U",
  #          sig.strength = 1)
  # 
  # 
  # U0.5 <- pmap(list(U_clim = lag_prcp[[n0.5]]$recent,
  #                   clim_sd = 1,
  #                   sig.strength = lag_prcp[[n0.5]]$sig.strength),
  #              matsU) %>% bind_rows %>%
  #   pivot_longer(cols = everything(), names_to = "cell", values_to = "value") %>%
  #   mutate(cell = as.integer(gsub("X", "", cell)),
  #          matrix = "U",
  #          sig.strength = 0.5)
  # 
  # F1 <- pmap(list(F_clim = lag_prcp[[n1]]$recent,
  #                 clim_sd = 1,
  #                 sig.strength = lag_prcp[[n1]]$sig.strength),
  #            matsF) %>% bind_rows %>%
  #   pivot_longer(cols = everything(), names_to = "cell", values_to = "value") %>%
  #   mutate(cell = as.integer(gsub("X", "", cell)),
  #          matrix = "F",
  #          sig.strength = 1)
  # 
  # F0.5 <- pmap(list(F_clim = lag_prcp[[n0.5]]$recent,
  #                   clim_sd = 1,
  #                   sig.strength = lag_prcp[[n0.5]]$sig.strength),
  #              matsF) %>% bind_rows %>%
  #   pivot_longer(cols = everything(), names_to = "cell", values_to = "value") %>%
  #   mutate(cell = as.integer(gsub("X", "", cell)),
  #          matrix = "F",
  #          sig.strength = 0.5)
  # 
  # actual_sd <- rbind(Ucell_values %>% select(sd) %>% mutate(matrix = "U") %>% rename(actual_sd = sd) %>%
  #                      tibble::rownames_to_column(var = "cell") %>% mutate(cell = as.integer(cell)),
  #                    Fcell_values %>% select(sd) %>% mutate(matrix = "F") %>% rename(actual_sd = sd) %>%
  #                      tibble::rownames_to_column(var = "cell") %>% mutate(cell = as.integer(cell)))
  # 
  # Uobs <- Umats %>% t %>% as.data.frame() %>%
  #   pivot_longer(cols = everything(), names_to = "cell", values_to = "value") %>%
  #   mutate(cell = as.integer(gsub("V", "", cell)))
  # 
  # Fobs <- Fmats  %>% t %>% as.data.frame() %>%
  #   pivot_longer(cols = everything(), names_to = "cell", values_to = "value") %>%
  #   mutate(cell = as.integer(gsub("V", "", cell)))
  # 
  # # plot sd of simulations and observed
  # sds <- rbind(U1,U0.5,F1,F0.5) %>%
  #   group_by(matrix, sig.strength, cell) %>%
  #   summarise(sd = sd(value, na.rm = T)) %>%
  #   left_join(., actual_sd) %>% ungroup() %>%
  #   pivot_longer(cols = c("sd", "actual_sd"), names_to = "sd_type", values_to = "sd") %>% as.data.frame %>%
  #   ggplot(.) +
  #   geom_bar(aes(x = matrix, y = sd, fill = factor(sd_type)), stat = "identity", position = "dodge") +
  #   facet_wrap(vars(cell), scales = "free", dir = "v") +
  #   theme(legend.position = "bottom") + scale_fill_discrete(name = "", labels = c("Observed sd", "Simulated sd"))
  # 
  # 
  # ## sd distributions with different sig.strengths
  # valuesUsig <- ggplot() +
  #   geom_histogram(data = U1, aes(x = value, y = stat(width*density), fill = factor(sig.strength))) +
  #   geom_histogram(data = U0.5, aes(x = value, y = -stat(width*density), fill = factor(sig.strength))) +
  #   facet_wrap(vars(cell), scales = "free", dir = "v") +
  #   scale_y_continuous(labels = percent_format()) + labs(title = paste(i, j), fill = "Climate signal strength") +
  #   ylab("in % of matrices") + xlab("cell value") + theme(legend.position = "bottom")
  # 
  # valuesFsig <- ggplot() +
  #   geom_histogram(data = F1, aes(x = value, y = stat(width*density), fill = factor(sig.strength))) +
  #   geom_histogram(data = F0.5, aes(x = value, y = -stat(width*density), fill = factor(sig.strength))) +
  #   facet_wrap(vars(cell), scales = "free", dir = "v") +
  #   scale_y_continuous(labels = percent_format()) + labs(title = paste(i, j), fill = "Climate signal strength") +
  #   ylab("in % of matrices") + xlab("cell value") + theme(legend.position = "bottom")
  # 
  # ## Compare simulated distribution with observed
  # valuesUobs <- ggplot() +
  #   geom_histogram(data = Uobs, aes(x = value, y = stat(width*density), fill = "Observed")) +
  #   geom_histogram(data = U0.5, aes(x = value, y = -stat(width*density), fill = "Simulated")) +
  #   facet_wrap(vars(cell), scales = "free", dir = "v") +
  #   scale_y_continuous(labels = percent_format()) + labs(title = "Observed vs. simulated matrix cell values" ,
  #                                                        subtitle = paste(i, j), fill = "U matrix cell values") +
  #   ylab("") + xlab("cell value") + theme(legend.position = "bottom")
  # 
  # valuesFobs <- ggplot() +
  #   geom_histogram(data = Fobs, aes(x = value, y = stat(width*density), fill = "Observed")) +
  #   geom_histogram(data = F0.5, aes(x = value, y = -stat(width*density), fill = "Simulated")) +
  #   facet_wrap(vars(cell), scales = "free", dir = "v") +
  #   scale_y_continuous(labels = percent_format()) + labs(title = "Observed vs. simulated F matrix cell values" ,
  #                                                        subtitle = paste(i, j), fill = "F matrix cell values") +
  #   ylab("") + xlab("cell value") + theme(legend.position = "bottom")
  # 
  # 
  # if(is.na(j)) {
  #   filename = file.path(output_dir, paste0("cell_values_", i, "_plots.pdf"))
  # } else{
  #   filename = file.path(output_dir, paste0("cell_values_", i, "_", j, "_plots.pdf"))
  # }
  # 
  # 
  # pdf(filename)
  # 
  # print(valuesUobs)
  # print(valuesFobs)
  # print(sds)
  # print(valuesUsig)
  # print(valuesFsig)
  # 
  # dev.off()
  
}



write.csv(df, file.path(output_dir, "species_information.csv"))





