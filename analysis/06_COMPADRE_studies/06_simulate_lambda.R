suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(popbio))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(faux))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(Rage))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyr))

start <- Sys.time()
start

taskID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste("TASKID:", taskID))

args <- commandArgs(TRUE)

if(length(args) == 0) {
  stop("Supply csv file with Species_authors", call. = FALSE)
}

# load most current compadre database
load('/data/lagged/COMPADRE_v.6.21.1.0.RData')

## Selection ids of Species_author
species <- read.csv(args[1])

i = species$SpeciesAuthor[taskID]
j = species$MatrixPopulation[taskID]

if(!is.na(j)){
id <- which(compadre$metadata$SpeciesAuthor == i & 
              compadre$metadata$MatrixPopulation == j &
              compadre$metadata$MatrixComposite == "Individual")
} else {
  id <- which(compadre$metadata$SpeciesAuthor == i & 
                compadre$metadata$MatrixComposite == "Individual")
  
}

#-----------------------------------------------------------
# Select mpm's 
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
# Simulate with the mpm's
#-----------------------------------------------------------

## create stochastic climate sensitive mpm's from the data
## I assume that when sig.strength is 1 the size of the climate effect is 1 SD. A sig.strength of 0.5 = 0.5 * SD for climate effect and random noise of Norm(0, 0.5*SD)

### Create environmental sequence ----------------------------
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

# "Stocastic" mpm -----------------------------------

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

clim_sd <- seq(from = 0.01, to = 2, length.out = 10)
clim_corr <- c(-0.9,0,0.9)

df <- expand.grid(clim_sd, clim_corr) 
df <- do.call("rbind", replicate(30, df, simplify = FALSE)) %>% `names<-`(., c("clim_sd", "clim_auto"))

lag_clim <- lapply(as.list(c(1:nrow(df))), function(x) 
  create_seq(n_it = 5000,clim_sd = df$clim_sd[x], clim_auto = df$clim_auto[x], lag = 1) %>% filter(complete.cases(.)))


# species specific mpm
mpm <- function(U_clim, F_clim, sig.strength = 1, clim_sd) {
  
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

  
  #### Lagged effect in U or F matrix
  
 # lag_u <- pblapply(lag_clim, 
 #                   function(x) st.lamb(env_U = x$lagged,
 #                                       env_F = x$recent,
 #                                       clim_sd = sd(x$recent, na.rm = T),
 #                                       clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
 #                                       sig.strength = 1)
 # ) 
  
 # lag_f <- pblapply(lag_clim, 
 #                   function(x) st.lamb(env_U = x$recent,
 #                                       env_F = x$lagged,
 #                                       clim_sd = sd(x$recent, na.rm = T),
 #                                       clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
 #                                       sig.strength = 1)
 # ) 
  
 # lag_n <- pblapply(lag_clim, 
 #                   function(x) st.lamb(env_U = x$recent,
 #                                       env_F = x$recent,
 #                                       clim_sd = sd(x$recent, na.rm = T),
 #                                       clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
 #                                       sig.strength = 1)
 # ) 
  
 # lag_uf <- list("Umatrix" = lag_u,
 #                "Fmatrix" = lag_f,
 #                "None" = lag_n)
  
#output_dir <- args[2]

#n_pop = length(unique(species$MatrixPopulation[which(species$SpeciesAuthor == i)]))

#if(is.na(j)) {
#  output_file <- paste0("mpm_", i, "_laguf.RDS")
#} else {
#  output_file <- paste("mpm", i, j, "laguf.RDS", sep = "_")
#}

#saveRDS(lag_uf, file.path(output_dir, output_file))


## Run mpm sequences again
# Randomly select a sequence with a low and high sd
low_sd <- sample(which(round(df$clim_sd, digits = 3) == 0.231), 1)
high_sd <- sample(which(round(df$clim_sd, digits = 3) == 2), 1)

cells <- list()

for(n in c(low_sd, high_sd)) {
  
  env_U <- lag_clim[[n]]$recent
  env_F <- lag_clim[[n]]$lagged
  
n_it = length(env_U)
  
  env <- list(U_clim = env_U,
              F_clim = env_F,
              clim_sd = c(0.231, 2)[n],
              sig.strength = 1)
  
  ### Get all mpm's
  mats <- pmap(env, mpm) 
  
  
  mats = sapply(mats, FUN = as.vector, USE.NAMES = T) %>% t %>% as.data.frame()


cells <- append(cells,
                mats)

}

if(is.na(j)) {
  saveRDS(cells, file.path(output_dir, paste0(i, "_cell_values_HnL.RDS")))
} else {
saveRDS(cells, file.path(output_dir, paste(i, j, "cell_values_HnL.RDS", sep= "_")))
}


Sys.time() - start
