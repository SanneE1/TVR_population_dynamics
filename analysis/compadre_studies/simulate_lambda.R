suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(popbio))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(faux))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(Rage))

start <- Sys.time()
start

taskID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste("TASKID:", taskID))

args <- commandArgs(TRUE)

if(length(args) == 0) {
  stop("Supply csv file with Species_authors", call. = FALSE)
}

# load most current compadre database
load('/data/lagged/COMPADRE_v.X.X.X.4.RData')

## Selection ids of Species_author
species <- read.csv(args[1])

i = species$SpeciesAuthor[taskID]
j = species$MatrixPopulation[taskID]

if(!is.na(j)){
id <- which(compadre$metadata$SpeciesAuthor == i & compadre$metadata$MatrixPopulation == j)
} else {
  id <- which(compadre$metadata$SpeciesAuthor == i)
  
}

# Required functions

logit <- function(x) log(x/(1-x))

inv_logit <- function(x) {
  return(
    1/(1 + exp(-(x)))
  )
}

#-----------------------------------------------------------
# Select mpm's 
#-----------------------------------------------------------

# select only individual matrices
id2 <- id[which(compadre$metadata$SpeciesAuthor[id] == i & compadre$metadata$MatrixComposite[id] == "Individual")]
# Retrieve all full matrices
Amats <- lapply(as.list(id2), function(x) as.vector(compadre$mat[x][[1]]$matA)) %>% bind_cols
## matrix cell position in vector -> c([1,1], [2,1], etc., [1,2], [2,2], etc.)
Umats <- lapply(as.list(id2), function(x) as.vector(compadre$mat[x][[1]]$matU)) %>% bind_cols 
Fmats <- lapply(as.list(id2), function(x) as.vector(compadre$mat[x][[1]]$matF)) %>% bind_cols

# get dimension of the current current study/species_author. All matrices should be of the same size!!
dim <- unique(compadre$metadata$MatrixDimension[id2])

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
## I assume that when signal_strength is 1 the size of the climate effect is 1 SD. A signal_strength of 0.5 = 0.5 * SD for climate effect and random noise of Norm(0, 0.5*SD)

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

st.lamb <- function(growth, reproduction, signal_strength){
  
  n_it = length(growth)
  
  env <- data.frame(growth = growth,
                    reproduction = reproduction)
  
  env <- env[complete.cases(env), ]
  
  env <- split(env, seq(nrow(env)))
  
  ### Get all mpm's
  mats <- lapply(env, function(x) mpm(U_clim = x$growth, F_clim = x$reproduction, signal_strength))
  
  lambda = stoch.growth.rate(mats, maxt = n_it, verbose = F)$sim
  
  return(data.frame(lambda = lambda,
                    actual_sd = sd(growth, na.rm = T),
                    actual_auto = round(acf(growth, plot = F, na.action = na.pass)$acf[2], digits = 3),
                    signal_strength = signal_strength))
}

clim_sd <- seq(from = 0.01, to = 2, length.out = 10)
clim_corr <- c(-0.9,0,0.9)
sig.strength <- c(1, 0.5, 0.1)

df <- expand.grid(clim_sd, clim_corr, sig.strength) 
df <- do.call("rbind", replicate(30, df, simplify = FALSE)) %>% `names<-`(., c("clim_sd", "clim_auto", "sig.strength"))

lag_clim <- lapply(as.list(c(1:nrow(df))), function(x) create_seq(clim_sd = df$clim_sd[x], clim_corr = df$clim_auto[x], sig.strength = df$sig.strength[x]) %>% filter(complete.cases(.)))


# species specific mpm
  mpm <- function(U_clim, F_clim, signal_strength) {
    
    Umat <- Ucell_values$mean + Ucell_values$sd * (U_clim * signal_strength) + 
      rnorm(n = length(Ucell_values$sd),
        mean = 0, sd = (Ucell_values$sd * (1-signal_strength))) 
    Fmat <- Fcell_values$mean + Fcell_values$sd * (F_clim * signal_strength) + 
      rnorm(n = length(Fcell_values$sd),
            mean = 0, sd = (Fcell_values$sd * (1-signal_strength))) 
    
    Umat <- case_when(Umat < 0 ~ 0,
                      Umat > 1 ~ 1,
                      between(Umat, 0, 1) ~ Umat)

    Fmat <- case_when(Fmat <= 0 ~ 0,
                      Fmat > 0 ~ Fmat)
    Amat <- Umat + Fmat
    mpm <- matrix(Amat, nrow = dim)
  
    return(mpm)  
  }

  
  #### Lagged effect in U or F matrix
  
  lag_u <- pblapply(lag_clim, 
                    function(x) st.lamb(growth = x$lagged,
                                        reproduction = x$recent,
                                        sig.strength = unique(x$sig.strength))
  ) %>% bind_rows
  
  lag_f <- pblapply(lag_clim, 
                    function(x) st.lamb(growth = x$recent,
                                        reproduction = x$lagged,
                                        sig.strength = unique(x$sig.strength))
  ) %>% bind_rows
  
  lag_n <- pblapply(lag_clim, 
                    function(x) st.lamb(growth = x$recent,
                                        reproduction = x$recent,
                                        sig.strength = unique(x$sig.strength))
  ) %>% bind_rows
  
  lag_uf <- list("Umatrix" = data.frame(lag_u,
                                        type = "Umatrix",
                                        clim_sd = clim_sd,
                                        clim_auto = clim_corr),
                 "Fmatrix" = data.frame(lag_f,
                                        type = "Fmatrix",
                                        clim_sd = clim_sd,
                                        clim_auto = clim_corr),
                 "None" = data.frame(lag_n,
                                     type = "None",
                                     clim_sd = clim_sd,
                                     clim_auto = clim_corr))
  
output_dir <- args[2]

n_pop = length(unique(species$MatrixPopulation[which(species$SpeciesAuthor == i)]))

if(length(n_pop == 1)) {
  output_file <- paste0("mpm_", i, "_laguf.RDS")
} else {
  output_file <- paste("mpm", i, j, "laguf.RDS", sep = "_")
}

saveRDS(lag_uf, file.path(output_dir, output_file))

Sys.time() - start
