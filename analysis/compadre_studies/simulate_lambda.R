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

id <- which(compadre$metadata$SpeciesAuthor == i & compadre$metadata$MatrixPopulation == j)


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


id2 <- id[which(compadre$metadata$SpeciesAuthor[id] == i & compadre$metadata$MatrixComposite[id] == "Individual")]
# Retrieve all full matrices
Amats <- lapply(as.list(id2), function(x) as.vector(compadre$mat[x][[1]]$matA)) %>% bind_cols
## matrix cell position in vector -> c([1,1], [2,1], etc., [1,2], [2,2], etc.)
Umats <- lapply(as.list(id2), function(x) as.vector(compadre$mat[x][[1]]$matU)) %>% bind_cols 
Fmats <- lapply(as.list(id2), function(x) as.vector(compadre$mat[x][[1]]$matF)) %>% bind_cols


dim <- unique(compadre$metadata$MatrixDimension[id2])

if(length(dim) != 1) {
  stop("different sized matrices in same study")
}

print(paste("dim =", dim))


## Get mean and standard deviation for each cell in the matrices ()
Acell_values <- Amats %>%
  summarise(mean = apply(., 1, mean),
            sd = apply(., 1, sd))%>%
  mutate(across(everything(), ~ replace(., is.na(.), 0)))

Ucell_values <- logit(Umats) %>%
  mutate(across(everything(), ~ case_when(. > 12 ~ 12,
                                          . < -12 ~ -12,
                                          between(., -12, 12) ~ .))) %>%
  summarise(mean = apply(., 1, mean),
            sd = apply(., 1, sd))%>%
  mutate(across(everything(), ~ replace(., is.na(.), 0)))

Fcell_values <- log(Fmats) %>%
  mutate(across(everything(), ~ case_when(. < -12 ~ -12,
                                          . >= -12 ~ .))) %>%
  summarise(mean = apply(., 1, mean),
            sd = apply(., 1, sd))%>%
  mutate(across(everything(), ~ replace(., is.na(.), 0)))



#-----------------------------------------------------------
# Simulate with the mpm's
#-----------------------------------------------------------

## create stochastic climate sensitive mpm's from the data
## I assume that the sd in the cell variables is the climate effects

### Create environmental sequence ----------------------------
create_seq <- function(n_it, clim_sd, clim_corr, lag) { 
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
                   lagged = lagged)
  return(df)
}

# "Stocastic" mpm -----------------------------------

st.lamb <- function(growth, reproduction){
  
  n_it = length(growth)
  
  env <- data.frame(growth = growth,
                    reproduction = reproduction)
  
  env <- env[complete.cases(env), ]
  
  env <- as.list(as.data.frame(t(env)))
  
  ### Get all mpm's
  mats <- lapply(env, function(x) mpm(U_clim = x[1],F_clim = x[2]))
  
  
  lambda = stoch.growth.rate(mats, maxt = n_it, verbose = F)$sim
  
  return(lambda)
}

clim_sd <- rep(seq(from = 0, to = 2, length.out = 10), 90)
clim_corr <- rep(rep(c(-0.9,0,0.9), each = 10), 30)

lag_clim <- lapply(as.list(c(1:900)), function(x) create_seq(10000, clim_sd = clim_sd[x], clim_corr = clim_corr[x], lag = 1) %>% filter(!is.na(.)))


# species specific mpm
  mpm <- function(U_clim, F_clim, signal_strength_U = 1, signal_strength_F = 1) {
    
    Umat <- inv_logit(Ucell_values$mean + Ucell_values$sd * (U_clim * signal_strength_U) + rnorm(1, sd = 1-signal_strength_U))
    Fmat <- exp(Fcell_values$mean + Fcell_values$sd * (F_clim * signal_strength_F))
    
    Amat <- Umat + Fmat
    
    mpm <- matrix(Amat, nrow = dim)
    
    return(mpm)  
  }

  
  #### Lagged effect in U or F matrix
  
  lag_u <- pblapply(lag_clim, 
                    function(x) st.lamb(growth = x$lagged,
                                        reproduction = x$recent)
  )
  
  lag_f <- pblapply(lag_clim, 
                    function(x) st.lamb(growth = x$recent,
                                        reproduction = x$lagged)
  )
  
  lag_n <- pblapply(lag_clim, 
                    function(x) st.lamb(growth = x$recent,
                                        reproduction = x$recent)
  )
  
  lag_uf <- list("Umatrix" = data.frame(lambda = lag_u %>% unlist,
                                        type = rep("Umatrix", length(clim_sd)),
                                        clim_sd = clim_sd,
                                        clim_auto = clim_corr),
                 "Fmatrix" = data.frame(lambda = lag_f %>% unlist,
                                        type = rep("Fmatrix", length(clim_sd)),
                                        clim_sd = clim_sd,
                                        clim_auto = clim_corr),
                 "None" = data.frame(lambda = lag_n %>% unlist,
                                     type = rep("None", length(clim_sd)),
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
