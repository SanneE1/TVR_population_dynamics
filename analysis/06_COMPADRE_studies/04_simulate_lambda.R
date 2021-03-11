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

st.lamb <- function(growth, reproduction, clim_sd, clim_auto, sig.strength){
  
  n_it = length(growth)
  
  env <- data.frame(growth = growth,
                    reproduction = reproduction)
  
  env <- env[complete.cases(env), ]
  
  env <- split(env, seq(nrow(env)))
  
  ### Get all mpm's
  mats <- lapply(env, function(x) mpm(U_clim = x$growth, F_clim = x$reproduction, sig.strength = sig.strength))
  
  df <- data.frame(lambda = stoch.growth.rate(mats, maxt = 1000, verbose = F)$sim,
                   clim_sd = clim_sd,
                   clim_auto = clim_auto)
  
  return(list(df = df,
              mats = data.frame(t(lapply(mats, as.vector) %>% bind_rows)) %>% 
                `colnames<-`(c("1,1", "2,1", "1,2", "2,2"))))
}

clim_sd <- seq(from = 0.01, to = 2, length.out = 10)
clim_corr <- c(-0.9,0,0.9)
sig.strength <- c(1, 0.5, 0.1)

df <- expand.grid(clim_sd, clim_corr, sig.strength) 
df <- do.call("rbind", replicate(30, df, simplify = FALSE)) %>% `names<-`(., c("clim_sd", "clim_auto", "sig.strength"))

lag_clim <- lapply(as.list(c(1:nrow(df))), function(x) create_seq(clim_sd = df$clim_sd[x], clim_corr = df$clim_auto[x], sig.strength = df$sig.strength[x]) %>% filter(complete.cases(.)))


# species specific mpm
  mpm <- function(U_clim, F_clim, sig.strength) {
    
    Umat <- Ucell_values$mean + Ucell_values$sd * (U_clim * sig.strength) + 
      rnorm(n = length(Ucell_values$sd),
        mean = 0, sd = (Ucell_values$sd * (1-sig.strength))) 
    Fmat <- Fcell_values$mean + Fcell_values$sd * (F_clim * sig.strength) + 
      rnorm(n = length(Fcell_values$sd),
            mean = 0, sd = (Fcell_values$sd * (1-sig.strength))) 
    
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
                                        clim_sd = sd(x$recent, na.rm = T),
                                        clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                        sig.strength = unique(x$sig.strength))
  ) 
  
  lag_f <- pblapply(lag_clim, 
                    function(x) st.lamb(growth = x$recent,
                                        reproduction = x$lagged,
                                        clim_sd = sd(x$recent, na.rm = T),
                                        clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                        sig.strength = unique(x$sig.strength))
  ) 
  
  lag_n <- pblapply(lag_clim, 
                    function(x) st.lamb(growth = x$recent,
                                        reproduction = x$recent,
                                        clim_sd = sd(x$recent, na.rm = T),
                                        clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                        sig.strength = unique(x$sig.strength))
  ) 
  
  lag_uf <- list("Umatrix" = lag_u,
                 "Fmatrix" = lag_f,
                 "None" = lag_n)
  
output_dir <- args[2]

n_pop = length(unique(species$MatrixPopulation[which(species$SpeciesAuthor == i)]))

if(length(n_pop == 1)) {
  output_file <- paste0("mpm_", i, "_laguf.RDS")
} else {
  output_file <- paste("mpm", i, j, "laguf.RDS", sep = "_")
}

saveRDS(lag_uf, file.path(output_dir, output_file))

Sys.time() - start
