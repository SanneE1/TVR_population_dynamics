suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(popbio))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(faux))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(Rage))

start <- Sys.time()
start

taskID <- as.integer(Sys.getenv("SGE_TASK_ID"))


### Species_author for this specific run 

args() = commandArgs(TRUE)  ## should be a location for the csv with Species_author

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

if(dim != 1) {
  stop("different sized matrices in same study")
}


## Get mean and standard deviation for each cell in the matrices ()
Acell_values <- Amats %>%
  summarise(mean = apply(., 1, mean),
            sd = apply(., 1, sd))%>%
  mutate(across(everything(), ~ replace(., is.na(.), 0)))

Ucell_values <- (Umats) %>%
  mutate(across(everything(), ~ case_when(. > 6 ~ 6,
                                          . < -6 ~ -6,
                                          between(., -6, 6) ~ .))) %>%
  summarise(mean = apply(., 1, mean),
            sd = apply(., 1, sd))%>%
  mutate(across(everything(), ~ replace(., is.na(.), 0)))

Fcell_values <- (Fmats) %>%
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
    
    Umat <- (Ucell_values$mean + Ucell_values$sd * (U_clim * signal_strength_U))
    Fmat <- (Fcell_values$mean + Fcell_values$sd * (F_clim * signal_strength_F))
    
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
  
  lag_uf <- list("Umatrix" = lag_u, "Fmatrix" = lag_f, "None" = lag_n)
  
  saveRDS(lag_uf, paste("work/evers/simulations/mpm/mpm_", i, "_laguf.RDS", sep = ""))




# #-----------------------------------------------------------
# # plot results
# #-----------------------------------------------------------
# 
#   laguf_df <- data.frame(clim_sd = clim_sd,
#                          clim_corr = clim_corr,
#                          lambda = unlist(lag_uf),
#                          type = rep(names(lag_uf), each = length(lag_uf[[1]])))
#   lagpf_p <- ggplot(laguf_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
#     labs(colour = "Lag type", title = "Lagged climate in P or F") +
#     facet_grid(cols = vars(clim_corr)) + 
#     theme(legend.position = "bottom",
#           plot.title = element_text(hjust = 0.5)) +
#     plot_annotation(title = i) 
#   
#   ggsave(lag, filename = paste0("work/evers/simulations/compadre_", i, "plots.png"))
#   
  




Sys.time() - start