suppressPackageStartupMessages(library(dplyr))
# suppressPackageStartupMessages(library(leaflet))
# suppressPackageStartupMessages(library(measurements))
suppressPackageStartupMessages(library(popbio))
suppressPackageStartupMessages(library(pbapply))
# suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(faux))
suppressPackageStartupMessages(library(patchwork))


start <- Sys.time()
start
#  ----------------------------------------------------------------------------------------------------------------------------
# get taskID
#  ----------------------------------------------------------------------------------------------------------------------------

taskID <- as.integer(Sys.getenv("SGE_TASK_ID"))

#-----------------------------------------------------------
# Retrieve data from actual species with 2x2 matrices
#-----------------------------------------------------------

# load most current compadre database
load('/data/lagged/COMPADRE_v.X.X.X.4.RData')

# Select mpm's 

## Subsection from Dalgleish
id <- which(compadre$metadata$MatrixDimension == 4)


## Get all matrices for the selected id's, grouped by species
species <- unique(compadre$metadata$SpeciesAuthor[id])

i = species[taskID]

id2 <- id[which(compadre$metadata$SpeciesAuthor[id] == i & compadre$metadata$MatrixComposite[id] == "Individual")]
# Retrieve all full matrices
Amats <- lapply(as.list(id2), function(x) as.vector(compadre$mat[x][[1]]$matA)) %>% bind_cols
## matrix cell position in vector -> c([1,1], [2,1], etc., [1,2], [2,2], etc.)
Umats <- lapply(as.list(id2), function(x) as.vector(compadre$mat[x][[1]]$matU)) %>% bind_cols 
Fmats <- lapply(as.list(id2), function(x) as.vector(compadre$mat[x][[1]]$matF)) %>% bind_cols

if(length(unique(compadre$metadata$MatrixDimension[id2])) != 1) {
  stop("different sized matrices in same study")
}

dim <- unique(compadre$metadata$MatrixDimension[id2])



## Get mean and standard deviation for each cell in the matrices ()
Acell_values <- Amats %>%
  summarise(mean = apply(., 1, mean),
            sd = apply(., 1, sd))

Ucell_values <- Umats %>%
  summarise(mean = apply(., 1, mean),
            sd = apply(., 1, sd))

Fcell_values <- Fmats %>%
  summarise(mean = apply(., 1, mean),
            sd = apply(., 1, sd))



#-----------------------------------------------------------
# Simulate with the mpm's
#-----------------------------------------------------------

## create stochastic climate sensitive mpm's from the data
## I assume that the sd in the cell variables is the climate effects

logit <- function(x) log(x/(1-x))

inv_logit <- function(x) {
  return(
    1/(1 + exp(-(x)))
  )
}

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
  mats <- lapply(env, function(x) mpm(x[1], x[2]))
  
  
  lambda = stoch.growth.rate(mats, maxt = n_it, verbose = F)$sim
  
  return(lambda)
}

clim_sd <- rep(seq(from = 0, to = 2, length.out = 10), 90)
clim_corr <- rep(rep(c(-0.9,0,0.9), each = 10), 30)

lag_clim <- lapply(as.list(c(1:900)), function(x) create_seq(10000, clim_sd = clim_sd[x], clim_corr = clim_corr[x], lag = 1) %>% filter(!is.na(.)))


# species specific mpm
  mpm <- function(growth, reproduction) {
    
    Umat <- Ucell_values$mean + Ucell_values$sd * growth
    Fmat <- Fcell_values$mean + Fcell_values$sd * reproduction
    
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
                    function(x) st.lamb(env_growth = x$recent,
                                        env_reproduction = x$lagged)
  )
  
  lag_n <- pblapply(lag_clim, 
                    function(x) st.lamb(env_growth = x$recent,
                                        env_reproduction = x$recent)
  )
  
  lag_uf <- list("Pkernel" = lag_p, "Fkernel" = lag_f, "none" = lag_n)
  
  saveRDS(lag_fp, paste("work/evers/simulations/mpm/mpm_", i, "_laguf.RDS", sep = ""))




#-----------------------------------------------------------
# plot results
#-----------------------------------------------------------

  laguf_df <- data.frame(clim_sd = clim_sd,
                         clim_corr = clim_corr,
                         lambda = unlist(lag_pf),
                         type = rep(names(lag_pf), each = length(lag_pf[[1]])))
  lagpf_p <- ggplot(lagpf_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    labs(colour = "Lag type", title = "Lagged climate in P or F") +
    facet_grid(cols = vars(clim_corr)) + 
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))
  


    plot_annotation(title = i) 
  
  ggsave(lag, filename = paste0("work/evers/simulations/compadre_", i, "plots.png"))
  
  
}



Sys.time() - start