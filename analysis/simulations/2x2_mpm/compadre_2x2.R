suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(measurements))
suppressPackageStartupMessages(library(popbio))
suppressPackageStartupMessages(library(pbapply))
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

# Select 2x2 mpm's  ---- I should probably get more than only the Dalgleish studies, but just taking all is not a good idea either

## Subsection from Dalgleish
id <- which(grepl("Dalgleish", compadre$metadata$Authors) & compadre$metadata$MatrixDimension == 2)
## All 2x2 mpm's 
# id <- which(compadre$metadata$MatrixDimension == 2)


## Get all matrices for the selected id's, grouped by species
species <- unique(compadre$metadata$SpeciesAuthor[id])

Amat_list <- as.list(rep(NA, length(species)))
Umat_list <- as.list(rep(NA, length(species)))
Fmat_list <- as.list(rep(NA, length(species)))
names(Amat_list) <- species
names(Umat_list) <- species
names(Fmat_list) <- species

i = species[taskID]
i

id2 <- id[which(compadre$metadata$SpeciesAuthor[id] == i)]
# Retrieve all full matrices
Amat_list[[i]] <- lapply(as.list(id2), function(x) as.vector(compadre$mat[x][[1]]$matA))  
## matrix cell position in vector -> c([1,1], [2,1], [1,2], [2,2])
Umat_list[[i]] <- lapply(as.list(id2), function(x) as.vector(compadre$mat[x][[1]]$matU))  
Fmat_list[[i]] <- lapply(as.list(id2), function(x) as.vector(compadre$mat[x][[1]]$matF))  



## Get mean and standard deviation for each cell in the 2x2 matrices (small)
Acell_values <- lapply(Amat_list, function(x) as.data.frame(matrix(unlist(x), ncol = 4, byrow = T)) %>%
                         `colnames<-`(c("S->S", "S->B", "B->S", "B->B")) %>%
                         summarise(mean = apply(., 2, mean),
                                   sd = apply(., 2, sd)) %>%
                         `rownames<-`(c("S->S", "S->B", "B->S", "B->B"))
) 

Ucell_values <- lapply(Umat_list, function(x) as.data.frame(matrix(unlist(x), ncol = 4, byrow = T)) %>%
                         `colnames<-`(c("S->S", "S->B", "B->S", "B->B")) %>%
                         summarise(mean = apply(., 2, mean),
                                   sd = apply(., 2, sd)) %>%
                         `rownames<-`(c("S->S", "S->B", "B->S", "B->B"))
) 

Fcell_values <- lapply(Fmat_list, function(x) as.data.frame(matrix(unlist(x), ncol = 4, byrow = T)) %>%
                         `colnames<-`(c("S->S", "S->B", "B->S", "B->B")) %>%
                         summarise(mean = apply(., 2, mean),
                                   sd = apply(., 2, sd)) %>%
                         `rownames<-`(c("S->S", "S->B", "B->S", "B->B"))
) 


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

st.lamb <- function(env_surv, env_growth, env_reproduction) {
  
  n_it = length(env_surv)
  
  env <- data.frame(survival = env_surv,
                    growth = env_growth,
                    reproduction = env_reproduction)
  
  env <- env[complete.cases(env), ]
  
  env <- as.list(as.data.frame(t(env)))
  
  ### Get all mpm's
  mats <- lapply(env, function(x) mpm(x[1], x[2], x[3]))
  
  
  df = stoch.growth.rate(mats, maxt = n_it, verbose = F)$sim
  
  return(df)
}

clim_sd <- rep(seq(from = 0, to = 2, length.out = 10), 90)
clim_corr <- rep(rep(c(-0.9,0,0.9), each = 10), 30)

print("creating lagged climate sequences")
lag_clim <- pblapply(as.list(c(1:900)), function(x) create_seq(10000, clim_sd = clim_sd[x], clim_corr = clim_corr[x], lag = 1) %>% filter(!is.na(.)))
print("done")

# species specific mpm
  
  mpm <- function(survival, growth, reproduction) {
    ## Basic mpm
    mpm <- matrix(0, nrow = 2, ncol = 2)
    
    # survival (but actually more like survival and staying in the same state)
    mpm[1,1] <- Acell_values[[i]][1,1] + survival * Ucell_values[[i]][1,2]         # + recruitment by small into small? 
    
    mpm[2,2] <- Acell_values[[i]][4,1] + survival * Ucell_values[[i]][4,2]                   # + reproduction by large into large?
    
    
    #growth (but actually survival and moving to next state)
    mpm[2,1] <- Acell_values[[i]][2,1] + growth * Ucell_values[[i]][2,2]                   # + recruitment by small into big 
    
    
    # reproduction 
    mpm[1,2] <- Acell_values[[i]][3,1] + reproduction * Fcell_values[[i]][3,2]                   # + shrinkage from big to small?
    
    
    return(mpm)  
  }
  
  
  #### Lagged effect within P "functions"
print("lagged effect in growth")
  lag_g <- pblapply(lag_clim,
                    function(x) st.lamb(env_surv = x$recent,
                                        env_growth = x$lagged,
                                        env_reproduction = rep(0,length(x$recent)))
  )
print("lagged effect in survival")  
  lag_s <- pblapply(lag_clim,
                    function(x) st.lamb(env_surv = x$lagged,
                                        env_growth = x$recent,
                                        env_reproduction = rep(0,length(x$recent)))
  )
print("lagged effect control")  
  lag_n <- pblapply(lag_clim,
                    function(x) st.lamb(env_surv = x$recent,
                                        env_growth = x$recent,
                                        env_reproduction = rep(0,length(x$recent)))
  )
  
  lag <- list("growth" = lag_g, "survival" = lag_s, "none" = lag_n)
  
  saveRDS(lag, paste("/work/evers/simulations/mpm_", i, "_lag.RDS", sep = ""))
  
  #### Lagged effect
 print("lagged effect in P kernel") 
  lag_p <- pblapply(lag_clim, 
                    function(x) st.lamb(env_surv = x$lagged,
                                        env_growth = x$lagged,
                                        env_reproduction = x$recent)
  )
 print("lagged effect in F kernel") 
  lag_f <- pblapply(lag_clim, 
                    function(x) st.lamb(env_surv = x$recent,
                                        env_growth = x$recent,
                                        env_reproduction = x$lagged)
  )
 print("lagged effect control") 
  lag_n <- pblapply(lag_clim, 
                    function(x) st.lamb(env_surv = x$recent,
                                        env_growth = x$recent,
                                        env_reproduction = x$recent)
  )
  
  lag_fp <- list("Pkernel" = lag_p, "Fkernel" = lag_f, "none" = lag_n)
  
  saveRDS(lag_fp, paste("/work/evers/simulations/mpm_", i, "_lagfp.RDS", sep = ""))
  





#-----------------------------------------------------------
# plot results
#-----------------------------------------------------------
print("start plotting")
# Summary plot function for easy looping
sum_plot <- function(lag, lag_pf){
  ### Climate variables
  clim_sd <- rep(seq(from = 0, to = 2, length.out = 10), 90)
  clim_corr <- rep(rep(c(-0.9,0,0.9), each = 10), 30)
  
  lag_df <- data.frame(clim_sd = clim_sd,
                       clim_corr = clim_corr,
                       lambda = unlist(lag),
                       type = rep(names(lag), each = length(lag[[1]])))
  lag_p <- ggplot(lag_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    labs(colour = "Lag type", title = "Lagged climate in growth or survival") +
    facet_grid(cols = vars(clim_corr)) + 
    scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + theme(legend.position = "bottom",
                                                                             plot.title = element_text(hjust = 0.5))
  
  lagpf_df <- data.frame(clim_sd = clim_sd,
                         clim_corr = clim_corr,
                         lambda = unlist(lag_pf),
                         type = rep(names(lag_pf), each = length(lag_pf[[1]])))
  lagpf_p <- ggplot(lagpf_df) + geom_smooth(aes(x = clim_sd, y = lambda, colour = as.factor(type)))+ 
    labs(colour = "Lag type", title = "Lagged climate in P or F") +
    facet_grid(cols = vars(clim_corr)) + 
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))
  
  return(list(lag = lag_p,
              lagpf = lagpf_p))
}

 
  ### Load results
  
  lag <- readRDS(paste0("/work/evers/simulations/mpm_", i, "_lag.RDS"))
  lagpf <- readRDS(paste0("/work/evers/simulations/mpm_", i, "_lagfp.RDS"))
  
  plots <- sum_plot(lag, lagpf)
  
  # mat <- matrix(nrow = 2, ncol = 2, dimnames = list(c("small", "big"), c("small", "big")))
  
  ## matrix cell position in vector -> c([1,1], [2,1], [1,2], [2,2])
  matss <- paste0(round(Acell_values[[i]]$mean[1], 2), " (", round(Ucell_values[[i]]$sd[1], 2), ")" )
  matsb <- paste0(round(Acell_values[[i]]$mean[2], 2), " (", round(Ucell_values[[i]]$sd[2], 2), ")" )
  matbs <- paste0(round(Acell_values[[i]]$mean[3], 2), " (", round(Ucell_values[[i]]$sd[3], 2), ")" )
  matbb <- paste0(round(Acell_values[[i]]$mean[4], 2), " (", round(Ucell_values[[i]]$sd[4], 2), ")" )
  
  ## Lagged effects
  design <- "AA
           BB
           CD"
  
  lag <- wrap_plots(A = plots$lag, B = plots$lagpf, D = guide_area(),
                    C = wrap_elements(grid::textGrob(
                      substitute(atop(paste(ss, "    ", bs,), paste(sb, "    ", bb)), 
                                 env = list(ss = matss, bs = matbs,
                                            sb = matsb, bb = matbb)))), 
                    design = design, guides = "collect") +
    plot_annotation(title = i) 
  
  ggsave(lag, filename = paste0("work/evers/simulations/compadre_", i, "plots.png"))




Sys.time() - start
