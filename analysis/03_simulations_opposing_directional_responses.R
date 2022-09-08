#
# Script to run the main simulations on effect of both lagged and recent responses to climate 
# in a matrix population model, where climate responses are in OPPOSING direction
#
# This script includes arguments so it can be submitted on the UFZ HPC "EVE"
# Also see "submit_simulation.sh" in the same folder for the job submission script
rm(list=ls())

start <- Sys.time()

set.seed(2)

library(dplyr)

# Get required arguments supplied during job submission
args = commandArgs(trailingOnly = T)

if(length(args)!=3) {  
stop("Provide (only) signal strength for the analysis", call.=F)
}

# Get task ID from array job -> determines for which lifehistory the script will run
taskID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Set the proportion of variance (p) explained by the climate driver in the temporal sequence.
# Here set to 1, 0.5, 0.25 or 0.05. See manuscript for more info
i = as.numeric(args[1])

# File location of the simulation script with SAME DIRECTIONAL responses to climate
output_dir <- args[2]

# location of lifehistories file
lh_location <- args[3]

# Number of iterations to run the simulations
n_it = 10000

# Print information for the log file
print(paste("sig.strength =", i))
print(paste("output dir =", output_dir))
print(paste("location of lifehistories file = ", lh_location))
print(paste("#iterations =", n_it))

# Load and prep lifehistories table  ------------------------------------------------------------------------------------

lh_df <- read.csv(lh_location)
lh_df <- lh_df[taskID,]

print("life history currently being used")
print(lh_df)

# Source script with R functions required

source("/gpfs0/home/evers/lagged_buffering/R/lambda_functions.R")

#---------------------------------------------------------------------------------------------------------
# Start of analyses
#---------------------------------------------------------------------------------------------------------

## set up sequences of climate standard deviation and autocorrelation values so that there are 30 duplicates 
## for each combination of sd & autocorrelation value

df_clim <- expand.grid(clim_sd = seq(from = 0.01, to = 1, length.out = 5),
                       clim_auto =  c(-0.9,0,0.9),
                       rep = c(1:30))

#### Create main temporal sequences
lag_clim <- lapply(as.list(c(1:nrow(df_clim))), 
                   function(x) create_seq(n_it, 
                                          clim_sd = df_clim$clim_sd[x], 
                                          clim_auto = df_clim$clim_auto[x], 
                                          lag = 1))

#### Lagged effect between U & F matrices
lag_p <- lapply(lag_clim,
                function(x) st.lamb_o(mpm_df = lh_df,
                                    env_surv = x$lagged,
                                    env_growth = x$lagged,
                                    env_reproduction = x$recent,
                                    clim_sd = sd(x$recent, na.rm = T),
                                    clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                    sig.strength = i)
) %>% bind_rows() %>%
  tibble::add_column(lag_type = "Umatrix")

print("done w/ P sim")


lag_n2 <- lapply(lag_clim,
                 function(x) st.lamb_o(mpm_df = lh_df,
                                     env_surv = x$recent,
                                     env_growth = x$recent,
                                     env_reproduction = x$recent,
                                     clim_sd = sd(x$recent, na.rm = T),
                                     clim_auto = acf(x$recent, plot = F, na.action = na.pass)$acf[2],
                                     sig.strength = i)
) %>% bind_rows() %>%
  tibble::add_column(lag_type = "none")

print("done w/ none sim")

lag_df <- rbind(lag_p, lag_n2)

write.csv(lag_df, file.path(output_dir, paste("mpm", i, taskID, "lagf_oposing.csv", sep = "_")), row.names = F)

Sys.time() - start
