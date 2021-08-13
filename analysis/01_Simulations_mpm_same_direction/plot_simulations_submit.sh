#!/bin/bash

#SBATCH -D /work/evers/

#SBATCH -o /work/%u/%x-%A.log

#Specify job name
#SBATCH -J plot_simulation

#Resources
# max running time
#SBATCH -t 0:15:00

# memory per core (hard limit)
#SBATCH --mem-per-cpu=20G


# create output direcotry per job

# Load modules
module load foss/2019b R/4.0.0-2


Rscript "$HOME"/lagged_buffering/analysis/01_Simulations_mpm_same_direction/plot_simulations.R \

