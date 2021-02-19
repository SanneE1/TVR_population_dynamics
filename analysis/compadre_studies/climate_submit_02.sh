#!/bin/bash

#SBATCH -D /work/evers/

#SBATCH -o /work/%u/%x-%A-%a.log

#Specify job name
#SBATCH -J mpms_climate

#Resources
# max running time
#SBATCH -t 30

# memory per core (hard limit)
#SBATCH --mem-per-cpu=8G

# Load modules
module load foss/2019b R/4.0.0-2

other_dir=$1

Rscript "$HOME"/lagged_buffering/analysis/compadre_studies/retrieve_climate_02.R \
"$other_dir"
