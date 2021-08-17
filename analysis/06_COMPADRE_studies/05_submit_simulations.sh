#!/bin/bash

#SBATCH -D /work/evers/

#SBATCH -o /work/%u/%x-%A.log

#Specify job name
#SBATCH -J simulate_diffdirection

#Resources
# max running time
#SBATCH -t 12:00:00

# memory per core (hard limit)
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=25

# Load modules
module load foss/2019b R/4.0.0-2

sigstrength=$1

export MC_CORES=${SLURM_CPUS_PER_TASK:-1}

Rscript "$HOME"/lagged_buffering/analysis/06_COMPADRE_studies/05_get_population_level_sd_auto.R \
"$sigstrength"