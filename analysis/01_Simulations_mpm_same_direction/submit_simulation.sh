#!/bin/bash

#SBATCH -D /work/evers/

#SBATCH -o /work/%u/%x-%A.log

#Specify job name
#SBATCH -J simulate_lag

#Resources
# max running time
#SBATCH -t 30:00:00

# memory per core (hard limit)
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=10
 
# create output direcotry per job
OUTPUT_PATH="/work/$USER/$SLURM_JOB_NAME-$SLURM_JOB_ID"
mkdir -p "$OUTPUT_PATH"

# Load modules
module load foss/2019b R/4.0.0-2

sigstrength=$1

export MC_CORES=${SLURM_CPUS_PER_TASK:-1}
 
Rscript "$HOME"/lagged_buffering/analysis/01_Simulations_mpm_same_direction/simulate_mpm.R \
"$sigstrength" \
"$OUTPUT_PATH"