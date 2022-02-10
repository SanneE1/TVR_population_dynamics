#!/bin/bash

#SBATCH -D /work/evers/

#SBATCH -o /work/%u/%x-%A.log

#Specify job name
#SBATCH -J compadre_sim

#Resources
# max running time
#SBATCH -t 24:00:00

# memory per core (hard limit)
#SBATCH --mem-per-cpu=18G
#SBATCH --cpus-per-task=21

# create output direcotry per job
OUTPUT_PATH="/work/$USER/$SLURM_JOB_NAME-$SLURM_JOB_ID"
mkdir -p "$OUTPUT_PATH"
mkdir -p "$OUTPUT_PATH/rds"

# Load modules
module load foss/2019b R/4.0.0-2

export MC_CORES=${SLURM_CPUS_PER_TASK:-1}

Rscript --vanilla "$HOME"/lagged_buffering/analysis/05_COMPADRE_studies/05_get_population_level_sd_auto.R \
"$OUTPUT_PATH"
