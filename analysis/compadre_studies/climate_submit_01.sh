#!/bin/bash

#SBATCH -D /work/evers/

#SBATCH -o /work/%u/%x-%A-%a.log

#Specify job name
#SBATCH -J mpms_climate

#Resources
# max running time
#SBATCH -t 15

# memory per core (hard limit)
#SBATCH --mem-per-cpu=8G

# Array numbers 
#SBATCH -a 1-4284

# create output direcotry per job
OUTPUT_PATH="/work/$USER/$SLURM_JOB_NAME-$SLURM_ARRAY_JOB_ID"
mkdir -p "$OUTPUT_PATH"

# Load modules
module load foss/2019b R/4.0.0-2


Rscript "$HOME"/lagged_buffering/analysis/compadre_studies/retrieve_climate_01.R \
"$OUTPUT_PATH"
