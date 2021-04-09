#!/bin/bash

#SBATCH -D /work/evers/

#SBATCH -o /work/%u/%x-%A-%a.log

#Specify job name
#SBATCH -J mpms_compadre

#Resources
# max running time


# memory per core (hard limit)
#SBATCH --mem-per-cpu=8G

# Array numbers 

# create output direcotry per job
OUTPUT_PATH="/work/$USER/$SLURM_JOB_NAME-$SLURM_ARRAY_JOB_ID"
mkdir -p "$OUTPUT_PATH"

# Load modules
module load foss/2019b R/4.0.0-2


Rscript "$HOME"/lagged_buffering/analysis/06_COMPADRE_studies/06_simulate_lambda.R \
"/data/lagged/species_authors.csv" \
"$OUTPUT_PATH"
