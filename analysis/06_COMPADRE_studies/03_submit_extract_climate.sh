#!/bin/bash

#SBATCH -D /work/evers/

#SBATCH -o /work/%u/%x-%A.log

#Specify job name
#SBATCH -J extract_CHELSA

#Resources
# max running time
#SBATCH -t 30:00:00

# memory per core (hard limit)
#SBATCH --mem-per-cpu=8G


# create output direcotry per job
OUTPUT_PATH="/work/$USER/$SLURM_JOB_NAME-$SLURM_ARRAY_JOB_ID"
mkdir -p "$OUTPUT_PATH"

# Load modules
module load foss/2019b R/4.0.0-2

Cdownloads=$1
LatLonfile=$2

Rscript "$HOME"/lagged_buffering/analysis/06_COMPADRE_studies/03_extract_climate.R \
"$Cdownloads" \
"$LatLonfile" \
"$OUTPUT_PATH"
