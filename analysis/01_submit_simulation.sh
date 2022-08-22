#!/bin/bash

#SBATCH -D /work/evers/

#SBATCH -o /work/%u/%x-%A-%a.log

#Specify job name
#SBATCH -J simulate_lag

#Resources
# max running time
#SBATCH -t 00:45:00

# memory per core (hard limit)
#SBATCH --mem-per-cpu=4G

# create output direcotry per job
OUTPUT_DIR="/work/$USER/$SLURM_JOB_NAME-$SLURM_JOB_ID"
mkdir -p $OUTPUT_DIR

# location lifehistories datafile 
lifehist="/data/lagged/data/life_histories_df.csv"
# Load modules
module load foss/2019b R/4.0.0-2

sigstrength=$1

export MC_CORES=${SLURM_CPUS_PER_TASK:-1}
 
Rscript --vanilla "$HOME"/lagged_buffering/analysis/01_simulations_same_directional_responses.R \
"$sigstrength" \
"$OUTPUT_DIR" \
"$lifehist"

