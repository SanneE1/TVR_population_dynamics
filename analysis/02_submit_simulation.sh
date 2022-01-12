#!/bin/bash

#SBATCH -D /work/evers/

#SBATCH -o /work/%u/%x-%A.log

#Specify job name
#SBATCH -J simulate_diffdirection

#Resources
# max running time
#SBATCH -t 24:00:00

# memory per core (hard limit)
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=25

# Load modules
module load foss/2019b R/4.0.0-2

sigstrength=$1
output_dir="/gpfs1/data/lagged/results/02_Simulations_mpm_opposing_directions"

export MC_CORES=${SLURM_CPUS_PER_TASK:-1}
 
Rscript "$HOME"/lagged_buffering/analysis/02_simulations_opposing_directional_responses \
"$sigstrength" \
"output_dir"
