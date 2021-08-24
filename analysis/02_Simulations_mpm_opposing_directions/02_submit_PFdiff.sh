#!/bin/bash

#SBATCH -D /work/evers/

#SBATCH -o /work/%u/%x-%A.log

#Specify job name
#SBATCH -J simulate_diffdirection

#Resources
# max running time
#SBATCH -t 01:00:00

# memory per core (hard limit)
#SBATCH --mem-per-cpu=10G

# Load modules
module load foss/2019b R/4.0.0-2

sigstrength=$1
sourcefile="$HOME/lagged_buffering/analysis/01_Simulations_mpm_same_direction/simulate_mpm.R"

export MC_CORES=${SLURM_CPUS_PER_TASK:-1}
 
Rscript "$HOME"/lagged_buffering/analysis/02_Simulations_mpm_opposing_directions/02_simulations_with_diff_in_PF.R \
"$sigstrength" \
"$sourcefile"
