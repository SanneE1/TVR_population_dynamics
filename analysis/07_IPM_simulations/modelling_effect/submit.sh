#!/bin/bash

#$ -S /bin/bash

#$ -wd /work/$USER

#$ -o /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.log
#$ -j y

#Specify job name
#$ -N  simulate_effect

#Resources
# max running time


# memory per core (hard limit)
#$ -l h_vmem=1G

# Array numbers 
#$ -t 1-100

#needed when submitting a non-parallel job
#$ -binding linear:1

#create a single output directory per job
output_dir="/work/$USER/$JOB_NAME-$JOB_ID"
mkdir -p "$output_dir"

module load foss/2018b R/3.5.1

output="$output_dir"/${JOB_NAME}_${JOB_ID}_$SGE_TASK_ID.rds


Rscript "$HOME"/lagged_buffering/analysis/simulations/modelling effect/simulate_modelling_effect_EVE.R \
"$output"
