#!/bin/bash

#$ -S /bin/bash

#$ -wd /work/$USER

#$ -o /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.log
#$ -j y

#Specify job name
#$ -N ipmr_lagged

#Resources
# max running time


# memory per core (hard limit)
#$ -l h_vmem=40G

# Array numbers 
#$ -t 1-960

#needed when submitting a non-parallel job
#$ -binding linear:1

#create a single output directory per job
output_dir="/work/$USER/$JOB_NAME-$JOB_ID"
mkdir -p "$output_dir"


module load foss/2018b R/3.5.1

function=$1
params=$2
climparams=$3
iterations=$4
meshpoints=$5
output="$output_dir"/${JOB_NAME}_${function}_${JOB_ID}_$SGE_TASK_ID.rds

Rscript "$HOME"/lagged_buffering/analysis/simulations/Full simulation/Eve_test_assumptions.R \
--function="$function" \
--iterations="$iterations" \
--meshpoints="$meshpoints" \
"$output" \
"$params" \
"$climparams"
