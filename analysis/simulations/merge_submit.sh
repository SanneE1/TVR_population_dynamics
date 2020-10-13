#!/bin/bash

#$ -S /bin/bash

#$ -wd /work/$USER

#$ -o /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.log
#$ -j y

#Specify job name
#$ -N merge_ipmr_lagged

#Resources
# max running time
#$ -l h_rt=900

# memory per core (hard limit)
#$ -l h_vmem=40G

# Array numbers 
#$ -t 1
#needed when submitting a non-parallel job
#$ -binding linear:1

module load foss/2018b R/3.5.1


output=$1

Rscript "$HOME"/lagged_buffering/analysis/simulations/merge_array_jobs.R \
"$output"
