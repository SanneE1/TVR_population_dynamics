#!/bin/bash

#$ -S /bin/bash

#$ -wd /work/$USER

#$ -o /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.log
#$ -j y

#Specify job name
#$ -N mpm_lagged

#Resources
# max running time


# memory per core (hard limit)
#$ -l h_vmem=8G

# Array numbers 
#$ -t 1-41

#needed when submitting a non-parallel job
#$ -binding linear:1

#create a single output directory per job

module load foss/2019b R

Rscript "$HOME"/lagged_buffering/analysis/simulations/2x2_mpm/compadre_2x2.R \
