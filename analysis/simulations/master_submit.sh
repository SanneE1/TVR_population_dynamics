#!/bin/bash

#$ -l h_rt=900

# memory per core (hard limit)
#$ -l h_vmem=8G

#needed when submitting a non-parallel job
#$ -binding linear:1

function=$1

$DEPENDENCY=$(qsub -terse analysis/simulations/submit_test.sh "$function")

qsub -hold_jid "$DEPENDENCY" analysis/simulations/merge_submit.sh /work/evers/ipmr_lagged-$DEPENDENCY/
