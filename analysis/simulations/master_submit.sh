#!/bin/bash

function=$1

$DEPENDENCY=$(qsub -terse analysis/simulations/submit_test.sh "$function")

qsub -hold_jid $DEPENDENCY analysis/simulations/merge_submit.sh /work/evers/ipmr_lagged-$DEPENDENCY/