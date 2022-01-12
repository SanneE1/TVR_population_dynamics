#!/bin/bash

#SBATCH -D /data/lagged/CHELSA/

#SBATCH -o /work/%u/%x-%A.log

#Specify job name
#SBATCH -J CHELSA_download

#Resources
# max running time

# memory per core (hard limit)
#SBATCH --mem-per-cpu=1G

wget --no-host-directories --force-directories --input-file="$HOME"/lagged_buffering/analysis/05_COMPADRE_studies/CHELSA_files.txt 
