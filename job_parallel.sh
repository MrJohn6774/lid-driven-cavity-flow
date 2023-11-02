#!/bin/bash
#PBS -P COMPENG
#PBS -l select=1:ncpus=24:mpiprocs=24:mem=4GB
#PBS -l walltime=00:10:00

cd "$PBS_O_WORKDIR"
module load python/3.6.5
module load gcc
module load openmpi-gcc

python3 ./mpirun_artemis.py
