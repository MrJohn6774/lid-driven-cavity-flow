#!/bin/bash
#PBS -P PANDORA
#PBS -l select=1:ncpus=6:mpiprocs=6:mem=2GB
#PBS -l walltime=00:01:00

cd "$PBS_O_WORKDIR"
module load python/3.6.5
module load gcc
module load openmpi-gcc

python3 ./mpirun.py
