#!/bin/bash
#PBS -P PANDORA
#PBS -l select=1:ncpus=1:mem=2GB
#PBS -l walltime=00:01:00

cd "$PBS_O_WORKDIR"
build/solver_serial
