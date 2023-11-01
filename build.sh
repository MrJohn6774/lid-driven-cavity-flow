#!/bin/sh

# load modules
module load gcc
module load openmpi-gcc

mkdir -p build

gfortran -o build/solver_serial -O3 app/main.f90
mpif90 -o build/solver_parallel -O3 app/major-domain.f90
