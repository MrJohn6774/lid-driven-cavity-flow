## Navier Stoke Solver

### Project setup
- a linux environment
- clone the repository
- install gfortran libopenmpi-dev
```bash
sudo apt update && apt install -y gfortran libopenmpi-dev
```

### Development
Instead of the traditional make/cmake build system, this project uses Fortran Package Manager (FPM). The executable is already included in the root of this repository so you don't need to install one.

To build the binary on Artemis
```bash
chmod u+x *.sh
./build.sh

# or
qsub ./job_build.sh
```

To build and run the executable (./fpm not available on Artemis because the hpc kernel is too old)
```bash
chmod u+x ./fpm

# run the serial solver
./fpm run solver_serial

# run the parallel solver
# npx * npy must equal to the number of max avail processors
./fpm run solver_parallel -- <npx> <npy> 21 21

# run the parallel solver with mpirun
mpirun -np 4 build/gfortran_6FB96CAE7088C0B9/app/solver_parallel 2 2 21 21

# run parallel solver in batch with python
python3 mpirun.py
```

To clean the build files and folders
```bash
./fpm clean -all
```

To generate graphs
```bash
pip3 install numpy pandas matplotlib

python3 mpirun.py
python3 plot_graphs.py    # graphs will be saved to outputs folder
```

More info on FPM can be found [here](https://github.com/fortran-lang/fpm)

Debugging is supported in VS Code. Install `gdb` and the recommended VS Code extensions

### Resources
- [Learn Fortran](https://fortran-lang.org/learn/)
- [Git User Guide](https://sydneyuni.atlassian.net/wiki/spaces/RC/pages/229277917/Git+-+What+you+need+to+know)
- [Artemis HPC Docs](https://sydneyuni.atlassian.net/wiki/spaces/RC/pages/1033929078/Artemis+HPC+documentation)