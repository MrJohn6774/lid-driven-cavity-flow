## Navier Stoke Solver

### Project setup
- a linux environment
- clone the repository
- install gfortran
```bash
sudo apt update && apt install -y gfortran
```

### Development
Instead of the traditional make/cmake build system, this project uses Fortran Package Manager (FPM). The executable is already included in the root of this repository so you don't need to install one.

To build and run the executable
```bash
# run the serial solver
./fpm run solver_serial

# run the parallel solver with maximum number of available processors
./fpm run solver_parallel
```

To clean the build files and folders
```bash
./fpm clean -all
```

More info on FPM can be found [here](https://github.com/fortran-lang/fpm)

Debugging is supported in VS Code. Install `gdb` and the recommended VS Code extensions

### Resources
- [Learn Fortran](https://fortran-lang.org/learn/)
- [Git User Guide](https://sydneyuni.atlassian.net/wiki/spaces/RC/pages/229277917/Git+-+What+you+need+to+know)
- [Artemis HPC Docs](https://sydneyuni.atlassian.net/wiki/spaces/RC/pages/1033929078/Artemis+HPC+documentation)