# n-body-solver in c++

A simple implementation of an n-body-solver written in c++ by using direct force calculation. Leapfrog is used as an integrating scheme.

This program was developed on CLion (with CMake) on macOS and uses [openMPI](https://www.open-mpi.org) to parallelise the calculation.
Run it using command line arguments for number of steps and path of input file:
```
mpirun -np 4 nbody 1000000 ../input/solar_jfc.dat
```
assuming that the binary nbody is located in a bin or cmake-build-debug folder.

This code has also been tested on the Piz Daint supercomputer at CSCS. There is a makefile in the `bin` folder.