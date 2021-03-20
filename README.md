# n-body-solver in c++

A simple implementation of an n-body-solver written in c++ by using direct force calculation. Leapfrog is used to propagate in time.

This program was developed on CLion on macOS and uses openMPI to parallelise the calculation.
Run it using command line arguments for number of steps and path of input file:
```
mpirun -np 4 nbody 1000000 ../input/solar_jfc.dat
```
