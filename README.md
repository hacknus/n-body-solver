# n-body-solver in c++

A simple implementation of an n-body-solver written in c++ by using direct force calculation. Leapfrog is used as an integrating scheme.

This program was developed on CLion (with CMake) on macOS and uses [openMPI](https://www.open-mpi.org) to parallelise the calculation.
Run it as follows:
```
mpirun -np 4 nbody
```
The program will read some configuration parameters from the file `input/input.conf` which looks as follows:
```
path=../input/solar_jfc_rev.dat
steps=100000
dt=86400
save_interval=10
ignore_bodies=62
G=6.67408e-11
```
If `dt=0` is specified, the program will calculate `dt` according to the internal function.
assuming that the binary nbody is located in a bin or cmake-build-debug folder.
The parameter `ignore_bodies` determines how many bodies are ignored when calculating the force on other bodies.
This is useful for small bodies like satellites, comets or asteroids that do not exert a force on large bodies such as planets.
For `ignore_bodies=62` the last 62 bodies of the dataset are ignored. 
This code has also been tested on the Piz Daint supercomputer at CSCS. There is a makefile in the `bin` folder.