//
// Created by Linus on 11.12.20.
//

#ifndef NBODY_MATH_UTILS_H
#define NBODY_MATH_UTILS_H

#include "body.h"
#include <vector>
#include <mpi.h>

using namespace std;

void calc_direct_force(vector<Body> &bodies, int a, int b);

void leapfrog(vector<Body> &bodies, double dt, int num_procs, int my_id, MPI_Datatype mpi_body_type);

double get_dt(vector<Body> &particles, int a, int b);


#endif //NBODY_MATH_UTILS_H
