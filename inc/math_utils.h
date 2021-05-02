//
// Created by Linus on 11.12.20.
//

#ifndef NBODY_MATH_UTILS_H
#define NBODY_MATH_UTILS_H

#include "body.h"
#include <vector>
#include <mpi.h>

using namespace std;

void calc_direct_force(vector<Body> &bodies, vector<Body>::size_type a, vector<Body>::size_type b,
                       unsigned long int ignore_bodies, float G, float softening);

void leapfrog(vector<Body> &bodies, double dt, int num_procs, int my_id, MPI_Datatype MPI_BODY_TYPE,
              unsigned long int ignore_bodies, float G, float softening);

void rk_accel(vector<Body> &bodies, vector<Body>::size_type a, vector<Body>::size_type b,
              unsigned long int ignore_bodies, float G, float softening, vector<double> &ax, vector<double> &ay,
              vector<double> &az);

void rk4(vector<Body> &bodies, double dt, int num_procs, int my_id, MPI_Datatype MPI_BODY_TYPE,
         unsigned long int ignore_bodies, float G, float softening);

void sync(vector<Body> &bodies, int num_procs, int myid, MPI_Datatype MPI_BODY_TYPE);

double get_dt(vector<Body> &particles, vector<Body>::size_type a, vector<Body>::size_type b, float softening);


#endif //NBODY_MATH_UTILS_H
