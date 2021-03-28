//
// Created by Linus on 11.12.20.
//

#ifndef NBODY_IO_H
#define NBODY_IO_H


#include <vector>
#include <string>
#include "body.h"

using namespace std;

void get_initial_values(string *path, unsigned long long int *steps, double *dt, unsigned long int *save_interval,
                        unsigned long int *ignore_bodies, float *G, float *softening);

void write_file(vector<Body> bodies, char filename[], double dt, double t);

vector<Body> read_initial(string path, float G);


#endif //NBODY_IO_H
