//
// Created by Linus on 11.12.20.
//

#ifndef NBODY_IO_H
#define NBODY_IO_H


#include <vector>
#include <string>
#include "body.h"

using namespace std;

void get_initial_values(string *path, uint64_t *steps, double *dt, uint32_t *save_interval, uint32_t *ignore_bodies,
                        float *G);

void write_file(vector<Body> bodies, char filename[], double dt, double t);

vector<Body> read_initial(string path, float G);


#endif //NBODY_IO_H
