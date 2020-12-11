//
// Created by Linus on 11.12.20.
//

#ifndef NBODY_MATH_UTILS_H
#define NBODY_MATH_UTILS_H

#include "body.h"
#include <vector>

using namespace std;

void calc_direct_force(vector<Body> &particles);

void leapfrog(vector<Body> &particles, double dt);

double get_dt(vector<Body> &particles);


#endif //NBODY_MATH_UTILS_H
