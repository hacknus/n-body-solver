//
// Created by Linus on 11.12.20.
//

#ifndef NBODY_IO_H
#define NBODY_IO_H
#include <iterator>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <math.h>
#include "body.h"

using namespace std;

void write_file(vector<Body> bodies, char filename[], double dt, double t);
vector<Body> read_initial(void);


#endif //NBODY_IO_H
