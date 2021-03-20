//
// Created by Linus on 11.12.20.
//

#ifndef NBODY_BODY_H
#define NBODY_BODY_H

#include <mpi.h>

struct Body {
    double m;
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
};

void init_body(Body *b, double _m, double _x, double _y, double _z, double _vx, double _vy, double _vz);

MPI_Datatype make_mpi_type();

#endif //NBODY_BODY_H
