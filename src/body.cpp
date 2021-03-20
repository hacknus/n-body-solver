//
// Created by Linus on 11.12.20.
//

#include "body.h"

void init_body(Body *b, double _m, double _x, double _y, double _z, double _vx, double _vy, double _vz) {
    b->m = _m;
    b->x = _x;
    b->y = _y;
    b->z = _z;
    b->vx = _vx;
    b->vy = _vy;
    b->vz = _vz;
    b->ax = 0.0;
    b->ay = 0.0;
    b->az = 0.0;
};

MPI_Datatype make_mpi_type() {
    MPI_Datatype mpi_body_type;
    const int nitems = 10;
    int blocklengths[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[10] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                              MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                              MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[10];

    offsets[0] = offsetof(Body, m);
    offsets[1] = offsetof(Body, x);
    offsets[2] = offsetof(Body, y);
    offsets[3] = offsetof(Body, z);
    offsets[4] = offsetof(Body, vx);
    offsets[5] = offsetof(Body, vy);
    offsets[6] = offsetof(Body, vz);
    offsets[7] = offsetof(Body, ax);
    offsets[8] = offsetof(Body, ay);
    offsets[9] = offsetof(Body, az);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_body_type);
    MPI_Type_commit(&mpi_body_type);
    return mpi_body_type;
}
