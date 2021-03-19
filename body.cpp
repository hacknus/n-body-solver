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
