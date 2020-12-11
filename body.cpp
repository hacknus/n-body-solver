//
// Created by Linus on 11.12.20.
//

#include "body.h"

Body::Body() {
};

void Body::init(double _m, double _x, double _y, double _z, double _vx, double _vy, double _vz) {
    m = _m;
    x = _x;
    y = _y;
    z = _z;
    vx = _vx;
    vy = _vy;
    vz = _vz;
};