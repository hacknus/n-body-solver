//
// Created by Linus on 11.12.20.
//

#ifndef NBODY_BODY_H
#define NBODY_BODY_H


class Body {
public:
    Body();

    void init(double, double, double, double, double, double, double);

    double m;
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;

    // calculated scalars
    double r;
    double ekin;
    double epot;

    double softening;
    double potential;

};


#endif //NBODY_BODY_H
