//
// Created by Linus on 11.12.20.
//

#include "math_utils.h"
#include "body.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <bits/stdc++.h>

using namespace std;

void calc_direct_force(vector<Body> &particles) {
    double G = 1;
    double softening = 0.01;
    double x, y, z;
    int counter = 0;

    for (int self = 0; self < particles.size(); self++) {
        //cout << "calculating " << counter << "... \n";
        counter++;
        particles[self].ax = 0;
        particles[self].ay = 0;
        particles[self].az = 0;

        for (int partner = 0; partner < particles.size(); partner++) {
            if (self != partner) {
                x = particles[self].x - particles[partner].x;
                y = particles[self].y - particles[partner].y;
                z = particles[self].z - particles[partner].z;
                particles[self].ax -=
                        G * particles[partner].m * x / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) + pow(softening, 2), 1.5);
                particles[self].ay -=
                        G * particles[partner].m * y / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) + pow(softening, 2), 1.5);
                particles[self].az -=
                        G * particles[partner].m * z / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) + pow(softening, 2), 1.5);
//                particles[self].epot +=
//                        G * particles[partner].m * particles[self].m /
//                        pow(pow(x, 2) + pow(y, 2) + pow(z, 2) + pow(softening, 2), 0.5);
            }
        }
    }
}

void leapfrog(vector<Body> &particles, double dt) {

    for (auto self = particles.begin(); self != particles.end(); ++self) {
        self->x = self->x + self->vx * 0.5 * dt;
        self->y = self->y + self->vy * 0.5 * dt;
        self->z = self->z + self->vz * 0.5 * dt;
    }

    calc_direct_force(particles);

    for (auto self = particles.begin(); self != particles.end(); ++self) {
        self->vx = self->vx + self->ax * dt;
        self->vy = self->vy + self->ay * dt;
        self->vz = self->vz + self->az * dt;
        // self->ekin = 0.5 * self->m * (pow(self->vx, 2) + pow(self->vy, 2) + pow(self->vz, 2));
        self->x = self->x + self->vx * 0.5 * dt;
        self->y = self->y + self->vy * 0.5 * dt;
        self->z = self->z + self->vz * 0.5 * dt;
    }

}

double get_dt(vector<Body> &particles) {

    unsigned int n_p = particles.size();
    double dt[n_p];
    int index = 0;
    double softening = 0.01;
    double min_dt = 0.001;
    double a_mag = 0;
    for (Body &p : particles) {
        a_mag = pow(pow(p.ax, 2) + pow(p.ay, 2) + pow(p.az, 2), 0.5);
        dt[index] = 0.1 * sqrt(softening / a_mag);
        index++;
    }
    int n_dt = sizeof(dt) / sizeof(dt[0]);
    min_dt = *min_element(dt, dt + n_dt);
    cout << "dt: " << min_dt << "\n";
    return min_dt;
}