//
// Created by Linus on 11.12.20.
//

#include "math_utils.h"
#include "body.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <mpi.h>
#include <algorithm> // Necessary for `std::copy`...

using namespace std;

void calc_direct_force(vector<Body> &bodies, int a, int b) {
    double G = 6.67408e-11;
    double softening = 0.0001;
    double x, y, z;

    for (int self = a; self < b; self++) {

        bodies[self].ax = 0;
        bodies[self].ay = 0;
        bodies[self].az = 0;

        // here we cheat: only the interactions with the first 9 bodies are calculated (planets)
        // comet to comet interactions are neglected.
        // for (int partner = 0; partner < bodies.size(); partner++) {
        for (int partner = 0; partner < 9; partner++) {
                if (self != partner) {
                    x = bodies[self].x - bodies[partner].x;
                    y = bodies[self].y - bodies[partner].y;
                    z = bodies[self].z - bodies[partner].z;
                    bodies[self].ax -=
                            G * bodies[partner].m * x / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) + pow(softening, 2), 1.5);
                    bodies[self].ay -=
                            G * bodies[partner].m * y / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) + pow(softening, 2), 1.5);
                    bodies[self].az -=
                            G * bodies[partner].m * z / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) + pow(softening, 2), 1.5);
                    //                particles[self].epot +=
                    //                        G * particles[partner].m * particles[self].m /
                    //                        pow(pow(x, 2) + pow(y, 2) + pow(z, 2) + pow(softening, 2), 0.5);
                }

        }
    }
}

void leapfrog(vector<Body> &bodies, double dt, int num_procs, int myid, MPI_Datatype mpi_body_type) {
    int a = bodies.size() / num_procs * myid;
    int b = bodies.size() / num_procs * (myid + 1);

    for (int i = a; i < b; i++) {
        bodies[i].x = bodies[i].x + bodies[i].vx * 0.5 * dt;
        bodies[i].y = bodies[i].y + bodies[i].vy * 0.5 * dt;
        bodies[i].z = bodies[i].z + bodies[i].vz * 0.5 * dt;
    }

    const int tag = 13;
    if (myid == 0) {
        MPI_Status status;
        for (int proc = 1; proc < num_procs; proc++) {
            vector<Body> recv;
            recv.resize(b - a);
            MPI_Recv(&recv.front(), recv.size(), mpi_body_type, proc, tag, MPI_COMM_WORLD, &status);
            copy(recv.begin(), recv.end(), bodies.begin() + bodies.size() / num_procs * proc);
        }
    } else {
        vector<Body> send;
        send = vector<Body>(bodies.begin() + a, bodies.begin() + b);
        MPI_Send(&send.front(), send.size(), mpi_body_type, 0, tag, MPI_COMM_WORLD);
    }

    MPI_Bcast(&bodies.front(), bodies.size(), mpi_body_type, 0, MPI_COMM_WORLD);

    calc_direct_force(bodies, a, b);

    for (int i = a; i < b; i++) {
        bodies[i].vx = bodies[i].vx + bodies[i].ax * dt;
        bodies[i].vy = bodies[i].vy + bodies[i].ay * dt;
        bodies[i].vz = bodies[i].vz + bodies[i].az * dt;
        // self->ekin = 0.5 * self->m * (pow(self->vx, 2) + pow(self->vy, 2) + pow(self->vz, 2));
        bodies[i].x = bodies[i].x + bodies[i].vx * 0.5 * dt;
        bodies[i].y = bodies[i].y + bodies[i].vy * 0.5 * dt;
        bodies[i].z = bodies[i].z + bodies[i].vz * 0.5 * dt;
    }

}

double get_dt(vector<Body> &bodies, int a, int b) {

    double dt[b-a];
    int index = 0;
    double softening = 0.01;
    double min_dt = 0.001;
    double min_dt_out = 0.001;
    double a_mag = 0;
    for (int i = a; i < b; i++) {
        a_mag = pow(pow(bodies[i].ax, 2) + pow(bodies[i].ay, 2) + pow(bodies[i].az, 2), 0.5);
        dt[i-a] = 0.1 * sqrt(softening / a_mag);
        index++;
    }
    int n_dt = sizeof(dt) / sizeof(dt[0]);
    min_dt = *min_element(dt, dt + n_dt);
    MPI_Allreduce(&min_dt, &min_dt_out, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return min_dt_out;
}