//
// Created by Linus on 11.12.20.
//

#include "math_utils.h"
#include <math.h>
#include <vector>
#include <mpi.h>
#include <algorithm> // Necessary for `std::copy`...

using namespace std;

void calc_direct_force(vector<Body> &bodies, vector<Body>::size_type a, vector<Body>::size_type b,
                       unsigned long int ignore_bodies, float G, float softening) {
    double x, y, z;

    for (vector<Body>::size_type self = a; self < b; self++) {

        bodies[self].ax = 0;
        bodies[self].ay = 0;
        bodies[self].az = 0;

        // here we cheat: only the interactions with the first 9 bodies are calculated (planets)
        // comet to comet interactions are neglected.
        // for (int partner = 0; partner < bodies.size(); partner++) {
        for (vector<Body>::size_type partner = 0; partner < (bodies.size() - ignore_bodies); partner++) {
            if (self != partner) {
                x = bodies[self].x - bodies[partner].x;
                y = bodies[self].y - bodies[partner].y;
                z = bodies[self].z - bodies[partner].z;
                bodies[self].ax -=
                        G * bodies[partner].m * x / pow(pow(x, 2) + pow(y, 2) + pow(z, 2)
                                                        + pow(softening, 2), 1.5);
                bodies[self].ay -=
                        G * bodies[partner].m * y / pow(pow(x, 2) + pow(y, 2) + pow(z, 2)
                                                        + pow(softening, 2), 1.5);
                bodies[self].az -=
                        G * bodies[partner].m * z / pow(pow(x, 2) + pow(y, 2) + pow(z, 2)
                                                        + pow(softening, 2), 1.5);
                //                particles[self].epot +=
                //                        G * particles[partner].m * particles[self].m /
                //                        pow(pow(x, 2) + pow(y, 2) + pow(z, 2) + pow(softening, 2), 0.5);
            }

        }
    }
}

void
leapfrog(vector<Body> &bodies, double dt, int num_procs, int myid, MPI_Datatype MPI_BODY_TYPE,
         unsigned long int ignore_bodies, float G, float softening) {
    vector<Body>::size_type a = floor((float) bodies.size() / (float) num_procs * myid);
    vector<Body>::size_type b = floor((float) bodies.size() / (float) num_procs * (myid + 1));

    for (vector<Body>::size_type i = a; i < b; i++) {
        bodies[i].x = bodies[i].x + bodies[i].vx * 0.5 * dt;
        bodies[i].y = bodies[i].y + bodies[i].vy * 0.5 * dt;
        bodies[i].z = bodies[i].z + bodies[i].vz * 0.5 * dt;
    }

    // gather all partial body slices from sub-processes on root-process
    const int tag = 13;
    if (myid == 0) {
        MPI_Status status;
        vector<Body>::size_type ai;
        vector<Body>::size_type bi;
        for (int proc = 1; proc < num_procs; proc++) {
            vector<Body> recv;
            ai = floor((float) bodies.size() / (float) num_procs * proc);
            bi = floor((float) bodies.size() / (float) num_procs * (proc + 1));
            recv.resize(bi - ai);
            MPI_Recv(&recv.front(), recv.size(), MPI_BODY_TYPE, proc, tag, MPI_COMM_WORLD, &status);
            copy(begin(recv), end(recv), begin(bodies) + ai);
        }
    } else {
        vector<Body> send;
        send = vector<Body>(bodies.begin() + a, bodies.begin() + b);
        MPI_Send(&send.front(), send.size(), MPI_BODY_TYPE, 0, tag, MPI_COMM_WORLD);
    }

    // distribute combined body slices to all sub-processes from root-processes
    MPI_Bcast(&bodies.front(), bodies.size(), MPI_BODY_TYPE, 0, MPI_COMM_WORLD);

    calc_direct_force(bodies, a, b, ignore_bodies, G, softening);

    for (vector<Body>::size_type i = a; i < b; i++) {
        bodies[i].vx = bodies[i].vx + bodies[i].ax * dt;
        bodies[i].vy = bodies[i].vy + bodies[i].ay * dt;
        bodies[i].vz = bodies[i].vz + bodies[i].az * dt;
        // self->ekin = 0.5 * self->m * (pow(self->vx, 2) + pow(self->vy, 2) + pow(self->vz, 2));
        bodies[i].x = bodies[i].x + bodies[i].vx * 0.5 * dt;
        bodies[i].y = bodies[i].y + bodies[i].vy * 0.5 * dt;
        bodies[i].z = bodies[i].z + bodies[i].vz * 0.5 * dt;
    }

}

double get_dt(vector<Body> &bodies, vector<Body>::size_type a, vector<Body>::size_type b, float softening) {

    double dt[b - a];
    double min_dt;
    double min_dt_out = 0.001;
    double a_mag;
    for (vector<Body>::size_type i = a; i < b; i++) {
        a_mag = pow(bodies[i].ax * bodies[i].ax + bodies[i].ay * bodies[i].ay + bodies[i].az * bodies[i].az, 0.5);
        dt[i - a] = sqrt(softening / a_mag);
    }
    int n_dt = sizeof(dt) / sizeof(dt[0]);
    min_dt = *min_element(dt, dt + n_dt);
    // gather minimum of min_dt and distribute to all processes in min_dt_out
    MPI_Allreduce(&min_dt, &min_dt_out, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return min_dt_out;
}