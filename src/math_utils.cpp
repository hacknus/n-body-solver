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
        bodies[self].epot = 0;

        // here we cheat: only the interactions with the first 9 bodies are calculated (planets)
        // comet to comet interactions are neglected.
        // for (int partner = 0; partner < bodies.size(); partner++) {
        for (vector<Body>::size_type partner = 0; partner < (bodies.size() - ignore_bodies); partner++) {
            if (self != partner) {
                x = bodies[self].x - bodies[partner].x;
                y = bodies[self].y - bodies[partner].y;
                z = bodies[self].z - bodies[partner].z;
                bodies[self].ax -= G * bodies[partner].m * x / pow(pow(x, 2) + pow(y, 2) + pow(z, 2)
                                                                   + pow(softening, 2), 1.5);
                bodies[self].ay -= G * bodies[partner].m * y / pow(pow(x, 2) + pow(y, 2) + pow(z, 2)
                                                                   + pow(softening, 2), 1.5);
                bodies[self].az -= G * bodies[partner].m * z / pow(pow(x, 2) + pow(y, 2) + pow(z, 2)
                                                                   + pow(softening, 2), 1.5);
                bodies[self].epot += G * bodies[partner].m * bodies[self].m /
                                     pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 0.5);
            }

        }
    }
}

void leapfrog(vector<Body> &bodies, double dt, int num_procs, int myid, MPI_Datatype MPI_BODY_TYPE,
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
        bodies[i].ekin = 0.5 * bodies[i].m * (pow(bodies[i].vx, 2) + pow(bodies[i].vy, 2) + pow(bodies[i].vz, 2));
        bodies[i].x = bodies[i].x + bodies[i].vx * 0.5 * dt;
        bodies[i].y = bodies[i].y + bodies[i].vy * 0.5 * dt;
        bodies[i].z = bodies[i].z + bodies[i].vz * 0.5 * dt;
    }

    // gather all partial body slices from sub-processes on root-process
    const int tag2 = 12;
    if (myid == 0) {
        MPI_Status status;
        vector<Body>::size_type ai;
        vector<Body>::size_type bi;
        for (int proc = 1; proc < num_procs; proc++) {
            vector<Body> recv;
            ai = floor((float) bodies.size() / (float) num_procs * proc);
            bi = floor((float) bodies.size() / (float) num_procs * (proc + 1));
            recv.resize(bi - ai);
            MPI_Recv(&recv.front(), recv.size(), MPI_BODY_TYPE, proc, tag2, MPI_COMM_WORLD, &status);
            copy(begin(recv), end(recv), begin(bodies) + ai);
        }
    } else {
        vector<Body> send;
        send = vector<Body>(bodies.begin() + a, bodies.begin() + b);
        MPI_Send(&send.front(), send.size(), MPI_BODY_TYPE, 0, tag2, MPI_COMM_WORLD);
    }

}

void rk_accel(vector<Body> &bodies, vector<Body>::size_type a, vector<Body>::size_type b,
              unsigned long int ignore_bodies, float G, float softening, vector<double> &ax, vector<double> &ay,
              vector<double> &az) {
    double x, y, z;

    for (vector<Body>::size_type self = a; self < b; self++) {

        bodies[self].ax = 0;
        bodies[self].ay = 0;
        bodies[self].az = 0;
        bodies[self].epot = 0;

        // here we cheat: only the interactions with the first 9 bodies are calculated (planets)
        // comet to comet interactions are neglected.
        // for (int partner = 0; partner < bodies.size(); partner++) {
        for (vector<Body>::size_type partner = 0; partner < (bodies.size() - ignore_bodies); partner++) {
            if (self != partner) {
                x = bodies[self].x - bodies[partner].x;
                y = bodies[self].y - bodies[partner].y;
                z = bodies[self].z - bodies[partner].z;
                ax[self] -= G * bodies[partner].m * x / pow(pow(x, 2) + pow(y, 2) + pow(z, 2)
                                                            + pow(softening, 2), 1.5);
                ay[self] -= G * bodies[partner].m * y / pow(pow(x, 2) + pow(y, 2) + pow(z, 2)
                                                            + pow(softening, 2), 1.5);
                az[self] -= G * bodies[partner].m * z / pow(pow(x, 2) + pow(y, 2) + pow(z, 2)
                                                            + pow(softening, 2), 1.5);
                bodies[self].epot += G * bodies[partner].m * bodies[self].m /
                                     pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 0.5);
            }

        }
    }
}


void rk4(vector<Body> &bodies, double dt, int num_procs, int myid, MPI_Datatype MPI_BODY_TYPE,
         unsigned long int ignore_bodies, float G, float softening) {
    vector<Body>::size_type a = floor((float) bodies.size() / (float) num_procs * myid);
    vector<Body>::size_type b = floor((float) bodies.size() / (float) num_procs * (myid + 1));

    sync(bodies, num_procs, myid, MPI_BODY_TYPE);

    vector<Body> bodies_temp;
    bodies_temp = bodies;

    vector<double> ax1(bodies.size(), 0.0);
    vector<double> ay1(bodies.size(), 0.0);
    vector<double> az1(bodies.size(), 0.0);

    vector<double> vx1(bodies.size(), 0.0);
    vector<double> vy1(bodies.size(), 0.0);
    vector<double> vz1(bodies.size(), 0.0);

    vector<double> ax2(bodies.size(), 0.0);
    vector<double> ay2(bodies.size(), 0.0);
    vector<double> az2(bodies.size(), 0.0);

    vector<double> vx2(bodies.size(), 0.0);
    vector<double> vy2(bodies.size(), 0.0);
    vector<double> vz2(bodies.size(), 0.0);

    vector<double> ax3(bodies.size(), 0.0);
    vector<double> ay3(bodies.size(), 0.0);
    vector<double> az3(bodies.size(), 0.0);

    vector<double> vx3(bodies.size(), 0.0);
    vector<double> vy3(bodies.size(), 0.0);
    vector<double> vz3(bodies.size(), 0.0);

    vector<double> ax4(bodies.size(), 0.0);
    vector<double> ay4(bodies.size(), 0.0);
    vector<double> az4(bodies.size(), 0.0);

    vector<double> vx4(bodies.size(), 0.0);
    vector<double> vy4(bodies.size(), 0.0);
    vector<double> vz4(bodies.size(), 0.0);

    rk_accel(bodies_temp, a, b, ignore_bodies, G, softening, ax1, ay1, az1);

    for (vector<Body>::size_type i = a; i < b; i++) {
        vx1[i] = bodies[i].vx;
        vy1[i] = bodies[i].vy;
        vz1[i] = bodies[i].vz;
        bodies_temp[i].x = bodies[i].x + vx1[i] * 0.5 * dt;
        bodies_temp[i].y = bodies[i].y + vy1[i] * 0.5 * dt;
        bodies_temp[i].z = bodies[i].z + vz1[i] * 0.5 * dt;
    }

    sync(bodies_temp, num_procs, myid, MPI_BODY_TYPE);
    rk_accel(bodies_temp, a, b, ignore_bodies, G, softening, ax2, ay2, az2);

    for (vector<Body>::size_type i = a; i < b; i++) {
        vx2[i] = bodies[i].vx + ax1[i] * dt;
        vy2[i] = bodies[i].vy + ay1[i] * dt;
        vz2[i] = bodies[i].vz + az1[i] * dt;
        bodies_temp[i].x = bodies[i].x + vx2[i] * 0.5 * dt;
        bodies_temp[i].y = bodies[i].y + vy2[i] * 0.5 * dt;
        bodies_temp[i].z = bodies[i].z + vz2[i] * 0.5 * dt;
    }

    sync(bodies_temp, num_procs, myid, MPI_BODY_TYPE);
    rk_accel(bodies_temp, a, b, ignore_bodies, G, softening, ax3, ay3, az3);

    for (vector<Body>::size_type i = a; i < b; i++) {
        vx3[i] = bodies[i].vx + ax2[i] * dt;
        vy3[i] = bodies[i].vy + ay2[i] * dt;
        vz3[i] = bodies[i].vz + az2[i] * dt;
        bodies_temp[i].x = bodies[i].x + vx3[i] * 0.5 * dt;
        bodies_temp[i].y = bodies[i].y + vy3[i] * 0.5 * dt;
        bodies_temp[i].z = bodies[i].z + vz3[i] * 0.5 * dt;
    }

    rk_accel(bodies_temp, a, b, ignore_bodies, G, softening, ax4, ay4, az4);


    sync(bodies, num_procs, myid, MPI_BODY_TYPE);

    for (vector<Body>::size_type i = a; i < b; i++) {
        vx4[i] = bodies[i].vx + ax4[i] * dt;
        vy4[i] = bodies[i].vy + ay4[i] * dt;
        vz4[i] = bodies[i].vz + az4[i] * dt;
        bodies[i].x = bodies[i].x + 1 / 6.0 * dt * (vx1[i] + 2 * vx2[i] + 2 * vx3[i] + vx4[i]);
        bodies[i].y = bodies[i].y + 1 / 6.0 * dt * (vy1[i] + 2 * vy2[i] + 2 * vy3[i] + vy4[i]);
        bodies[i].z = bodies[i].z + 1 / 6.0 * dt * (vz1[i] + 2 * vz2[i] + 2 * vz3[i] + vz4[i]);
        bodies[i].vx = bodies[i].vx + 1 / 6.0 * dt * (ax1[i] + 2 * ax2[i] + 2 * ax3[i] + ax4[i]);
        bodies[i].vy = bodies[i].vy + 1 / 6.0 * dt * (ay1[i] + 2 * ay2[i] + 2 * ay3[i] + ay4[i]);
        bodies[i].vz = bodies[i].vz + 1 / 6.0 * dt * (az1[i] + 2 * az2[i] + 2 * az3[i] + az4[i]);
        bodies[i].ekin = 0.5 * bodies[i].m * (pow(bodies[i].vx, 2) + pow(bodies[i].vy, 2) + pow(bodies[i].vz, 2));
    }

    // gather all partial body slices from sub-processes on root-process
    const int tag2 = 12;
    if (myid == 0) {
        MPI_Status status;
        vector<Body>::size_type ai;
        vector<Body>::size_type bi;
        for (int proc = 1; proc < num_procs; proc++) {
            vector<Body> recv;
            ai = floor((float) bodies.size() / (float) num_procs * proc);
            bi = floor((float) bodies.size() / (float) num_procs * (proc + 1));
            recv.resize(bi - ai);
            MPI_Recv(&recv.front(), recv.size(), MPI_BODY_TYPE, proc, tag2, MPI_COMM_WORLD, &status);
            copy(begin(recv), end(recv), begin(bodies) + ai);
        }
    } else {
        vector<Body> send;
        send = vector<Body>(bodies.begin() + a, bodies.begin() + b);
        MPI_Send(&send.front(), send.size(), MPI_BODY_TYPE, 0, tag2, MPI_COMM_WORLD);
    }

}

void sync(vector<Body> &bodies, int num_procs, int myid, MPI_Datatype MPI_BODY_TYPE) {
    vector<Body>::size_type a = floor((float) bodies.size() / (float) num_procs * myid);
    vector<Body>::size_type b = floor((float) bodies.size() / (float) num_procs * (myid + 1));
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

}


double get_dt(vector<Body> &bodies, vector<Body>::size_type a, vector<Body>::size_type b, float softening) {

    double dt[b - a];
    double min_dt;
    double min_dt_out = 0.001;
    double a_mag;
    for (vector<Body>::size_type i = a; i < b; i++) {
        a_mag = sqrt(bodies[i].ax * bodies[i].ax + bodies[i].ay * bodies[i].ay + bodies[i].az * bodies[i].az);
        dt[i - a] = 1 * sqrt(softening / a_mag);
    }
    int n_dt = sizeof(dt) / sizeof(dt[0]);
    min_dt = *min_element(dt, dt + n_dt);
    // gather minimum of min_dt and distribute to all processes in min_dt_out
    MPI_Allreduce(&min_dt, &min_dt_out, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return min_dt_out;
}