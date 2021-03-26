#include <iostream>
#include <vector>
#include <string>
#include "body.h"
#include "math_utils.h"
#include "io.h"
#include <mpi.h>

using namespace std;


int main(int argc, char **argv) {

    int num_procs, myid;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (myid == 0) {
        cout << "|-------------------------------------------|\n";
        cout << "|--------GRAVITATIONAL N-BODY SOLVER--------|\n";
        cout << "|-------------------------------------------|\n";
    }

    vector<Body> bodies;

    // create MPI_Datatype to broadcast vector of custom structs
    MPI_Datatype mpi_body_type = make_mpi_type();

    vector<Body>::size_type a, b;
    double t = 0;
    char filename[128]; // make sure it's big enough
    unsigned long long int num_steps = 0;
    double dt = 0;
    unsigned long int save_interval = 0;
    unsigned long int ignore_bodies = 0;
    float G = 1;

    if (myid == 0) {
        // only root process reads the input file
        string path;
        get_initial_values(&path, &num_steps, &dt, &save_interval, &ignore_bodies, &G);
        cout << "[OK] path for initial conditions is: " << path << "\n";
        cout << "[OK] simulation steps: " << num_steps << "\n";
        cout << "[OK] dt (internal calculation if 0): " << dt << "\n";
        cout << "[OK] save interval is: " << save_interval << "\n";
        cout << "[OK] ignoring bodies: " << ignore_bodies << "\n";
        cout << "[OK] G is: " << G << "\n";
        // read initial file
        bodies = read_initial(path, G);
    }

    unsigned long int size = bodies.size();

    // broadcast size of bodies vector of root-process to all sub-processes
    MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (myid != 0) bodies.resize(size); // sub-processes resize their vector
    // broadcast bodies vector with initial conditions from root-process to sub-processes
    MPI_Bcast(&bodies.front(), size, mpi_body_type, 0, MPI_COMM_WORLD);

    // broadcast configuration parameters
    MPI_Bcast(&num_steps, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&save_interval, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ignore_bodies, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&G, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    a = bodies.size() / num_procs * myid;
    b = bodies.size() / num_procs * (myid + 1);


    if (myid == 0) {
        cout << "[OK] found " << size << " bodies" << "\n";
        cout << "[OK] starting simulation: \n"
                "     for " << num_steps << " steps \n "
                "     on " << num_procs << " cores. \n";
        cout << "\n|-------------------------------------------|\n\n";
    }
    // calculate forces (accelerations) once in order to determine initial time-step
    calc_direct_force(bodies, 0, bodies.size(), ignore_bodies, G);

    // begin simulation
    for (int step = 0; step < num_steps; step++) {
        if (dt == 0) dt = get_dt(bodies, a, b);
        t += dt;
        leapfrog(bodies, dt, num_procs, myid, mpi_body_type, ignore_bodies, G);

        if ((myid == 0) && (step % save_interval == 0)) {
            // only root process saves all the data
            // important:   the slice of the bodies that are distributed to the root process are already dt/2
            //              propagated further than all the others. This deviation is minimal and
            //              not visible when plotting the orbits.
            cout << "[OK] step " << step << "/" << num_steps << " completed.\n";
            snprintf(filename, sizeof(filename), "../output/out_%07d.dat", step);
            write_file(bodies, filename, dt, t);
        }
    }

    if (myid == 0) cout << "[OK] simulation completed after " << num_steps << " steps \n";
    if (myid == 0) cout << "\n|-------------------------------------------|\n\n";

    MPI_Type_free(&mpi_body_type);
    MPI_Finalize();

    return 0;

}