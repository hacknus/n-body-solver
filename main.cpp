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
        cout << "\n";
        cout << "GRAVITATIONAL N-BODY SOLVER \n";
        cout << "\n";
    }

    vector<Body>::size_type index = 0;
    vector<Body> bodies;


    MPI_Datatype mpi_body_type = make_mpi_type();


    uint64_t num_steps = 0;

    if (myid == 0) {

        if (argc < 3){
            cout << "[ERROR] no input file/step number specified! \n";
            cout << "        correct usage: mpirun -np $num_cores nBody $path_to_input_file $num_steps \n";
            exit(EXIT_FAILURE);
        }

        cout << "[OK] root is reading input data from " << argv[2] << "\n";
        char *pCh;
        num_steps = strtoul (argv[1], &pCh, 10);

        bodies = read_initial(argv[2]);
    }

    uint32_t size = bodies.size();

    MPI_Bcast(&size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    if (myid != 0) bodies.resize(size);
    MPI_Bcast(&bodies.front(), bodies.size(), mpi_body_type, 0, MPI_COMM_WORLD);

    MPI_Bcast(&num_steps, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

    if (myid == 0){
        cout << "[OK] found " << size << " bodies found" << "\n";
        cout << "[OK] starting simulation for " << num_steps << " steps \n";
        cout << "\n-------------------------------------------------------------------\n\n";
    }

    calc_direct_force(bodies, 0, bodies.size());

    int a = bodies.size() / num_procs * myid;
    int b = bodies.size() / num_procs * (myid + 1);
    double dt = get_dt(bodies, a, b);
    double t = 0;

    char filename[32]; // make sure it's big enough

    for (int step = 0; step < num_steps; step++) {
        dt = get_dt(bodies, a, b);
        dt = 24*60*60;
        t += dt;
        leapfrog(bodies, dt, num_procs, myid, mpi_body_type);

        if ((myid == 0) && (step % 10 == 0)){
            if (myid == 0) cout << "dt: " << dt << "\n";
            snprintf(filename, sizeof(filename), "../output/out_%07d.dat", step);
            write_file(bodies, filename, dt, t);
        }
    }

    MPI_Type_free(&mpi_body_type);
    MPI_Finalize();

    return 0;

}