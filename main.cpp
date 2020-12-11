#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <math.h>
#include "body.h"
#include "math_utils.h"
#include "io.h"

using namespace std;


int main() {
    cout << "\n";
    cout << "GRAVITATIONAL N-BODY SOLVER \n";
    cout << " \n";


    cout << "reading data... \n";

    ifstream file("cdata.csv");

    vector<Body>::size_type index = 0;
    vector<Body> bodies;

    bodies = read_initial();

    int n = bodies.size();
    cout << "[OK] found " << n << " particles\n";

    calc_direct_force(bodies);

    cout << "ax: " << bodies[0].ax << "\n";
    cout << "ay: " << bodies[0].ay << "\n";
    cout << "az: " << bodies[0].az << "\n";

    cout << "ax: " << bodies[1].ax << "\n";
    cout << "ay: " << bodies[1].ay << "\n";
    cout << "az: " << bodies[1].az << "\n";

    double dt = get_dt(bodies);
    double t = 0;

    char filename[32]; // make sure it's big enough

    for (int step = 0; step < 5000; step++) {

        dt = get_dt(bodies);
        t += dt;
        leapfrog(bodies, dt);

        snprintf(filename, sizeof(filename), "../output/out_%05d.bin", step);
        write_file(bodies, filename, dt, t);

    }

    return 0;


}