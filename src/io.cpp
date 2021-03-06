//
// Created by Linus on 11.12.20.
//

#include "io.h"
#include <iterator>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

void get_initial_values(string *path, unsigned long long int *steps, double *dt, unsigned long int *save_interval,
                        unsigned long int *ignore_bodies, float *G, float *softening) {
    string input_path = "../input/input.conf";

    // std::ifstream is RAII, i.e. no need to call close
    ifstream cFile(input_path);
    if (cFile.is_open()) {
        string line;
        while (getline(cFile, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }
            auto delimiterPos = line.find("=");
            auto name = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);
            char *pCh;
            if (name == "path") *path = value;
            if (name == "steps") *steps = strtoul(value.c_str(), &pCh, 10);
            if (name == "dt") *dt = strtod(value.c_str(), &pCh);
            if (name == "save_interval") *save_interval = stoul(value.c_str());
            if (name == "ignore_bodies") *ignore_bodies = stoul(value.c_str());
            if (name == "G") *G = stof(value.c_str());
            if (name == "softening") *softening = stof(value.c_str());
        }
    } else {
        cout << "[ERROR] no input file in : '" << input_path << "'!\n";
        exit(EXIT_FAILURE);
    }
}


void write_file(vector<Body> bodies, char filename[], double dt, double t) {
    fstream outFile;
    outFile.open(filename, ios::out | ios::binary);
    for (vector<Body>::size_type index = 0; index < bodies.size(); index++) {
        // TODO: save name as well
        outFile.write((char *) &bodies[index].m, sizeof(double));
        outFile.write((char *) &bodies[index].x, sizeof(double));
        outFile.write((char *) &bodies[index].y, sizeof(double));
        outFile.write((char *) &bodies[index].z, sizeof(double));
        outFile.write((char *) &bodies[index].vx, sizeof(double));
        outFile.write((char *) &bodies[index].vy, sizeof(double));
        outFile.write((char *) &bodies[index].vz, sizeof(double));
        outFile.write((char *) &bodies[index].ax, sizeof(double));
        outFile.write((char *) &bodies[index].ay, sizeof(double));
        outFile.write((char *) &bodies[index].az, sizeof(double));
        outFile.write((char *) &dt, sizeof(double));
        outFile.write((char *) &t, sizeof(double));
        outFile.write((char *) &bodies[index].epot, sizeof(double));
        outFile.write((char *) &bodies[index].ekin, sizeof(double));
    }
    outFile.close();
}

vector<Body> read_initial(string path, float G) {
    ifstream csvFile;
    csvFile.open(path);

    if (!csvFile.is_open()) {
        cout << "[ERROR] no such file: '" << path << "'!\n";
        exit(EXIT_FAILURE);
    }

    vector<Body> bodies;
    string line;
    vector<string> vec;
    getline(csvFile, line); // skip the 1st line
    while (getline(csvFile, line)) {

        if (line.empty()) // skip empty lines
        {
            continue;
        }

        istringstream iss(line);
        string lineStream;
        string::size_type sz;

        vector<double> row;

        uint32_t counter = 0;
        while (getline(iss, lineStream, ',')) {
            counter++;
            if (counter == 1) continue; //skip first entry (name of object, not a double)
            // TODO: save name in bodies vector and eventually in binary .dat output.
            row.push_back(stod(lineStream, &sz)); // convert to double
        }
        float au, m_sol, day;
        if (G != 1) {
            au = 1.5e11;
            m_sol = 2e30;
            day = 24.0 * 60.0 * 60.0;
        } else {
            au = 1;
            m_sol = 1;
            day = 1;
        }

        Body b = Body();
        init_body(&b, row[0] * m_sol,
                  row[1] * au, row[2] * au, row[3] * au,
                  row[4] * au / day, row[5] * au / day, row[6] * au / day);
        bodies.push_back(b);
    }

    return bodies;
}