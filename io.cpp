//
// Created by Linus on 11.12.20.
//

#include "io.h"
#include "body.h"

void write_file(vector<Body> bodies, char filename[], double dt, double t){

    cout << "writing " << filename << "\n";
    int n = bodies.size();
    int index = 0;
    fstream outFile;
    outFile.open(filename, ios::out | ios::binary);
    for (index = 0; index < n; index++) {
        outFile.write((char *) &bodies[index].m, sizeof(double));
        outFile.write((char *) &bodies[index].x, sizeof(double));
        outFile.write((char *) &bodies[index].y, sizeof(double));
        outFile.write((char *) &bodies[index].z, sizeof(double));
        outFile.write((char *) &bodies[index].vx, sizeof(double));
        outFile.write((char *) &bodies[index].vy, sizeof(double));
        outFile.write((char *) &bodies[index].vz, sizeof(double));
        outFile.write((char *) &dt, sizeof(double));
        outFile.write((char *) &t, sizeof(double));
    }
    outFile.close();
}

vector<Body> read_initial(void) {
    ifstream csvFile;
    csvFile.open("../cdata.csv");

    if (!csvFile.is_open()) {
        cout << "no such file." << endl;
        exit(EXIT_FAILURE);
    }

    vector<Body> bodies;
    string line;
    vector<string> vec;
    getline(csvFile, line); // skip the 1st line

    while (getline(csvFile, line)) {
        if (line.empty()) // skip empty lines:
        {
            //cout << "empty line!" << endl;
            continue;
        }

        istringstream iss(line);
        string lineStream;
        string::size_type sz;

        vector<double> row;

        while (getline(iss, lineStream, ',')) {
            row.push_back(stod(lineStream, &sz)); // convert to double
        }
        Body b = Body();
        init_body(&b, row[0], row[1], row[2], row[3], row[4], row[5], row[6]);
        bodies.push_back(b);
    }

    //cout << "size ts = " << timeStampIMU.size() << endl;
    for (size_t i = 0; i < bodies.size(); i++) {
        cout << "mass = " << bodies[i].m << endl;
        cout << "--------------------------------" << endl;
    }

    return bodies;
}