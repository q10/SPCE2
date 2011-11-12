#include "common.h"

ostream & operator<<(ostream & out, Simulation * simulation) {
    // The policy is that the config file output only contains the configuration of the box and its particles;
    // all other system parameters will be set in code and compiled before running;
    // otherwise, the complexity of managing the config file will grow exponentially
    out << "BOX_LENGTH\t" << setprecision(10) << simulation->BOX_LENGTH << endl
            << "BOX_Z_LENGTH\t" << simulation->BOX_Z_LENGTH << endl
            << "EWALD_ALPHA\t" << simulation->EWALD_ALPHA << endl
            << "EWALD_NXY\t" << simulation->EWALD_NXY << endl
            << "EWALD_NZ\t" << simulation->EWALD_NZ << endl;

    for (unsigned int i = 0; i < simulation->WATERS.size(); i++)
        out << "WATER\t" << simulation->WATERS[i] << endl;

    for (unsigned int i = 0; i < simulation->IONS.size(); i++)
        out << "ION\t" << simulation->IONS[i] << endl;

    return out;
}

string Simulation::to_lammpstrj(int time_step) {
    stringstream lammpstrj_string;
    int atom_count = 0, ion_id;
    double *coords;
    lammpstrj_string << setprecision(10) << "ITEM: TIMESTEP" << endl << time_step << endl
            << "ITEM: NUMBER OF ATOMS" << endl << (3 * WATERS.size() + IONS.size()) << endl
            << "ITEM: BOX BOUNDS" << endl
            << "0 " << BOX_LENGTH << endl
            << "0 " << BOX_LENGTH << endl
            << "0 " << BOX_Z_LENGTH << endl
            << "ITEM: ATOMS id type x y z diameter q" << endl;
    for (unsigned int i = 0; i < WATERS.size(); i++) {
        coords = WATERS[i]->coords;
        lammpstrj_string << ++atom_count << " 1 " << coords[0] << " " << coords[1] << " " << coords[2] << " " << Water::SIGMA_O << " " << Water::Q_O << endl;
        lammpstrj_string << ++atom_count << " 2 " << coords[3] << " " << coords[4] << " " << coords[5] << " " << Water::SIGMA_H << " " << Water::Q_H << endl;
        lammpstrj_string << ++atom_count << " 2 " << coords[6] << " " << coords[7] << " " << coords[8] << " " << Water::SIGMA_H << " " << Water::Q_H << endl;
    }

    for (unsigned int i = 0; i < IONS.size(); i++) {
        coords = IONS[i]->coords;
        ion_id = (IONS[i]->charge < 0.0) ? 3 : 4;
        lammpstrj_string << ++atom_count << " " << ion_id << " "
                << coords[0] << " " << coords[1] << " " << coords[2] << " " << IONS[i]->SIGMA << " " << IONS[i]->charge << endl;
    }
    return lammpstrj_string.str();
}
