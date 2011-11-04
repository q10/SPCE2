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

string Simulation::to_vmd(int time_step) {
    stringstream vmd_string;
    int atom_count = 0, ion_id;
    double *coords;
    vmd_string << setprecision(10) << "ITEM: TIMESTEP" << endl << time_step << endl
            << "ITEM: NUMBER OF ATOMS" << endl << (3 * WATERS.size() + IONS.size()) << endl
            << "ITEM: BOX BOUNDS" << endl
            << "0 " << BOX_LENGTH << endl
            << "0 " << BOX_LENGTH << endl
            << "0 " << BOX_Z_LENGTH << endl
            << "ITEM: ATOMS id type x y z diameter q" << endl;
    for (unsigned int i = 0; i < WATERS.size(); i++) {
        coords = WATERS[i]->coords;
        vmd_string << ++atom_count << " 1 " << coords[0] << " " << coords[1] << " " << coords[2] << " " << Water::SIGMA_O << " " << Water::Q_O << endl;
        vmd_string << ++atom_count << " 2 " << coords[3] << " " << coords[4] << " " << coords[5] << " " << Water::SIGMA_H << " " << Water::Q_H << endl;
        vmd_string << ++atom_count << " 2 " << coords[6] << " " << coords[7] << " " << coords[8] << " " << Water::SIGMA_H << " " << Water::Q_H << endl;
    }

    for (unsigned int i = 0; i < IONS.size(); i++) {
        coords = IONS[i]->coords;
        ion_id = (IONS[i]->charge < 0.0) ? 3 : 4;
        vmd_string << ++atom_count << " " << ion_id << " "
                << coords[0] << " " << coords[1] << " " << coords[2] << " " << IONS[i]->SIGMA << " " << IONS[i]->charge << endl;
    }
    return vmd_string.str();
}

void test_vmd_output() {
    cerr << "---- BEGIN TEST - LAMMPSTRJ (VMD) FILE OUTPUT ----" << endl;
    Simulation * s = new Simulation();
    s->DATA_SAMPLING_RATE = 2;
    s->NUM_MC_SWEEPS = 10;
    s->run_mc();
    cerr << "---- END TEST - LAMMPSTRJ (VMD) FILE OUTPUT ----" << endl;
    return;
}

void test_config_output() {
    cerr << "---- BEGIN TEST - CONFIG FILE OUTPUT ----" << endl;
    Simulation * s = new Simulation();
    s->NUM_MC_SWEEPS = 100;
    s->run_mc();
    s->SAMPLER_SET->write_config_snapshot();
    cerr << "---- END TEST - CONFIG FILE OUTPUT ----" << endl;
    return;
}
