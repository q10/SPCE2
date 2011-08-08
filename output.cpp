#include "common.h"

ostream & operator<<(ostream & out, Simulation * simulation) {
    out << "TEMPERATURE\t" << simulation->TEMPERATURE << endl;

    for (int i = 0; i < simulation->WATERS.size(); i++)
        out << "WATER\t" << simulation->WATERS[i] << endl;

    for (int i = 0; i < simulation->IONS.size(); i++)
        out << "ION\t" << simulation->IONS[i] << endl;

    return out;
}

string Simulation::to_vmd(int time_step) {
    stringstream vmd_string;
    int atom_count = 0, charge;
    double *coords;
    vmd_string << "ITEM: TIMESTEP" << endl << time_step << endl
            << "ITEM: NUMBER OF ATOMS" << endl << (3 * WATERS.size() + IONS.size()) << endl
            << "ITEM: BOX BOUNDS" << endl
            << "0 " << BOX_LENGTH << endl
            << "0 " << BOX_LENGTH << endl
            << "0 " << BOX_Z_LENGTH << endl
            << "ITEM: ATOMS" << endl;
    for (int i = 0; i < WATERS.size(); i++) {
        coords = WATERS[i]->coords;
        vmd_string << ++atom_count << " 1 " << coords[0] / BOX_LENGTH << " " << coords[1] / BOX_LENGTH << " " << coords[2] / BOX_Z_LENGTH << endl;
        vmd_string << ++atom_count << " 2 " << coords[3] / BOX_LENGTH << " " << coords[4] / BOX_LENGTH << " " << coords[5] / BOX_Z_LENGTH << endl;
        vmd_string << ++atom_count << " 2 " << coords[6] / BOX_LENGTH << " " << coords[7] / BOX_LENGTH << " " << coords[8] / BOX_Z_LENGTH << endl;
    }

    for (int i = 0; i < IONS.size(); i++) {
        coords = IONS[i]->coords;
        charge = (IONS[i]->charge < 0.0) ? 3 : 4;
        vmd_string << ++atom_count << " " << charge << " "
                << coords[0] / BOX_LENGTH << " " << coords[1] / BOX_LENGTH << " " << coords[2] / BOX_Z_LENGTH << endl;
    }
    return vmd_string.str();
}

void test_vmd_output() {
    Simulation * s = new Simulation();
    s->NUM_EQUILIBRATION_SWEEPS = 100000;
    s->equilibrate();
    cout << s->to_vmd(42) << endl;
    return;
}
