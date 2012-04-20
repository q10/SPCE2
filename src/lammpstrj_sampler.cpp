#include "common.h"

LAMMPSTRJSampler::LAMMPSTRJSampler(WaterSystem * s) {
    system = s;
}

LAMMPSTRJSampler::~LAMMPSTRJSampler() {
    if (LAMMPSTRJ_FILE.is_open())
        LAMMPSTRJ_FILE.close();
}

void LAMMPSTRJSampler::start() {
    timestep = lammpstrj_snapshot_counter = 0;
    if (!LAMMPSTRJ_FILE.is_open()) {
        string filename = system->NAME + ".lammpstrj";
        LAMMPSTRJ_FILE.open(filename.c_str());
        ASSERT(LAMMPSTRJ_FILE.is_open(), "Could not open LAMMPSTRJ output file.");
    }
    return;
}

void LAMMPSTRJSampler::sample() {
    if (++timestep % RELATIVE_LAMMPSTRJ_SNAPSHOT_RATE == 0) {
        int atom_count = 0, ion_type;
        double *coords;

        LAMMPSTRJ_FILE << setprecision(10) << "ITEM: TIMESTEP" << endl << ++lammpstrj_snapshot_counter << endl
                << "ITEM: NUMBER OF ATOMS" << endl << (3 * system->WATERS.size() + system->IONS.size()) << endl
                << "ITEM: BOX BOUNDS" << endl
                << "0 " << system->BOX_LENGTH << endl
                << "0 " << system->BOX_LENGTH << endl
                << "0 " << system->BOX_Z_LENGTH << endl
                << "ITEM: ATOMS id type x y z diameter q" << endl;

        for (unsigned int i = 0; i < system->WATERS.size(); i++) {
            coords = system->WATERS[i]->coords;
            LAMMPSTRJ_FILE << setprecision(10) << ++atom_count << " 1 " << coords[0] << " " << coords[1] << " " << coords[2] << " " << Water::SIGMA_O << " " << Water::Q_O << endl
                    << ++atom_count << " 2 " << coords[3] << " " << coords[4] << " " << coords[5] << " " << Water::SIGMA_H << " " << Water::Q_H << endl
                    << ++atom_count << " 2 " << coords[6] << " " << coords[7] << " " << coords[8] << " " << Water::SIGMA_H << " " << Water::Q_H << endl;
        }

        for (unsigned int i = 0; i < system->IONS.size(); i++) {
            coords = system->IONS[i]->coords;
            ion_type = (system->IONS[i]->charge < 0.0) ? 3 : 4;
            LAMMPSTRJ_FILE << ++atom_count << " " << ion_type << " "
                    << coords[0] << " " << coords[1] << " " << coords[2] << " " << system->IONS[i]->SIGMA << " " << system->IONS[i]->charge << endl;
        }
    }
    return;
}

void LAMMPSTRJSampler::finish() {
    if (!LAMMPSTRJ_FILE.is_open())
        LAMMPSTRJ_FILE.close();
    return;
}

string LAMMPSTRJSampler::results() {
    stringstream results;
    results << "LAMMPS trajectory data is available in " << system->NAME << ".lammpstrj" << endl;
    return results.str();
}
