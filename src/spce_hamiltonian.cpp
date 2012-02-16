#include "common.h"

SPCEHamiltonian::SPCEHamiltonian(WaterSystem &sys)
: EWALD_ALPHA(sys.EWALD_ALPHA), EWALD_NXY(sys.EWALD_NXY),
EWALD_NZ(sys.EWALD_NZ), TARGET_WATER_DENSITY(sys.TARGET_WATER_DENSITY),
BOX_LENGTH(sys.BOX_LENGTH), HALF_BOX_LENGTH(sys.HALF_BOX_LENGTH),
BOX_Z_LENGTH(sys.BOX_Z_LENGTH), HALF_BOX_Z_LENGTH(sys.HALF_BOX_Z_LENGTH),
BOX_VOLUME(sys.BOX_VOLUME), WATERS(sys.WATERS), IONS(sys.IONS),
TOTAL_ENERGY(sys.TOTAL_ENERGY), TEMP_INDEX(sys.TEMP_INDEX) {
}

SPCEHamiltonian::~SPCEHamiltonian() {
}

void SPCEHamiltonian::initialize_calculations() {
    cerr << "SPCEHamiltonian setting up Ewald tables and initializing energy calcultions......";
    initialize_all_ewald_tables();
    TOTAL_ENERGY = total_real_space_energy() + total_ewald_energy();
    cerr << "done." << endl;
}

double SPCEHamiltonian::total_energy_difference() {
    TEMP_ENERGY_DIFF = ewald_diff(TEMP_INDEX) + energy_diff_of_particle_with_index(TEMP_INDEX);
    TOTAL_ENERGY += TEMP_ENERGY_DIFF;
    return TEMP_ENERGY_DIFF;
}

void SPCEHamiltonian::undo_calculations() {
    TOTAL_ENERGY -= TEMP_ENERGY_DIFF;

    // reset partial rho_k's
    dcomplex *column;
    int NUM_TOTAL_PARTICLES = WATERS.size() + IONS.size();
    for (unsigned int k = 0; k < K_VECTORS.size(); k++) {
        column = RHO_K_VALUES[k];
        column[NUM_TOTAL_PARTICLES] += column[NUM_TOTAL_PARTICLES + 1] - column[TEMP_INDEX];
        column[TEMP_INDEX] = column[NUM_TOTAL_PARTICLES + 1];
    }
}
