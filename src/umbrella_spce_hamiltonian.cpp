#include "common.h"

UmbrellaSPCEHamiltonian::UmbrellaSPCEHamiltonian(System &sys)
: SPCEHamiltonian(sys), WINDOW_LOWER_BOUND(sys.WINDOW_LOWER_BOUND),
WINDOW_UPPER_BOUND(sys.WINDOW_UPPER_BOUND) {
    OVER_THE_WINDOW = false;
}

UmbrellaSPCEHamiltonian::~UmbrellaSPCEHamiltonian() {
}

void UmbrellaSPCEHamiltonian::initialize_calculations() {
    cerr << "UmbrellaSPCEHamiltonian setting up......" << endl;
    ASSERT(IONS.size() >= 2, "Not enough ions to start window sampling.");
    // ADD SOME CHECKS TO ENFORCE PARTICLES INTO A WINDOW ?
    SPCEHamiltonian::initialize_calculations();
    cerr << "UmbrellaSPCEHamiltonian setup done." << endl;
}

double UmbrellaSPCEHamiltonian::total_energy_difference(int index) {
    if (index >= (int) WATERS.size()) {
        double ion_dist = IONS[0]->distance_from(IONS[1]);
        if (ion_dist < WINDOW_LOWER_BOUND or ion_dist >= WINDOW_UPPER_BOUND) {
            OVER_THE_WINDOW = true;
            return D_INFINITY;
        }
    }
    return SPCEHamiltonian::total_energy_difference(index);
}

void UmbrellaSPCEHamiltonian::undo_calculations(int index) {
    if (OVER_THE_WINDOW and index >= (int) WATERS.size()) {
        OVER_THE_WINDOW = false;
        return;
    } else
        SPCEHamiltonian::undo_calculations(index);
}

double UmbrellaSPCEHamiltonian::energy_between_two_ions(int i, int j, WHICH_TYPE typ) {
    double r = (typ == OLD) ? IONS[i]->old_distance_from(IONS[j]) : IONS[i]->distance_from(IONS[j]);
    double rb = Ion::SIGMA / r;
    return -pow(r - HALF_BOX_LENGTH, 2.0) +
            4.0 * Ion::EPSILON * (pow(rb, 12) - pow(rb, 6)) +
            PCONSTANTS::ELECTROSTATIC_K * IONS[i]->charge * IONS[j]->charge * ERFC_TABLE[int(r * 1000.0)] / r;
}
