#include "common.h"

const double Simulation::BOLTZMANN_K = 0.00831447122,
        Simulation::ELECTROSTATIC_K = 1389.354325379097;
const dcomplex Simulation::COMPLEX_ONE(1.0, 0.0);

Simulation::Simulation(int num_waters, int num_ions) {
    default_initialize_system_parameters(num_waters, num_ions);
    default_initialize_waters(num_waters);
    default_initialize_ions(num_ions);
    initialize_all_ewald_tables(0.0784, 5, 5);
    calculate_and_init_energy();
    sampler = new Sampler(this);
}

Simulation::~Simulation() {
    delete sampler;
}

void Simulation::default_initialize_system_parameters(int num_waters, int num_ions) {
    TEMPERATURE = 300.0;
    BETA = 1.0 / (BOLTZMANN_K * TEMPERATURE);
    TARGET_WATER_DENSITY = Water::STD_DENSITY;
    BOX_VOLUME = (num_waters + num_ions) / TARGET_WATER_DENSITY;
    BOX_LENGTH = pow(BOX_VOLUME, 1.0 / 3.0);

    BOX_Z_LENGTH = BOX_LENGTH;
    HALF_BOX_LENGTH = BOX_LENGTH / 2.0;
    HALF_BOX_Z_LENGTH = HALF_BOX_LENGTH;

    DISPLACEMENT_DISTANCE = 0.2;
    DISPLACEMENT_ROTATION = 0.17 * M_PI;
    NUM_EQUILIBRATION_SWEEPS = 5000;
    NUM_MC_SWEEPS = 1000000;
    NUM_MC_ATTEMPTS_PER_SWEEP = 1000;
    total_attempted_mc_translations = total_attempted_mc_rotations = num_successful_mc_translations = num_successful_mc_rotations = 0;
    return;
}

void Simulation::default_initialize_waters(int num_waters) {
    WATERS.clear();
    Water *w;
    double HOH_ANGLE_RAD = DEG2RAD(Water::HOH_ANGLE_DEG), *coords = new double[9];
    double tmp_r = Water::OH_LENGTH * sin(HOH_ANGLE_RAD), rand_angle_rad;

    for (int i = 0; i < num_waters; i++) {
        // Oxygen
        for (int j = 0; j < 3; j++)
            coords[j] = RAN3() * BOX_LENGTH;

        // First Hydrogen
        coords[3] = coords[0];
        coords[4] = coords[1];
        coords[5] = coords[2] + Water::OH_LENGTH;

        // Second Hydrogen
        rand_angle_rad = 2.0 * M_PI * RAN3();
        coords[6] = coords[0] + tmp_r * cos(rand_angle_rad);
        coords[7] = coords[1] + tmp_r * sin(rand_angle_rad);
        coords[8] = coords[2] + Water::OH_LENGTH * cos(HOH_ANGLE_RAD);

        w = new Water(coords, DISPLACEMENT_DISTANCE, DISPLACEMENT_ROTATION, BOX_LENGTH, BOX_Z_LENGTH);
        WATERS.push_back(w);
    }
    delete [] coords;
    return;
}

void Simulation::default_initialize_ions(int num_ions) {
    IONS.clear();
    Ion *ion;
    double *coords = new double[3], charge;

    for (int i = 0; i < num_ions; i++) {
        charge = (RAN3() < 0.5) ? -1.0 : 1.0;
        for (int j = 0; j < 3; j++)
            coords[j] = RAN3() * BOX_LENGTH;
        ion = new Ion(coords, charge, DISPLACEMENT_DISTANCE, BOX_LENGTH, BOX_Z_LENGTH);
        IONS.push_back(ion);
    }
    delete [] coords;
    return;
}

void Simulation::set_temperature(double new_temp) {
    TEMPERATURE = 300.0;
    BETA = 1.0 / (BOLTZMANN_K * TEMPERATURE);
    return;
}

void Simulation::expand_box_z_direction(double new_len) {
    // first increase the box z length, then move all particles up by displacement dz, where dz is such that after the move, water slab will be in the center
    // set water's box z length when expanding box
    if (new_len == 0.0)
        new_len = 2 * BOX_LENGTH;
    ASSERT(new_len > BOX_LENGTH, "INVALID BOX Z LENGTH - MUST BE GREATER THAN BOX LENGTH");
    BOX_Z_LENGTH = new_len;
    HALF_BOX_Z_LENGTH = BOX_Z_LENGTH / 2.0;
    BOX_VOLUME = BOX_LENGTH * BOX_LENGTH * BOX_Z_LENGTH;

    for (int i = 0; i < WATERS.size(); i++) {
        for (int j = 2; j < 9; j += 3)
            WATERS[i]->coords[j] += (BOX_Z_LENGTH - BOX_LENGTH) / 2.0;
        WATERS[i]->BOX_Z_LENGTH = BOX_Z_LENGTH;
    }

    for (int i = 0; i < IONS.size(); i++) {
        IONS[i]->coords[2] += (BOX_Z_LENGTH - BOX_LENGTH) / 2.0;
        IONS[i]->BOX_Z_LENGTH = BOX_Z_LENGTH;
    }

    // redo all Ewald tables and recalculate energies
    int times = int(ceil(BOX_Z_LENGTH / BOX_LENGTH));
    initialize_all_ewald_tables(EWALD_ALPHA, EWALD_NXY, EWALD_NZ * times);
    calculate_and_init_energy();
    return;
}
