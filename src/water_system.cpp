#include "common.h"

WaterSystem::WaterSystem(int num_waters, int num_ions) {
    default_initialize_system_parameters(num_waters, num_ions);
    default_initialize_ewald_parameters(0.0784, 5, 5);
    default_initialize_waters(num_waters);
    default_initialize_ions(num_ions);
}

WaterSystem::~WaterSystem() {
    // should use pointers and call delete WATERS for auto deallocation of waters and ions, or maybe use &WATERS
    WATERS.clear();
    IONS.clear();
    SAMPLERS.clear();
}

void WaterSystem::default_initialize_system_parameters(int num_waters, int num_ions) {
    set_temperature(300.0);
    NAME = "SPCE_" + TIMESTAMP();

    TARGET_WATER_DENSITY = Water::STD_DENSITY;
    BOX_VOLUME = (num_waters + num_ions) / TARGET_WATER_DENSITY;
    BOX_LENGTH = BOX_Z_LENGTH = pow(BOX_VOLUME, 1.0 / 3.0);
    HALF_BOX_LENGTH = HALF_BOX_Z_LENGTH = BOX_LENGTH / 2.0;

    ION_PROBABILITY_WEIGHT = 3;
    DISPLACEMENT_DISTANCE = 0.2;
    DISPLACEMENT_ROTATION = 0.17 * M_PI;
    NUM_EQUILIBRATION_SWEEPS = 5000;
    NUM_MC_SWEEPS = 100000;
    NUM_MC_ATTEMPTS_PER_SWEEP = 1000;

    WINDOW_LOWER_BOUND = WINDOW_UPPER_BOUND = -1;
    DATA_SAMPLING_RATE = 20;

    num_ion_disp_attempts = num_ion_disp_successes = 0;
}

void WaterSystem::default_initialize_waters(int num_waters) {
    WATERS.clear();
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

        WATERS.push_back(new Water(this, coords));
    }
    delete [] coords;
}

void WaterSystem::default_initialize_ions(int num_ions) {
    IONS.clear();
    double *coords = new double[3], charge;

    for (int i = 0; i < num_ions; i++) {
        charge = (RAN3() < 0.5) ? -1.0 : 1.0;
        for (int j = 0; j < 3; j++)
            coords[j] = RAN3() * BOX_LENGTH;
        IONS.push_back(new Ion(this, coords, charge));
    }
    delete [] coords;
}

void WaterSystem::default_initialize_ewald_parameters(double alpha, int nxy, int nz) {
    EWALD_ALPHA = alpha;
    EWALD_NXY = nxy;
    EWALD_NZ = nz;
}

void WaterSystem::set_temperature(double new_temp) {
    TEMPERATURE = new_temp;
    BETA = 1.0 / (PCONSTANTS::BOLTZMANN_K * TEMPERATURE);
}

void WaterSystem::expand_box_z_direction(double new_len) {
    // first increase the box z length, then move all particles up by displacement dz, where dz is such that after the move, water slab will be in the center
    // set water's box z length when expanding box
    if (new_len == 0.0)
        new_len = 2 * BOX_LENGTH;
    ASSERT(new_len > BOX_LENGTH, "INVALID BOX Z LENGTH - MUST BE GREATER THAN BOX LENGTH");
    BOX_Z_LENGTH = new_len;
    HALF_BOX_Z_LENGTH = BOX_Z_LENGTH / 2.0;
    BOX_VOLUME = BOX_LENGTH * BOX_LENGTH * BOX_Z_LENGTH;

    double shift = (BOX_Z_LENGTH - BOX_LENGTH) / 2.0;
    for (unsigned int i = 0; i < WATERS.size(); i++) {
        for (int j = 2; j < 9; j += 3)
            WATERS[i]->coords[j] += shift;
    }

    for (unsigned int i = 0; i < IONS.size(); i++)
        IONS[i]->coords[2] += shift;

    // redo all Ewald tables and recalculate energies - the appropriate update calls will be made later in the runtime
    EWALD_NZ *= int(ceil(BOX_Z_LENGTH / BOX_LENGTH));
}

void WaterSystem::mc_move() {
    if (RAN3() < 0.5) {
        num_ion_disp_attempts++;
        num_ion_disp_successes++;
        TEMP_INDEX = WATERS.size() + RANDINT(0, IONS.size());
        IONS[TEMP_INDEX - WATERS.size()]->mc_translate();
    } else {
        TEMP_INDEX = RANDINT(0, WATERS.size());
        (RAN3() < 0.5) ? WATERS[TEMP_INDEX]->mc_translate() : WATERS[TEMP_INDEX]->mc_rotate();
    }

    /*
    TEMP_INDEX = RANDINT(0, WATERS.size() + IONS.size() * ION_PROBABILITY_WEIGHT);
    if (TEMP_INDEX >= int(WATERS.size())) {
        TEMP_INDEX = WATERS.size() + ((TEMP_INDEX - WATERS.size()) / ION_PROBABILITY_WEIGHT);
        IONS[TEMP_INDEX - WATERS.size()]->mc_translate();
    } else
        (RAN3() < 0.5) ? WATERS[TEMP_INDEX]->mc_translate() : WATERS[TEMP_INDEX]->mc_rotate();
     */

    /*    
     if (RAN3() < 0.5) {
        TEMP_INDEX = RANDINT(0, WATERS.size() + IONS.size() * ION_PROBABILITY_WEIGHT);
        if (TEMP_INDEX >= (int) WATERS.size())
            TEMP_INDEX = WATERS.size() + ((TEMP_INDEX - WATERS.size()) / ION_PROBABILITY_WEIGHT);
        (TEMP_INDEX < (int) WATERS.size()) ? WATERS[TEMP_INDEX]->mc_translate() : IONS[TEMP_INDEX - WATERS.size()]->mc_translate();
     } else {
         TEMP_INDEX = RANDINT(0, WATERS.size());
         WATERS[TEMP_INDEX]->mc_rotate();
     }
     */
}

void WaterSystem::undo_mc_move() {
    if (TEMP_INDEX >= int(WATERS.size()))
        num_ion_disp_successes--;
    (TEMP_INDEX < int(WATERS.size())) ? WATERS[TEMP_INDEX]->undo_move() : IONS[TEMP_INDEX - WATERS.size()]->undo_move();
}

void WaterSystem::add_rdf_sampler() {
    RDFSampler * sampler = new RDFSampler(this);
    SAMPLERS.push_back(dynamic_cast<Sampler *> (sampler));
}

void WaterSystem::add_lammpstrj_sampler() {
    LAMMPSTRJSampler * sampler = new LAMMPSTRJSampler(this);
    SAMPLERS.push_back(dynamic_cast<Sampler *> (sampler));
}

void WaterSystem::add_ion_pair_distance_sampler() {
    IonPairDistanceSampler * sampler = new IonPairDistanceSampler(this);
    SAMPLERS.push_back(dynamic_cast<Sampler *> (sampler));
}

void WaterSystem::write_config_snapshot() {
    ofstream config_file;
    string config_filename = NAME + ".config";
    config_file.open(config_filename.c_str());
    ASSERT(config_file.is_open(), "Could not open config output file.");
    config_file << this;
    config_file.close();
}

void WaterSystem::initialize_sampling() {
    for (unsigned int i = 0; i < SAMPLERS.size(); i++)
        SAMPLERS[i]->start();
}

void WaterSystem::sample_data() {
    for (unsigned int i = 0; i < SAMPLERS.size(); i++)
        SAMPLERS[i]->sample();
}

void WaterSystem::finish_sampling() {
    for (unsigned int i = 0; i < SAMPLERS.size(); i++)
        SAMPLERS[i]->finish();
}

void WaterSystem::print_individual_sampler_results() {
    cerr << "\nSIMULATION STATISTICS:" << setprecision(10) << endl
            << "Acceptance ratio of ion displacement moves" << double(num_ion_disp_successes) / double(num_ion_disp_attempts) << endl;
    for (unsigned int i = 0; i < SAMPLERS.size(); i++)
        cout << SAMPLERS[i]->results() << endl;
}

ostream & operator<<(ostream & out, WaterSystem * system) {
    out << "BOX_LENGTH\t" << setprecision(10) << system->BOX_LENGTH << endl
            << "BOX_Z_LENGTH\t" << system->BOX_Z_LENGTH << endl
            << "EWALD_ALPHA\t" << system->EWALD_ALPHA << endl
            << "EWALD_NXY\t" << system->EWALD_NXY << endl
            << "EWALD_NZ\t" << system->EWALD_NZ << endl;

    out << "ION_PAIR_DISTANCE_WINDOW\t" << system->WINDOW_LOWER_BOUND << "\t" << system->WINDOW_UPPER_BOUND << endl;

    for (unsigned int i = 0; i < system->WATERS.size(); i++)
        out << "WATER\t" << system->WATERS[i] << endl;

    for (unsigned int i = 0; i < system->IONS.size(); i++)
        out << "ION\t" << system->IONS[i] << endl;

    return out;
}
