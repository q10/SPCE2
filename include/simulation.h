/* 
 * File:   simulation.h
 * Author: BENSON J MA
 *
 * Created on August 3, 2011, 4:58 PM
 */

#ifndef SIMULATION_H
#define	SIMULATION_H

template <class EF> class Simulation;
template <class EF> std::ostream & operator<<(std::ostream & out, Simulation<EF> * simulation);

template <class EF> class Simulation {
private:
    EF * ENERGY_FUNCTION;

    struct timeval start_time, end_time;

    int total_attempted_mc_translations, total_attempted_mc_rotations, num_successful_mc_translations, num_successful_mc_rotations;
    void mc_sweep();
    void mc_translate();
    void mc_rotate();
    bool mc_accept(int index);

    int lammpstrj_timestep, lammpstrj_snapshot_counter;
    std::ofstream LAMMPSTRJ_FILE;

public:
    System SYSTEM;

    Simulation(int num_waters = 200, int num_ions = 2);
    ~Simulation();

    void default_initialize_system_parameters(int num_waters, int num_ions);
    void default_initialize_sampling_parameters();
    void default_initialize_waters(int num_waters);
    void default_initialize_ions(int num_ions);
    void default_initialize_ewald_parameters(double alpha, int nxy, int nz);

    void equilibrate();
    void run_mc();

    void set_temperature(double new_temp);
    void expand_box_z_direction(double new_len = 0.0);


    std::vector <Sampler *> SAMPLERS;
    int DATA_SAMPLING_RATE;
    int RELATIVE_LAMMPSTRJ_SNAPSHOT_RATE;

    void initialize_sampling();
    void sample_data();
    void finish_sampling();
    void print_individual_sampler_results();

    //void turn_on_lammpstrj_sampler();
    void add_rdf_sampler();
    void add_ion_pair_distance_sampler();

    std::string config_filename, lammpstrj_filename;
    void write_lammpstrj_snapshot();
    void write_config_snapshot();

    bool IS_LAMMPSTRJ_SAMPLING;

    friend std::ostream & operator<< <> (std::ostream & out, Simulation<EF> * simulation);
    std::string to_lammpstrj(int time_step);
};

template <class EF> Simulation<EF>::Simulation(int num_waters, int num_ions) {
    default_initialize_system_parameters(num_waters, num_ions);
    default_initialize_ewald_parameters(0.0784, 5, 5);
    default_initialize_waters(num_waters);
    default_initialize_ions(num_ions);
    default_initialize_sampling_parameters();
    ENERGY_FUNCTION = new EF(SYSTEM);
}

template <class EF> Simulation<EF>::~Simulation() {
    delete ENERGY_FUNCTION;
    SAMPLERS.clear();
    if (LAMMPSTRJ_FILE.is_open())
        LAMMPSTRJ_FILE.close();
}

template <class EF> void Simulation<EF>::default_initialize_system_parameters(int num_waters, int num_ions) {
    set_temperature(300.0);
    SYSTEM.NAME = "SPCE_" + TIMESTAMP();

    SYSTEM.TARGET_WATER_DENSITY = Water::STD_DENSITY;
    SYSTEM.BOX_VOLUME = (num_waters + num_ions) / SYSTEM.TARGET_WATER_DENSITY;
    SYSTEM.BOX_LENGTH = SYSTEM.BOX_Z_LENGTH = pow(SYSTEM.BOX_VOLUME, 1.0 / 3.0);
    SYSTEM.HALF_BOX_LENGTH = SYSTEM.HALF_BOX_Z_LENGTH = SYSTEM.BOX_LENGTH / 2.0;

    SYSTEM.ION_PROBABILITY_WEIGHT = 3;
    SYSTEM.DISPLACEMENT_DISTANCE = 0.2;
    SYSTEM.DISPLACEMENT_ROTATION = 0.17 * M_PI;
    SYSTEM.NUM_EQUILIBRATION_SWEEPS = 5000;
    SYSTEM.NUM_MC_SWEEPS = 1000000;
    SYSTEM.NUM_MC_ATTEMPTS_PER_SWEEP = 1000;
    total_attempted_mc_translations = total_attempted_mc_rotations = num_successful_mc_translations = num_successful_mc_rotations = 0;

    SYSTEM.WINDOW_LOWER_BOUND = SYSTEM.WINDOW_UPPER_BOUND = -1;
}

template <class EF> void Simulation<EF>::default_initialize_waters(int num_waters) {
    SYSTEM.WATERS.clear();
    double HOH_ANGLE_RAD = DEG2RAD(Water::HOH_ANGLE_DEG), *coords = new double[9];
    double tmp_r = Water::OH_LENGTH * sin(HOH_ANGLE_RAD), rand_angle_rad;

    for (int i = 0; i < num_waters; i++) {
        // Oxygen
        for (int j = 0; j < 3; j++)
            coords[j] = RAN3() * SYSTEM.BOX_LENGTH;

        // First Hydrogen
        coords[3] = coords[0];
        coords[4] = coords[1];
        coords[5] = coords[2] + Water::OH_LENGTH;

        // Second Hydrogen
        rand_angle_rad = 2.0 * M_PI * RAN3();
        coords[6] = coords[0] + tmp_r * cos(rand_angle_rad);
        coords[7] = coords[1] + tmp_r * sin(rand_angle_rad);
        coords[8] = coords[2] + Water::OH_LENGTH * cos(HOH_ANGLE_RAD);

        SYSTEM.WATERS.push_back(new Water(&SYSTEM, coords));
    }
    delete [] coords;
}

template <class EF> void Simulation<EF>::default_initialize_ions(int num_ions) {
    SYSTEM.IONS.clear();
    double *coords = new double[3], charge;

    for (int i = 0; i < num_ions; i++) {
        charge = (RAN3() < 0.5) ? -1.0 : 1.0;
        for (int j = 0; j < 3; j++)
            coords[j] = RAN3() * SYSTEM.BOX_LENGTH;
        SYSTEM.IONS.push_back(new Ion(&SYSTEM, coords, charge));
    }
    delete [] coords;
}

template <class EF> void Simulation<EF>::default_initialize_ewald_parameters(double alpha, int nxy, int nz) {
    SYSTEM.EWALD_ALPHA = alpha;
    SYSTEM.EWALD_NXY = nxy;
    SYSTEM.EWALD_NZ = nz;
}

template <class EF> void Simulation<EF>::set_temperature(double new_temp) {
    SYSTEM.TEMPERATURE = new_temp;
    SYSTEM.BETA = 1.0 / (PCONSTANTS::BOLTZMANN_K * SYSTEM.TEMPERATURE);
}

template <class EF> void Simulation<EF>::expand_box_z_direction(double new_len) {
    // first increase the box z length, then move all particles up by displacement dz, where dz is such that after the move, water slab will be in the center
    // set water's box z length when expanding box
    if (new_len == 0.0)
        new_len = 2 * SYSTEM.BOX_LENGTH;
    ASSERT(new_len > SYSTEM.BOX_LENGTH, "INVALID BOX Z LENGTH - MUST BE GREATER THAN BOX LENGTH");
    SYSTEM.BOX_Z_LENGTH = new_len;
    SYSTEM.HALF_BOX_Z_LENGTH = SYSTEM.BOX_Z_LENGTH / 2.0;
    SYSTEM.BOX_VOLUME = SYSTEM.BOX_LENGTH * SYSTEM.BOX_LENGTH * SYSTEM.BOX_Z_LENGTH;

    double shift = (SYSTEM.BOX_Z_LENGTH - SYSTEM.BOX_LENGTH) / 2.0;
    for (unsigned int i = 0; i < SYSTEM.WATERS.size(); i++) {
        for (int j = 2; j < 9; j += 3)
            SYSTEM.WATERS[i]->coords[j] += shift;
    }

    for (unsigned int i = 0; i < SYSTEM.IONS.size(); i++)
        SYSTEM.IONS[i]->coords[2] += shift;

    // redo all Ewald tables and recalculate energies
    SYSTEM.EWALD_NZ *= int(ceil(SYSTEM.BOX_Z_LENGTH / SYSTEM.BOX_LENGTH));
    ENERGY_FUNCTION->initialize_calculations();
}

template <class EF> void Simulation<EF>::default_initialize_sampling_parameters() {
    DATA_SAMPLING_RATE = 20;
    IS_LAMMPSTRJ_SAMPLING = false;
    RELATIVE_LAMMPSTRJ_SNAPSHOT_RATE = 10;
    lammpstrj_timestep = 0;
    lammpstrj_snapshot_counter = 0;
}

template <class EF> void Simulation<EF>::equilibrate() {
    ENERGY_FUNCTION->initialize_calculations();
    gettimeofday(&start_time, NULL);
    for (int h = 0; h < SYSTEM.NUM_EQUILIBRATION_SWEEPS; h++) {
        mc_sweep();
        std::cerr << "Equilibration MC sweep " << h + 1 << " of " << SYSTEM.NUM_EQUILIBRATION_SWEEPS << " complete." << std::endl;
    }
    gettimeofday(&end_time, NULL);
    std::cerr << "Equilibration completed in " << std::setprecision(10) << timeval_diff(&end_time, &start_time) / 1000000.0 / 3600.0 << " hours." << std::endl;
    return;
}

template <class EF> void Simulation<EF>::run_mc() {
    initialize_sampling();
    ENERGY_FUNCTION->initialize_calculations();
    gettimeofday(&start_time, NULL);
    for (int h = 0; h < SYSTEM.NUM_MC_SWEEPS; h++) {
        mc_sweep();
        if (h % DATA_SAMPLING_RATE == 0)
            sample_data();
        std::cerr << "MC sweep " << h + 1 << " of " << SYSTEM.NUM_MC_SWEEPS << " complete." << std::endl;
    }
    finish_sampling();
    gettimeofday(&end_time, NULL);
    std::cerr << "MC run completed in " << std::setprecision(10) << timeval_diff(&end_time, &start_time) / 1000000.0 / 3600.0 << " hours." << std::endl;
    return;
}

template <class EF> void Simulation<EF>::mc_sweep() {
    static struct timeval sweep_start, sweep_end;
    gettimeofday(&sweep_start, NULL);
    for (int i = 0; i < SYSTEM.NUM_MC_ATTEMPTS_PER_SWEEP; i++) {
        if (RAN3() < 0.5) {
            total_attempted_mc_translations++;
            mc_translate();
        } else {
            total_attempted_mc_rotations++;
            mc_rotate();
        }
    }
    gettimeofday(&sweep_end, NULL);
    std::cerr << std::setprecision(10) << timeval_diff(&sweep_end, &sweep_start) / 1000000.0 << std::endl;
    return;
}

template <class EF> void Simulation<EF>::mc_translate() {
    int rand_i = RANDINT(0, SYSTEM.WATERS.size() + SYSTEM.IONS.size() * SYSTEM.ION_PROBABILITY_WEIGHT);
    if (rand_i >= (int) SYSTEM.WATERS.size())
        rand_i = SYSTEM.WATERS.size() + ((rand_i - SYSTEM.WATERS.size()) / SYSTEM.ION_PROBABILITY_WEIGHT);

    (rand_i < (int) SYSTEM.WATERS.size()) ? SYSTEM.WATERS[rand_i]->mc_translate() : SYSTEM.IONS[rand_i - SYSTEM.WATERS.size()]->mc_translate();
    if (mc_accept(rand_i))
        num_successful_mc_translations++;
    return;
}

template <class EF> void Simulation<EF>::mc_rotate() {
    int rand_i = RANDINT(0, SYSTEM.WATERS.size());
    SYSTEM.WATERS[rand_i]->mc_rotate();
    if (mc_accept(rand_i))
        num_successful_mc_rotations++;
    return;
}

template <class EF> bool Simulation<EF>::mc_accept(int index) {
    // total_energy_difference also updates TOTAL_ENERGY
    if (RAN3() < exp(-SYSTEM.BETA * ENERGY_FUNCTION->total_energy_difference(index)))
        return true;
    else {
        // undo the move if move not accepted
        (index < (int) SYSTEM.WATERS.size()) ? SYSTEM.WATERS[index]->undo_move() : SYSTEM.IONS[index - SYSTEM.WATERS.size()]->undo_move();
        ENERGY_FUNCTION->undo_calculations(index);
        return false;
    }
}

template <class EF> void Simulation<EF>::initialize_sampling() {
    lammpstrj_timestep = 0;
    for (unsigned int i = 0; i < SAMPLERS.size(); i++)
        SAMPLERS[i]->start();
}

template <class EF> void Simulation<EF>::sample_data() {
    for (unsigned int i = 0; i < SAMPLERS.size(); i++)
        SAMPLERS[i]->sample();
    if (IS_LAMMPSTRJ_SAMPLING and ++lammpstrj_snapshot_counter % RELATIVE_LAMMPSTRJ_SNAPSHOT_RATE == 0)
        write_lammpstrj_snapshot();
}

template <class EF> void Simulation<EF>::finish_sampling() {
    for (unsigned int i = 0; i < SAMPLERS.size(); i++)
        SAMPLERS[i]->finish();
}

template <class EF> void Simulation<EF>::print_individual_sampler_results() {
    for (unsigned int i = 0; i < SAMPLERS.size(); i++)
        std::cout << SAMPLERS[i]->results() << std::endl;
}

template <class EF> void Simulation<EF>::add_rdf_sampler() {
    RDFSampler * sampler = new RDFSampler(&SYSTEM);
    SAMPLERS.push_back(dynamic_cast<Sampler *> (sampler));
}

template <class EF> void Simulation<EF>::add_ion_pair_distance_sampler() {
    IonPairDistanceSampler * sampler = new IonPairDistanceSampler(&SYSTEM);
    SAMPLERS.push_back(dynamic_cast<Sampler *> (sampler));
}

template <class EF> void Simulation<EF>::write_lammpstrj_snapshot() {
    if (!LAMMPSTRJ_FILE.is_open()) {
        if (lammpstrj_filename.compare("") == 0)
            lammpstrj_filename = SYSTEM.NAME + ".lammpstrj";
        LAMMPSTRJ_FILE.open(lammpstrj_filename.c_str());
        ASSERT(LAMMPSTRJ_FILE.is_open(), "Could not open LAMMPSTRJ output file.");
    }
    LAMMPSTRJ_FILE << to_lammpstrj(++lammpstrj_timestep);
}

template <class EF> void Simulation<EF>::write_config_snapshot() {
    std::ofstream config_file;
    if (config_filename.compare("") == 0)
        config_filename = SYSTEM.NAME + ".config";
    config_file.open(config_filename.c_str());
    ASSERT(config_file.is_open(), "Could not open config output file.");
    config_file << this;
    config_file.close();
}

template <class EF> std::ostream & operator<<(std::ostream & out, Simulation<EF> * simulation) {
    // The policy is that the config file output only contains the configuration of the box and its particles;
    // all other system parameters will be set in code and compiled before running;
    // otherwise, the complexity of managing the config file will grow exponentially
    out << "BOX_LENGTH\t" << std::setprecision(10) << simulation->SYSTEM.BOX_LENGTH << std::endl
            << "BOX_Z_LENGTH\t" << simulation->SYSTEM.BOX_Z_LENGTH << std::endl
            << "EWALD_ALPHA\t" << simulation->SYSTEM.EWALD_ALPHA << std::endl
            << "EWALD_NXY\t" << simulation->SYSTEM.EWALD_NXY << std::endl
            << "EWALD_NZ\t" << simulation->SYSTEM.EWALD_NZ << std::endl;

    out << "ION_PAIR_DISTANCE_WINDOW\t" << simulation->SYSTEM.WINDOW_LOWER_BOUND << "\t" << simulation->SYSTEM.WINDOW_UPPER_BOUND << std::endl;

    for (unsigned int i = 0; i < simulation->SYSTEM.WATERS.size(); i++)
        out << "WATER\t" << simulation->SYSTEM.WATERS[i] << std::endl;

    for (unsigned int i = 0; i < simulation->SYSTEM.IONS.size(); i++)
        out << "ION\t" << simulation->SYSTEM.IONS[i] << std::endl;

    return out;
}

template <class EF> std::string Simulation<EF>::to_lammpstrj(int time_step) {
    std::stringstream lammpstrj_string;
    int atom_count = 0, ion_id;
    double *coords;
    lammpstrj_string << std::setprecision(10) << "ITEM: TIMESTEP" << std::endl << time_step << std::endl
            << "ITEM: NUMBER OF ATOMS" << std::endl << (3 * SYSTEM.WATERS.size() + SYSTEM.IONS.size()) << std::endl
            << "ITEM: BOX BOUNDS" << std::endl
            << "0 " << SYSTEM.BOX_LENGTH << std::endl
            << "0 " << SYSTEM.BOX_LENGTH << std::endl
            << "0 " << SYSTEM.BOX_Z_LENGTH << std::endl
            << "ITEM: ATOMS id type x y z diameter q" << std::endl;
    for (unsigned int i = 0; i < SYSTEM.WATERS.size(); i++) {
        coords = SYSTEM.WATERS[i]->coords;
        lammpstrj_string << ++atom_count << " 1 " << coords[0] << " " << coords[1] << " " << coords[2] << " " << Water::SIGMA_O << " " << Water::Q_O << std::endl;
        lammpstrj_string << ++atom_count << " 2 " << coords[3] << " " << coords[4] << " " << coords[5] << " " << Water::SIGMA_H << " " << Water::Q_H << std::endl;
        lammpstrj_string << ++atom_count << " 2 " << coords[6] << " " << coords[7] << " " << coords[8] << " " << Water::SIGMA_H << " " << Water::Q_H << std::endl;
    }

    for (unsigned int i = 0; i < SYSTEM.IONS.size(); i++) {
        coords = SYSTEM.IONS[i]->coords;
        ion_id = (SYSTEM.IONS[i]->charge < 0.0) ? 3 : 4;
        lammpstrj_string << ++atom_count << " " << ion_id << " "
                << coords[0] << " " << coords[1] << " " << coords[2] << " " << SYSTEM.IONS[i]->SIGMA << " " << SYSTEM.IONS[i]->charge << std::endl;
    }
    return lammpstrj_string.str();
}

#endif	/* SIMULATION_H */
