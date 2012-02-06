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

    int total_attempted_mc_translations, total_attempted_mc_rotations, num_successful_mc_translations, num_successful_mc_rotations;

    struct timeval start_time, end_time;

    void default_initialize_sampling_parameters();

    void mc_sweep();
    void mc_translate();
    void mc_rotate();
    bool mc_accept(int index);

    int lammpstrj_timestep, lammpstrj_snapshot_counter;
    std::ofstream LAMMPSTRJ_FILE;

public:
    WaterSystem SYSTEM;

    Simulation(int num_waters = 200, int num_ions = 2);
    ~Simulation();

    void equilibrate();
    void run_mc();

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
};

template <class EF> Simulation<EF>::Simulation(int num_waters, int num_ions) {
    total_attempted_mc_translations = total_attempted_mc_rotations = num_successful_mc_translations = num_successful_mc_rotations = 0;
    ENERGY_FUNCTION = new EF(SYSTEM);
}

template <class EF> Simulation<EF>::~Simulation() {
    delete ENERGY_FUNCTION;
    SAMPLERS.clear();
    if (LAMMPSTRJ_FILE.is_open())
        LAMMPSTRJ_FILE.close();
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
    LAMMPSTRJ_FILE << SYSTEM.to_lammpstrj(++lammpstrj_timestep);
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
    return out << &(simulation->SYSTEM);
}

#endif	/* SIMULATION_H */
