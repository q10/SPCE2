/* 
 * File:   simulation.h
 * Author: BENSON J MA
 *
 * Created on August 3, 2011, 4:58 PM
 */

#ifndef SIMULATION_H
#define	SIMULATION_H

template <class EF, class SYS, class SAM> class Simulation;
template <class EF, class SYS, class SAM> std::ostream & operator<<(std::ostream & out, Simulation<EF, SYS, SAM> * simulation);

template <class EF, class SYS, class SAM> class Simulation {
private:
    EF * ENERGY_FUNCTION;
    struct timeval start_time, end_time, sweep_start, sweep_end;
    int lammpstrj_timestep, lammpstrj_snapshot_counter;
    std::ofstream LAMMPSTRJ_FILE;

    void default_initialize_sampling_parameters();
    void mc_sweep();

public:
    SYS SYSTEM;

    Simulation();
    ~Simulation();

    void equilibrate();
    void run_mc();

    std::vector <SAM *> SAMPLERS;
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

    friend std::ostream & operator<< <> (std::ostream & out, Simulation<EF, SYS, SAM> * simulation);
};

template <class EF, class SYS, class SAM> Simulation<EF, SYS, SAM>::Simulation() {
    default_initialize_sampling_parameters();
    ENERGY_FUNCTION = new EF(SYSTEM);
}

template <class EF, class SYS, class SAM> Simulation<EF, SYS, SAM>::~Simulation() {
    delete ENERGY_FUNCTION;
    SAMPLERS.clear();
    if (LAMMPSTRJ_FILE.is_open())
        LAMMPSTRJ_FILE.close();
}

template <class EF, class SYS, class SAM> void Simulation<EF, SYS, SAM>::default_initialize_sampling_parameters() {
    DATA_SAMPLING_RATE = 20;
    IS_LAMMPSTRJ_SAMPLING = false;
    RELATIVE_LAMMPSTRJ_SNAPSHOT_RATE = 10;
    lammpstrj_timestep = 0;
    lammpstrj_snapshot_counter = 0;
}

template <class EF, class SYS, class SAM> void Simulation<EF, SYS, SAM>::equilibrate() {
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

template <class EF, class SYS, class SAM> void Simulation<EF, SYS, SAM>::run_mc() {
    ENERGY_FUNCTION->initialize_calculations();
    initialize_sampling();
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

template <class EF, class SYS, class SAM> void Simulation<EF, SYS, SAM>::mc_sweep() {
    gettimeofday(&sweep_start, NULL);
    for (int i = 0; i < SYSTEM.NUM_MC_ATTEMPTS_PER_SWEEP; i++) {
        SYSTEM.mc_move();
        if (RAN3() >= exp(-SYSTEM.BETA * ENERGY_FUNCTION->total_energy_difference())) {
            SYSTEM.undo_mc_move();
            ENERGY_FUNCTION->undo_calculations();
        }
    }
    gettimeofday(&sweep_end, NULL);
    std::cerr << std::setprecision(10) << timeval_diff(&sweep_end, &sweep_start) / 1000000.0 << std::endl;
}

template <class EF, class SYS, class SAM> void Simulation<EF, SYS, SAM>::initialize_sampling() {
    lammpstrj_timestep = 0;
    for (unsigned int i = 0; i < SAMPLERS.size(); i++)
        SAMPLERS[i]->start();
}

template <class EF, class SYS, class SAM> void Simulation<EF, SYS, SAM>::sample_data() {
    for (unsigned int i = 0; i < SAMPLERS.size(); i++)
        SAMPLERS[i]->sample();
    if (IS_LAMMPSTRJ_SAMPLING and ++lammpstrj_snapshot_counter % RELATIVE_LAMMPSTRJ_SNAPSHOT_RATE == 0)
        write_lammpstrj_snapshot();
}

template <class EF, class SYS, class SAM> void Simulation<EF, SYS, SAM>::finish_sampling() {
    for (unsigned int i = 0; i < SAMPLERS.size(); i++)
        SAMPLERS[i]->finish();
}

template <class EF, class SYS, class SAM> void Simulation<EF, SYS, SAM>::print_individual_sampler_results() {
    for (unsigned int i = 0; i < SAMPLERS.size(); i++)
        std::cout << SAMPLERS[i]->results() << std::endl;
}

template <class EF, class SYS, class SAM> void Simulation<EF, SYS, SAM>::add_rdf_sampler() {
    RDFSampler * sampler = new RDFSampler(&SYSTEM);
    SAMPLERS.push_back(dynamic_cast<SAM *> (sampler));
}

template <class EF, class SYS, class SAM> void Simulation<EF, SYS, SAM>::add_ion_pair_distance_sampler() {
    IonPairDistanceSampler * sampler = new IonPairDistanceSampler(&SYSTEM);
    SAMPLERS.push_back(dynamic_cast<SAM *> (sampler));
}

template <class EF, class SYS, class SAM> void Simulation<EF, SYS, SAM>::write_lammpstrj_snapshot() {
    if (!LAMMPSTRJ_FILE.is_open()) {
        if (lammpstrj_filename.compare("") == 0)
            lammpstrj_filename = SYSTEM.NAME + ".lammpstrj";
        LAMMPSTRJ_FILE.open(lammpstrj_filename.c_str());
        ASSERT(LAMMPSTRJ_FILE.is_open(), "Could not open LAMMPSTRJ output file.");
    }
    LAMMPSTRJ_FILE << SYSTEM.to_lammpstrj(++lammpstrj_timestep);
}

template <class EF, class SYS, class SAM> void Simulation<EF, SYS, SAM>::write_config_snapshot() {
    std::ofstream config_file;
    if (config_filename.compare("") == 0)
        config_filename = SYSTEM.NAME + ".config";
    config_file.open(config_filename.c_str());
    ASSERT(config_file.is_open(), "Could not open config output file.");
    config_file << this;
    config_file.close();
}

template <class EF, class SYS, class SAM> std::ostream & operator<<(std::ostream & out, Simulation<EF, SYS, SAM> * simulation) {
    return out << &(simulation->SYSTEM);
}

#endif	/* SIMULATION_H */
