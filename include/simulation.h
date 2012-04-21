/* 
 * File:   simulation.h
 * Author: BENSON J MA
 *
 * Created on August 3, 2011, 4:58 PM
 */

#ifndef SIMULATION_H
#define	SIMULATION_H

template <class EF, class SYS> class Simulation;
template <class EF, class SYS> std::ostream & operator<<(std::ostream & out, Simulation<EF, SYS> * simulation);

template <class EF, class SYS> class Simulation {
private:
    EF * ENERGY_FUNCTION;
    struct timeval start_time, end_time, sweep_start, sweep_end;
    void mc_sweep();

public:
    SYS * SYSTEM;

    Simulation(SYS * another_system = NULL);
    ~Simulation();

    void equilibrate();
    void run_mc();

    void print_individual_sampler_results();
    friend std::ostream & operator<< <> (std::ostream & out, Simulation<EF, SYS> * simulation);
};

template <class EF, class SYS> Simulation<EF, SYS>::Simulation(SYS * another_system) {
    SYSTEM = (another_system) ? new SYS(*another_system) : new SYS();
    ENERGY_FUNCTION = new EF(*SYSTEM); // need to fix this copy method to make a true copy of stuff
}

template <class EF, class SYS> Simulation<EF, SYS>::~Simulation() {
    delete SYSTEM;
    delete ENERGY_FUNCTION;
}

template <class EF, class SYS> void Simulation<EF, SYS>::equilibrate() {
    ENERGY_FUNCTION->initialize_calculations();
    gettimeofday(&start_time, NULL);
    for (int h = 0; h < SYSTEM->NUM_EQUILIBRATION_SWEEPS; h++) {
        mc_sweep();
        std::cerr << "Equilibration MC sweep " << h + 1 << " of " << SYSTEM->NUM_EQUILIBRATION_SWEEPS << " complete." << std::endl;
    }
    gettimeofday(&end_time, NULL);
    std::cerr << "Equilibration completed in " << std::setprecision(10) << timeval_diff(&end_time, &start_time) / 1000000.0 / 3600.0 << " hours." << std::endl;
    return;
}

template <class EF, class SYS> void Simulation<EF, SYS>::run_mc() {
    ENERGY_FUNCTION->initialize_calculations();
    SYSTEM->initialize_sampling();
    gettimeofday(&start_time, NULL);
    for (int h = 0; h < SYSTEM->NUM_MC_SWEEPS; h++) {
        mc_sweep();
        if (h % SYSTEM->DATA_SAMPLING_RATE == 0)
            SYSTEM->sample_data();
        std::cerr << "MC sweep " << h + 1 << " of " << SYSTEM->NUM_MC_SWEEPS << " complete." << std::endl;
    }
    SYSTEM->finish_sampling();
    gettimeofday(&end_time, NULL);
    std::cerr << "MC run completed in " << std::setprecision(10) << timeval_diff(&end_time, &start_time) / 1000000.0 / 3600.0 << " hours." << std::endl;
    return;
}

template <class EF, class SYS> void Simulation<EF, SYS>::mc_sweep() {
    gettimeofday(&sweep_start, NULL);
    for (int i = 0; i < SYSTEM->NUM_MC_ATTEMPTS_PER_SWEEP; i++) {
        SYSTEM->mc_move();
        if (RAN3() >= exp(-SYSTEM->BETA * ENERGY_FUNCTION->total_energy_difference())) {
            SYSTEM->undo_mc_move();
            ENERGY_FUNCTION->undo_calculations();
        }
    }
    gettimeofday(&sweep_end, NULL);
    std::cerr << std::setprecision(10) << timeval_diff(&sweep_end, &sweep_start) / 1000000.0 << std::endl;
}

template <class EF, class SYS> void Simulation<EF, SYS>::print_individual_sampler_results() {
    for (unsigned int i = 0; i < SYSTEM->SAMPLERS.size(); i++)
        std::cout << SYSTEM->SAMPLERS[i]->results() << std::endl;
}

template <class EF, class SYS> std::ostream & operator<<(std::ostream & out, Simulation<EF, SYS> * simulation) {
    return out << simulation->SYSTEM << std::endl << simulation->ENERGY_FUNCTION << std::endl;
}

#endif	/* SIMULATION_H */
