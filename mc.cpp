#include "common.h"

void Simulation::equilibrate() {
    static struct timeval equilibrate_start, equilibrate_end;
    gettimeofday(&equilibrate_start, NULL);
    for (int h = 0; h < NUM_EQUILIBRATION_SWEEPS; h++) {
        mc_sweep();
        cerr << "Equilibration MC sweep " << h + 1 << " of " << NUM_EQUILIBRATION_SWEEPS << " complete." << endl;
    }
    gettimeofday(&equilibrate_end, NULL);
    cerr << "Equilibration completed in " << setprecision(10) << timeval_diff(&equilibrate_end, &equilibrate_start) / 1000000.0 / 3600.0 << " hours." << endl;
    return;
}

void Simulation::run_mc() {
    static struct timeval run_start, run_end;
    sampler->start();
    gettimeofday(&run_start, NULL);
    for (int h = 0; h < NUM_MC_SWEEPS; h++) {
        mc_sweep();
        if (h % sampler->DATA_SAMPLING_RATE == 0)
            sampler->sample_data();
        cerr << "MC sweep " << h + 1 << " of " << NUM_MC_SWEEPS << " complete." << endl;
    }
    sampler->finish();
    gettimeofday(&run_end, NULL);
    cerr << "MC run completed in " << setprecision(10) << timeval_diff(&run_end, &run_start) / 1000000.0 / 3600.0 << " hours." << endl;
    return;
}

void Simulation::mc_sweep() {
    static struct timeval sweep_start, sweep_end;
    gettimeofday(&sweep_start, NULL);
    for (int i = 0; i < NUM_MC_ATTEMPTS_PER_SWEEP; i++) {
        if (RAN3() < 0.5) {
            total_attempted_mc_translations++;
            mc_translate();
        } else {
            total_attempted_mc_rotations++;
            mc_rotate();
        }
    }
    gettimeofday(&sweep_end, NULL);
    cerr << setprecision(10) << timeval_diff(&sweep_end, &sweep_start) / 1000000.0 << endl;
    return;
}

void Simulation::mc_translate() {
    int rand_i = RANDINT(0, WATERS.size() + IONS.size() * ION_PROBABILITY_WEIGHT);
    double old_energy_particle_i = energy_of_particle_with_index(rand_i);
    if (rand_i < WATERS.size())
        WATERS[rand_i]->mc_translate();
    else {
        rand_i = WATERS.size() + ((rand_i - WATERS.size()) % ION_PROBABILITY_WEIGHT);
        IONS[rand_i]->mc_translate();
    }
    if (mc_accept(rand_i, old_energy_particle_i))
        num_successful_mc_translations++;
    return;
}

void Simulation::mc_rotate() {
    int rand_i = RANDINT(0, WATERS.size());
    double old_energy_particle_i = energy_of_particle_with_index(rand_i);
    WATERS[rand_i]->mc_rotate();
    if (mc_accept(rand_i, old_energy_particle_i))
        num_successful_mc_rotations++;
    return;
}

bool Simulation::mc_accept(int index, double old_energy_particle_i) {
    double total_energy_diff = ewald_diff(index) + energy_of_particle_with_index(index) - old_energy_particle_i;
    if (RAN3() < exp(-BETA * total_energy_diff))
        TOTAL_ENERGY += total_energy_diff;
    else {
        // undo the move if move not accepted
        (index < WATERS.size()) ? WATERS[index]->undo_move() : IONS[index - WATERS.size()]->undo_move();

        // reset partial rho_k's
        dcomplex *column;
        int NUM_TOTAL_PARTICLES = WATERS.size() + IONS.size();
        for (int k = 0; k < K_VECTORS.size(); k++) {
            column = RHO_K_VALUES[k];
            column[NUM_TOTAL_PARTICLES] += column[NUM_TOTAL_PARTICLES + 1] - column[index];
            column[index] = column[NUM_TOTAL_PARTICLES + 1];
        }
        column = NULL;
        return false;
    }
    return true;
}
