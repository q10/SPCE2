#include "common.h"

void Simulation::equilibrate() {
    int equilibrate_start = clock();
    for (int h = 0; h < NUM_EQUILIBRATION_SWEEPS; h++) {
        mc_sweep();
        cerr << "Equilibration MC sweep " << h + 1 << " of " << NUM_EQUILIBRATION_SWEEPS << " complete." << endl;
    }
    cerr << "Equilibration completed in " << setprecision(10) << ((double) (clock() - equilibrate_start) / (double) CLOCKS_PER_SEC) / 3600 << " hours." << endl;
    return;
}

void Simulation::run_mc() {
    int run_start = clock();
    for (int h = 0; h < NUM_MC_SWEEPS; h++) {
        mc_sweep();
        if (h % Sampler::DATA_SAMPLING_RATE == 0)
            sampler->sample_data();
        cerr << "MC sweep " << h + 1 << " of " << NUM_MC_SWEEPS << " complete." << endl;
    }
    cerr << "MC run completed in " << setprecision(10) << ((double) (clock() - run_start) / (double) CLOCKS_PER_SEC) / 3600 << " hours." << endl;
    return;
}

void Simulation::mc_sweep() {
    int sweep_start = clock();
    for (int i = 0; i < NUM_MC_ATTEMPTS_PER_SWEEP; i++) {
        if (RAN3() < 0.5) {
            total_attempted_mc_translations++;
            mc_translate();
        } else {
            total_attempted_mc_rotations++;
            mc_rotate();
        }
    }
    cerr << (double) (clock() - sweep_start) / (double) CLOCKS_PER_SEC << endl;
    return;
}

void Simulation::mc_translate() {
    int rand_i = RANDINT(0, WATERS.size() + IONS.size());
    double old_energy_particle_i = energy_of_particle_with_index(rand_i);
    (rand_i < WATERS.size()) ? WATERS[rand_i]->mc_translate() : IONS[rand_i - WATERS.size()]->mc_translate();
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
        return false;
    }
    return true;
}
