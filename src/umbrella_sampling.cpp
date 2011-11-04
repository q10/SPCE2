#include "common.h"

void UmbrellaSampling::run_umbrella_system() {
    cerr << "---- BEGIN - UMBRELLA SAMPLING ----" << endl;
    Simulation * simulation = new Simulation();
    simulation->IONS[0]->charge = -1.0;
    simulation->IONS[1]->charge = 1.0;
    // set ion positions
    simulation->NUM_EQUILIBRATION_SWEEPS = 100000;
    simulation->equilibrate();

    // set the new potential energy wall
    // set up how the sampler works, customize output
    simulation->DATA_SAMPLING_RATE = 20;
    simulation->NUM_MC_SWEEPS = 100000;
    simulation->run_mc();

    simulation->SAMPLER_SET->write_config_snapshot();
    cerr << "\n---- END - UMBRELLA SAMPLING ----\n" << endl;
    return;
}
