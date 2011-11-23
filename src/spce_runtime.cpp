#include "common.h"

void SPCERuntime::run_umbrella_system() {
    cerr << "---- BEGIN - UMBRELLA SAMPLING ----" << endl;
    double window_lower_bound = 0.0, window_upper_bound = 2.0;
    Simulation * simulation = new Simulation();

    Ion * anion = simulation->IONS[0];
    Ion * cation = simulation->IONS[1];
    anion->charge = -1.0;
    cation->charge = 1.0;

    cerr << "Initializing anion-cation distance to be inside window [" << window_lower_bound << ", " << window_upper_bound << "] Angstroms......";
    while (anion->distance_from(cation) < window_lower_bound or anion->distance_from(cation) >= window_upper_bound)
        anion->set_random_coords();
    cerr << "done.\n" << endl;

    simulation->NUM_EQUILIBRATION_SWEEPS = 100000;
    simulation->turn_on_window_sampling_mc(window_lower_bound, window_upper_bound);
    simulation->equilibrate();

    simulation->SAMPLER_SET->turn_on_lammpstrj_sampler();
    simulation->SAMPLER_SET->add_ion_pair_distance_sampler();
    simulation->DATA_SAMPLING_RATE = 20;
    simulation->NUM_MC_SWEEPS = 500000;
    simulation->run_mc();

    simulation->SAMPLER_SET->write_config_snapshot();
    cerr << "\n---- END - UMBRELLA SAMPLING ----\n" << endl;
    return;
}

void SPCERuntime::run_all_tests(int argc, char** argv) {
    run_umbrella_system();
    //test_water_rotation();
    //test_lammpstrj_output();
    //test_config_output();
    //test_config_input();
    //test_radial_dist();
    //test_ion_pair_dist();
    return;
}

void SPCERuntime::test_config_input() {
    cerr << "---- BEGIN TEST - CONFIG FILE INPUT ----" << endl
            << "Reading from input file sample.config...\n" << endl;
    Simulation * simulation = ConfigReader::new_simulation_with_config("sample.config");
    cout << simulation << endl;
    cerr << "\n---- END TEST - CONFIG FILE INPUT ----\n" << endl;
    return;
}

void SPCERuntime::test_lammpstrj_output() {
    cerr << "---- BEGIN TEST - LAMMPSTRJ (VMD) FILE OUTPUT ----" << endl;
    Simulation * s = new Simulation();
    s->SAMPLER_SET->turn_on_lammpstrj_sampler();
    s->DATA_SAMPLING_RATE = 2;
    s->NUM_MC_SWEEPS = 10;
    s->run_mc();
    cerr << "\n---- END TEST - LAMMPSTRJ (VMD) FILE OUTPUT ----\n" << endl;
    return;
}

void SPCERuntime::test_config_output() {
    cerr << "---- BEGIN TEST - CONFIG FILE OUTPUT ----" << endl;
    Simulation * s = new Simulation();
    s->IONS[0]->charge = -1.0;
    s->IONS[1]->charge = 1.0;
    s->NUM_MC_SWEEPS = 10;
    s->turn_on_window_sampling_mc(9.12, 12.56);
    s->run_mc();
    s->SAMPLER_SET->write_config_snapshot();
    cerr << "\n---- END TEST - CONFIG FILE OUTPUT ----" << endl;
    return;
}

void SPCERuntime::test_radial_dist() {
    cerr << "---- BEGIN TEST - RADIAL DISTRIBUTION SAMPLER ----" << endl;
    Simulation * simulation = new Simulation();
    simulation->IONS[0]->charge = -1.0;
    simulation->IONS[1]->charge = 1.0;
    simulation->NUM_MC_SWEEPS = 50000;
    simulation->SAMPLER_SET->add_rdf_sampler();
    simulation->run_mc();
    simulation->SAMPLER_SET->print_individual_sampler_results();
    cerr << "\n---- END TEST - RADIAL DISTRIBUTION SAMPLER ----\n" << endl;
    return;
}

void SPCERuntime::test_ion_pair_dist() {
    cerr << "---- BEGIN TEST - ION PAIR DISTANCE SAMPLER ----" << endl;
    Simulation * simulation = new Simulation();
    simulation->IONS[0]->charge = -1.0;
    simulation->IONS[1]->charge = 1.0;
    simulation->NUM_MC_SWEEPS = 50000;
    simulation->SAMPLER_SET->add_ion_pair_distance_sampler();
    simulation->run_mc();
    simulation->SAMPLER_SET->print_individual_sampler_results();
    cerr << "---- END TEST - ION PAIR DISTANCE SAMPLER ----\n" << endl;
    return;
}

void SPCERuntime::test_water_rotation() {
    cerr << "---- BEGIN TEST - WATER ROTATION ----" << endl;
    double * coords = new double [9];
    coords[0] = 1.34;
    coords[1] = 2.0;
    coords[2] = 3.2;
    for (int i = 3; i < 9; i++)
        coords[i] = RAN3()*5.0;
    Water * w = new Water(coords, 0.2, 0.17, 20.0, 20.0);
    delete [] coords;
    for (int k = 0; k < 10000; k++) {
        w->mc_rotate();
        cout << w << endl;
    }
    cerr << "---- END TEST - WATER ROTATION ----\n" << endl;
    return;
}
