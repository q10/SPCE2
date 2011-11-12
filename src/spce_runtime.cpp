#include "common.h"

void SPCERuntime::run_umbrella_system() {
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
    simulation->NUM_MC_SWEEPS = 500000;
    simulation->run_mc();

    simulation->SAMPLER_SET->write_config_snapshot();
    cerr << "\n---- END - UMBRELLA SAMPLING ----\n" << endl;
    return;
}

void SPCERuntime::run_all_tests(int argc, char** argv) {
    //test_water_rotation();
    //test_vmd_output();
    //test_config_output();
    //test_config_input();
    //test_radial_dist();
    test_ion_pair_dist();
    return;
}

void SPCERuntime::test_config_input() {
    cerr << "---- BEGIN TEST - CONFIG FILE INPUT ----" << endl
            << "reading from input file sample.config" << endl;
    Simulation * simulation = ConfigReader::new_simulation_with_config("sample.config");
    cout << simulation << endl;
    cerr << "---- END TEST - CONFIG FILE INPUT ----" << endl;
    return;
}

void SPCERuntime::test_vmd_output() {
    cerr << "---- BEGIN TEST - LAMMPSTRJ (VMD) FILE OUTPUT ----" << endl;
    Simulation * s = new Simulation();
    s->SAMPLER_SET->turn_on_lammpstrj_sampler();
    s->DATA_SAMPLING_RATE = 2;
    s->NUM_MC_SWEEPS = 10;
    s->run_mc();
    cerr << "---- END TEST - LAMMPSTRJ (VMD) FILE OUTPUT ----" << endl;
    return;
}

void SPCERuntime::test_config_output() {
    cerr << "---- BEGIN TEST - CONFIG FILE OUTPUT ----" << endl;
    Simulation * s = new Simulation();
    s->NUM_MC_SWEEPS = 100;
    s->run_mc();
    s->SAMPLER_SET->write_config_snapshot();
    cerr << "---- END TEST - CONFIG FILE OUTPUT ----" << endl;
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
    simulation->SAMPLER_SET->add_ion_pair_distance_sampler();
    simulation->run_mc();
    simulation->SAMPLER_SET->print_individual_sampler_results();
    cerr << "---- END TEST - ION PAIR DISTANCE SAMPLER ----" << endl;
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
