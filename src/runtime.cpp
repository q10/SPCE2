#include "common.h"

void SPCERuntime::run_umbrella_system(int argc, char** argv) {
    cerr << "---- BEGIN - UMBRELLA SAMPLING ----" << endl;
    ASSERT(argc >= 3, "need program parameters - window_lower_bound, window_upper_bound, simulation_name!");

    double window_lower_bound = atof(argv[1]), window_upper_bound = atof(argv[2]);
    WaterSystem * system = new WaterSystem();
    system->NAME = argv[3];
    system->WINDOW_LOWER_BOUND = window_lower_bound;
    system->WINDOW_UPPER_BOUND = window_upper_bound;    
    system->DATA_SAMPLING_RATE = 20;
    system->NUM_EQUILIBRATION_SWEEPS = 100000;
    system->NUM_MC_SWEEPS = 1000000;

    Ion * anion = system->IONS[0];
    Ion * cation = system->IONS[1];
    anion->charge = -1.0;
    cation->charge = 1.0;
    cerr << "Initializing anion-cation distance to be inside window [" << window_lower_bound << ", " << window_upper_bound << "] Angstroms......";
    while (anion->distance_from(cation) < window_lower_bound or anion->distance_from(cation) >= window_upper_bound)
        anion->set_random_coords();
    cerr << "done.\n" << endl;
    
    system->add_rdf_sampler();
    system->add_lammpstrj_sampler();
    system->add_ion_pair_distance_sampler();
    cout << system << endl;
    
    Simulation<UmbrellaSPCEHamiltonian, WaterSystem> simulation(system); // simulation makes a new copy of system and everything in it to make things easier
    simulation.equilibrate();
    simulation.run_mc();
    
    simulation.SYSTEM->write_config_snapshot(); // working with the copy of system that's inside this simulation
    cerr << "\n---- END - UMBRELLA SAMPLING ----\n" << endl;
    return;
}

void SPCERuntime::run_all_tests(int argc, char** argv) {
    run_umbrella_system(argc, argv);
    //test_water_rotation();
    //test_lammpstrj_output();
    //test_config_output();
    //test_config_input();
    //test_radial_dist();
    //test_ion_pair_dist();
    return;
}

/*
void SPCERuntime::test_config_input() {
    cerr << "---- BEGIN TEST - CONFIG FILE INPUT ----" << endl
            << "Reading from input file sample.config...\n" << endl;
    Simulation * simulation = ConfigReader::new_simulation_with_config("sample.config");
    cout << simulation << endl;
    cerr << "\n---- END TEST - CONFIG FILE INPUT ----\n" << endl;
    return;
}
 **/
void SPCERuntime::test_lammpstrj_output() {
    cerr << "---- BEGIN TEST - LAMMPSTRJ (VMD) FILE OUTPUT ----" << endl;
    Simulation<SPCEHamiltonian, WaterSystem> s;
    s.SYSTEM->add_lammpstrj_sampler();
    s.SYSTEM->DATA_SAMPLING_RATE = 2;
    s.SYSTEM->NUM_MC_SWEEPS = 10;
    s.run_mc();
    cerr << "\n---- END TEST - LAMMPSTRJ (VMD) FILE OUTPUT ----\n" << endl;
    return;
}

void SPCERuntime::test_config_output() {
    cerr << "---- BEGIN TEST - CONFIG FILE OUTPUT ----" << endl;
    Simulation<SPCEHamiltonian, WaterSystem> s;
    s.SYSTEM->IONS[0]->charge = -1.0;
    s.SYSTEM->IONS[1]->charge = 1.0;
    s.SYSTEM->NUM_MC_SWEEPS = 10;
    s.SYSTEM->WINDOW_LOWER_BOUND = 9.12;
    s.SYSTEM->WINDOW_UPPER_BOUND = 12.56;
    s.run_mc();
    s.SYSTEM->write_config_snapshot();
    cerr << "\n---- END TEST - CONFIG FILE OUTPUT ----" << endl;
    return;
}

void SPCERuntime::test_radial_dist() {
    cerr << "---- BEGIN TEST - RADIAL DISTRIBUTION SAMPLER ----" << endl;
    Simulation<SPCEHamiltonian, WaterSystem> s;
    s.SYSTEM->IONS[0]->charge = -1.0;
    s.SYSTEM->IONS[1]->charge = 1.0;
    s.SYSTEM->NUM_MC_SWEEPS = 50000;
    s.SYSTEM->add_rdf_sampler();
    s.run_mc();
    s.print_individual_sampler_results();
    cerr << "\n---- END TEST - RADIAL DISTRIBUTION SAMPLER ----\n" << endl;
    return;
}

void SPCERuntime::test_ion_pair_dist() {
    cerr << "---- BEGIN TEST - ION PAIR DISTANCE SAMPLER ----" << endl;
    Simulation<SPCEHamiltonian, WaterSystem> s;
    s.SYSTEM->IONS[0]->charge = -1.0;
    s.SYSTEM->IONS[1]->charge = 1.0;
    s.SYSTEM->NUM_MC_SWEEPS = 50000;
    s.SYSTEM->add_ion_pair_distance_sampler();
    s.run_mc();
    s.print_individual_sampler_results();
    cerr << "---- END TEST - ION PAIR DISTANCE SAMPLER ----\n" << endl;
    return;
}
/*
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
 **/
