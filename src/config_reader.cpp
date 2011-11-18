#include "common.h"

const program_flags_t ConfigReader::PROGRAM_FLAGS[] = {
    { "input", required_argument, NULL, 'r'},
    { "output_config", required_argument, NULL, 'w'},
    { "output_lammpstrj", required_argument, NULL, 'v'},
    { NULL, no_argument, NULL, 0}
};

Simulation * ConfigReader::new_simulation(int argc, char** argv) {
    Simulation * simulation = new Simulation();
    read_program_flags(argc, argv, simulation);
    return simulation;
}

Simulation * ConfigReader::new_simulation_with_config(string input_config_filename) {
    Simulation * simulation = new Simulation();
    load_configuration_file(input_config_filename, simulation);
    return simulation;
}

void ConfigReader::read_program_flags(int argc, char** argv, Simulation * simulation) {
    int option;

    while (1) {
        // getopt_long stores the option index here.
        int options_i = 0;
        string argstr = "r:w:v:ABCDEFGHIJKLMNNOPQRSTUVWXYZabcdefghijklmopqsuxyz";

        // "abc:d:f:" means that a and b don't have args, while c, d, f do
        option = getopt_long_only(argc, argv, argstr.c_str(), PROGRAM_FLAGS, &options_i);

        // Detect the end of the options.
        if (option == -1 or option == '?')
            break;
        else if (option == 0) {
            if (PROGRAM_FLAGS[options_i].flag != 0)
                break;
            printf("option %s", PROGRAM_FLAGS[options_i].name);
            if (optarg)
                printf(" with arg %s", optarg);
            printf("\n");
            break;
        } else if (option == 'r') {
            string input_config_filename = optarg;
            if (input_config_filename[0] == '=')
                input_config_filename = input_config_filename.substr(1);
            load_configuration_file(input_config_filename, simulation);
        } else if (option == 'w') {
            string config_filename = optarg;
            if (config_filename[0] == '=')
                config_filename = config_filename.substr(1);
            simulation->SAMPLER_SET->config_filename = config_filename;
        } else if (option == 'v') {
            string lammpstrj_filename = optarg;
            if (lammpstrj_filename[0] == '=')
                lammpstrj_filename = lammpstrj_filename.substr(1);
            simulation->SAMPLER_SET->lammpstrj_filename = lammpstrj_filename;
        } else
            ASSERT(false, "Invalid program flag, invalid flag parameter, or missing a flag parameter.");
    }
    return;
}

void ConfigReader::load_configuration_file(string input_config_filename, Simulation * simulation) {
    // For code simplicity purposes, this function works as intended ONLY when the config file format is 
    // exactly the same as that produced by the simulation's << operator

    ifstream input_filestream(input_config_filename.c_str());
    ASSERT(input_filestream.is_open(), "Could not open input configuration file " + STRING(input_config_filename));

    bool window_sampling_mode = false;
    double window_lower_bound = 0, window_upper_bound = 0;
    unsigned int line_num = 0, num_ions = 0, num_waters = 0, ewald_nxy = 0, ewald_nz = 0;
    string line, line_key;
    vector< double > coords;
    double ewald_alpha = 0.0, val;

    while (getline(input_filestream, line)) {
        istringstream iss(line);
        line_num++;
        iss >> line_key;

        if (line_key.compare("BOX_LENGTH") == 0) {
            double len;
            iss >> len;
            if (abs(len) > 0.0)
                simulation->BOX_LENGTH = abs(len);
            else
                cerr << "WARNING: Bad box length at line " << line_num << " in config file - Using defaults instead." << endl;

        } else if (line_key.compare("BOX_Z_LENGTH") == 0) {
            double zlen;
            iss >> zlen;
            if (abs(zlen) > 0.0)
                simulation->BOX_LENGTH = abs(zlen);
            else
                cerr << "WARNING: Bad box length at line " << line_num << " in config file - Using defaults instead." << endl;

        } else if (line_key.compare("EWALD_ALPHA") == 0) {
            iss >> ewald_alpha;
        } else if (line_key.compare("EWALD_NXY") == 0) {
            iss >> ewald_nxy;
        } else if (line_key.compare("EWALD_NZ") == 0) {
            iss >> ewald_nz;
        } else if (line_key.compare("ION_PAIR_DISTANCE_WINDOW") == 0) {
            window_sampling_mode = true;
            iss >> window_lower_bound;
            iss >> window_upper_bound;
        } else if (line_key.compare("ION") == 0) {
            coords.clear();
            while (iss >> val)
                coords.push_back(val);
            ASSERT((int) coords.size() == 4, "Not enough parameters on line " + STRING(line_num) + " of config file.");
            ASSERT(coords[3] != 0.0, "Charge is zero for ion on line " + STRING(line_num) + " of config file.");

            if (simulation->IONS.size() < ++num_ions) {
                Ion * ion = new Ion(&coords[0], coords[3], simulation->DISPLACEMENT_DISTANCE,
                        simulation->BOX_LENGTH, simulation->BOX_Z_LENGTH);
                simulation->IONS.push_back(ion);
            } else {
                simulation->IONS[num_ions - 1]->set_coords(&coords[0]);
                simulation->IONS[num_ions - 1]->charge = coords[3];
            }
        } else if (line_key.compare("WATER") == 0) {
            coords.clear();
            while (iss >> val)
                coords.push_back(val);
            ASSERT((int) coords.size() == 9, "Not enough coordinates on line " + STRING(line_num) + " of config file.");

            if (simulation->WATERS.size() < ++num_waters) {
                Water * water = new Water(&coords[0], simulation->DISPLACEMENT_DISTANCE,
                        simulation->DISPLACEMENT_ROTATION, simulation->BOX_LENGTH, simulation->BOX_Z_LENGTH);
                simulation->WATERS.push_back(water);
            } else
                simulation->WATERS[num_waters - 1]->set_coords(&coords[0]);

        } else
            cerr << "WARNING: Malformed line " << line_num << " in config file - Ignoring this line." << endl;
    }

    // removes excess particles if any
    while (simulation->WATERS.size() > num_waters)
        simulation->WATERS.pop_back();
    while (simulation->IONS.size() > num_ions)
        simulation->IONS.pop_back();

    // recalculate energies and ewald tables
    simulation->initialize_all_ewald_tables(ewald_alpha, ewald_nxy, ewald_nz);
    simulation->calculate_and_init_energy();

    // sets window sampling mode on if selected
    if (window_sampling_mode)
        simulation->turn_on_window_sampling_mc(window_lower_bound, window_upper_bound);
    return;
}
