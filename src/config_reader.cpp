#include "common.h"

const program_flags_t ConfigReader::PROGRAM_FLAGS[] = {
    { "input", required_argument, NULL, 'r'},
    { "output_config", required_argument, NULL, 'w'},
    { "output_lammpstrj", required_argument, NULL, 'v'},
    { NULL, no_argument, NULL, 0}
};

/*
void ConfigReader::read_program_flags(int argc, char** argv, system * system) {
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
            load_configuration_file(input_config_filename, system);
        } else if (option == 'w') {
            string config_filename = optarg;
            if (config_filename[0] == '=')
                config_filename = config_filename.substr(1);
            system->SAMPLER_SET->config_filename = config_filename;
        } else if (option == 'v') {
            string lammpstrj_filename = optarg;
            if (lammpstrj_filename[0] == '=')
                lammpstrj_filename = lammpstrj_filename.substr(1);
            system->SAMPLER_SET->lammpstrj_filename = lammpstrj_filename;
        } else
            ASSERT(false, "Invalid program flag, invalid flag parameter, or missing a flag parameter.");
    }
    return;
}
 */

// For code simplicity purposes, this function works as intended ONLY when the config file format is
// exactly the same as that produced by the system's << operator

WaterSystem * ConfigReader::new_water_system(string input_config_filename) {
    WaterSystem * system = new WaterSystem();


    ifstream input_filestream(input_config_filename.c_str());
    ASSERT(input_filestream.is_open(), "Could not open input configuration file " + STRING(input_config_filename));

    unsigned int line_num = 0, num_ions = 0, num_waters = 0;
    string line, line_key;
    while (getline(input_filestream, line)) {
        istringstream iss(line);
        line_num++;
        iss >> line_key;

        if (line_key.compare("BOX_LENGTH") == 0 or line_key.compare("BOX_Z_LENGTH") == 0)
            handle_box_length();
        else if (line_key.compare("EWALD_ALPHA") == 0)
            iss >> system->EWALD_ALPHA;
        else if (line_key.compare("EWALD_NXY") == 0)
            iss >> system->EWALD_NXY;
        else if (line_key.compare("EWALD_NZ") == 0)
            iss >> system->EWALD_NZ;
        else if (line_key.compare("ION_PAIR_DISTANCE_WINDOW") == 0)
            handle_ion_pair_distance_window();
        else if (line_key.compare("ION") == 0)
            handle_ion();
        else if (line_key.compare("WATER") == 0)
            handle_water();
        else
            cerr << "WARNING: Malformed line " << line_num << " in config file - Ignoring this line." << endl;
    }

    // removes excess particles if any
    while (system->WATERS.size() > num_waters)
        system->WATERS.pop_back();
    while (system->IONS.size() > num_ions)
        system->IONS.pop_back();

    return system;
}

void ConfigReader::handle_box_length(WaterSystem * system, int & line_num) {
    double len;
    iss >> len;
    if (abs(len) > 0.0)
        system->BOX_LENGTH = abs(len);
    else
        cerr << "WARNING: Bad box length at line " << line_num << " in config file - Using defaults instead." << endl;


    //box z

    else if () {
        double zlen;
        iss >> zlen;
        if (abs(zlen) > 0.0)
            system->BOX_LENGTH = abs(zlen);
        else
            cerr << "WARNING: Bad box length at line " << line_num << " in config file - Using defaults instead." << endl;

    }

}

void ConfigReader::handle_ion_pair_distance_window(WaterSystem * system, int & line_num) {
    window_sampling_mode = true;
    iss >> window_lower_bound;
    iss >> window_upper_bound;

}

void ConfigReader::handle_ion(WaterSystem * system, int & line_num) {
    vector< double > coords;
    coords.clear();
    while (iss >> val)
        coords.push_back(val);
    ASSERT((int) coords.size() == 4, "Not enough parameters on line " + STRING(line_num) + " of config file.");
    ASSERT(coords[3] != 0.0, "Charge is zero for ion on line " + STRING(line_num) + " of config file.");

    if (system->IONS.size() < ++num_ions) {
        Ion * ion = new Ion(system, &coords[0], coords[3]);
        system->IONS.push_back(ion);
    } else {
        system->IONS[num_ions - 1]->set_coords(&coords[0]);
        system->IONS[num_ions - 1]->charge = coords[3];
    }


}

void ConfigReader::handle_water(WaterSystem * system, int & line_num) {
    vector< double > coords;
    coords.clear();
    while (iss >> val)
        coords.push_back(val);
    ASSERT((int) coords.size() == 9, "Not enough coordinates on line " + STRING(line_num) + " of config file.");

    if (system->WATERS.size() < ++num_waters) {
        Water * water = new Water(system, &coords[0]);
        system->WATERS.push_back(water);
    } else
        system->WATERS[num_waters - 1]->set_coords(&coords[0]);


}