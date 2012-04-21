#include "common.h"

// For code simplicity purposes, this function works as intended ONLY when the config file format is
// exactly the same as that produced by the system's << operator

ConfigReader::ConfigReader(std::string input_config_filename) {
    ifstream input_filestream(input_config_filename.c_str());
    ASSERT(input_filestream.is_open(), "Could not open input configuration file " + STRING(input_config_filename));

    line_num = num_ions = num_waters = 0;
    system = new WaterSystem();

    while (getline(input_filestream, line)) {
        istringstream iss(line);
        line_num++;
        iss >> line_key;

        if (line_key.compare("BOX_LENGTH") == 0 or line_key.compare("BOX_Z_LENGTH") == 0)
            handle_box_length(iss);
        else if (line_key.compare("EWALD_ALPHA") == 0 or line_key.compare("EWALD_NXY") == 0 or line_key.compare("EWALD_NZ") == 0)
            handle_ewald_parameters(iss);
        else if (line_key.compare("ION_PAIR_DISTANCE_WINDOW") == 0)
            handle_ion_pair_distance_window(iss);
        else if (line_key.compare("ION") == 0)
            handle_ion(iss);
        else if (line_key.compare("WATER") == 0)
            handle_water(iss);
        else
            cerr << "WARNING: Malformed line " << line_num << " in config file - Ignoring this line." << endl;
    }

    // removes excess particles if any
    while (system->WATERS.size() > num_waters)
        system->WATERS.pop_back();
    while (system->IONS.size() > num_ions)
        system->IONS.pop_back();
}

WaterSystem * ConfigReader::new_water_system() {
    WaterSystem * copy = new WaterSystem(*system);
    return copy;
}

void ConfigReader::handle_box_length(istringstream & iss) {
    double len;
    iss >> len;
    len = abs(len);
    if (len > 0.0) {
        if (line_key.compare("BOX_LENGTH") == 0)
            system->BOX_LENGTH = len;
        else
            system->BOX_Z_LENGTH = len;
    } else
        cerr << "WARNING: Bad box length at line " << line_num << " in config file - Using defaults instead." << endl;
}

void ConfigReader::handle_ewald_parameters(std::istringstream & iss) {
    double x = 0;
    iss >> x;
    x = abs(x);
    int y = int(x);
    if (line_key.compare("EWALD_ALPHA") == 0 and x > 0.0)
        system->EWALD_ALPHA = x;
    else if (line_key.compare("EWALD_NXY") == 0 and y > 0)
        system->EWALD_NXY = y;
    else if (line_key.compare("EWALD_NZ") == 0 and y > 0)
        system->EWALD_NZ = y;
    else
        cerr << "WARNING: Bad Ewald parameter at line " << line_num << " in config file - Using defaults instead." << endl;
}

void ConfigReader::handle_ion_pair_distance_window(istringstream & iss) {
    tmp_coords.clear();
    double val;
    while (iss >> val)
        tmp_coords.push_back(val);
    ASSERT((int) tmp_coords.size() >= 2, "Not enough parameters on line " + STRING(line_num) + " of config file.");
    system->WINDOW_LOWER_BOUND = tmp_coords[0];
    system->WINDOW_UPPER_BOUND = tmp_coords[1];
}

void ConfigReader::handle_ion(istringstream & iss) {
    tmp_coords.clear();
    double val;
    while (iss >> val)
        tmp_coords.push_back(val);
    ASSERT((int) tmp_coords.size() >= 4, "Not enough parameters on line " + STRING(line_num) + " of config file.");
    ASSERT(tmp_coords[3] != 0.0, "Charge is zero for ion on line " + STRING(line_num) + " of config file.");

    if (system->IONS.size() < ++num_ions) {
        Ion * ion = new Ion(system, &tmp_coords[0], tmp_coords[3]);
        system->IONS.push_back(ion);
    } else {
        system->IONS[num_ions - 1]->set_coords(&tmp_coords[0]);
        system->IONS[num_ions - 1]->charge = tmp_coords[3];
    }
}

void ConfigReader::handle_water(istringstream & iss) {
    tmp_coords.clear();
    double val;
    while (iss >> val)
        tmp_coords.push_back(val);
    ASSERT((int) tmp_coords.size() >= 9, "Not enough coordinates on line " + STRING(line_num) + " of config file.");

    if (system->WATERS.size() < ++num_waters) {
        Water * water = new Water(system, &tmp_coords[0]);
        system->WATERS.push_back(water);
    } else
        system->WATERS[num_waters - 1]->set_coords(&tmp_coords[0]);
}
