/* 
 * File:   config_reader.h
 * Author: BENSON J MA
 *
 * Created on August 4, 2011, 1:35 AM
 */

#ifndef CONFIG_READER_H
#define	CONFIG_READER_H

typedef struct option program_flags_t;

class ConfigReader {
private:
    unsigned int line_num, num_ions, num_waters;
    std::string line, line_key;
    std::vector< double > tmp_coords;

    WaterSystem * system;

    void handle_water(std::istringstream & iss);
    void handle_ion(std::istringstream & iss);
    void handle_box_length(std::istringstream & iss);
    void handle_ewald_parameters(std::istringstream & iss);
    void handle_ion_pair_distance_window(std::istringstream & iss);

public:
    ConfigReader(std::string input_config_filename);
    ~ConfigReader();
    WaterSystem * new_water_system();
};

#endif	/* CONFIG_READER_H */
