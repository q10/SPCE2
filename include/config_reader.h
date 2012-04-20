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
    static const program_flags_t PROGRAM_FLAGS[];

    //static void read_program_flags(int argc, char** argv, Simulation * simulation);
    static void handle_water(WaterSystem * system, int & line_num);
    static void handle_ion(WaterSystem * system, int & line_num, int &num_ione);
    static void handle_box_length(WaterSystem * system, int & line_num);
    static void handle_ion_pair_distance_window(WaterSystem * system, int & line_num);

public:
    static WaterSystem * new_water_system(std::string input_config_filename);
};

#endif	/* CONFIG_READER_H */
