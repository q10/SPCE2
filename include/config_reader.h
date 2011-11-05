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

    static void read_program_flags(int argc, char** argv, Simulation * simulation);
    static void load_configuration_file(std::string input_config_filename, Simulation * simulation);

public:
    static Simulation * new_simulation(int argc, char** argv);
    static Simulation * new_simulation_with_config(std::string input_config_filename);
};

#endif	/* CONFIG_READER_H */
