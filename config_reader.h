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
    static void load_configuration_file(Simulation * simulation);

public:
    static Simulation * get_new_simulation(int argc, char** argv);
};

#endif	/* CONFIG_READER_H */
