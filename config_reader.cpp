#include "common.h"

const program_flags_t ConfigReader::PROGRAM_FLAGS[] = {
    { "input", required_argument, NULL, 'r'},
    { "output_config", required_argument, NULL, 'w'},
    { "output_vmd", required_argument, NULL, 'v'},
    { NULL, no_argument, NULL, 0}
};
/*
static Simulation * ConfigReader::get_new_simulation(int argc, char** argv) {
    Simulation * simulation = new Simulation();
    read_program_flags(argc, argv, simulation);
    return simulation;
}

static void ConfigReader::read_program_flags(int argc, char** argv, Simulation * simulation) {
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
            //using_input_config_file = true;
            input_config_filename = optarg;
            if (input_config_filename[0] == '=')
                input_config_filename = input_config_filename.substr(1);
            load_configuration_file(input_config_filename, simulation);
        } else if (option == 'w') {
            use_custom_output_config_filename = true;
            output_config_filename = optarg;
            if (output_config_filename[0] == '=')
                output_config_filename = output_config_filename.substr(1);
        } else if (option == 'v') {
            use_custom_output_vmd_filename = true;
            output_vmd_filename = optarg;
            if (output_vmd_filename[0] == '=')
                output_vmd_filename = output_vmd_filename.substr(1);
        } else
            ASSERT(false, "Invalid program flag, invalid flag parameter, or missing a flag parameter.");
    }
    return;
}

static void ConfigReader::load_configuration_file(Simulation * simulation) {
}
 */