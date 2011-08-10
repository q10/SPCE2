/* 
 * File:   main.cpp
 * Author: BENSON J MA
 *
 * Created on August 3, 2011, 4:48 PM
 */

#include "common.h"

void run_tests(int argc, char** argv);

int main(int argc, char** argv) {
    run_tests(argc, argv);
    return 0;
}

void run_tests(int argc, char** argv) {
    //test_water_rotation();
    //test_vmd_output();
    //test_config_output();
    //test_config_input();
    test_radial_dist();
    return;
}
