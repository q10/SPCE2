/* 
 * File:   spce_runtime.h
 * Author: BENSON J MA
 *
 * Created on November 5, 2011, 1:07 AM
 */

#ifndef RUNTIME_H
#define	RUNTIME_H

namespace SPCERuntime {
    void run_umbrella_system(int argc, char** argv);

    void run_all_tests(int argc, char** argv);
    void test_water_rotation();
    void test_config_input();
    void test_lammpstrj_output();
    void test_config_output();
    void test_radial_dist(int argc, char** argv);
    void test_ion_pair_dist();
};

#endif	/* SPCE_RUNTIME_H */
