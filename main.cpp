/* 
 * File:   main.cpp
 * Author: BENSON J MA
 *
 * Created on August 3, 2011, 4:48 PM
 */

#include "common.h"

void run_tests(int argc, char** argv);

int main(int argc, char** argv) {
    vector<vector<float> > m_data;
    m_data = vector<vector<float> >(10);
    for (int i = 0; i < 10; i++) {
        
        m_data[i] = vector<float>(3);
    }
    for (int i = 1; i < 10; i++)
        cout << m_data[i][0] << endl;
    //run_tests(argc, argv);
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
