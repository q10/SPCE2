/* 
 * File:   simulation.h
 * Author: BENSON J MA
 *
 * Created on August 3, 2011, 4:58 PM
 */

#ifndef SIMULATION_H
#define	SIMULATION_H

class SamplerSet;

class Simulation {
private:
    int total_attempted_mc_translations, total_attempted_mc_rotations, num_successful_mc_translations, num_successful_mc_rotations;
    void mc_sweep();
    void mc_translate();
    void mc_rotate();
    bool mc_accept(int index, double old_energy_particle_i);


    double energy_between_two_waters(int i, int j);
    double energy_between_ion_and_water(int i, int j);
    double energy_between_two_ions(int i, int j);


    std::vector <double> ERFC_TABLE;
    std::vector <double *> K_VECTORS;
    std::vector <dcomplex *> RHO_K_VALUES;
    dcomplex *exp_kr_O[3], *exp_kr_H1[3], *exp_kr_H2[3], *exp_kr_ion[3];
    void initialize_erfc_table();
    void initialize_k_vectors_table();
    void initialize_rho_k_values_table();
    void initialize_exp_kr_cache_tables();
    dcomplex partial_rho(int index, double * k_coords);
    double ewald_diff(int index);
    double ewald_diff_water(int water_index);
    double ewald_diff_ion(int index);
    void set_exp_kr_table_for_water(int water_index);
    void set_exp_kr_table_for_ion(int ion_index);


public:
    static const double BOLTZMANN_K, ELECTROSTATIC_K;
    static const dcomplex COMPLEX_ONE;


    double TEMPERATURE, BETA, TARGET_WATER_DENSITY, BOX_LENGTH, HALF_BOX_LENGTH, BOX_Z_LENGTH, HALF_BOX_Z_LENGTH, BOX_VOLUME;
    std::vector <Water *> WATERS;
    std::vector <Ion *> IONS;
    SamplerSet * SAMPLER_SET;

    Simulation(int num_waters = 200, int num_ions = 2);
    ~Simulation();

    void default_initialize_system_parameters(int num_waters, int num_ions);
    void default_initialize_waters(int num_waters);
    void default_initialize_ions(int num_ions);
    void default_initialize_sampling_parameters();

    void set_temperature(double new_temp);
    void expand_box_z_direction(double new_len = 0.0);


    double DISPLACEMENT_DISTANCE, DISPLACEMENT_ROTATION;
    int NUM_EQUILIBRATION_SWEEPS, NUM_MC_SWEEPS, NUM_MC_ATTEMPTS_PER_SWEEP, ION_PROBABILITY_WEIGHT;
    void equilibrate();
    void run_mc();


    double TOTAL_ENERGY;
    void calculate_and_init_energy();
    double energy_of_particle_with_index(int index);
    double total_real_space_energy();


    double EWALD_ALPHA;
    int EWALD_NXY, EWALD_NZ, NUM_K_VECTORS, K_111_INDEX;
    void initialize_all_ewald_tables(double ewald_alpha, int ewald_nxy, int ewald_nz);
    double total_ewald_energy();



    bool WINDOW_SAMPLING_MODE;
    int DATA_SAMPLING_RATE;

    friend std::ostream & operator<<(std::ostream & out, Simulation * simulation);
    std::string to_vmd(int time_step);
};


void test_vmd_output();
void test_config_output();

#endif	/* SIMULATION_H */
