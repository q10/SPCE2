/* 
 * File:   spce_hamiltonian.h
 * Author: BENSON J MA
 *
 * Created on February 1, 2012, 3:55 PM
 */

#ifndef SPCE_HAMILTONIAN_H
#define	SPCE_HAMILTONIAN_H

class SPCEHamiltonian {
protected:
    double &TEMPERATURE, &BETA, &TARGET_WATER_DENSITY, &BOX_LENGTH, &HALF_BOX_LENGTH, &BOX_Z_LENGTH, &HALF_BOX_Z_LENGTH, &BOX_VOLUME;
    std::vector <Water *> &WATERS;
    std::vector <Ion *> &IONS;

    double &EWALD_ALPHA;
    int &EWALD_NXY, &EWALD_NZ, NUM_K_VECTORS, K_111_INDEX;
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
    void initialize_all_ewald_tables();
    double total_ewald_energy();


    double &TOTAL_ENERGY, TEMP_ENERGY_DIFF;
    enum WHICH_TYPE {OLD, CURRENT, WATER, ION};
    double energy_of_particle_with_index(int index);
    double energy_diff_of_particle_with_index(int index);
    double energy_between_two_waters(int i, int j, WHICH_TYPE typ);
    double energy_between_ion_and_water(int i, int j, WHICH_TYPE typ);
    virtual double energy_between_two_ions(int i, int j, WHICH_TYPE typ);
    double total_real_space_energy();



public:
    SPCEHamiltonian(System &sys);
    virtual ~SPCEHamiltonian();
    virtual void initialize_calculations();
    virtual double total_energy_difference(int index);
    virtual void undo_calculations(int index);
};

#endif
