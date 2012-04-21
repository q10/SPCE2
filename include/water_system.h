/* 
 * File:   system.h
 * Author: BENSON J MA
 *
 * Created on February 1, 2012, 8:59 PM
 */

#ifndef WATER_SYSTEM_H
#define	WATER_SYSTEM_H

class Water;
class Ion;
class Sampler;

class WaterSystem {
private:
    void default_initialize_system_parameters(int num_waters, int num_ions);
    void default_initialize_sampling_parameters();
    void default_initialize_waters(int num_waters);
    void default_initialize_ions(int num_ions);
    void default_initialize_ewald_parameters(double alpha, int nxy, int nz);

public:
    std::string NAME;

    double EWALD_ALPHA;
    int EWALD_NXY, EWALD_NZ;

    double WINDOW_LOWER_BOUND, WINDOW_UPPER_BOUND;

    double TEMPERATURE, BETA, TARGET_WATER_DENSITY, BOX_LENGTH, HALF_BOX_LENGTH, BOX_Z_LENGTH, HALF_BOX_Z_LENGTH, BOX_VOLUME;
    std::vector <Water *> WATERS;
    std::vector <Ion *> IONS;

    double DISPLACEMENT_DISTANCE, DISPLACEMENT_ROTATION;
    int NUM_EQUILIBRATION_SWEEPS, NUM_MC_SWEEPS, NUM_MC_ATTEMPTS_PER_SWEEP, ION_PROBABILITY_WEIGHT;

    double TOTAL_ENERGY;
    int TEMP_INDEX;

    std::vector <Sampler *> SAMPLERS;
    int DATA_SAMPLING_RATE;


    WaterSystem(int num_waters = 200, int num_ions = 2);
    ~WaterSystem();

    void set_temperature(double new_temp);
    void expand_box_z_direction(double new_len = 0.0);

    void mc_move();
    void undo_mc_move();

    void add_rdf_sampler();
    void add_lammpstrj_sampler();
    void add_ion_pair_distance_sampler();

    void write_config_snapshot();
    void initialize_sampling();
    void sample_data();
    void finish_sampling();

    friend std::ostream & operator<<(std::ostream & out, WaterSystem * system);
};

#endif	/* SYSTEM_H */
