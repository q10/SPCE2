/* 
 * File:   sampler.h
 * Author: BENSON J MA
 *
 * Created on August 4, 2011, 1:07 AM
 */

#ifndef SAMPLER_SET_H
#define	SAMPLER_SET_H

class SamplerSet {
private:
    bool is_lammpstr_sampling;
    int lammpstrj_timestep, lammpstrj_snapshot_counter;
    std::ofstream LAMMPSTRJ_FILE;
    std::vector <Sampler *> samplers;

public:
    int RELATIVE_LAMMPSTRJ_SNAPSHOT_RATE;
    Simulation * simulation;

    SamplerSet(Simulation * s);
    ~SamplerSet();

    void start();
    void sample_data();
    void finish();
    void print_individual_sampler_results();
    
    void turn_on_lammpstrj_sampler();
    void add_rdf_sampler();
    void add_ion_pair_distance_sampler();

    std::string config_filename, lammpstrj_filename;
    void write_lammpstrj_snapshot();
    void write_config_snapshot();
};

#endif	/* SAMPLER_SET_H */
