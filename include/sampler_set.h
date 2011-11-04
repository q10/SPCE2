/* 
 * File:   sampler.h
 * Author: BENSON J MA
 *
 * Created on August 4, 2011, 1:07 AM
 */

#ifndef SAMPLER_H
#define	SAMPLER_H

class SamplerSet {
private:
    int vmd_timestep, vmd_snapshot_counter;
    std::ofstream VMD_FILE, CONFIG_FILE;

    std::vector <Sampler *> samplers;
public:
    int RELATIVE_VMD_SNAPSHOT_RATE;

    Simulation * simulation;
    std::string SIMULATION_TIME_STAMP;

    SamplerSet(Simulation * s);
    ~SamplerSet();

    void start();
    void sample_data();
    void finish();

    std::string config_filename, vmd_filename;
    void write_vmd_snapshot();
    void write_config_snapshot();

};

void test_radial_dist();

#endif	/* SAMPLER_H */
