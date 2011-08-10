/* 
 * File:   sampler.h
 * Author: BENSON J MA
 *
 * Created on August 4, 2011, 1:07 AM
 */

#ifndef SAMPLER_H
#define	SAMPLER_H

class Sampler {
private:
    int num_gr, radial_dist_num_his_bars;
    double *radial_dist_distance, *water_water_RDF, *ion_water_RDF, *ion_ion_RDF, delg;

    void initialize_radial_dist_sampler();
    void radial_dist_sample();
    void compute_radial_dist_results();

    int vmd_timestep;
    std::ofstream VMD_FILE, CONFIG_FILE;
    void close_vmd_file();

public:
    int DATA_SAMPLING_RATE;

    Simulation * simulation;
    std::string SIMULATION_TIME_STAMP;

    Sampler(Simulation * s);
    ~Sampler();

    void sample_data();
    void start();
    void finish();


    std::string config_filename, vmd_filename;
    void write_vmd_snapshot();
    void write_config_snapshot();

    std::string radial_dist_results();
};

void test_radial_dist();

#endif	/* SAMPLER_H */
