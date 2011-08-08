/* 
 * File:   sampler.h
 * Author: BENSON J MA
 *
 * Created on August 4, 2011, 1:07 AM
 */

#ifndef SAMPLER_H
#define	SAMPLER_H

class Sampler {
public:
    static const int DATA_SAMPLING_RATE = 5;

    Simulation * simulation;

    Sampler(Simulation * s);
    ~Sampler();

    void sample_data();

    
    int num_gr, radial_dist_num_his_bars;
    double *radial_dist_distance, *water_water_RDF, *ion_water_RDF, *ion_ion_RDF, delg;

    void initialize_radial_dist_sampler();
    void radial_dist_sample();
    void compute_radial_dist_results();
    void print_radial_dist_results();
};

void test_radial_dist_sampler();

#endif	/* SAMPLER_H */
