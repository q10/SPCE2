/* 
 * File:   rdf_sampler.h
 * Author: BENSON J MA
 *
 * Created on November 4, 2011, 2:57 PM
 */

#ifndef RDF_SAMPLER_H
#define	RDF_SAMPLER_H

class RDFSampler : public Sampler {
private:
    double delg;
    int num_gr, radial_dist_num_his_bars;
    std::vector <double> radial_dist_distance, water_water_RDF, anion_water_RDF, cation_water_RDF, ion_ion_RDF;
    Simulation * simulation;

public:
    RDFSampler(Simulation * s);
    ~RDFSampler();

    void start();
    void sample();
    void finish();

    std::string radial_dist_results();
};


#endif	/* RDF_SAMPLER_H */
