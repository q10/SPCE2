/* 
 * File:   lammpstrj_sampler.h
 * Author: benson
 *
 * Created on April 20, 2012, 12:22 PM
 */

#ifndef LAMMPSTRJ_SAMPLER_H
#define	LAMMPSTRJ_SAMPLER_H

class LAMMPSTRJSampler : public Sampler {
private:
    int timestep, lammpstrj_snapshot_counter;
    std::ofstream LAMMPSTRJ_FILE;

public:
    LAMMPSTRJSampler(WaterSystem * s);
    ~LAMMPSTRJSampler();

    int RELATIVE_LAMMPSTRJ_SNAPSHOT_RATE;

    void start();
    void sample();
    void finish();
    std::string results();
};

#endif	/* LAMMPSTRJ_SAMPLER_H */
