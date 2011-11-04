/* 
 * File:   ipair_distance_sampler.h
 * Author: BENSON J MA
 *
 * Created on November 4, 2011, 3:29 PM
 */

#ifndef IPAIR_DISTANCE_SAMPLER_H
#define	IPAIR_DISTANCE_SAMPLER_H

class IonPairDistanceSampler : public Sampler {
private:
    void initialize_time_series_sampler();
    void time_series_sample();
    void finish_time_series_sampler();

public:
    IonPairDistanceSampler();
    ~IonPairDistanceSampler();

    void start();
    void sample();
    void finish();
};

#endif	/* IPAIR_DISTANCE_SAMPLER_H */
