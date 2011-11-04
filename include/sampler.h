/* 
 * File:   sampler.h
 * Author: BENSON J MA
 *
 * Created on November 4, 2011, 2:58 PM
 */

#ifndef SAMPLER_H
#define	SAMPLER_H

class Sampler {
public:
    Sampler() {};
    virtual ~Sampler() {};
    virtual void start() = 0;
    virtual void sample() = 0;
    virtual void finish() = 0;
};

#endif	/* SAMPLER_H */
