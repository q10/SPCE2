/* 
 * File:   system.h
 * Author: BENSON J MA
 *
 * Created on February 1, 2012, 8:59 PM
 */

#ifndef SYSTEM_H
#define	SYSTEM_H

class Water;
class Ion;

class System {
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

    System() {}
    ~System() {}
};

#endif	/* SYSTEM_H */
