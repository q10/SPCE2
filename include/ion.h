/* 
 * File:   ion.h
 * Author: BENSON J MA
 *
 * Created on August 3, 2011, 4:57 PM
 */

#ifndef ION_H
#define	ION_H

class System;

class Ion {
public:
    static double SIGMA, EPSILON;
    double *coords, *old_coords;
    double charge, &DISPLACEMENT_DISTANCE, &BOX_LENGTH, &BOX_Z_LENGTH;

    Ion(System * sys, double * tmp_coords, double tmp_charge);
    ~Ion();

    // Must be friend to access private members.
    friend std::ostream & operator<<(std::ostream & out, Ion * ion);

    void set_coords(double * tmp_coords);
    void set_random_coords();

    void mc_translate();
    void keep_inside_box();
    void undo_move();

    double distance_from(Water * other_water);
    double distance_from(Ion * other_ion);
    double old_distance_from(Water * other_water);
    double old_distance_from(Ion * other_ion);
    double squared_distance_from(Water * other_water);
    double squared_distance_from(Ion * other_ion);
};

#endif	/* ION_H */
