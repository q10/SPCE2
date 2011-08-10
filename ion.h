/* 
 * File:   ion.h
 * Author: BENSON J MA
 *
 * Created on August 3, 2011, 4:57 PM
 */

#ifndef ION_H
#define	ION_H

class Ion {
public:
    static const double SIGMA, EPSILON;
    double *coords, *old_coords;
    double charge, DISPLACEMENT_DISTANCE, BOX_LENGTH, BOX_Z_LENGTH;

    Ion(double * tmp_coords, double tmp_charge, double tmp_disp_dist, double box_length, double box_z_length);
    ~Ion();

    // Must be friend to access private members.
    friend std::ostream & operator<<(std::ostream & out, Ion * ion);

    void set_coords(double * tmp_coords);
    
    void mc_translate();
    void keep_inside_box();
    void undo_move();

    double distance_from(Water * other_water);
    double distance_from(Ion * other_ion);
};

#endif	/* ION_H */
