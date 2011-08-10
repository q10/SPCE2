/* 
 * File:   water.h
 * Author: BENSON J MA
 *
 * Created on August 3, 2011, 4:48 PM
 */

#ifndef WATER_H
#define	WATER_H

class Ion;

class Water {
private:
    void set_center_of_mass_of_water();
    void set_rotation_matrix(double * rand_unit_vector, double theta_rad);

public:
    static const double SIGMA,
    EPSILON,
    STD_DENSITY,
    Q_H,
    Q_O,
    H_MASS,
    O_MASS,
    WATER_MASS,
    OH_LENGTH,
    HOH_ANGLE_DEG;

    double *TMP_CENTER_OF_MASS, **ROTATION_MATRIX;
    double *coords, *old_coords, DISPLACEMENT_DISTANCE, DISPLACEMENT_ROTATION, BOX_LENGTH, BOX_Z_LENGTH;

    Water(double * tmp_coords, double tmp_disp_dist, double tmp_disp_rot, double box_length, double box_z_length);
    ~Water();

    // Must be friend to access private members.
    friend std::ostream & operator<<(std::ostream & out, Water * water);

    void set_coords(double * tmp_coords);

    void mc_translate();
    void mc_rotate();
    void keep_inside_box();
    void undo_move();

    double distance_from(Water * other_water);
    double distance_from(Ion * other_ion);
};


void test_water_rotation();

#endif	/* WATER_H */
