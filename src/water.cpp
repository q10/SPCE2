#include "common.h"

const double Water::SIGMA = 3.166,
        Water::EPSILON = 0.650,
        Water::STD_DENSITY = 0.0334225755,
        Water::Q_H = 0.4238,
        Water::Q_O = -0.8476,
        Water::SIGMA_H = 2.4,
        Water::SIGMA_O = 3.04,
        Water::H_MASS = 1.0,
        Water::O_MASS = 16.0,
        Water::WATER_MASS = 18.0,
        Water::OH_LENGTH = 1.0000,
        Water::HOH_ANGLE_DEG = 109.47;

ostream & operator<<(ostream & out, Water * water) {
    out << setprecision(10) << water->coords[0];
    for (int i = 1; i < 9; i++)
        out << "\t" << water->coords[i];
    return out;
}

Water::Water(double * tmp_coords, double tmp_disp_dist, double tmp_disp_rot, double box_length, double box_z_length) {
    DISPLACEMENT_DISTANCE = tmp_disp_dist;
    DISPLACEMENT_ROTATION = tmp_disp_rot;
    BOX_LENGTH = box_length;
    BOX_Z_LENGTH = box_z_length;

    TMP_CENTER_OF_MASS = new double [3];
    ROTATION_MATRIX = new double * [3];
    for (int i = 0; i < 3; i++)
        ROTATION_MATRIX[i] = new double [3];
    old_coords = new double [9];
    coords = new double [9];
    set_coords(tmp_coords);
}

Water::~Water() {
    delete [] coords;
    delete [] old_coords;
    delete [] TMP_CENTER_OF_MASS;
    for (int i = 0; i < 3; i++)
        delete [] ROTATION_MATRIX[i];
    delete [] ROTATION_MATRIX;
}

void Water::set_coords(double * tmp_coords) {
    for (int i = 0; i < 9; i++)
        coords[i] = tmp_coords[i];
    return;
}

void Water::mc_translate() {
    double rand_displacement;
    for (int g = 0; g < 9; g++)
        old_coords[g] = coords[g];

    for (int j = 0; j < 3; j++) {
        rand_displacement = DISPLACEMENT_DISTANCE * (2.0 * RAN3() - 1.0);
        coords[j] += rand_displacement;
        coords[j + 3] += rand_displacement;
        coords[j + 6] += rand_displacement;
    }
    keep_inside_box();
    return;
}

void Water::mc_rotate() {
    set_center_of_mass_of_water();
    set_rotation_matrix(RANDUNITVECTOR(), DISPLACEMENT_ROTATION * (2.0 * RAN3() - 1.0));

    // save old position 
    // shift water such that its center of mass is now the origin (use the old_coords set of coords)
    for (int g = 0; g < 9; g++) {
        old_coords[g] = coords[g];
        old_coords[g] -= TMP_CENTER_OF_MASS[g % 3];
    }

    // apply rotation matrix to all 9 coordinates of water
    for (int g = 0; g < 9; g += 3) {
        coords[g] = ROTATION_MATRIX[0][0] * old_coords[g] + ROTATION_MATRIX[0][1] * old_coords[g + 1] + ROTATION_MATRIX[0][2] * old_coords[g + 2];
        coords[g + 1] = ROTATION_MATRIX[1][0] * old_coords[g] + ROTATION_MATRIX[1][1] * old_coords[g + 1] + ROTATION_MATRIX[1][2] * old_coords[g + 2];
        coords[g + 2] = ROTATION_MATRIX[2][0] * old_coords[g] + ROTATION_MATRIX[2][1] * old_coords[g + 1] + ROTATION_MATRIX[2][2] * old_coords[g + 2];
    }

    // un-shift water (use the coords), and restore old_coords to proper original coords
    for (int g = 0; g < 9; g++) {
        coords[g] += TMP_CENTER_OF_MASS[g % 3];
        old_coords[g] += TMP_CENTER_OF_MASS[g % 3];
    }

    // make sure the oxygen stays inside boundaries; otherwise shift appropriately
    keep_inside_box();
    return;
}

void Water::keep_inside_box() {
    for (int j = 0; j < 2; j++) {
        if (coords[j] > BOX_LENGTH) {
            coords[j] -= BOX_LENGTH;
            coords[j + 3] -= BOX_LENGTH;
            coords[j + 6] -= BOX_LENGTH;
        }
        if (coords[j] < 0.0) {
            coords[j] += BOX_LENGTH;
            coords[j + 3] += BOX_LENGTH;
            coords[j + 6] += BOX_LENGTH;
        }
    }
    if (coords[2] > BOX_Z_LENGTH) {
        coords[2] -= BOX_Z_LENGTH;
        coords[5] -= BOX_Z_LENGTH;
        coords[8] -= BOX_Z_LENGTH;
    }
    if (coords[2] < 0.0) {
        coords[2] += BOX_Z_LENGTH;
        coords[5] += BOX_Z_LENGTH;
        coords[8] += BOX_Z_LENGTH;
    }
    return;
}

void Water::undo_move() {
    for (int j = 0; j < 9; j++)
        coords[j] = old_coords[j];
    return;
}

double Water::distance_from(Water * other_water) {
    return sqrt(squared_distance_from(other_water));
}

double Water::distance_from(Ion * other_ion) {
    return sqrt(squared_distance_from(other_ion));
}

double Water::squared_distance_from(Water * other_water) {
    double * other_coords = other_water->coords;
    double dx = coords[0] - other_coords[0];
    double dy = coords[1] - other_coords[1];
    double dz = coords[2] - other_coords[2];
    dx -= BOX_LENGTH * ROUND(dx / BOX_LENGTH);
    dy -= BOX_LENGTH * ROUND(dy / BOX_LENGTH);
    dz -= BOX_Z_LENGTH * ROUND(dz / BOX_Z_LENGTH);
    return dx * dx + dy * dy + dz * dz;
}

double Water::squared_distance_from(Ion * other_ion) {
    double * other_coords = other_ion->coords;
    double dx = coords[0] - other_coords[0];
    double dy = coords[1] - other_coords[1];
    double dz = coords[2] - other_coords[2];
    dx -= BOX_LENGTH * ROUND(dx / BOX_LENGTH);
    dy -= BOX_LENGTH * ROUND(dy / BOX_LENGTH);
    dz -= BOX_Z_LENGTH * ROUND(dz / BOX_Z_LENGTH);
    return dx * dx + dy * dy + dz * dz;
}
