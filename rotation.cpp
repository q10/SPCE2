#include "common.h"

void Water::set_center_of_mass_of_water() {
    for (int j = 0; j < 3; j++)
        TMP_CENTER_OF_MASS[j] = (O_MASS * coords[j] + H_MASS * coords[j + 3] + H_MASS * coords[j + 6]) / WATER_MASS;
    return;
}

void Water::set_rotation_matrix(double * rand_unit_vector, double theta_rad) {
    double rx = rand_unit_vector[0], ry = rand_unit_vector[1], rz = rand_unit_vector[2], sin_t = sin(theta_rad);
    double cos_t = sqrt(1.0 - sin_t * sin_t);
    double l_cos_t = 1.0 - cos_t, rx2 = rx*rx, ry2 = ry*ry, rz2 = rz*rz, rxry = rx*ry, ryrz = ry*rz, rxrz = rx*rz;

    ROTATION_MATRIX[0][0] = cos_t + rx2 * l_cos_t;
    ROTATION_MATRIX[0][1] = rxry * l_cos_t - rz * sin_t;
    ROTATION_MATRIX[0][2] = rxrz * l_cos_t + ry * sin_t;

    ROTATION_MATRIX[1][0] = rxry * l_cos_t + rz * sin_t;
    ROTATION_MATRIX[1][1] = cos_t + ry2 * l_cos_t;
    ROTATION_MATRIX[1][2] = ryrz * l_cos_t - rx * sin_t;

    ROTATION_MATRIX[2][0] = rxrz * l_cos_t - ry * sin_t;
    ROTATION_MATRIX[2][1] = ryrz * l_cos_t + rx * sin_t;
    ROTATION_MATRIX[2][2] = cos_t + rz2 * l_cos_t;

    delete [] rand_unit_vector;
    return;
}

void test_water_rotation() {
    cout << "---- BEGIN TEST - ROTATION ----" << endl;

    double * coords = new double [9];
    coords[0] = 1.34;
    coords[1] = 2.0;
    coords[2] = 3.2;
    for (int i = 3; i < 9; i++)
        coords[i] = RAN3()*5.0;

    Water * w = new Water(coords, 0.2, 0.17, 20.0, 20.0);
    delete [] coords;
    for (int k = 0; k < 10000; k++) {
        w->mc_rotate();
        cout << w << endl;
    }

    cout << "---- END TEST - ROTATION ----\n" << endl;
    return;
}
