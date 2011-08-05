#include "common.h"

const double Ion::SIGMA = 5.0,
        Ion::EPSILON = 1.0;

ostream & operator<<(ostream & out, Ion * ion) {
    return out << setprecision(10) << ion->coords[0] << "\t" << ion->coords[1] << "\t" << ion->coords[2] << "\t" << ion->charge;
}

Ion::Ion(double * tmp_coords, double tmp_charge, double tmp_disp_dist, double box_length, double box_z_length) {
    charge = tmp_charge;
    DISPLACEMENT_DISTANCE = tmp_disp_dist;
    BOX_LENGTH = box_length;
    BOX_Z_LENGTH = box_z_length;
    for (int i = 0; i < 3; i++)
        coords[i] = tmp_coords[i];
}

Ion::~Ion() {
}

void Ion::mc_translate() {
    for (int g = 0; g < 3; g++) {
        old_coords[g] = coords[g];
        coords[g] += DISPLACEMENT_DISTANCE * (2.0 * RAN3() - 1.0);
    }
    keep_inside_box();
    return;
}

void Ion::keep_inside_box() {
    for (int j = 0; j < 2; j++) {
        if (coords[j] > BOX_LENGTH)
            coords[j] -= BOX_LENGTH;
        if (coords[j] < 0.0)
            coords[j] += BOX_LENGTH;
    }
    if (coords[2] > BOX_Z_LENGTH)
        coords[2] -= BOX_Z_LENGTH;
    if (coords[2] < 0.0)
        coords[2] += BOX_Z_LENGTH;
    return;
}

void Ion::undo_move() {
    for (int j = 0; j < 3; j++)
        coords[j] = old_coords[j];
    return;
}

double Ion::distance_from(Water * other_water) {
    double * other_coords = other_water->coords;
    double dx = coords[0] - other_coords[0];
    double dy = coords[1] - other_coords[1];
    double dz = coords[2] - other_coords[2];
    dx -= BOX_LENGTH * ROUND(dx / BOX_LENGTH);
    dy -= BOX_LENGTH * ROUND(dy / BOX_LENGTH);
    dz -= BOX_Z_LENGTH * ROUND(dz / BOX_Z_LENGTH);
    return sqrt(dx * dx + dy * dy + dz * dz);
}

double Ion::distance_from(Ion * other_ion) {
    double * other_coords = other_ion->coords;
    double dx = coords[0] - other_coords[0];
    double dy = coords[1] - other_coords[1];
    double dz = coords[2] - other_coords[2];
    dx -= BOX_LENGTH * ROUND(dx / BOX_LENGTH);
    dy -= BOX_LENGTH * ROUND(dy / BOX_LENGTH);
    dz -= BOX_Z_LENGTH * ROUND(dz / BOX_Z_LENGTH);
    return sqrt(dx * dx + dy * dy + dz * dz);
}
