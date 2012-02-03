#include "common.h"

double Ion::SIGMA = 5.0,
        Ion::EPSILON = 1.0;

ostream & operator<<(ostream & out, Ion * ion) {
    return out << setprecision(10) << ion->coords[0] << "\t" << ion->coords[1] << "\t" << ion->coords[2] << "\t" << ion->charge;
}

Ion::Ion(System * sys, double * tmp_coords, double tmp_charge)
: DISPLACEMENT_DISTANCE(sys->DISPLACEMENT_DISTANCE), BOX_LENGTH(sys->BOX_LENGTH), BOX_Z_LENGTH(sys->BOX_Z_LENGTH) {
    charge = tmp_charge;
    old_coords = new double [3];
    coords = new double [3];
    set_coords(tmp_coords);
}

Ion::~Ion() {
    delete [] old_coords;
    delete [] coords;
}

void Ion::set_coords(double * tmp_coords) {
    for (int i = 0; i < 3; i++)
        coords[i] = tmp_coords[i];
    return;
}

void Ion::set_random_coords() {
    coords[0] = RAN3() * BOX_LENGTH;
    coords[1] = RAN3() * BOX_LENGTH;
    coords[2] = RAN3() * BOX_Z_LENGTH;
    return;
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
    return sqrt(squared_distance_from(other_water));
}

double Ion::distance_from(Ion * other_ion) {
    return sqrt(squared_distance_from(other_ion));
}

double Ion::old_distance_from(Water * other_water) {
    double * other_coords = other_water->old_coords;
    double dx = coords[0] - other_coords[0];
    double dy = coords[1] - other_coords[1];
    double dz = coords[2] - other_coords[2];
    dx -= BOX_LENGTH * ROUND(dx / BOX_LENGTH);
    dy -= BOX_LENGTH * ROUND(dy / BOX_LENGTH);
    dz -= BOX_Z_LENGTH * ROUND(dz / BOX_Z_LENGTH);
    return dx * dx + dy * dy + dz * dz;
}

double Ion::old_distance_from(Ion * other_ion) {
    double * other_coords = other_ion->old_coords;
    double dx = coords[0] - other_coords[0];
    double dy = coords[1] - other_coords[1];
    double dz = coords[2] - other_coords[2];
    dx -= BOX_LENGTH * ROUND(dx / BOX_LENGTH);
    dy -= BOX_LENGTH * ROUND(dy / BOX_LENGTH);
    dz -= BOX_Z_LENGTH * ROUND(dz / BOX_Z_LENGTH);
    return dx * dx + dy * dy + dz * dz;
}

double Ion::squared_distance_from(Water * other_water) {
    double * other_coords = other_water->coords;
    double dx = coords[0] - other_coords[0];
    double dy = coords[1] - other_coords[1];
    double dz = coords[2] - other_coords[2];
    dx -= BOX_LENGTH * ROUND(dx / BOX_LENGTH);
    dy -= BOX_LENGTH * ROUND(dy / BOX_LENGTH);
    dz -= BOX_Z_LENGTH * ROUND(dz / BOX_Z_LENGTH);
    return dx * dx + dy * dy + dz * dz;
}

double Ion::squared_distance_from(Ion * other_ion) {
    double * other_coords = other_ion->coords;
    double dx = coords[0] - other_coords[0];
    double dy = coords[1] - other_coords[1];
    double dz = coords[2] - other_coords[2];
    dx -= BOX_LENGTH * ROUND(dx / BOX_LENGTH);
    dy -= BOX_LENGTH * ROUND(dy / BOX_LENGTH);
    dz -= BOX_Z_LENGTH * ROUND(dz / BOX_Z_LENGTH);
    return dx * dx + dy * dy + dz * dz;
}
