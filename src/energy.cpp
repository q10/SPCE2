#include "common.h"

double SPCEHamiltonian::energy_diff_of_particle_with_index(int index) {
    double energy_diff = 0.0;
    int NUM_WATERS = WATERS.size(), NUM_IONS = IONS.size();
    if (index < NUM_WATERS) {
        for (int i = 0; i < NUM_WATERS; i++) {
            if (i != index)
                energy_diff += energy_between_two_waters(index, i, CURRENT) - energy_between_two_waters(index, i, OLD);
        }
        for (int i = 0; i < NUM_IONS; i++)
            energy_diff += energy_between_ion_and_water(i, index, CURRENT) - energy_between_ion_and_water(i, index, WATER);
    } else {
        index -= NUM_WATERS;
        for (int i = 0; i < NUM_IONS; i++) {
            if (i != index)
                energy_diff += energy_between_two_ions(index, i, CURRENT) - energy_between_two_ions(index, i, OLD);
        }
        for (int i = 0; i < NUM_WATERS; i++)
            energy_diff += energy_between_ion_and_water(index, i, CURRENT) - energy_between_ion_and_water(index, i, ION);
    }
    return energy_diff;
}

double SPCEHamiltonian::energy_of_particle_with_index(int index) {
    double energy = 0.0;
    int NUM_WATERS = WATERS.size(), NUM_IONS = IONS.size();
    if (index < NUM_WATERS) {
        for (int i = 0; i < NUM_WATERS; i++) {
            if (i != index)
                energy += energy_between_two_waters(index, i, CURRENT);
        }
        for (int i = 0; i < NUM_IONS; i++)
            energy += energy_between_ion_and_water(i, index, CURRENT);
    } else {
        index -= NUM_WATERS;
        for (int i = 0; i < NUM_IONS; i++) {
            if (i != index)
                energy += energy_between_two_ions(index, i, CURRENT);
        }
        for (int i = 0; i < NUM_WATERS; i++)
            energy += energy_between_ion_and_water(index, i, CURRENT);
    }
    return energy;
}

double SPCEHamiltonian::total_real_space_energy() {
    double real_space_energy = 0.0;
    for (unsigned int i = 0; i < WATERS.size() + IONS.size(); i++)
        real_space_energy += energy_of_particle_with_index(i);
    return real_space_energy / 2.0;
}

double SPCEHamiltonian::energy_between_ion_and_water(int i, int j, WHICH_TYPE typ) {
    double *ion_i = IONS[i]->coords, *water_j = WATERS[j]->coords;
    if (typ == WATER)
        water_j = WATERS[j]->old_coords;
    else if (typ == ION)
        ion_i = IONS[i]->old_coords;

    double dx, dy, dz, old_dx, old_dy, old_dz, r, r2, tmp_energy = 0.0;
    bool use_same_x, use_same_y, use_same_z;

    for (int atom = 0; atom < 9; atom += 3) {
        dx = old_dx = abs(ion_i[0] - water_j[atom]);
        dy = old_dy = abs(ion_i[1] - water_j[atom + 1]);
        dz = old_dz = abs(ion_i[2] - water_j[atom + 2]);

        // Fix the distances to use the same nearest mirror image for each molecule based on O-O distance
        if (atom == 0) {
            dx -= BOX_LENGTH * ROUND(dx / BOX_LENGTH);
            dy -= BOX_LENGTH * ROUND(dy / BOX_LENGTH);
            dz -= BOX_Z_LENGTH * ROUND(dz / BOX_Z_LENGTH);
            if (dx == old_dx)
                use_same_x = true;
            if (dy == old_dy)
                use_same_y = true;
            if (dz == old_dz)
                use_same_z = true;
        } else {
            if (!use_same_x)
                dx = BOX_LENGTH - old_dx;
            if (!use_same_y)
                dy = BOX_LENGTH - old_dy;
            if (!use_same_z)
                dz = BOX_Z_LENGTH - old_dz;
        }

        r = sqrt(dx * dx + dy * dy + dz * dz);

        if (atom == 0) {
            r2 = (Water::SIGMA + Ion::SIGMA) / (2 * r);
            tmp_energy += 4.0 * Water::EPSILON * (pow(r2, 12) - pow(r2, 6)) +
                    PCONSTANTS::ELECTROSTATIC_K * IONS[i]->charge * Water::Q_O * ERFC_TABLE[int(r * 1000.0)] / r;
        } else
            tmp_energy += PCONSTANTS::ELECTROSTATIC_K * IONS[i]->charge * Water::Q_H * ERFC_TABLE[int(r * 1000.0)] / r;
    }
    ion_i = water_j = NULL;
    return tmp_energy;
}

double SPCEHamiltonian::energy_between_two_ions(int i, int j, WHICH_TYPE typ) {
    double r = (typ == OLD) ? IONS[i]->old_distance_from(IONS[j]) : IONS[i]->distance_from(IONS[j]);
    double rb = Ion::SIGMA / r;
    return 4.0 * Ion::EPSILON * (pow(rb, 12) - pow(rb, 6)) +
            PCONSTANTS::ELECTROSTATIC_K * IONS[i]->charge * IONS[j]->charge * ERFC_TABLE[int(r * 1000.0)] / r;
}

double SPCEHamiltonian::energy_between_two_waters(int i, int j, WHICH_TYPE typ) {
    double *water_i = (typ == OLD) ? WATERS[i]->old_coords : WATERS[i]->coords;
    double dx, dy, dz, old_dx, old_dy, old_dz, r, r2, *water_j = WATERS[j]->coords, tmp_energy = 0.0;
    bool use_same_x, use_same_y, use_same_z;

    use_same_x = use_same_y = use_same_z = false;

    // Loop over atoms within the two waters
    for (int atom = 0; atom < 9; atom += 3) {
        for (int atom2 = 0; atom2 < 9; atom2 += 3) {
            dx = old_dx = abs(water_i[atom] - water_j[atom2]);
            dy = old_dy = abs(water_i[atom + 1] - water_j[atom2 + 1]);
            dz = old_dz = abs(water_i[atom + 2] - water_j[atom2 + 2]);

            // Fix the distances to use the same nearest mirror image for each molecule based on O-O distance
            if (atom == 0 and atom2 == 0) {
                dx -= BOX_LENGTH * ROUND(dx / BOX_LENGTH);
                dy -= BOX_LENGTH * ROUND(dy / BOX_LENGTH);
                dz -= BOX_Z_LENGTH * ROUND(dz / BOX_Z_LENGTH);
                if (dx == old_dx)
                    use_same_x = true;
                if (dy == old_dy)
                    use_same_y = true;
                if (dz == old_dz)
                    use_same_z = true;
            } else {
                if (!use_same_x)
                    dx = BOX_LENGTH - old_dx;
                if (!use_same_y)
                    dy = BOX_LENGTH - old_dy;
                if (!use_same_z)
                    dz = BOX_Z_LENGTH - old_dz;
            }

            r = sqrt(dx * dx + dy * dy + dz * dz);

            // Covers the O-O case, the H-O/O-H case, and the H-H case for electrostatics
            if (atom == 0 and atom2 == 0) {
                r2 = Water::SIGMA / r;
                tmp_energy += 4.0 * Water::EPSILON * (pow(r2, 12) - pow(r2, 6)) +
                        PCONSTANTS::ELECTROSTATIC_K * Water::Q_O * Water::Q_O * ERFC_TABLE[int(r * 1000.0)] / r;
            } else if (atom == 0 or atom2 == 0)
                tmp_energy += PCONSTANTS::ELECTROSTATIC_K * Water::Q_O * Water::Q_H * ERFC_TABLE[int(r * 1000.0)] / r;
            else
                tmp_energy += PCONSTANTS::ELECTROSTATIC_K * Water::Q_H * Water::Q_H * ERFC_TABLE[int(r * 1000.0)] / r;
        }
    }
    water_i = water_j = NULL;
    return tmp_energy;
}
