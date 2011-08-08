#include "common.h"

void Simulation::initialize_all_ewald_tables(double ewald_alpha, int ewald_nxy, int ewald_nz) {
    ASSERT(ewald_alpha > 0.0 and ewald_nxy > 0 and ewald_nz > 0, "INVALID EWALD PARAMETERS ALPHA, NXY, OR NZ");
    EWALD_ALPHA = ewald_alpha;
    EWALD_NXY = ewald_nxy;
    EWALD_NZ = ewald_nz;
    K_111_INDEX = (3 * EWALD_NXY + 2) * (2 * EWALD_NZ + 1) + EWALD_NZ;
    NUM_K_VECTORS = (EWALD_NXY + 1) * (2 * EWALD_NXY + 1) * (2 * EWALD_NZ + 1) - 1;

    initialize_erfc_table();
    initialize_k_vectors_table();
    initialize_rho_k_values_table();
    initialize_exp_kr_cache_tables();
    return;
}

void Simulation::initialize_erfc_table() {
    ERFC_TABLE.clear();
    double sqrt_alpha_1000 = sqrt(EWALD_ALPHA) / 1000.0, key;
    for (double r = 0.0; r < max(HALF_BOX_LENGTH, HALF_BOX_Z_LENGTH); r += 0.001) {
        // to optimize map search by using one less division, we keep the keys as x10^4
        key = floor(r * 1000.0);
        ERFC_TABLE[key] = erfc(key * sqrt_alpha_1000);
    }
    return;
}

void Simulation::initialize_k_vectors_table() {
    K_VECTORS.clear();

    // K vector order as follows (go through all ny and nzvalues descending, then ascending):    
    double *arr, tmp_k2, half_factor, alpha_inv_4 = -1.0 / (4.0 * EWALD_ALPHA), four_pi_volume_ek = 4.0 * ELECTROSTATIC_K * M_PI / BOX_VOLUME;

    for (int nx = 0; nx <= EWALD_NXY; nx++) {
        for (int ny = -EWALD_NXY; ny <= EWALD_NXY; ny++) {
            for (int nz = -EWALD_NZ; nz <= EWALD_NZ; nz++) {
                if (nx != 0 || ny != 0 || nz != 0) { // removes the K=0 case
                    arr = new double[4];
                    arr[0] = 2.0 * M_PI * nx / BOX_LENGTH;
                    arr[1] = 2.0 * M_PI * ny / BOX_LENGTH;
                    arr[2] = 2.0 * M_PI * nz / BOX_Z_LENGTH;

                    // k_entry[3] will be the Fourier coefficient exp(-K^2 / (4*alpha)) / K^2, where K^2, or kx^2 + ky^2 + kz^2
                    tmp_k2 = arr[0] * arr[0] + arr[1] * arr[1] + arr[2] * arr[2];
                    arr[3] = exp(tmp_k2 * alpha_inv_4) / tmp_k2;

                    // placing the 4*pi/V here allows for smaller number of multiplications and divisions later on
                    // half_factor accounts for double weighting of Fourier coefficients along the nx=0 plane
                    half_factor = (nx == 0) ? 0.5 : 1.0;
                    arr[3] *= half_factor * four_pi_volume_ek;
                    K_VECTORS.push_back(arr);
                }
            }
        }
    }
    ASSERT(NUM_K_VECTORS == K_VECTORS.size(), "NUMBER OF K VECTORS NOT MATCHING.");
    return;
}

void Simulation::initialize_rho_k_values_table() {
    RHO_K_VALUES.clear();
    dcomplex *column, column_sum;
    int NUM_TOTAL_PARTICLES = WATERS.size() + IONS.size();

    for (int k = 0; k < K_VECTORS.size(); k++) {
        column = new dcomplex [NUM_TOTAL_PARTICLES + 2];
        column_sum = dcomplex(0.0, 0.0);

        for (int w = 0; w < NUM_TOTAL_PARTICLES; w++) {
            column[w] = partial_rho(w, K_VECTORS[k]);
            column_sum += column[w];
        }

        column[NUM_TOTAL_PARTICLES] = column_sum;
        RHO_K_VALUES.push_back(column);
    }
    return;
}

void Simulation::initialize_exp_kr_cache_tables() {
    for (int i = 0; i < 3; i++) {
        if (exp_kr_O[i] != NULL)
            delete [] exp_kr_O[i];
        if (exp_kr_H1[i] != NULL)
            delete [] exp_kr_H1[i];
        if (exp_kr_H2[i] != NULL)
            delete [] exp_kr_H2[i];
        if (exp_kr_ion[i] != NULL)
            delete [] exp_kr_ion[i];
    }

    int tmp_n;
    for (int i = 0; i < 3; i++) {
        tmp_n = (i < 2) ? EWALD_NXY : EWALD_NZ;
        exp_kr_O[i] = new dcomplex [2 * tmp_n + 1];
        exp_kr_H1[i] = new dcomplex [2 * tmp_n + 1];
        exp_kr_H2[i] = new dcomplex [2 * tmp_n + 1];
        exp_kr_ion[i] = new dcomplex [2 * tmp_n + 1];
    }

    for (int i = 0; i < 3; i++) {
        tmp_n = (i < 2) ? EWALD_NXY : EWALD_NZ;
        exp_kr_O[i][tmp_n] = exp_kr_H1[i][tmp_n] = exp_kr_H2[i][tmp_n] = exp_kr_ion[i][tmp_n] = COMPLEX_ONE;
    }

    return;
}

dcomplex Simulation::partial_rho(int index, double * k_coords) {
    double q, *coords;
    dcomplex part_rho(0.0, 0.0);

    if (index < WATERS.size()) {
        for (int atom = 0; atom < 9; atom += 3) {
            coords = WATERS[index]->coords;
            q = (atom == 0) ? Water::Q_O : Water::Q_H;
            part_rho += q * exp(dcomplex(0.0, k_coords[0] * coords[atom] + k_coords[1] * coords[atom + 1] + k_coords[2] * coords[atom + 2]));
        }
    } else {
        index -= WATERS.size();
        coords = IONS[index]->coords;
        part_rho = IONS[index]->charge * exp(dcomplex(0.0, k_coords[0] * coords[0] + k_coords[1] * coords[1] + k_coords[2] * coords[2]));
    }
    return part_rho;
}

double Simulation::total_ewald_energy() {
    double ewald_energy = 0.0;
    int NUM_TOTAL_PARTICLES = WATERS.size() + IONS.size();
    for (int k = 0; k < 725; k++)
        ewald_energy += norm(RHO_K_VALUES[k][NUM_TOTAL_PARTICLES]) * K_VECTORS[k][3];
    return ewald_energy;
}

double Simulation::ewald_diff(int index) {
    return (index < WATERS.size()) ? ewald_diff_water(index) : ewald_diff_ion(index);
}

double Simulation::ewald_diff_water(int water_index) {
    double sum_of_ewald_diffs = 0.0, old_pk2;
    dcomplex *column, tmp_x_O, tmp_y_O, tmp_x_H1, tmp_y_H1, tmp_x_H2, tmp_y_H2;
    int k = 0, NUM_TOTAL_PARTICLES = WATERS.size() + IONS.size();
    set_exp_kr_table_for_water(water_index);

    for (int nx = 0; nx <= EWALD_NXY; nx++) {
        tmp_x_O = Water::Q_O * exp_kr_O[0][nx];
        tmp_x_H1 = Water::Q_H * exp_kr_H1[0][nx];
        tmp_x_H2 = Water::Q_H * exp_kr_H2[0][nx];

        for (int ny = 0; ny <= 2 * EWALD_NXY; ny++) {
            tmp_y_O = tmp_x_O * exp_kr_O[1][ny];
            tmp_y_H1 = tmp_x_H1 * exp_kr_H1[1][ny];
            tmp_y_H2 = tmp_x_H2 * exp_kr_H2[1][ny];

            for (int nz = 0; nz <= 2 * EWALD_NZ; nz++) {
                if (nx != 0 || ny != EWALD_NXY || nz != EWALD_NZ) {
                    column = RHO_K_VALUES[k];
                    old_pk2 = norm(column[NUM_TOTAL_PARTICLES]);

                    // save old rho(K, R)
                    // calculate and save new rho(K, R)
                    // update total rho to rho(K, R)_new - rho(K, R)_old
                    column[NUM_TOTAL_PARTICLES + 1] = column[water_index];
                    column[water_index] = tmp_y_O * exp_kr_O[2][nz] + tmp_y_H1 * exp_kr_H1[2][nz] + tmp_y_H2 * exp_kr_H2[2][nz];
                    column[NUM_TOTAL_PARTICLES] += column[water_index] - column[NUM_TOTAL_PARTICLES + 1];

                    sum_of_ewald_diffs += (norm(column[NUM_TOTAL_PARTICLES]) - old_pk2) * K_VECTORS[k][3];
                    k++;
                }
            }
        }
    }
    return sum_of_ewald_diffs;
}

double Simulation::ewald_diff_ion(int index) {
    int k = 0, ion_index = index - WATERS.size(), NUM_TOTAL_PARTICLES = WATERS.size() + IONS.size();
    double sum_of_ewald_diffs = 0.0, old_pk2;
    dcomplex *column, tmp_x_ion, tmp_y_ion;
    set_exp_kr_table_for_ion(ion_index);

    for (int nx = 0; nx <= EWALD_NXY; nx++) {
        tmp_x_ion = IONS[ion_index]->charge * exp_kr_ion[0][nx];

        for (int ny = 0; ny <= 2 * EWALD_NXY; ny++) {
            tmp_y_ion = tmp_x_ion * exp_kr_ion[1][ny];

            for (int nz = 0; nz <= 2 * EWALD_NZ; nz++) {
                if (nx != 0 || ny != EWALD_NXY || nz != EWALD_NZ) {
                    column = RHO_K_VALUES[k];
                    old_pk2 = norm(column[NUM_TOTAL_PARTICLES]);

                    // save old rho(K, R)
                    // calculate and save new rho(K, R)
                    // update total rho to rho(K, R)_new - rho(K, R)_old
                    column[NUM_TOTAL_PARTICLES + 1] = column[index];
                    column[index] = tmp_y_ion * exp_kr_ion[2][nz];
                    column[NUM_TOTAL_PARTICLES] += column[index] - column[NUM_TOTAL_PARTICLES + 1];

                    sum_of_ewald_diffs += (norm(column[NUM_TOTAL_PARTICLES]) - old_pk2) * K_VECTORS[k][3];
                    k++;
                }
            }
        }
    }
    return sum_of_ewald_diffs;
}

void Simulation::set_exp_kr_table_for_water(int water_index) {
    double *coords = WATERS[water_index]->coords;
    int tmp_n, tmp_m, tmp_o;
    
    for (int i = 0; i < 3; i++) {
        tmp_o = (i < 2) ? EWALD_NXY : EWALD_NZ;
        tmp_n = tmp_o + 1;
        tmp_m = tmp_o - 1;

        exp_kr_O[i][tmp_n] = exp(dcomplex(0.0, K_VECTORS[K_111_INDEX][i] * coords[i]));
        exp_kr_H1[i][tmp_n] = exp(dcomplex(0.0, K_VECTORS[K_111_INDEX][i] * coords[i + 3]));
        exp_kr_H2[i][tmp_n] = exp(dcomplex(0.0, K_VECTORS[K_111_INDEX][i] * coords[i + 6]));

        exp_kr_O[i][tmp_m] = COMPLEX_ONE / exp_kr_O[i][tmp_n];
        exp_kr_H1[i][tmp_m] = COMPLEX_ONE / exp_kr_H1[i][tmp_n];
        exp_kr_H2[i][tmp_m] = COMPLEX_ONE / exp_kr_H2[i][tmp_n];

        for (int j = 2; j <= tmp_o; j++) {
            exp_kr_O[i][tmp_o + j] = exp_kr_O[i][tmp_m + j] * exp_kr_O[i][tmp_n];
            exp_kr_O[i][tmp_o - j] = exp_kr_O[i][tmp_n - j] * exp_kr_O[i][tmp_m];

            exp_kr_H1[i][tmp_o + j] = exp_kr_H1[i][tmp_m + j] * exp_kr_H1[i][tmp_n];
            exp_kr_H1[i][tmp_o - j] = exp_kr_H1[i][tmp_n - j] * exp_kr_H1[i][tmp_m];

            exp_kr_H2[i][tmp_o + j] = exp_kr_H2[i][tmp_m + j] * exp_kr_H2[i][tmp_n];
            exp_kr_H2[i][tmp_o - j] = exp_kr_H2[i][tmp_n - j] * exp_kr_H2[i][tmp_m];
        }
    }
    return;
}

void Simulation::set_exp_kr_table_for_ion(int ion_index) {
    int tmp_n, tmp_m, tmp_o;
    double *coords = IONS[ion_index]->coords;
    for (int i = 0; i < 3; i++) {
        tmp_o = (i < 2) ? EWALD_NXY : EWALD_NZ;
        tmp_n = tmp_o + 1;
        tmp_m = tmp_o - 1;

        exp_kr_ion[i][tmp_n] = exp(dcomplex(0.0, K_VECTORS[K_111_INDEX][i] * coords[i]));
        exp_kr_ion[i][tmp_m] = COMPLEX_ONE / exp_kr_ion[i][tmp_n];

        for (int j = 2; j <= tmp_o; j++) {
            exp_kr_ion[i][tmp_o + j] = exp_kr_ion[i][tmp_m + j] * exp_kr_ion[i][tmp_n];
            exp_kr_ion[i][tmp_o - j] = exp_kr_ion[i][tmp_n - j] * exp_kr_ion[i][tmp_m];
        }
    }
    return;
}
