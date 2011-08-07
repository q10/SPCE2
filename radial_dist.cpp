#include "common.h"

void Sampler::initialize_radial_dist_sampler() {
    radial_dist_num_his_bars = 100;
    radial_dist_distance = new double [radial_dist_num_his_bars];
    radial_dist_data = new double [radial_dist_num_his_bars];
    ion_radial_dist_data = new double [radial_dist_num_his_bars];

    num_gr = 0;
    delg = simulation->BOX_LENGTH / (2 * radial_dist_num_his_bars);
    for (int i = 0; i < radial_dist_num_his_bars; i++)
        radial_dist_data[i] = ion_radial_dist_data[i] = 0;
    return;
}

void Sampler::radial_dist_sample() {
    num_gr++;
    double dx, dy, dz, dr, *coords, *other_coords;

    // water-water
    for (int i = 0; i < simulation->WATERS.size() - 1; i++) {
        for (int j = i + 1; j < simulation->WATERS.size(); j++) {
            coords = simulation->WATERS[i]->coords;
            other_coords = simulation->WATERS[j]->coords;
            dx = abs(coords[0] - other_coords[0]);
            dy = abs(coords[1] - other_coords[1]);
            dz = abs(coords[2] - other_coords[2]);
            dx -= simulation->BOX_LENGTH * ROUND(dx / simulation->BOX_LENGTH);
            dy -= simulation->BOX_LENGTH * ROUND(dy / simulation->BOX_LENGTH);
            dz -= simulation->BOX_Z_LENGTH * ROUND(dz / simulation->BOX_Z_LENGTH);
            dr = sqrt(dx * dx + dy * dy + dz * dz);
            if (dx < simulation->HALF_BOX_LENGTH and dy < simulation->HALF_BOX_LENGTH and dz < simulation->HALF_BOX_Z_LENGTH) {
                int ig = int(dr / delg);
                radial_dist_data[ig] += 2;
            }
        }
    }

    // ion-water
    for (int i = 0; i < simulation->IONS.size(); i++) {
        for (int j = 0; j < simulation->WATERS.size(); j++) {
            coords = simulation->IONS[i]->coords;
            other_coords = simulation->WATERS[j]->coords;
            dx = abs(coords[0] - other_coords[0]);
            dy = abs(coords[1] - other_coords[1]);
            dz = abs(coords[2] - other_coords[2]);
            dx -= simulation->BOX_LENGTH * ROUND(dx / simulation->BOX_LENGTH);
            dy -= simulation->BOX_LENGTH * ROUND(dy / simulation->BOX_LENGTH);
            dz -= simulation->BOX_Z_LENGTH * ROUND(dz / simulation->BOX_Z_LENGTH);
            if (dx < simulation->HALF_BOX_LENGTH and dy < simulation->HALF_BOX_LENGTH and dz < simulation->HALF_BOX_Z_LENGTH) {
                dr = sqrt(dx * dx + dy * dy + dz * dz);
                int ig = int(dr / delg);
                ion_radial_dist_data[ig] += 2;
            }
        }
    }
    return;
}

void Sampler::compute_radial_dist_results() {
    double r, vb, nid;
    for (int i = 0; i < radial_dist_num_his_bars; i++) {
        r = delg * (i + 0.5);
        radial_dist_distance[i] = r;
        vb = (pow(i + 1, 3.0) - pow(i, 3.0)) * pow(delg, 3.0);
        nid = (4 / 3) * M_PI * vb * Water::STD_DENSITY;
        radial_dist_data[i] /= num_gr * simulation->WATERS.size() * nid;

        //nid = (4 / 3) * M_PI * vb * ION_DENSITY;
        ion_radial_dist_data[i] /= num_gr * simulation->IONS.size() * nid;
    }
    return;
}

void Sampler::print_radial_dist_results() {
    cout << "\nr(Angstroms)\tg(r)[O-O]\tg(r)[ion-O]" << endl;
    for (int k = 0; k < radial_dist_num_his_bars; k++)
        cout << setprecision(10) << radial_dist_distance[k] << "\t" << radial_dist_data[k] << "\t" << ion_radial_dist_data[k] << endl;
    return;
}

/*
void test_radial_dist_sampler() {
    cout << "---- BEGIN TEST - RADIAL DISTRIBUTION SAMPLER ----" << endl;

    NUM_EQUILIBRATION_SWEEPS = 0;
    NUM_MC_SWEEPS = 100000;
    initialize();
    run_mc();
    compute_radial_dist_results();
    print_radial_dist_results();

    cout << "\n---- END TEST - RADIAL DISTRIBUTION SAMPLER ----\n" << endl;
    return;
}
 */
