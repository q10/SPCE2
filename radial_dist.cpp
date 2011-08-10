#include "common.h"

void Sampler::initialize_radial_dist_sampler() {
    if (radial_dist_distance != NULL)
        delete [] radial_dist_distance;
    if (water_water_RDF != NULL)
        delete [] water_water_RDF;
    if (ion_water_RDF != NULL)
        delete [] ion_water_RDF;
    if (ion_ion_RDF != NULL)
        delete [] ion_ion_RDF;

    radial_dist_num_his_bars = 100;
    radial_dist_distance = new double [radial_dist_num_his_bars];
    water_water_RDF = new double [radial_dist_num_his_bars];
    ion_water_RDF = new double [radial_dist_num_his_bars];
    ion_ion_RDF = new double [radial_dist_num_his_bars];

    num_gr = 0;
    delg = simulation->BOX_LENGTH / (2 * radial_dist_num_his_bars);
    for (int i = 0; i < radial_dist_num_his_bars; i++)
        water_water_RDF[i] = ion_water_RDF[i] = ion_ion_RDF[i] = 0;
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
            if (dx < simulation->HALF_BOX_LENGTH and dy < simulation->HALF_BOX_LENGTH and dz < simulation->HALF_BOX_Z_LENGTH) {
                dr = sqrt(dx * dx + dy * dy + dz * dz);
                int ig = int(dr / delg);
                water_water_RDF[ig] += 2;
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
                ion_water_RDF[ig] += 2;
            }
        }
    }

    // ion-ion
    for (int i = 0; i < simulation->IONS.size() - 1; i++) {
        for (int j = i + 1; j < simulation->IONS.size(); j++) {
            coords = simulation->IONS[i]->coords;
            other_coords = simulation->IONS[j]->coords;
            dx = abs(coords[0] - other_coords[0]);
            dy = abs(coords[1] - other_coords[1]);
            dz = abs(coords[2] - other_coords[2]);
            dx -= simulation->BOX_LENGTH * ROUND(dx / simulation->BOX_LENGTH);
            dy -= simulation->BOX_LENGTH * ROUND(dy / simulation->BOX_LENGTH);
            dz -= simulation->BOX_Z_LENGTH * ROUND(dz / simulation->BOX_Z_LENGTH);
            if (dx < simulation->HALF_BOX_LENGTH and dy < simulation->HALF_BOX_LENGTH and dz < simulation->HALF_BOX_Z_LENGTH) {
                dr = sqrt(dx * dx + dy * dy + dz * dz);
                int ig = int(dr / delg);
                ion_ion_RDF[ig] += 2;
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
        water_water_RDF[i] /= num_gr * simulation->WATERS.size() * nid;

        //nid = (4 / 3) * M_PI * vb * ION_DENSITY;
        ion_water_RDF[i] /= num_gr * simulation->IONS.size() * nid;
        ion_ion_RDF[i] /= num_gr * simulation->IONS.size() * nid;
    }
    return;
}

string Sampler::radial_dist_results() {
    stringstream rad_dist_results;
    rad_dist_results << "\nr(Angstroms)\tg(r)[O-O]\tg(r)[ion-O]\tg(r)[ion-ion]" << endl;
    for (int k = 0; k < radial_dist_num_his_bars; k++)
        rad_dist_results << setprecision(10) << radial_dist_distance[k] << "\t"
            << water_water_RDF[k] << "\t" << ion_water_RDF[k] << "\t" << ion_ion_RDF[k] << endl;
    return rad_dist_results.str();
}

void test_radial_dist() {
    //cout << "---- BEGIN TEST - RADIAL DISTRIBUTION SAMPLER ----" << endl;

    Simulation * simulation = new Simulation();
    simulation->NUM_MC_SWEEPS = 10;
    simulation->run_mc();
    cout << simulation->sampler->radial_dist_results();
    //cout << "\n---- END TEST - RADIAL DISTRIBUTION SAMPLER ----\n" << endl;
    return;
}
