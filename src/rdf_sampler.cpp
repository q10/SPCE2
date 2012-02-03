#include "common.h"

RDFSampler::RDFSampler(System * s) {
    system = s;
    radial_dist_num_his_bars = 200;
}

RDFSampler::~RDFSampler() {
}

void RDFSampler::start() {
    radial_dist_distance.clear();
    water_water_RDF.clear();
    anion_water_RDF.clear();
    cation_water_RDF.clear();
    ion_ion_RDF.clear();

    num_gr = 0;
    delg = system->BOX_Z_LENGTH / radial_dist_num_his_bars;
    for (int i = 0; i < radial_dist_num_his_bars; i++) {
        radial_dist_distance.push_back(0.0);
        water_water_RDF.push_back(0.0);
        anion_water_RDF.push_back(0.0);
        cation_water_RDF.push_back(0.0);
        ion_ion_RDF.push_back(0.0);
    }
    return;
}

void RDFSampler::sample() {
    num_gr++;
    double dx, dy, dz, dr, *coords, *other_coords;

    // water-water
    for (unsigned int i = 0; i < system->WATERS.size() - 1; i++) {
        for (unsigned int j = i + 1; j < system->WATERS.size(); j++) {
            coords = system->WATERS[i]->coords;
            other_coords = system->WATERS[j]->coords;
            dx = abs(coords[0] - other_coords[0]);
            dy = abs(coords[1] - other_coords[1]);
            dz = abs(coords[2] - other_coords[2]);
            dx -= system->BOX_LENGTH * ROUND(dx / system->BOX_LENGTH);
            dy -= system->BOX_LENGTH * ROUND(dy / system->BOX_LENGTH);
            dz -= system->BOX_Z_LENGTH * ROUND(dz / system->BOX_Z_LENGTH);
            if (dx < system->HALF_BOX_LENGTH and dy < system->HALF_BOX_LENGTH and dz < system->HALF_BOX_Z_LENGTH) {
                dr = sqrt(dx * dx + dy * dy + dz * dz);
                int ig = int(dr / delg);
                water_water_RDF[ig] += 2;
            }
        }
    }

    // ion-water
    for (unsigned int i = 0; i < system->IONS.size(); i++) {
        for (unsigned int j = i; j < system->WATERS.size(); j++) {
            coords = system->IONS[i]->coords;
            other_coords = system->WATERS[j]->coords;
            dx = abs(coords[0] - other_coords[0]);
            dy = abs(coords[1] - other_coords[1]);
            dz = abs(coords[2] - other_coords[2]);
            dx -= system->BOX_LENGTH * ROUND(dx / system->BOX_LENGTH);
            dy -= system->BOX_LENGTH * ROUND(dy / system->BOX_LENGTH);
            dz -= system->BOX_Z_LENGTH * ROUND(dz / system->BOX_Z_LENGTH);
            if (dx < system->HALF_BOX_LENGTH and dy < system->HALF_BOX_LENGTH and dz < system->HALF_BOX_Z_LENGTH) {
                dr = sqrt(dx * dx + dy * dy + dz * dz);
                int ig = int(dr / delg);
                if (system->IONS[i]->charge < 0.0)
                    anion_water_RDF[ig] += 2;
                else
                    cation_water_RDF[ig] += 2;
            }
        }
    }

    // ion-ion
    for (unsigned int i = 0; i < system->IONS.size() - 1; i++) {
        for (unsigned int j = i + 1; j < system->IONS.size(); j++) {
            coords = system->IONS[i]->coords;
            other_coords = system->IONS[j]->coords;
            dx = abs(coords[0] - other_coords[0]);
            dy = abs(coords[1] - other_coords[1]);
            dz = abs(coords[2] - other_coords[2]);
            dx -= system->BOX_LENGTH * ROUND(dx / system->BOX_LENGTH);
            dy -= system->BOX_LENGTH * ROUND(dy / system->BOX_LENGTH);
            dz -= system->BOX_Z_LENGTH * ROUND(dz / system->BOX_Z_LENGTH);
            if (dx < system->HALF_BOX_LENGTH and dy < system->HALF_BOX_LENGTH and dz < system->HALF_BOX_Z_LENGTH) {
                dr = sqrt(dx * dx + dy * dy + dz * dz);
                int ig = int(dr / delg);
                ion_ion_RDF[ig] += 2;
            }
        }
    }
    return;
}

void RDFSampler::finish() {
    // computes radial distribution results
    double r, vb, nid;
    for (int i = 0; i < radial_dist_num_his_bars; i++) {
        r = delg * (i + 0.5);
        radial_dist_distance[i] = r;
        vb = (pow(i + 1, 3.0) - pow(i, 3.0)) * pow(delg, 3.0);
        nid = (4 / 3) * M_PI * vb * Water::STD_DENSITY;
        water_water_RDF[i] /= num_gr * system->WATERS.size() * nid;

        //nid = (4 / 3) * M_PI * vb * ION_DENSITY;
        int num_anions = 0;
        for (unsigned int k = 0; k < system->IONS.size(); k++)
            if (system->IONS[k]->charge < 0.0)
                num_anions++;
        anion_water_RDF[i] /= num_gr * num_anions * nid;
        cation_water_RDF[i] /= num_gr * (system->IONS.size() - num_anions) * nid;
        ion_ion_RDF[i] /= num_gr * system->IONS.size() * nid;
    }
    return;
}

string RDFSampler::results() {
    stringstream rad_dist_results;
    // Format of results is as such:
    // r(Angstroms)     g(r)[O-O]       g(r)[anion-O]   g(r)[cation-O]  g(r)[ion-ion]
    //rad_dist_results << "r(Angstroms)\tg(r)[O-O]\tg(r)[anion-O]\tg(r)[cation-O]\tg(r)[ion-ion]" << endl;
    for (int k = 0; k < radial_dist_num_his_bars; k++)
        rad_dist_results << setprecision(10) << radial_dist_distance[k] << "\t"
            << water_water_RDF[k] << "\t" << anion_water_RDF[k] << "\t" << cation_water_RDF[k] << "\t" << ion_ion_RDF[k] << endl;
    return rad_dist_results.str();
}
