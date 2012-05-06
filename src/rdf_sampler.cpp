#include "common.h"

RDFSampler::RDFSampler(WaterSystem * s) {
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
    // water-water
    for (unsigned int i = 0; i < system->WATERS.size() - 1; i++) {
        for (unsigned int j = i + 1; j < system->WATERS.size(); j++) {
            double dr = system->WATERS[i]->distance_from(system->WATERS[j]);
            int ig = int(dr / delg);
            if (ig < radial_dist_num_his_bars)
                water_water_RDF[ig] += 2;
        }
    }
    // ion-water
    for (unsigned int i = 0; i < system->IONS.size(); i++) {
        for (unsigned int j = i; j < system->WATERS.size(); j++) {
            double dr = system->IONS[i]->distance_from(system->WATERS[j]);
            int ig = int(dr / delg);
            if (ig < radial_dist_num_his_bars) {
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
            double dr = system->IONS[i]->distance_from(system->IONS[j]);
            int ig = int(dr / delg);
            if (ig < radial_dist_num_his_bars)
                ion_ion_RDF[ig] += 2;
        }
    }
    return;
}

void RDFSampler::finish() {
    // computes radial distribution results
    // The RDF is defined as the ratio between the average number density 
    // rho(r) observed at a distance r from an atom and the density at a 
    // distance r from an atom in an ideal gas at the same overall density
    int num_anions = 0;
    for (unsigned int k = 0; k < system->IONS.size(); k++)
        if (system->IONS[k]->charge < 0.0)
            num_anions++;
    int num_cations = system->IONS.size() - num_anions;
    double anion_density = num_anions / system->BOX_VOLUME;
    double cation_density = num_cations / system->BOX_VOLUME;
    double water_density = system->WATERS.size() / system->BOX_VOLUME;

    for (int i = 0; i < radial_dist_num_his_bars; i++) {
        radial_dist_distance[i] = delg * (i + 0.5);
        double shell_volume = (4 / 3) * M_PI * (pow((i + 1) * delg, 3.0) - pow(i*delg, 3.0));

        // need to divide by num_gr and num_waters to get time and number average
        water_water_RDF[i] /= (num_gr * system->WATERS.size() * shell_volume * water_density);
        anion_water_RDF[i] /= (num_gr * (num_anions + system->WATERS.size()) * shell_volume * anion_density);
        cation_water_RDF[i] /= (num_gr * (num_cations + system->WATERS.size()) * shell_volume * cation_density);
        ion_ion_RDF[i] /= (num_gr * system->IONS.size() * shell_volume * anion_density);
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
