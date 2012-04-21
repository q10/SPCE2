#include "common.h"

IonPairDistanceSampler::IonPairDistanceSampler(WaterSystem * s) {
    system = s;
}

IonPairDistanceSampler::~IonPairDistanceSampler() {
    if (ION_PAIR_DISTANCE_FILE.is_open())
        ION_PAIR_DISTANCE_FILE.close();
}

void IonPairDistanceSampler::start() {
    CATIONS.clear();
    ANIONS.clear();

    for (unsigned int i = 0; i < system->IONS.size(); i++) {
        if (system->IONS[i]->charge > 0)
            CATIONS.push_back(system->IONS[i]);
        else if (system->IONS[i]->charge < 0)
            ANIONS.push_back(system->IONS[i]);
    }
    ASSERT(!CATIONS.empty() and !ANIONS.empty(), "IonPairDistanceSampler: There needs to be at least one cation and one anion in the system before IonPairDistanceSampler can be initialized.");


    if (!ION_PAIR_DISTANCE_FILE.is_open()) {
        string filename = system->NAME + ".ipair_dist";
        ION_PAIR_DISTANCE_FILE.open(filename.c_str());
        ASSERT(ION_PAIR_DISTANCE_FILE.is_open(), "Could not open config output file.");
    }
    return;
}

void IonPairDistanceSampler::sample() {
    for (unsigned int i = 0; i < ANIONS.size(); i++)
        for (unsigned int j = 0; j < ANIONS.size(); j++)
            ION_PAIR_DISTANCE_FILE << ANIONS[i]->distance_from(CATIONS[j]) << endl;
    return;
}

void IonPairDistanceSampler::finish() {
    if (ION_PAIR_DISTANCE_FILE.is_open())
        ION_PAIR_DISTANCE_FILE.close();
    return;
}

string IonPairDistanceSampler::results() {
    stringstream results;
    results << "Ion pair distance sampling results are available in " << system->NAME << ".ipair_dist" << endl;
    return results.str();
}
