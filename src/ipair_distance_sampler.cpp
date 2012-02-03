#include "common.h"

IonPairDistanceSampler::IonPairDistanceSampler(System * s) {
    system = s;
    CATION = ANION = NULL;
    for (unsigned int i = 0; i < system->IONS.size(); i++) {
        Ion * ion = system->IONS[i];
        if (ion->charge > 0)
            CATION = ion;
        else if (ion->charge < 0)
            ANION = ion;
    }
    ASSERT((CATION != NULL) && (ANION != NULL), "There needs to be at least one cation and one anion in the system.");
}

IonPairDistanceSampler::~IonPairDistanceSampler() {
    if (ION_PAIR_DISTANCE_FILE.is_open())
        ION_PAIR_DISTANCE_FILE.close();
}

void IonPairDistanceSampler::start() {
    if (!ION_PAIR_DISTANCE_FILE.is_open()) {
        string filename = system->NAME + ".ipair_dist";
        ION_PAIR_DISTANCE_FILE.open(filename.c_str());
        ASSERT(ION_PAIR_DISTANCE_FILE.is_open(), "Could not open config output file.");
    }
    return;
}

void IonPairDistanceSampler::sample() {
    ION_PAIR_DISTANCE_FILE << CATION->distance_from(ANION) << endl;
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
