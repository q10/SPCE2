#include "common.h"

Sampler::Sampler(Simulation * s) {
    simulation = s;
    initialize_radial_dist_sampler();
}

Sampler::~Sampler() {
}

void Sampler::sample_data() {
    radial_dist_sample();
    return;
}
