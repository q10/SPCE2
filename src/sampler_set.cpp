#include "common.h"

SamplerSet::SamplerSet(Simulation * s) {
    simulation = s;
    is_lammpstr_sampling = false;
    RELATIVE_LAMMPSTRJ_SNAPSHOT_RATE = 10;
    lammpstrj_timestep = 0;
    lammpstrj_snapshot_counter = 0;
}

SamplerSet::~SamplerSet() {
    samplers.clear();
    if (LAMMPSTRJ_FILE.is_open())
        LAMMPSTRJ_FILE.close();
}

void SamplerSet::start() {
    for (unsigned int i = 0; i < samplers.size(); i++)
        samplers[i]->start();
    return;
}

void SamplerSet::sample_data() {
    for (unsigned int i = 0; i < samplers.size(); i++)
        samplers[i]->sample();
    if (is_lammpstr_sampling and ++lammpstrj_snapshot_counter % RELATIVE_LAMMPSTRJ_SNAPSHOT_RATE == 0)
        write_lammpstrj_snapshot();
    return;
}

void SamplerSet::finish() {
    for (unsigned int i = 0; i < samplers.size(); i++)
        samplers[i]->finish();
    return;
}

void SamplerSet::print_individual_sampler_results() {
    for (unsigned int i = 0; i < samplers.size(); i++)
        cout << samplers[i]->results() << endl;
    return;
}

void SamplerSet::turn_on_lammpstrj_sampler() {
    is_lammpstr_sampling = true;
    return;
}

void SamplerSet::add_rdf_sampler() {
    RDFSampler * sampler = new RDFSampler(simulation);
    samplers.push_back(dynamic_cast<Sampler *> (sampler));
    return;
}

void SamplerSet::add_ion_pair_distance_sampler() {
    IonPairDistanceSampler * sampler = new IonPairDistanceSampler(simulation);
    samplers.push_back(dynamic_cast<Sampler *> (sampler));
    return;
}

void SamplerSet::write_lammpstrj_snapshot() {
    if (!LAMMPSTRJ_FILE.is_open()) {
        if (lammpstrj_filename.compare("") == 0)
            lammpstrj_filename = simulation->NAME + ".lammpstrj";
        LAMMPSTRJ_FILE.open(lammpstrj_filename.c_str());
        ASSERT(LAMMPSTRJ_FILE.is_open(), "Could not open LAMMPSTRJ output file.");
    }
    LAMMPSTRJ_FILE << simulation->to_lammpstrj(++lammpstrj_timestep);
    return;
}

void SamplerSet::write_config_snapshot() {
    ofstream config_file;
    if (config_filename.compare("") == 0)
        config_filename = simulation->NAME + ".config";
    config_file.open(config_filename.c_str());
    ASSERT(config_file.is_open(), "Could not open config output file.");
    config_file << simulation;
    config_file.close();
    return;
}
