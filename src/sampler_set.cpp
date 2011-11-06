#include "common.h"

SamplerSet::SamplerSet(Simulation * s) {
    simulation = s;
    RELATIVE_VMD_SNAPSHOT_RATE = 100;
    vmd_timestep = 0;
    vmd_snapshot_counter = 0;
}

SamplerSet::~SamplerSet() {
    samplers.clear();
    if (VMD_FILE.is_open())
        VMD_FILE.close();
}

void SamplerSet::start() {
    for (unsigned int i = 0; i < samplers.size(); i++)
        samplers[i]->start();
    return;
}

void SamplerSet::sample_data() {
    for (unsigned int i = 0; i < samplers.size(); i++)
        samplers[i]->sample();
    if (++vmd_snapshot_counter % RELATIVE_VMD_SNAPSHOT_RATE == 0)
        write_vmd_snapshot();
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

void SamplerSet::write_vmd_snapshot() {
    if (!VMD_FILE.is_open()) {
        if (vmd_filename.compare("") == 0)
            vmd_filename = simulation->NAME + ".lammpstrj";
        VMD_FILE.open(vmd_filename.c_str());
        ASSERT(VMD_FILE.is_open(), "Could not open VMD output file.");
    }
    VMD_FILE << simulation->to_vmd(++vmd_timestep);
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
