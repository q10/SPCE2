#include "common.h"
#include "sampler.h"

SamplerSet::SamplerSet(Simulation * s) {
    simulation = s;

    SIMULATION_TIME_STAMP = TIMESTAMP();
    RELATIVE_VMD_SNAPSHOT_RATE = 100;

    vmd_timestep = 0;
    vmd_snapshot_counter = 0;
}

SamplerSet::~SamplerSet() {
    if (VMD_FILE.is_open())
        VMD_FILE.close();
    if (CONFIG_FILE.is_open())
        CONFIG_FILE.close();
}

void SamplerSet::sample_data() {
    for (unsigned int i = 0; i < samplers.size(); i++)
        samplers[i]->sample();
    if (++vmd_snapshot_counter % RELATIVE_VMD_SNAPSHOT_RATE == 0)
        write_vmd_snapshot();
    return;
}

void SamplerSet::start() {
    for (unsigned int i = 0; i < samplers.size(); i++)
        samplers[i]->start();
    return;
}

void SamplerSet::finish() {
    for (unsigned int i = 0; i < samplers.size(); i++)
        samplers[i]->finish();
    return;
}

void SamplerSet::write_vmd_snapshot() {
    if (!VMD_FILE.is_open()) {
        if (vmd_filename.compare("") == 0)
            vmd_filename = "SPCE_" + SIMULATION_TIME_STAMP + ".lammpstrj";
        VMD_FILE.open(vmd_filename.c_str());
        ASSERT(VMD_FILE.is_open(), "Could not open VMD output file.");
    }
    VMD_FILE << simulation->to_vmd(++vmd_timestep);
    return;
}

void SamplerSet::write_config_snapshot() {
    if (config_filename.compare("") == 0) {
        if (config_filename.compare("") == 0)
            config_filename = "SPCE_" + SIMULATION_TIME_STAMP + ".config";
        CONFIG_FILE.open(config_filename.c_str());
        ASSERT(CONFIG_FILE.is_open(), "Could not open config output file.");
    }
    CONFIG_FILE << simulation;
    CONFIG_FILE.close();
    return;
}
