#include "common.h"

Sampler::Sampler(Simulation * s) {
    simulation = s;

    SIMULATION_TIME_STAMP = TIMESTAMP();
    DATA_SAMPLING_RATE = 10;
    RELATIVE_VMD_SNAPSHOT_RATE = 100;

    vmd_timestep = 0;
    vmd_snapshot_counter = 0;
}

Sampler::~Sampler() {
    if (VMD_FILE.is_open())
        VMD_FILE.close();
    if (CONFIG_FILE.is_open())
        CONFIG_FILE.close();
}

void Sampler::sample_data() {
    radial_dist_sample();
    if (++vmd_snapshot_counter % RELATIVE_VMD_SNAPSHOT_RATE == 0)
        write_vmd_snapshot();
    return;
}

void Sampler::start() {
    initialize_radial_dist_sampler();
    return;
}

void Sampler::finish() {
    compute_radial_dist_results();
    return;
}

void Sampler::write_vmd_snapshot() {
    if (!VMD_FILE.is_open()) {
        if (vmd_filename.compare("") == 0)
            vmd_filename = "SPCE_" + SIMULATION_TIME_STAMP + ".lammpstrj";
        VMD_FILE.open(vmd_filename.c_str());
        ASSERT(VMD_FILE.is_open(), "Could not open VMD output file.");
    }
    VMD_FILE << simulation->to_vmd(++vmd_timestep);
    return;
}

void Sampler::write_config_snapshot() {
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
