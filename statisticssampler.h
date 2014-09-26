#pragma once
#include <system.h>

class StatisticsSampler
{
public:
    StatisticsSampler();
    void sample(System *system);
    void sampleKineticEnergy(System *system);
};
