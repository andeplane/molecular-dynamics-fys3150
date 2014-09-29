#pragma once
#include <system.h>

class StatisticsSampler
{
public:
    StatisticsSampler();
    ~StatisticsSampler();
    void sample(System *system);
    void sampleKineticEnergy(System *system);
};
